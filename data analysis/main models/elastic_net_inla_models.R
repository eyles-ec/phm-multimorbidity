# Load packages
library(glmnet)
library(dplyr)
library(sf)
library(sfdep)
library(spdep)
library(FNN)
library(INLA)
library(tmap)


#inla has to be installed like this:
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#if it doesn't have Rgraphviz and graph then use this:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)

#Function to impute one variable using k nearest neighbour, distance weighting, error catch if no missing
impute_knn <- function(var, 
                       bnssg, 
                       coords, 
                       k = 5) {
  missing_idx <- which(is.na(bnssg[[var]]))
  present_idx <- which(!is.na(bnssg[[var]]))
  
  if (length(missing_idx) == 0) return(bnssg[[var]])
  
  nn <- get.knnx(coords[present_idx, ], coords[missing_idx, ], k = k)
  
  for (i in seq_along(missing_idx)) {
    neighbours <- present_idx[nn$nn.index[i, ]]
    distances <- nn$nn.dist[i, ]
    weights <- 1 / (distances + 1e-6)
    bnssg[[var]][missing_idx[i]] <- sum(weights * bnssg[[var]][neighbours]) / sum(weights)
  }
  
  bnssg[[var]]
}

#function for enet models
enet <- function(df,
                 outcome,
                 lw,                
                 exclude_vars = NULL,
                 alpha = 0.5,
                 nfolds = 10,
                 seed = 123) {
  
  df <- as.data.frame(df)
  
  #remove existing lag columns if they exist
  df <- df |> dplyr::select(-dplyr::matches("^lag_"))
  
  #by outcome spatial lag calculation 
  lag_name <- paste0("lag_", outcome)
  df[[lag_name]] <- as.numeric(spdep::lag.listw(lw, df[[outcome]]))
  
  #exclude outcomes from predictors
  predictors_only <- df |> dplyr::select(-dplyr::all_of(outcome))
  
  #build model matrix, specifying dummy variables (year, local authority)
  X <- model.matrix(
    ~ . + factor(year) + factor(LAD22NM),
    data = predictors_only
  )
  
  #sanitise column names for inla
  colnames(X) <- make.names(colnames(X))
  
  #fit the elastic net model
  set.seed(seed)
  enet <- glmnet::cv.glmnet(
    x = X,
    y = df[[outcome]],
    alpha = alpha,
    nfolds = nfolds,
    family = "gaussian"
  )
  
  #calculate residuals
  fitted <- as.numeric(predict(enet, newx = X, s = "lambda.min"))
  resid_enet <- paste0("resid_", outcome)
  df[[resid_enet]] <- df[[outcome]] - fitted
  
  #pull out coefficients for inla, remove intercept
  coef_mat <- coef(enet, s = "lambda.min")
  selected <- rownames(coef_mat)[coef_mat[,1] != 0]
  selected <- selected[selected != "(Intercept)" & selected != "Intercept"]
  
  #sanitise and remove intercept in case it persists 
  selected <- make.names(selected)
  selected <- selected[selected != "X.Intercept."]
  
  #remove dummy-coded factor variables (LAD22NM and year)
  selected <- selected[!grepl("^LAD22NM", selected)] 
  selected <- selected[!grepl("^factor\\.", selected)]
  selected <- selected[!grepl("^year", selected)]
  
                        
  #return everything in a list
  return(list(
    df = df,
    X = X,
    enet_model = enet,
    selected_vars = selected,
    lambda.min = enet$lambda.min,
    lambda.1se = enet$lambda.1se,
    lag_col = lag_name,
    resid = resid_enet
  ))
}

#function for maps. call using sf_df (spatial dataframe), the outcome (map_col), the local gi (gi+col)
#give the first one (outcome map) title with map_title =, defaults to map_col if none provided
#resid_type allows for dynamic titling of residuals map (map 2), can set it to "enet" for elastic net, "inla" for inla
#if nothing, it is set to just 'residuals.' Put a file path and name in save_png to export the map
#makes a three panel map with the outcome, its modelled residuals, and then the gi of residuals (hot/cold spots)
#the gi of residuals represents unmodelled geographic autocorrelation/clustering remaining post-model
gi_tmap <- function(sf_df, 
                    map_col, 
                    resid_col, 
                    gi_col, 
                    map_title = NULL, 
                    resid_type = NULL, 
                    save_png = NULL) {
  
  if (is.null(map_title)) {
    map_title <- map_col
  }
  
  resid_title <- dplyr::case_when(
    resid_type == "enet" ~ "Elastic Net residuals",
    resid_type == "inla" ~ "INLA residuals",
    TRUE ~ "Residuals"
  )
  
  tmap_mode("plot")
 
  tm1 <- tm_shape(sf_df) +
    tm_fill(map_col,
            palette = "viridis",
            style = "quantile",
            title = map_title) +
    tm_borders()
  
  tm2 <- tm_shape(sf_df) +
    tm_fill(resid_col,
            palette = "-RdBu",
            style = "quantile",
            title = resid_title) +
    tm_borders()
  
  tm3 <- tm_shape(sf_df) +
    tm_fill(gi_col,
            palette = "-RdBu",
            style = "quantile",
            title = "Gi* (hot/cold spots)") +
    tm_borders()
  
  map <- tmap_arrange(tm1, tm2, tm3, ncol = 3)
  
  if (!is.null(save_png)) {
    tmap_save(map, filename = save_png, dpi = 300, width = 10, height = 4)
    message("Saved map to: ", save_png)
  }
  
  map
}

#inla function, needs sf_df (spatial dataframe), outcome as a string, selected variables vector from enet (ENET_OUTPUT$selected_vars)
#nb, the neighbour list we generated earlier
inla_model <- function(sf_df, 
                       outcome, 
                       selected, 
                       nb){
  
  #adjacency matrix/weight for inla, makes an inla 'graph'
  inlaw <- nb2mat(nb, style = "B")
  inla.write.graph(inlaw, "bnssg.graph")
  
  #area index
  sf_df$id <- 1:nrow(sf_df)
  
  #ensure selected variables are sanitary for use in inla formula
  selected <- unique(selected)
  selected <- setdiff(selected, c("(Intercept)", "Intercept", "X.Intercept."))
  selected <- selected[!grepl("^factor\\.", selected)]
  selected <- selected[!grepl("^LAD22NM", selected)]
  selected <- setdiff(selected, "year")
  
  #drop any missing variables
  selected <- intersect(selected, names(sf_df))

  #add BYM2 spatial random effect to model and build formula, collapse it into a string
  #BYM2 decomposes latent effects into weighted sums of independent AND spatial effects (structured and unstructured)
  #paper: https://arxiv.org/pdf/1601.1180
  #other paper which is useful: https://journals.sagepub.com/doi/10.1177/09622802241293776
  predictors <- paste(c(selected, 
                      "factor(year)", 
                      "factor(LAD22NM)", 
                      "f(id, model='bym2', graph = 'bnssg.graph')"),
                      collapse = " + ")
  
  #build inla formula
  formula_inla <- as.formula(sprintf("%s ~ %s", outcome, predictors))
  
  cat("INLA formula:\n", deparse(formula_inla), "\n")
  
  model_inla <- inla(
    formula_inla,
    family = "gaussian",
    data = as.data.frame(sf_df),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  sf_df$fitted_inla <- model_inla$summary.fitted.values$mean
  sf_df$resid_inla <- sf_df[[outcome]] - sf_df$fitted_inla
  
  #Bayesian p values for fixed effects (for paper: two-sided posterior tail probability)
  marginals <- model_inla$marginals.fixed
  fixed_sum <- model_inla$summary.fixed
  
  bayes_p <- sapply(marginals, function(m) {
    p_pos <- INLA::inla.pmarginal(0, m)
    2 * min(p_pos, 1 - p_pos)
  })
  
  fixed_sum$bayes_p_value <- bayes_p

  return(list(
    data       = sf_df,
    inla_model = model_inla,
    formula    = formula_inla,
    summary_fixed   = fixed_sum
  ))
  
}

#inla postestimation function 
inla_postestimation <- function(inla_out, outcome, outcome_type = c("score","change"), threshold = 10, nb) {
    
    model <- inla_out$inla_model
    sf_df <- inla_out$data
    
    #fixed effects
    fixed_effects <- inla_out$summary_fixed
    
    #Hyperparameters (BYM2 diagnostics)
    #phi (mixing parameter - proportion of random effect variance that is spatially structured) 0 - spatial noise, 0.50 noise/structured equally, 1 entirely spatially structured
    #precision - (1/variance) high -low residual variability, e.g how much remaining variation is there
    #high phi, high precision - strong/smooth spatial paterning, good model fit
    #high phi, low precision - strong spatial pattern but high heterogeneity
    #phi is not causal! or a goodness of fit!
    hyperparameters <- model$summary.hyperpar
    
    #Fitted values & uncertainty of these
    fitted <- model$summary.fitted.values |>
      dplyr::mutate(
        LSOA21CD = sf_df$LSOA21CD,
        year = sf_df$year,
        observed = sf_df[[outcome]],
        
        residual = observed - mean,
        
        #stability/precision of mean estimate
        ci_width = `0.975quant` - `0.025quant`
      ) |>
      dplyr::rename(
        fitted_mean = mean, 
        fitted_sd   = sd,
        fitted_lci  = `0.025quant`,
        fitted_uci  = `0.975quant`
      )
    if (outcome_type == "change") {
      #Given the data and the model, how likely is it that year-on-year change is positive or negative?
      fitted <- fitted |>
        dplyr::mutate(
          p_positive_change = 1 - pnorm(0, fitted_mean, fitted_sd), # if this is close to 1, evidence of increase (e.g. increase in seg4-5)
          p_negative_change = pnorm(0, fitted_mean, fitted_sd)#if this is close to 1, evidence of decrease (e.g. decrease in seg 4-5)
        )
      
    } else if (outcome_type == "score") {
      #Given the data and the model, how likely is it that the score is above or below the threshold (e.g. 10%)?
      fitted <- fitted |>
        dplyr::mutate(
          p_above_threshold = 1 - pnorm(threshold, fitted_mean, fitted_sd), #close to 1, evidence CMS 4-5 is higher than 10%
          p_below_threshold = pnorm(threshold, fitted_mean, fitted_sd)#close to 1, evidence CMS 4-5 is lower than 10%
        )
    
    } else {
      stop("outcome type incorrectly specified (must be change or score")
    }
    
  #Residual spatial autocorrelation using moran's i
  #looks at whether there's unexplained spatial autocorrelation
  #close to 0 and p > 0.05, resids are spatially uncorrelated
  #if not, some structure may be unexplained due to missing covariates
  moran_res <- moran.test(fitted$residual, nb2listw(nb, style = "W"))
  
  return(list(
    fixed_effects = fixed_effects,
    hyperparameters = hyperparameters,
    fitted = fitted,
    moran_residuals = moran_res
  ))
}

#inla postestimation mapping function
inla_post_map <- function(
    sf_df,
    fitted_mean,
    outcome,
    exceed_prob,
    save_png = NULL,
    title_prefix = "INLA"
) {
  
  tmap_mode("plot")
  
  
  tm1 <- tm_shape(sf_df) +
    tm_polygons(
      fill = fitted_mean,
      col  = NULL,
      fill.scale = tm_scale_intervals(
        style = "quantile",
        n = 5,
        values = "viridis"
      ),
      fill.legend = tm_legend(
        title = "Posterior mean"
      )
    )
  
  tm2 <- tm_shape(sf_df) +
    tm_polygons(
      fill = outcome,
      col  = NULL,
      fill.scale = tm_scale_intervals(
        style = "quantile",
        n = 5,
        values = "viridis"
      ),
      fill.legend = tm_legend(
        title = "Observed values"
      )
    )
  
  tm3 <- tm_shape(sf_df) +
    tm_polygons(
      fill = exceed_prob,
      col  = NULL,
      fill.scale = tm_scale_intervals(
        style = "fixed",
        breaks = c(0, 0.2, 0.5, 0.8, 1),
        values = "plasma"
      ),
      fill.legend = tm_legend(
        title = "Exceedance probability"
      )
    )
  
  map <- tmap_arrange(tm1, tm2, tm3, ncol = 3)
  
  if (!is.null(save_png)) {
    tmap_save(map, save_png, dpi = 300, width = 10, height = 4)
  }
  
  map
}

#results export function
export_model_results <- function(
    enet_model,
    inla_out,
    inla_post,
    outcome_name,
    output_dir = "."
) {
  
  #enet coefficients
  enet_coef <- coef(enet_model$enet_model, s = "lambda.min")
  
  enet_df <- data.frame(
    term = rownames(enet_coef),
    enet_coef = as.numeric(enet_coef)
  ) |>
    dplyr::filter(enet_coef != 0)
  
  #fixed effects from inla
  inla_df <- inla_out$summary_fixed |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::select(
      term,
      mean,
      sd,
      `0.025quant`,
      `0.975quant`,
      bayes_p_value
    )
  
  #join up enet and inla model results
  model_results <- dplyr::full_join(enet_df, inla_df, by = "term")
  
  #fitted values from postestimation
  fitted_df <- inla_post$fitted
  
  #hyperparameters from postestimation
  hyper_df <- inla_post$hyperparameters |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter")
  
  #moran's I on residuals from inla
  moran <- inla_post$moran_residuals
  
  moran_df <- data.frame(
    moran_I     = as.numeric(moran$estimate[["Moran I statistic"]]),
    expected_I  = as.numeric(moran$estimate[["Expectation"]]),
    variance    = as.numeric(moran$estimate[["Variance"]]),
    z_value     = as.numeric(moran$statistic),
    p_value     = as.numeric(moran$p.value),
    alternative = moran$alternative,
    method      = moran$method,
    data_name   = moran$data.name,
    stringsAsFactors = FALSE
  )
  
  #write CSV
  write.csv(
    model_results,
    file = file.path(output_dir, paste0("enet_inla_fixed_", outcome_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    fitted_df,
    file = file.path(output_dir, paste0("inla_fitted_", outcome_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    hyper_df,
    file = file.path(output_dir, paste0("inla_hyperparameters_", outcome_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    moran_df,
    file = file.path(output_dir, paste0("moran_residuals_", outcome_name, ".csv")),
    row.names = FALSE
  )
  
  #save all the models 
  saveRDS(
    enet_model,
    file = file.path(output_dir, paste0("enet_model_", outcome_name, ".rds"))
  )
  
  saveRDS(
    inla_out,
    file = file.path(output_dir, paste0("inla_model_", outcome_name, ".rds"))
  )

  invisible(list(
    fixed_effects = model_results,
    fitted        = fitted_df,
    hyperparameters = hyper_df
  ))
}


#pull wd from paths.R (put in .gitignore)
source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset
bnssg_csv <- read.csv("./BNSSG/linked/bnssg_long.csv")

#pull out LSOAS
bnssg_lsoas <- bnssg_csv$LSOA21CD

#link in LSOAs
#load lsoa data and then restrict to BNSSG LSOAS by filtering out LSOA codes not contained in bnssg

england <- st_read(file.path("./LSOA 2021/LSOA/LSOA_2021_EW_BGC.shp"))
#Projected CRS: OSGB36 / British National Grid

bnssg_map <- england[england$LSOA21CD %in% bnssg_lsoas,]

#save memory by removing english sf data
rm(england)

#join bnssg csv to map
bnssg <- bnssg_map %>% left_join(bnssg_csv, by = "LSOA21CD")

#Variables to exclude from the model matrix
exclude_vars <- c(
  "LSOA21CD",
  "LSOA21NM.x",
  "GlobalID",
  "LSOA21NM.y",
  "CHGIND",
  "population",
  "good_health",
  "fair_health",
  "poor_health",
  "Registered_pop_total",
  "total_road_m",
  "n_intersections",
  "area_km2",
  "IMD25",
  "swd_median_age",
  "swd_mean_cms", #rm for now
  "swd_median_cms", #rm for now
  "swd_pct_seg4_5_abs_change",
  "swd_mean_cms_abs_change" ,     
  "swd_mean_cms_yoy_change" 
)

#impute pollution with k nearest neighbours
pollution_vars <- c(
  "benzene_2024_mean_ugm3",
  "nox_2024_mean_ugm3",
  "ozone_2024days_over_120ugm3",
  "pm10_2024_mean_ugm3",
  "pm2_5_2024_mean_ugm3",
  "so2_2024_mean_ugm3"
)

#calculate centroids (for knn)
coords <- st_coordinates(st_centroid(bnssg))

#call the function over all the pollution variables
for (v in pollution_vars) {
  bnssg[[v]] <- impute_knn(v, bnssg, coords, k = 5)
}

# Remove excluded predictor variables
bnssg_clean <- bnssg %>%
  select(-any_of(exclude_vars))

#fix any remaining NAs (land use) by setting to 0 (NA means no presence, i.e. 0)
bnssg_clean <- bnssg_clean %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x)))

#neighbourlist, binary spatial weights, bnssg is functionally the same as bnssg_clean here
nb <- poly2nb(bnssg, queen = TRUE)
lw <- nb2listw(nb, style = "B")

#drop geometry for modelling, only necessary if running as a pipeline
bnssg_clean <- st_drop_geometry(bnssg_clean)

#remove other outcome
bnssg_score <- bnssg_clean |>
  select(
    -swd_pct_seg4_5_yoy_change
  )

bnssg_change <- bnssg_clean |>
  select(
    -swd_pct_seg4_5
  )

#elastic net model
enet_model_score <- enet(bnssg_score, "swd_pct_seg4_5", lw)
enet_model_change <- enet(bnssg_change, "swd_pct_seg4_5_yoy_change", lw)

#pull out enet names
resid_col_score <- enet_model_score$resid
lag_col_score   <- enet_model_score$lag_col

resid_col_change <- enet_model_change$resid
lag_col_change <- enet_model_change$lag_col

#add residual and lag back to sf object
bnssg[[resid_col_score]] <- enet_model_score$df[[resid_col_score]]
bnssg[[lag_col_score]]   <- enet_model_score$df[[lag_col_score]]
bnssg[[resid_col_change]] <- enet_model_change$df[[resid_col_change]]
bnssg[[lag_col_change]]   <- enet_model_change$df[[lag_col_change]]

#create container for gi 
gi_col_score <- paste0("gi_", resid_col_score)
gi_col_change <- paste0("gi_", resid_col_change)

#calculate local gi
bnssg[[gi_col_score]] <- as.numeric(localG(bnssg[[resid_col_score]], lw))
bnssg[[gi_col_change]] <- as.numeric(localG(bnssg[[resid_col_change]], lw))

#generate three panel map, save using function
enet_map_score <- gi_tmap(bnssg, "swd_pct_seg4_5", "resid_swd_pct_seg4_5", "gi_resid_swd_pct_seg4_5", map_title = "Proportion in CMS 4-5", resid_type = "enet", save_png = "enet_score.png")
enet_map_change <- gi_tmap(bnssg, "swd_pct_seg4_5_yoy_change", "resid_swd_pct_seg4_5_yoy_change", "gi_resid_swd_pct_seg4_5_yoy_change", map_title = "Yearly Change in CMS 4-5", resid_type = "enet", save_png = "enet_change.png")

#create inla models as residuals still show clustering

#subset out outcomes
bnssg_score <- bnssg|>
  select(
    -swd_pct_seg4_5_yoy_change
  )

bnssg_change <- bnssg |>
  select(
    -swd_pct_seg4_5
  )

#run inla function
inla_score <- inla_model(bnssg_score, "swd_pct_seg4_5", enet_model_score$selected_vars, nb)
inla_change <- inla_model(bnssg_change, "swd_pct_seg4_5_yoy_change", enet_model_change$selected_vars, nb)

#pull out inla names
resid_col_inla_score <- "resid_inla_score"
gi_col_inla_score    <- "gi_resid_inla_score"
resid_col_inla_change <- "resid_inla_change"
gi_col_inla_change <- "gi_resid_inla_change"

#add residual and local gi
bnssg[[resid_col_inla_score]] <- inla_score$data$resid_inla
bnssg[[gi_col_inla_score]] <- as.numeric(localG(bnssg[[resid_col_inla_score]], lw))
bnssg[[resid_col_inla_change]] <- inla_change$data$resid_inla
bnssg[[gi_col_inla_change]] <- as.numeric(localG(bnssg[[resid_col_inla_change]], lw))


#generate three panel map, save
inla_map_score <- gi_tmap(bnssg, "swd_pct_seg4_5", resid_col_inla_score, gi_col_inla_score, map_title = "Proportion in CMS 4-5", resid_type = "inla", save_png = "inla_score.png")
inla_map_change <- gi_tmap(bnssg, "swd_pct_seg4_5_yoy_change", resid_col_inla_change, gi_col_inla_change, map_title = "Change in CMS 4-5", resid_type = "inla", save_png = "inla_change.png")

#other inla postestimation
inla_post_score <- inla_postestimation(inla_score, "swd_pct_seg4_5", outcome_type = "score",threshold=10, nb)
inla_post_change<- inla_postestimation(inla_change, "swd_pct_seg4_5_yoy_change", outcome_type = "change", threshold = 0, nb)

#join fitted values to inla data by LSOA and year
sf_map_score <- inla_score$data |>
  dplyr::left_join(
    inla_post_score$fitted,
    by = c("LSOA21CD", "year")
  )

sf_map_change <- inla_change$data |>
  dplyr::left_join(
    inla_post_change$fitted,
    by = c("LSOA21CD", "year")
  )

inla_post_score_map <-inla_post_map(
                      sf_df = sf_map_score,    
                      fitted_mean  = "fitted_mean",
                      outcome   = "observed",
                      exceed_prob  = "p_above_threshold",
                      save_png = "inla_post_score.png")

inla_post_score_map <-inla_post_map(
                        sf_df = sf_map_change,    
                        fitted_mean  = "fitted_mean",
                        outcome   = "observed",
                        exceed_prob  = "p_positive_change",
                        save_png = "inla_post_change.png")

#export everything using function 
dir.create("results", showWarnings = FALSE)
export_model_results(enet_model_score, inla_score, inla_post_score, outcome_name = "swd_pct_seg4_5", output_dir = "./results")
export_model_results(enet_model_change, inla_change, inla_post_change, outcome_name = "swd_pct_seg4_5_yoy_change", output_dir = "./results")

