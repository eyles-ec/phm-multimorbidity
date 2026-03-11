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
impute_knn <- function(var, bnssg, coords, k = 5) {
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

#function for maps

#inla function


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
  "Resident_pop_total",
  "chn2024",
  "IMD25",
  "swd_median_age"
)

outcomes <- c("swd_pct_seg4_5", "swd_mean_cms", "swd_median_cms")

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

#calculate a spatial lag for each outcome
bnssg_clean$lag_swd_pct_seg4_5 <- lag.listw(lw, bnssg_clean$swd_pct_seg4_5)
bnssg_clean$lag_swd_mean_cms <- lag.listw(lw, bnssg_clean$swd_mean_cms)
bnssg_clean$lag_swd_median_cms <- lag.listw(lw, bnssg_clean$swd_median_cms)

#drop geometry for modelling, only necessary if running as a pipeline
bnssg_clean <- st_drop_geometry(bnssg_clean)

predictors_only <- bnssg_clean %>% 
  select(-all_of(outcomes))

#build model matrix, specifying dummy variables (year, local authority)
X <- model.matrix(
  ~ . + factor(year) + factor(LAD22NM),
  data = predictors_only
)


#sanitise names for later use in INLA in model matrix
colnames(X) <- make.names(colnames(X))

#vectors for each outcome 
Y1 <- bnssg_clean$swd_pct_seg4_5
Y2 <- bnssg_clean$swd_mean_cms
Y3 <- bnssg_clean$swd_median_cms

#fit LASSO models (alpha = 1), data are too multicollinear
# lasso_1 <- cv.glmnet(X, Y1, alpha = 1)
# lasso_2 <- cv.glmnet(X, Y2, alpha = 1)
# lasso_3 <- cv.glmnet(X, Y3, alpha = 1)
# 
#fit Ridge models (alpha = 0), doesn't select variables
# ridge_1 <- cv.glmnet(X, Y1, alpha = 0)
# ridge_2 <- cv.glmnet(X, Y2, alpha = 0)
# ridge_3 <- cv.glmnet(X, Y3, alpha = 0)

#fit Elastic Net models (alpha = 0.5), more appropriate as variables are likely multicollinear in a grouped fashion
#can perform variable selection here 
enet_1 <- cv.glmnet(X, Y1, alpha = 0.5)
enet_2 <- cv.glmnet(X, Y2, alpha = 0.5)
enet_3 <- cv.glmnet(X, Y3, alpha = 0.5)

#calc resids
fitted_enet1  <- predict(enet_1,  newx = X, s = "lambda.min")
fitted_enet2  <- predict(enet_2,  newx = X, s = "lambda.min")
fitted_enet3  <- predict(enet_3,  newx = X, s = "lambda.min")

bnssg$resid_enet1 <- as.numeric(Y1 - fitted_enet1)
bnssg$resid_enet2 <- as.numeric(Y2 - fitted_enet2)
bnssg$resid_enet3  <- as.numeric(Y1 - fitted_enet3)

#neighbourlist, binary spatial weights
nb <- poly2nb(bnssg, queen = TRUE)
lw <- nb2listw(nb, style = "B")

#calc local getis ord on residuals using lw and map 
bnssg$gi_enet1 <- as.numeric(localG(bnssg$resid_enet1, lw))
bnssg$gi_enet2 <- as.numeric(localG(bnssg$resid_enet2, lw))
bnssg$gi_enet3  <- as.numeric(localG(bnssg$resid_enet3,  lw))


#tmap resids and local gi 
tmap_mode("plot")

tm1 <- tm_shape(bnssg) +
  tm_fill("swd_pct_seg4_5",
          palette = "viridis",
          style = "quantile",
          title = "% in CMS 4-5") +
  tm_borders()

tm2 <- tm_shape(bnssg) +
  tm_fill("resid_enet1",
          palette = "-RdBu",
          style = "quantile",
          title = "Elastic Net residuals") +
  tm_borders()

tm3 <- tm_shape(bnssg) +
  tm_fill("gi_enet1",
          palette = "-RdBu",
          style = "quantile",
          title = "Gi* (hot/cold spots)") +
  tm_borders()

tmap_arrange(tm1, tm2, tm3, ncol = 3)

#inla models as residuals still show clustering

#calc new matrix weights (not list) for inla
inlaw  <- nb2mat(nb, style = "B")
inla.write.graph(inlaw, "bnssg.graph")

#index for areas
bnssg$id <- 1:nrow(bnssg)

#pull out enet covars
coef_enet1 <- coef(enet_1, s = "lambda.min")
selected_vars1 <- rownames(coef_enet1)[coef_enet1[,1] != 0]
selected_vars1 <- setdiff(selected_vars1, "(Intercept)")

#sanitise selected_vars to match model matrix X from enet
selected_vars1 <- make.names(selected_vars1)


#collapse variable names into a single string for selected variables from enet
rhs <- paste(selected_vars1, collapse = " + ")

#add the BYM2 spatial random effect to model
rhs <- paste(rhs, "+ f(id, model = 'bym2', graph = 'bnssg.graph')")

#build the inla formula
formula_inla <- as.formula(paste("swd_pct_seg4_5 ~", rhs))

X_inla <- model.matrix(
  ~ . + factor(year) + factor(LAD22NM),
  data = bnssg_clean
)

colnames(X_inla) <- make.names(colnames(X_inla))

bnssg <- cbind(bnssg, X_inla)

rhs <- paste(selected_vars1, collapse = " + ")
rhs <- paste(rhs, "+ f(id, model = 'bym2', graph = 'bnssg.graph')")
formula_inla <- as.formula(paste("swd_pct_seg4_5 ~", rhs))

#formula with selected covariates from elastic net
formula_inla <- swd_pct_seg4_5 ~ 
  var1 + var2 + var3 +   # <- replace with your selected vars
  f(id, model = "bym2", graph = "bnssg.graph")


res_inla <- inla(
  formula_inla,
  family = "gaussian",
  data   = as.data.frame(bnssg),
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

bnssg$fitted_inla  <- res_inla$summary.fitted.values$mean
bnssg$resid_inla   <- bnssg$swd_pct_seg4_5 - bnssg$fitted_inla
