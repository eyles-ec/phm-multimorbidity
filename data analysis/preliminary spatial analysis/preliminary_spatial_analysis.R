library(sf)
library(sfdep)
library(spdep)
library(dplyr)
library(tidyr)
library(ggplot2)

#Getis-ord analysis partially following this tutorial: https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut 

getis_ord <- function(
    sf_data,
    outcome,
    nsim = 5000,
    title = NULL
) {
  
  outcome_name <- rlang::as_name(rlang::ensym(outcome))
  
  #create a list of neighbours, queen is a type of contiguity criteria which is when a common edge/vertex is shared
  neighbour_list <- poly2nb(sf_data, queen = TRUE)
  
  #check for empty neighbour sets (e.g. LSOAs without any neighbours)
  #analysis may not run correctly if there are any. in this case, none
  if (any(card(neighbour_list) == 0)) {
    stop("Some areas have no neighbours; analysis may not run correctly.")
  }
  
  #assign binary weights, 1 to all neighbouring features, 0 to other features
  neighbour_weights <- nb2listw(neighbour_list, style = "B")
  
  #calculate global G statistic (Getis-Ord global G)
  global_go <- globalG.test(
    sf_data[[outcome_name]],
    neighbour_weights
  )
  
  ##calculate spatial lag and get local data
  sf_local <- sf_data %>%
    mutate(
      nbours  = st_contiguity(geometry),
      weights = st_weights(nbours),
      local_lag = st_lag({{ outcome }}, nbours, weights)
    )
  
  #local getis ord (gi)
  sf_gi <- sf_local %>%
    mutate(
      Gi = local_g_perm(
        {{ outcome }},
        nbours,
        weights,
        nsim = nsim #number of simulations from argument above, 5000 is default
      )
    ) %>%
    unnest(Gi)
  
  #hotspot classification following the guide referenced above
  #with the columns 'gi' and 'p_folded_sim"
  #'p_folded_sim' is the p-value of a folded permutation test
  #adding a new column called classification
  sf_classified <- sf_gi %>%
    mutate(
      classification = case_when(
        gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
        gi > 0 & p_folded_sim <= 0.1  ~ "Somewhat hot",
        gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
        gi < 0 & p_folded_sim <= 0.1  ~ "Somewhat cold",
        TRUE                         ~ "Insignificant"
      ),
      classification = factor(
        classification,
        levels = c(
          "Very hot", "Hot", "Somewhat hot",
          "Insignificant",
          "Somewhat cold", "Cold", "Very cold"
        )
      )
    )
  
  #map it, +ve values indicate neighbours with HIGH shared values, and -ve indicate neighbours with LOW shared values
  map <- ggplot(sf_classified, aes(fill = classification)) +
    geom_sf(color = "black", lwd = 0.1) +
    scale_fill_brewer(type = "div", palette = 5) +
    theme_void() +
    labs(
      fill = "Hot Spot Classification",
      title = title %||% paste("Getis-Ord hotspots for", outcome_name)
    )
  
  #return the sf with classes, global gi results, and the map
  list(
    sf = sf_classified,
    global_go = global_go,
    map = map
  )
}

#if there's negative values, such as in the change in score etc we need moran's I instead
moran_i <- function(
    sf_data,
    outcome,
    nsim = 5000,
    title = NULL
) {
  
  outcome_name <- rlang::as_name(rlang::ensym(outcome))
  
  #create neighbours (queen contiguity)
  neighbour_list <- poly2nb(sf_data, queen = TRUE)
  
  #check for empty neighbour sets
  if (any(card(neighbour_list) == 0)) {
    stop("Some areas have no neighbours; analysis may not run correctly.")
  }
  
  #create spatial weights (binary)
  neighbour_weights <- nb2listw(neighbour_list, style = "B")
  
  #global moran's i using two sided to test for any +/- spatial autocorrelation
  global_mi <- moran.test(
    sf_data[[outcome_name]],
    neighbour_weights,
    alternative = "two.sided"
  )
  
  #calculate spatial lag
  sf_local <- sf_data %>%
    mutate(
      local_lag = lag.listw(neighbour_weights, sf_data[[outcome_name]])
    )
  
  #local moran's i (same type as local gi)
  local_mi <- localmoran_perm(
    sf_data[[outcome_name]],
    neighbour_weights,
    nsim = nsim
  )
  
  #bind local mi to sf
  sf_lisa <- bind_cols(sf_local, as_tibble(local_mi)) 
  
  #rename p value so that it isn't a nightmare to map
  sf_lisa <- sf_lisa %>%
    rename(p_sim = `Pr(folded) Sim`)
  
  
  # ---- LISA classification ----
  #spatial high clusters/low clusters (kind of like hot/cold but not exactly)
  #spatial outliers e.g. not matching neighbours kinda vibe (highlow/lowhigh)
  sf_classified <- sf_lisa %>%
    mutate(
      classification = case_when(
        {{ outcome }} > 0 & local_lag > 0 & p_sim <= 0.01 ~ "Very high–high",
        {{ outcome }} > 0 & local_lag > 0 & p_sim <= 0.05 ~ "High–high",
        {{ outcome }} > 0 & local_lag > 0 & p_sim <= 0.1  ~ "Somewhat high–high",

        {{ outcome }} < 0 & local_lag < 0 & p_sim <= 0.01 ~ "Very low–low",
        {{ outcome }} < 0 & local_lag < 0 & p_sim <= 0.05 ~ "Low–low",
        {{ outcome }} < 0 & local_lag < 0 & p_sim <= 0.1  ~ "Somewhat low–low",

        {{ outcome }} > 0 & local_lag < 0 & p_sim <= 0.05 ~ "High–low",
        {{ outcome }} < 0 & local_lag > 0 & p_sim <= 0.05 ~ "Low–high",

        TRUE ~ "Insignificant"
      ),
      classification = factor(
        classification,
        levels = c(
          "Very high–high", "High–high", "Somewhat high–high",
          "High–low",
          "Insignificant",
          "Low–high",
          "Somewhat low–low", "Low–low", "Very low–low"
        )
      )
    )

  #map it
  map <- ggplot(sf_classified, aes(fill = classification)) +
    geom_sf(color = "black", lwd = 0.1) +
    scale_fill_brewer(type = "div", palette = 5) +
    theme_void() +
    labs(
      fill = "LISA Classification",
      title = title %||% paste("Local Moran's I for", outcome_name)

    )

  #return results
  list(
    sf = sf_classified,
    global_mi = global_mi,
    map = map
  )
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

#join bnssg csv to map
bnssg <- bnssg_map %>% left_join(bnssg_csv, by = "LSOA21CD")

#save memory by removing english sf data and extra bnssg files
rm(england, bnssg_csv, bnssg_map)

#run
gi_results <- getis_ord(
  sf_data = bnssg,
  outcome = swd_pct_seg4_5,
  title = "Cambridge Multimorbidity Score, % in Segments 4–5 hotspots"
)

#save the global getis ord results, stick it into a df first
gi_df <- data.frame(
  global_G     = unname(results$global_go$statistic),
  p_value      = results$global_go$p.value,
  expectation  = results$global_go$estimate[1],
  variance     = results$global_go$estimate[2]
)

write.csv(
  gi_df,
  "bnssg_global_getis_ord_swd_pct_seg4_5.csv",
  row.names = FALSE
)

#save the map
ggsave(
  filename = "bnssg_getis_ord_swd_pct_seg4_5_hotspots.png",
  plot = results$map,
  width = 8,
  height = 10,
  dpi = 300
)

#run moran's i for cms score segments
moran_results_score <- moran_i(sf_data = bnssg,  outcome = swd_pct_seg4_5,
                          title = "Cambridge Multimorbidity Score,  % Segments 4–5 hotspots")

#run moran's i for change
moran_results_change <- moran_i(sf_data = bnssg,  outcome = swd_pct_seg4_5_yoy_change,
                       title = "Cambridge Multimorbidity Score, Change in % Segments 4–5 Clustering")

#save the global Moran’s I results (score)
mi_score_df <- data.frame(
  moran_I     = unname(moran_results_score$global_mi$estimate["Moran I statistic"]),
  p_value     = moran_results_score$global_mi$p.value,
  expectation = unname(moran_results_score$global_mi$estimate["Expectation"]),
  variance    = unname(moran_results_score$global_mi$estimate["Variance"])
)

write.csv(
  mi_score_df,
  "bnssg_global_moran_i_swd_pct_seg4_5.csv",
  row.names = FALSE
)

#save the global Moran’s I results (change)
mi_change_df <- data.frame(
  moran_I     = unname(moran_results_change$global_mi$estimate["Moran I statistic"]),
  p_value     = moran_results_change$global_mi$p.value,
  expectation = unname(moran_results_change$global_mi$estimate["Expectation"]),
  variance    = unname(moran_results_change$global_mi$estimate["Variance"])
)

write.csv(
  mi_change_df,
  "bnssg_global_moran_i_swd_pct_seg4_5_yoy_change.csv",
  row.names = FALSE
)

#save Moran’s I map — score
ggsave(
  filename = "bnssg_moran_i_swd_pct_seg4_5_lisa.png",
  plot = moran_results_score$map,
  width = 8,
  height = 10,
  dpi = 300
)

#save Moran’s I map — change
ggsave(
  filename = "bnssg_moran_i_swd_pct_seg4_5_yoy_change_lisa.png",
  plot = moran_results_change$map,
  width = 8,
  height = 10,
  dpi = 300
)


