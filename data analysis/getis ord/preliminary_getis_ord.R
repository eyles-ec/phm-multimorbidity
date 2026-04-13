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
  title = "Cambridge Multimorbidity Score, % Segments 4–5 hotspots"
)

gi_results$global_go
gi_results$map

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
