library(dplyr)
library(tidyr)

#cribbed table 1 function from my own here: https://github.com/eyles-ec/common-ambition/blob/main/tables/table_1.R
#removed categorgical block, gt table (interesting, but not useful currently), then no need for bind
#added overall summary, as it is relevant to this study 
generate_table1 <- function(data, group_var, continuous_vars,
                            variable_labels, row_groups,
                            title = "Table 1: Descriptive Summary by Group") {
  
  #grouped summary 
  cont_by_group <- data %>%
    dplyr::select(all_of(c(group_var, continuous_vars))) %>%
    pivot_longer(
      cols = -all_of(group_var),
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable, .data[[group_var]]) %>%
    summarise(
      mean = round(mean(value, na.rm = TRUE), 1),
      sd   = round(sd(value, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(mean, " (", sd, ")")) %>%
    dplyr::select(variable, all_of(group_var), label) %>%
    pivot_wider(names_from = all_of(group_var), values_from = label)
  
  #overall summary
  cont_overall <- data %>%
    dplyr::select(all_of(continuous_vars)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable) %>%
    summarise(
      mean = round(mean(value, na.rm = TRUE), 1),
      sd   = round(sd(value, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(Overall = paste0(mean, " (", sd, ")")) %>%
    dplyr::select(variable, Overall)
  
  #merge group summary + overall summary
  table1 <- cont_by_group %>%
    left_join(cont_overall, by = "variable")
  
  #apply variable labels
  table1$variable <- variable_labels[table1$variable]
  
  #reorder rows according to row_groups
  ordered_vars <- unlist(row_groups)
  ordered_labels <- variable_labels[ordered_vars]
  
  #mutate into correct order
  table1 <- table1 %>%
    mutate(variable = factor(variable, levels = ordered_labels)) %>%
    arrange(variable)
  
  return(table1)
}


#pull wd from paths.R (put in .gitignore)
source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset
bnssg <- read.csv("./BNSSG/linked/bnssg.csv")

#variable labels for the table (nicer than just using the name) 
variable_labels <- c(
  
  #demographics
  swd_pct_female_24 = "Percent Female",
  swd_median_age_24 = "Median Age",
  swd_mean_age_24 = "Mean Age",
  popden = "Population Density",
  unemp = "Percent unemployed",
  poor_health = "Percent Bad/Very Bad Health",
  asian = "Percent Asian",
  black = "Percent Black",
  mixed = "Percent Mixed",
  white = "Percent White",
  other = "Percent Other",
  IMD25_decile = "IMD Decile",
  
  #CMS 2021
  swd_pct_seg4_5_21 = "Percent Segment 4-5, 2021",
  # swd_median_cms_21 = "Median CMS, 2021",
  # swd_mean_cms_21 = "Mean CMS, 2021",
  
  #CMS 2022
  swd_pct_seg4_5_22 = "Percent Segment 4-5, 2022",
  # swd_median_cms_22 = "Median CMS, 2022",
  # swd_mean_cms_22 = "Mean CMS, 2022",
  
  #CMS 2023
  swd_pct_seg4_5_23 = "Percent Segment 4-5, 2023",
  # swd_median_cms_23 = "Median CMS, 2023",
  # swd_mean_cms_23 = "Mean CMS, 2023",
  
  #CMS 2024
  swd_pct_seg4_5_24 = "Percent Segment 4-5, 2024",
  # swd_median_cms_24 = "Median CMS, 2024",
  # swd_mean_cms_24 = "Mean CMS, 2024",
  
  #Pollution, using unicode to represent ugm3 in its proper form [µg/m³]
  benzene_2024_mean_ugm3 = "Benzene mean \u00B5g/m\u00B3",
  nox_2024_mean_ugm3 = "NOx mean \u00B5g/m\u00B3",
  ozone_2024days_over_120ugm3 = "Ozone days over 120 \u00B5g/m\u00B3",
  pm10_2024_mean_ugm3 = "PM10 mean \u00B5g/m\u00B3",
  pm2_5_2024_mean_ugm3 = "PM2.5 mean \u00B5g/m\u00B3",
  so2_2024_mean_ugm3 = "SO2 mean \u00B5g/m\u00B3",
  
  #Built Environment, same for km2 
  road_density_km_per_km2 = "Road density, km/km\u00B2",
  intersection_density_per_km2 = "Intersection density/km\u00B2",
  
  #Land use
  Agriculture = "Percent Agricultural land",
  Forest..open.land.and.water = "Percent Forest, Open Land, and Water",
  Outdoor.recreation = "Percent Outdoor Recreational land",
  Residential.gardens = "Percent Residential gardens "
)

#pull out variables from label vector
continuous_vars <- names(variable_labels)

#ordering list to order table, can be changed 
row_groups <- list(
  "Demographics" = c(
    "swd_pct_female_24",
    "swd_median_age_24",
    "swd_mean_age_24",
    "popden",
    "unemp",
    "poor_health",
    "asian",
    "black",
    "mixed",
    "white",
    "other",
    "IMD25_decile"
  ),
  
  "Cambridge Multimorbidity Score 2021" = c(
    "swd_pct_seg4_5_21"
    # "swd_median_cms_21",
    # "swd_mean_cms_21"
  ),
  
  "Cambridge Multimorbidity Score 2022" = c(
    "swd_pct_seg4_5_22"
    # "swd_median_cms_22",
    # "swd_mean_cms_22"
  ),
  
  "Cambridge Multimorbidity Score 2023" = c(
    "swd_pct_seg4_5_23"
    # "swd_median_cms_23",
    # "swd_mean_cms_23"
  ),
  
  "Cambridge Multimorbidity Score 2024" = c(
    "swd_pct_seg4_5_24"
    # "swd_median_cms_24",
    # "swd_mean_cms_24"
  ),
  
  "Air Pollution" = c(
    "benzene_2024_mean_ugm3",
    "nox_2024_mean_ugm3",
    "ozone_2024days_over_120ugm3",
    "pm10_2024_mean_ugm3",
    "pm2_5_2024_mean_ugm3",
    "so2_2024_mean_ugm3"
  ),
  
  "Built Environment" = c(
    "road_density_km_per_km2",
    "intersection_density_per_km2"
  ),
  
  "Land Use" = c(
    "Agriculture",
    "Forest..open.land.and.water",
    "Outdoor.recreation",
    "Residential.gardens"
  )
)

#grouping variable
split_group <- "LAD22NM"

#run table 1
t1 <- generate_table1(bnssg, split_group, continuous_vars, variable_labels, row_groups)

#save table file for output, create directory if an issue
dir.create("./BNSSG/output", recursive = TRUE, showWarnings = FALSE)

write.csv(t1, "./BNSSG/output/table1.csv", row.names = FALSE)