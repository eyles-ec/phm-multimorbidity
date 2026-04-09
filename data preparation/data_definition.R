library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(rlang)

#function to process CSVs labelled by year, then add year to the columns for easy joining
process_csv <- function(path) {
  
  #extract the year from each file name
  year <- str_extract(basename(path), "\\d{2}(?=\\.csv$)")
  
  #read in
  df <- read_csv(path, show_col_types = FALSE)
  
  #identify columns add year to
  rename_to <- setdiff(names(df), keep)

  #apply renaming, ensuring if there's any missing columns the function doesn't break (rename_with)
  df %>%
    rename_with(~ paste0(.x, "_", year), all_of(rename_to))
  
}

#calculate a change outcome in df, using 'outcome.' this calculates absolute change from 2021 and change through time year on year. 
change_outcome <- function(df, outcome) {
  
  outcome_name <- as_name(enquo(outcome))
  
  #add in a change variable to the long format data (change from 2021), and a relative change which is change from one previous lag
  df <- df %>%
    group_by(LSOA21CD) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      !!paste0(outcome_name, "_abs_change") := {{outcome}} - {{outcome}}[year == 2021],
      !!paste0(outcome_name, "_yoy_change") := ifelse(year == 2021, 0, {{outcome}}- lag({{outcome}}))
    ) %>%
    ungroup()
  
  return(df)
}

source("paths.R")

#set working directory 
setwd(wd)

#load all linked data (see https://github.com/eyles-ec/epri-england/data_linkages)
all <- read.csv("./linked/all_england_combined_CSV.csv")

#remove uneeded columns by keeping project-relevant covariates
all <- all %>%
  select(LSOA21CD, population, households, popden, unemp, good_health, fair_health, poor_health, asian, black, mixed, white, 
         other, IMD25, IMD25_decile, benzene_2024_mean_ugm3, nox_2024_mean_ugm3, ozone_2024days_over_120ugm3,
         pm10_2024_mean_ugm3, pm2_5_2024_mean_ugm3, so2_2024_mean_ugm3, total_road_m, n_intersections, area_km2, 
         road_density_km_per_km2, intersection_density_per_km2, chn2024, Agriculture, Forest..open.land.and.water, 
         Outdoor.recreation, Residential.gardens)

#load in BNSSG data

#set some parameters for the join
keep <- c("LSOA21CD", "LSOA21NM", "CHGIND", "LAD22NM")
files <- list.files("../BNSSG", pattern = "\\.csv$", full.names = TRUE)

#load and process the files, appending a year to each relevant variable
dfs <- map(files, process_csv)

#merge into one, but you get the weird duplication patterning
bnssg <- reduce(dfs, full_join, by = "LSOA21CD")

#loop through to make just one column for each of the codes
remove_extra <- c("LSOA21NM", "CHGIND", "LAD22NM")

for (c in remove_extra) {
  #collect all matches
  matches <- grep(paste0("^", c), names(bnssg), value = TRUE)
  #keep the first match, and rename it to remove the . if needed
  keep_col <- matches[1]
  names(bnssg)[names(bnssg) == keep_col] <- c
  
  #drop the duplicates
  if (length(matches) > 1) {
    dup <- matches[-1]
    bnssg <- bnssg[, !names(bnssg) %in% dup]
  }
}

#add in exposome data, clipped to BNSSG footprint, and reorder so that the outcomes are first
bnssg_cols <- names(bnssg)

bnssg <- all %>%
  inner_join(bnssg, by = "LSOA21CD") %>%
  select(LSOA21CD,
         bnssg_cols[bnssg_cols != "LSOA21CD"],
         everything())

#make a long format dataset for analysis

bnssg_long <- bnssg %>%
  pivot_longer(
    cols = matches("_(21|22|23|24)$"),
    names_to = c(".value", "year"),
    names_pattern = "(.*)_(\\d{2})$"
  ) %>%
  mutate(
    year = as.numeric(paste0("20", year))
  )

#create change for mean CMS and for pct of seg 4-5
bnssg_long <- change_outcome(bnssg_long, swd_mean_cms)
bnssg_long <- change_outcome(bnssg_long, swd_pct_seg4_5)

#save complete long and wide files for analysis, create directory if an issue
dir.create("./BNSSG/linked", recursive = TRUE, showWarnings = FALSE)

write.csv(bnssg, "./BNSSG/linked/bnssg_wide.csv", row.names = FALSE)
write.csv(bnssg_long, "./BNSSG/linked/bnssg_long.csv", row.names = FALSE)

