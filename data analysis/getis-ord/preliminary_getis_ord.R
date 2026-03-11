library(sf)
library(sfdep)
library(spdep)
library(dplyr)
library(tidyr)
library(ggplot2)

#Getis-ord analysis partially following this tutorial: https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut 

#pull wd from paths.R (put in .gitignore)
source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset
bnssg_csv <- read.csv("./BNSSG/linked/bnssg_long.csv")

#pull out LSOAS
bnssg_lsoas <- bnssg$LSOA21CD

#link in LSOAs
#load lsoa data and then restrict to BNSSG LSOAS by filtering out LSOA codes not contained in bnssg

england <- st_read(file.path("./LSOA 2021/LSOA/LSOA_2021_EW_BGC.shp"))
#Projected CRS: OSGB36 / British National Grid

bnssg_map <- england[england$LSOA21CD %in% bnssg_lsoas,]

#save memory by removing english sf data
rm(england)

#join bnssg csv to map
bnssg <- bnssg_map %>% left_join(bnssg_csv, by = "LSOA21CD")

#create a list of neighbours, queen is a type of contiguity criteria which is when a common edge/vertex is shared
neighbour_list <- poly2nb(bnssg, queen=TRUE)

#check for empty neighbour sets (e.g. LSOAs without any neighbours)
#analysis may not run correctly if there are any. in this case, none
neighbour_empty <- which(card(neighbour_list) == 0)
neighbour_empty

#assign binary weights, 1 to all neighbouring features, 0 to other features
neighbour_weights <- nb2listw(neighbour_list, style = "B")

#calculate spatial lag of % in segments 4_5
cms_lag <- lag.listw(neighbour_weights, bnssg$swd_pct_seg4_5)

#calculate global G statistic (Getis-Ord global G)
globalG.test(bnssg$swd_pct_seg4_5, neighbour_weights)

#test local spatial autocorrelation

#neighbour identification, create weights, calculate spatial lag
bnssg_nbours <- bnssg %>%
  mutate(
    nbours = st_contiguity(geometry),
    weights = st_weights(nbours),
    cms_local_lag = st_lag(bnssg$swd_pct_seg4_5, nbours, weights)
  )

#calculate local getis-ord using local_g_perm
bnssg_hotspots <- bnssg_nbours %>%
  mutate(
    Gi = local_g_perm(bnssg$swd_pct_seg4_5, nbours, weights, nsim = 5000) 
    #use 5000 MonteCarlo simulations
  ) %>%
  unnest(Gi) #data structure step to remove gi from dataframe

#map it, +ve values indicate neighbours with HIGH shared values, and -ve indicate neighbours with LOW shared values
bnssg_hotspots %>%
  ggplot((aes(fill=gi)))+
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2() #this sets the midpoint to 0, which is 'truly random'

#code from the tutorial to make a hot/cold spot map based on the getis ord p values for each LSOA
bnssg_hotspots |> 
  # with the columns 'gi' and 'p_folded_sim"
  # 'p_folded_sim' is the p-value of a folded permutation test
  select(gi, p_folded_sim) |> 
  mutate(
    # Add a new column called "classification"
    classification = case_when(
      # Classify based on the following criteria:
      gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(fill = classification)) +
  geom_sf(color = "black", lwd = 0.1) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "Cambridge Multimorbidity Score, % Segments 4_5 hotspots"
  )

#next step is to run the same analysis on the residuals, wrap this into a function so we can do it for each year
