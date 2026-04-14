library(sf)
library(sfdep)
library(spdep)
library(dplyr)
library(tidyr)
library(ggplot2)

#map function taken from epri-england repo: https://github.com/eyles-ec/epri-england/tree/main/outputs
descriptive_maps <- function(
    bnssg,
    var,
    palette = "viridis",
    n_classes = 5,
    style = "quantile",
    overlay = NULL,
    overlay_name = NULL,   
    overlay_col = "black",
    legend_title = NULL
) {
  
  #Base map
  tm <- tm_shape(bnssg) +
    tm_polygons(
      fill = {{ var }},
      col = NULL,   # removes borders in tmap v4
      fill.scale = tm_scale_intervals(
        style = style,
        n = n_classes,
        values = palette
      ),
      fill.legend = tm_legend(
        title = legend_title %||% rlang::as_string(rlang::ensym(var))
      )
    )
  
  #Overlay (optional, no legend)
  if (!is.null(overlay)) {
    
    #coerce coordinate system first
    overlay <- sf::st_transform(overlay, sf::st_crs(bnssg))
    # Clip overlay to bnssg
    overlay_clipped <- sf::st_intersection(overlay, bnssg)
    
    tm <- tm +
      tm_shape(overlay_clipped) +
      tm_lines(col = overlay_col, lwd = 1)
  }
  
  # Layout and map decorations
  tm +
    tm_scalebar(position = c("center", "bottom")) +
    tm_compass(type = "arrow", position = c("right", "top")) +
    tm_layout(
      frame = TRUE,
      frame.lwd = 1,
      frame.col = "black",
      legend.outside = FALSE,
      legend.position = c("left", "center")
    )
}

#export function also taken from there 
export_maps <- function(
    maps,
    filenames,
    width = 1920,
    height = 1080,
    dpi = 300
) {
  
  # Error avoidance 
  if (length(maps) != length(filenames)) {
    stop("Length of 'maps' and 'filenames' must match.")
  }
  
  # Loop through and export each map
  for (i in seq_along(maps)) {
    tmap_save(
      maps[[i]],
      filename = filenames[[i]],
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  message("All maps exported successfully.")
}


#pull wd from paths.R (put in .gitignore)
source("../paths.R")
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

#tbd functions and saving

#calculate above or below 10% for seg 4/5 (expected is 10% of population)
bnssg$above10 <- bnssg$swd_pct_seg4_5 > 10

#mean centre the % 4-5 
bnssg$recentred <- bnssg$swd_pct_seg4_5 - 10

tmap_mode("plot")

tm_shape(bnssg) +
  tm_fill("above10",
          palette = c("lightblue", "pink"),
          title = "> 10%") +
  tm_borders() +
  tm_layout(legend.outside = TRUE)

tmap_mode("plot")
tm_shape(bnssg) +
  tm_fill("recentred",
          palette = "-RdBu",
          style   = "quantile",
          title   = "Recentred (0 at 10%)") +
  tm_borders()