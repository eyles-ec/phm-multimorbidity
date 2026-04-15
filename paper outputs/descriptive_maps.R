library(sf)
library(tmap)
library(dplyr)
library(rlang)
library(purrr)
library(tibble)

#map function adapted from epri-england repo: https://github.com/eyles-ec/epri-england/tree/main/outputs
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
  
  #coercion of var type
  var <- rlang::ensym(var)
  
  #Base map
  tm <- tm_shape(bnssg) 
  
 #binary/categorical variables
  if (is.null(style)) {
      
      tm <- tm +
        tm_fill(
          col     = rlang::as_string(var),
          palette = palette,
          title   = legend_title %||%
            rlang::as_string(var)
        ) +
        tm_borders(col = "grey40", lwd = 0.5)
      
    #continuous variables
    } else {
      
      tm <- tm +
        tm_polygons(
          fill = rlang::as_string(var),
          col  = NULL,
          fill.scale = tm_scale_intervals(
            style  = style,
            n      = n_classes,
            values = palette
          ),
          fill.legend = tm_legend(
            title = legend_title %||%
              rlang::as_string(var)
          )
        )
    }
  
  #Optional overlay (e.g. roads, boundaries), artefact of epri repo
  if (!is.null(overlay)) {
    
    overlay <- sf::st_transform(overlay, sf::st_crs(bnssg))
    overlay <- sf::st_intersection(overlay, bnssg)
    
    tm <- tm +
      tm_shape(overlay) +
      tm_lines(col = overlay_col, lwd = 1)
  }
  
  #Layout and map decorations
  tm +
    tm_scalebar(position = c("center", "bottom")) +
    tm_compass(type = "arrow", position = c("right", "top")) +
    tm_layout(
      frame           = TRUE,
      frame.lwd       = 1,
      frame.col       = "black",
      legend.outside  = FALSE,
      legend.position = c("left", "top")
    )
}


#map export function taken from above repo

export_maps <- function(
    maps,
    filenames,
    width  = 1920,
    height = 1080,
    dpi    = 300
) {
  
  #error avoidance - both maps and filenames have to match in length
  if (length(maps) != length(filenames)) {
    stop("Length of 'maps' and 'filenames' must match.")
  }
  
  for (i in seq_along(maps)) {
    tmap_save(
      maps[[i]],
      filename = filenames[[i]],
      width    = width,
      height   = height,
      dpi      = dpi
    )
  }
  
  message("All maps exported successfully.")
}


#pull wd from paths.R (gitignored)
source("../paths.R")
setwd(wd)


#load analysis dataset
bnssg_csv <- read.csv("./BNSSG/linked/bnssg_long.csv")

#extract LSOA codes
bnssg_lsoas <- bnssg_csv$LSOA21CD

#load England LSOA boundaries
england <- sf::st_read("./LSOA 2021/LSOA/LSOA_2021_EW_BGC.shp")

#restrict to BNSSG only
bnssg_map <- england[england$LSOA21CD %in% bnssg_lsoas, ]

#join in csv data
bnssg <- bnssg_map |>
  dplyr::left_join(bnssg_csv, by = "LSOA21CD")

#tidy up memory to save time/space
rm(england, bnssg_csv, bnssg_map)


#derivations
# >10% in seg 4/5, expected is 10% of population (7% in 4, 3% in 5)
bnssg$above10 <- bnssg$swd_pct_seg4_5 > 10

# recentred at 10%
bnssg$recentred <- bnssg$swd_pct_seg4_5 - 10


#map specifications in tibble, which allows us to call the mapping function over this list of maps

map_specs <- tribble(
  ~var,                  ~palette,        ~style,      ~n_classes, ~legend_title,                                 ~filename,
  
  # Segment 4–5
  "swd_pct_seg4_5",      "viridis",        "quantile",  5,          "% in segment 4–5",                            "seg45_pct.png",
  "above10",             c("lightblue",
                           "pink"),        NULL,        2,          ">10% in seg 4–5",                             "seg45_above10.png",
  "recentred",           "-RdBu",          "quantile",  5,          "Recentred (10% = 0)",                          "seg45_recentred.png",
  
  # Deprivation
  "IMD25",        "-viridis",       "quantile",     10,         "IMD 2025 decile (1 = most deprived)",         "imd25_decile.png",
  
  # Environment
  "pm10_2024_mean_ugm3", "inferno",         "quantile",  5,          "PM10 (µg/m³, 2024 mean)",           "pm10_2024_mean.png",
  
  # Demographics
  "swd_mean_age",        "plasma",          "quantile",  5,          "Mean age (years)",                           "mean_age.png",
  "white",               "viridis",         "quantile",  5,          "% White",                                    "pct_white.png",
  
  # Built environment
  "Residential.gardens", "viridis",         "quantile",  5,          "Residential gardens (%)",                   "residential_gardens.png",
  "road_density",        "magma",            "quantile",  5,          "Road density (km/km²)",                     "road_density.png",
  "popden",              "turbo",            "quantile",  5,          "Population density (per km²)",              "population_density.png"
)


#generate maps, using purrr's pmap, pmap applies a function once per row, using multiple inputs
#for each row of map_specs, call the function descriptive_maps, using each column as an argument 
#for descriptive_maps, returning a list of maps, which we can then export
tmap_mode("plot")

maps <- pmap(
  map_specs,
  function(var, palette, style, n_classes, legend_title, filename) {
    descriptive_maps(
      bnssg        = bnssg,
      var          = !!rlang::sym(var),
      palette      = palette,
      style        = style,
      n_classes    = n_classes,
      legend_title = legend_title
    )
  }
)

#export the list of maps from the above, as pmap will return a list of maps, use the filenames from map_specs
export_maps(
  maps      = maps,
  filenames = file.path("./outputs/maps", map_specs$filename),
  width     = 1920,
  height    = 1080,
  dpi       = 300
)
