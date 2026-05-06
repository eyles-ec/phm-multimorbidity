library(dplyr)
library(ggplot2)
library(rlang)
library(dplyr)
library(ggplot2)
library(rlang)

#plot outcomes (or anything really) over time by LAD and overall
#supplying a filename saves it with default quality values
plot_outcomes <- function(data,
                                 outcome,
                                 y_label,
                                 title,
                                 subtitle = NULL,
                                 group_var = LAD22NM,
                                 filename = NULL,
                                 width = 8,
                                 height = 6,
                                 dpi = 300) {
  
  outcome <- enquo(outcome)
  group_var <- enquo(group_var)
  
  #LAD-level means by year
  lad_means <- data %>%
    group_by(year, !!group_var) %>%
    summarise(
      mean_value = mean(!!outcome, na.rm = TRUE),
      .groups = "drop"
    )
  
  #Overall mean by year
  overall_means <- data %>%
    group_by(year) %>%
    summarise(
      mean_value = mean(!!outcome, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot() +
    #LAD trends
    geom_line(
      data = lad_means,
      aes(x = year,
          y = mean_value,
          colour = !!group_var,
          group = !!group_var),
      alpha = 0.6,
      linewidth = 0.8
    ) +
    #Overall trend
    geom_line(
      data = overall_means,
      aes(x = year,
          y = mean_value),
      colour = "black",
      linewidth = 1,
      linetype = "dashed"
    ) +
    geom_point(
      data = overall_means,
      aes(x = year,
          y = mean_value),
      colour = "black",
      size = 2
    ) +
    labs(
      x = "Year",
      y = y_label,
      title = title,
      subtitle = subtitle,
      colour = "LAD"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right"
    ) +
  scale_x_continuous(
    breaks = sort(unique(data$year)),
    labels = function(x) as.integer(x)
    )
  
  
  #export plot if a filename is supplied
  if (!is.null(filename)) {
    dir.create("./BNSSG/output", recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = file.path("./BNSSG/output", filename),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(p)
}

#pull wd from paths.R (put in .gitignore)
source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset
bnssg <- read.csv("./BNSSG/linked/bnssg_long.csv")

#percent in segs 4-5
plot_outcomes(
  data = bnssg,
  outcome = swd_pct_seg4_5,
  y_label = "Percent in Segment 4–5",
  title = "Percent of Population in Cambridge Multimorbidity Segments 4–5",
  subtitle = "Overall and LAD-specific trends, 2021–2024",
  filename = "percent_seg4_5_overall_lad.png"
)

#percent in segs 4-5 yoy change, excluding 2021 which is baseline year (so change of 0)
plot_outcomes(
  data = bnssg %>% filter(year != "2021"),
  outcome = swd_pct_seg4_5_yoy_change,
  y_label = "Year on Year Change in Segment 4–5 (%)",
  title = "Year on Year Change in Percent in CMS 4–5",
  subtitle = "Overall and LAD-specific trends, 2021–2024",
  filename = "percent_seg4_5_yoy_overall_lad.png"
)

#mean cms
plot_outcomes(
  data = bnssg,
  outcome = swd_mean_cms,
  y_label = "Mean Cambridge Multimorbidity Score",
  title = "Mean Cambridge Multimorbidity Score Over Time",
  subtitle = "Overall and LAD-specific trends, 2021–2024",
  filename = "mean_cms_overall_lad.png"
)

#mean cms yoy change, excluding 2021 which is baseline year (so change of 0)
plot_outcomes(
  data = bnssg %>% filter(year != 2021),
  outcome = swd_mean_cms_yoy_change,
  y_label = "Change in Mean Cambridge Multimorbidity Score",
  title = "Year on Year Change in Mean Cambridge Multimorbidity Score",
  subtitle = "Overall and LAD-specific trends, 2021–2024",
  filename = "mean_cms_yoy_overall_lad.png"
)