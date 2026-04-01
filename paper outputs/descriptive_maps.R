bnssg$above10 <- bnssg$swd_pct_seg4_5 > 10
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