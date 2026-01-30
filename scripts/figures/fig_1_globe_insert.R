# Overview map (globe insert) for Figure 1
# Jakob J. Assmann jakob.assmann@uzh.ch 30 January 2026

# Dependencies
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(cowplot)
library(rnaturalearth)
library(magick)

# Load plot data
plot_locations <- read_sf("data/geometries/plot_locations.gpkg")

## Generate pan-Arctic overview map

# Set crs by rotating NSIDC projection onto Canada
custom_crs <- st_crs(3413)
cat(custom_crs$wkt)
custom_crs$wkt <- gsub('"Longitude of origin",-45,', '"Longitude of origin",-100,', custom_crs$wkt)
cat(custom_crs$wkt)

# Define Polygon with cut off at 50 deg
bounds_poly <- data.frame(
  x = -100, # 0,
  y = 77.5 # 90
) %>%
  st_as_sf(coords = c("x", "y"),
           crs = 4326) %>%
  st_transform(custom_crs) %>%
  st_buffer(2750000)#3314693)

# Get outline of countries
countries_poly <- ne_countries(scale = 50, returnclass = "sf") %>%
  st_transform(custom_crs) %>%
  st_intersection(bounds_poly) %>%
  st_union() 

# Bring it all together into an overview map
pan_arctic_overview <- ggplot() +
  geom_sf(data = st_difference(st_as_sfc(st_bbox(bounds_poly)), bounds_poly),
          colour = NA,
          fill = "#000000",
          size = 1.5) +
  # geom_sf(data = bounds_poly,
  #         fill = "#012240", 
  #         colour = "NA", 
  #         size = 1,
  #         alpha = 0.5) +
  geom_sf(data = countries_poly,
          size = 0.01, 
          fill = "#e5e5e5", 
          colour = "#000000CC") +
  geom_sf(data = st_centroid(summarise(plot_locations)),
          shape = 21,
          fill = "#FF0000CC",
          colour = "#000000CC",
          size = 6) +
  geom_sf(data = st_buffer(bounds_poly,
                        -(st_bbox(bounds_poly)[3] - st_bbox(bounds_poly)[1]) * 0.005),
          fill = "NA",
          colour = "#FFFFFFCC",
          linewidth = 1.5) +
  coord_sf(xlim = st_bbox(bounds_poly)[c(1,3)],
           ylim = st_bbox(bounds_poly)[c(2,4)],
           crs = st_crs(bounds_poly),
           expand = FALSE) +
  theme_map() +
  theme(panel.grid.major = element_line(color = "#FFFFFF88",
                                        linewidth = 0.5,
                                        linetype = 1),
        plot.margin = margin(0, 0, 0, 0, unit = "cm"))
save_plot("figures/overview_map.png", pan_arctic_overview,
          base_asp = 1, bg = "#012541CC")

image_read("figures/overview_map.png") %>%
image_transparent("#000000", fuzz = 10) %>%
image_write("figures/overview_map.png")
