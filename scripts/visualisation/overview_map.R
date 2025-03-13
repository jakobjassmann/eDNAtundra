# Make a quick overview map

# Dependencies
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(cowplot)
library(rnaturalearth)

# Load plot data
plot_locations <- read_sf("data/geometries/plot_locations.gpkg")



## Generate pan-Arctic overview map

# Set crs by rotating NSIDC projection onto Canada
custom_crs <- st_crs(3413)
cat(custom_crs$wkt)
custom_crs$wkt <- gsub('"Longitude of origin",-45,', '"Longitude of origin",-90,', custom_crs$wkt)
cat(custom_crs$wkt)

# Define Polygon with cut off at 50 deg
bounds_poly <- data.frame(
  x = 0,
  y = 90
) %>%
  st_as_sf(coords = c("x", "y"),
           crs = 4326) %>%
  st_transform(custom_crs) %>%
  st_buffer(3314693)

# Get outline of countries
countries_poly <- ne_countries(scale = 50, returnclass = "sf") %>%
  st_transform(custom_crs) %>%
  st_intersection(bounds_poly) %>%
  st_union() 

# Bring it all together into an overview map
pan_arctic_overview <- ggplot() +
  geom_sf(data = bounds_poly,
          fill = "#82C4F5CC", 
          colour = "NA", 
          size = 1) +
  geom_sf(data = countries_poly,
          size = 0.01, 
          fill = "#e5e5e5CC", 
          colour = "#000000CC") +
  geom_sf(data = st_centroid(summarise(plot_locations)),
          shape = 21,
          fill = "#FF0000CC",
          colour = "#000000CC",
          size = 6) +
  geom_sf(data = bounds_poly,
          fill = "NA",
          colour = "#000000CC",
          size = 1.5) +
  coord_sf(xlim = st_bbox(bounds_poly)[c(1,3)],
           ylim = st_bbox(bounds_poly)[c(2,4)],
           crs = st_crs(bounds_poly)) +
  # scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0)) +
  theme_nothing()
save_plot("figures/overview_map.png", pan_arctic_overview,
          base_asp = 1)
