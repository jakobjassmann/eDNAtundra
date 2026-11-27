# Quick drone maps

# Dependencies
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(cowplot)
library(rnaturalearth)

# Load plot locations
plot_locations <- read_sf("data/geometries/plot_locations.gpkg")

# Load site boundaries
site_bounds <- read_sf("data/geometries/site_boundaries.gpkg") 

# Buffer site boundaries
site_bounds_buffered <- st_buffer(site_bounds, 20)


# Create mask for outside the boundaries
mask_site_1 <- st_difference(st_buffer(site_bounds[1,], 1000),
                             site_bounds_buffered[1,])
mask_site_2 <- st_difference(st_buffer(site_bounds[2,], 1000),
                             site_bounds_buffered[2,])
mask_site_3 <- st_difference(st_buffer(site_bounds[3,], 1000),
                             site_bounds_buffered[3,])

# Load drone data
site_1 <- rast("C:/Users/jakob/Desktop/eDNA drone mosaics/eDNA_site_1_rgb_transparent_mosaic_group1.tif")
site_2 <- rast("C:/Users/jakob/Desktop/eDNA drone mosaics/eDNA_site_2_rgb_transparent_mosaic_group1.tif")
site_3 <- rast("C:/Users/jakob/Desktop/eDNA drone mosaics/eDNA_site_3_rgb_transparent_mosaic_group1.tif")


# Plot rasters

plot_site_1 <- 
  ggplot() +
  geom_spatraster_rgb(data = site_1, alpha = 0.9) +
  geom_sf(data = filter(site_bounds, site == "Site 1"),
          colour = "#F2D64B", fill = "NA",
          linewidth = 1,
          linetype = "dotted") +
  geom_sf(data = filter(plot_locations, site == "Site 1") %>%
            mutate(plot_type = case_when(substr(Name, 3, 3) == "R" ~ "stratified random",
                                         substr(Name, 3, 3) == "S" ~ "systematic")),
          mapping = aes(shape = plot_type, size = plot_type),
          colour = "white",
          fill = "#F2D64B") +
  scale_shape_manual(values = c(24, 21)) +
  scale_size_manual(values = c(2, 3)) +
  coord_sf(xlim = c(ext(site_bounds_buffered[1,])[1], ext(site_bounds_buffered[1,])[2]),
           ylim = c(ext(site_bounds_buffered[1,])[3], ext(site_bounds_buffered[1,])[4]),
           clip = "off") +
  theme_nothing() +
  theme(plot.background = element_rect(fill = "black"))
save_plot("figures/not_included_in_manuscript/site_1_rgb.png", plot_site_1,
          base_asp = 1)


plot_site_2 <- 
  ggplot() +
  geom_spatraster_rgb(data = site_2, alpha = 0.9) +
  geom_sf(data = filter(site_bounds, site == "Site 2"),
          colour = "#5F1773", fill = "NA",
          linewidth = 1,
          linetype = "dotted") +
  geom_sf(data = filter(plot_locations, site == "Site 2") %>%
            mutate(plot_type = case_when(substr(Name, 3, 3) == "R" ~ "stratified random",
                                         substr(Name, 3, 3) == "S" ~ "systematic")),
          mapping = aes(shape = plot_type, size = plot_type),
          colour = "white",
          fill = "#5F1773") +
  scale_shape_manual(values = c(24, 21)) +
  scale_size_manual(values = c(2, 3)) +
  coord_sf(xlim = c(ext(site_bounds_buffered[2,])[1], ext(site_bounds_buffered[2,])[2]),
           ylim = c(ext(site_bounds_buffered[2,])[3], ext(site_bounds_buffered[2,])[4]),
           clip = "off") +
  theme_nothing() +
  theme(plot.background = element_rect(fill = "white"))
save_plot("figures/not_included_in_manuscript/site_2_rgb.png", plot_site_2,
          base_asp = 1)

plot_site_3 <- 
  ggplot() +
  geom_spatraster_rgb(data = site_3, alpha = 0.8) +
  geom_sf(data = filter(site_bounds, site == "Site 3"),
          colour = "#F2507B", fill = "NA",
          linewidth = 1,
          linetype = "dotted") +
  geom_sf(data = filter(plot_locations, site == "Site 3") %>%
            mutate(plot_type = case_when(substr(Name, 3, 3) == "R" ~ "stratified random",
                                         substr(Name, 3, 3) == "S" ~ "systematic")),
          mapping = aes(shape = plot_type, size = plot_type),
          colour = "white",
          fill = "#F2507B") +
  scale_shape_manual(values = c(24, 21)) +
  scale_size_manual(values = c(2, 3)) +
  coord_sf(xlim = c(ext(site_bounds_buffered[3,])[1], ext(site_bounds_buffered[3,])[2]),
           ylim = c(ext(site_bounds_buffered[3,])[3], ext(site_bounds_buffered[3,])[4]),
           clip = "off") +
  theme_nothing() +
  theme(plot.background = element_rect(fill = "black"))
save_plot("figures/not_included_in_manuscript/site_3_rgb.png", plot_site_3,
          base_asp = 1)
