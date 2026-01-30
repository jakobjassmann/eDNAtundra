# Figure 1 - Map and sampling overview
# Jakob J. Assmann jakob.assmann@uzh.ch 30 January 2026

# Dependencies
library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)

# Load Sentinel-2 imagery
sentinel_2 <- rast("data/sentinel/2023-08-08-Sentinel-2_L2A_true_color_-_Cambridge_Bay.tiff")

# Load site locations
site_locations <- read_sf("data/geometries/site_boundaries.gpkg")

# Cambridge Bay
cambridge_bay <- matrix(c(
    496808.045, 7668534.366,
    496808.045, 7666812.132,
    498572.031, 7666801.694,
    498572.031, 7668534.366,
    496808.045, 7668534.366
), byrow = T, ncol = 2) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(crs = st_crs(sentinel_2))

# Plot
site_map <- ggplot() +
    # Sentinel-2 background layers
    geom_spatraster_rgb(data = sentinel_2) +
    # Site locations
    geom_sf(
        data = site_locations,
        colour = "white",
        linewidth = 2,
    ) +
    geom_sf_text(
        data = st_centroid(site_locations),
        nudge_x = 1300,
        nudge_y = 1100,
        aes(label = site),
        colour = "white",
        size = 22 / .pt
    ) +
    # Cambridge Bay
    geom_sf(
        data = cambridge_bay,
        fill = NA,
        linewidth = 1,
        colour = "white"
    ) +
    annotate("text",
        x = st_coordinates(st_centroid(cambridge_bay))[1],
        y = st_coordinates(st_centroid(cambridge_bay))[2] + 1750,
        hjust = 0.85,
        label = "Cambridge Bay",
        colour = "white",
        size = 18 / .pt
    ) +
    # Scale
    annotate("segment",
        x = 480000, xend = 485000,
        y = 7659500, yend = 7659500,
        colour = "white",
        linewidth = 2
    ) +
    annotate("text",
        x = (485000 + 480000) / 2, y = 7660500,
        label = "5 km",
        color = "white",
        size = 18 / .pt
    ) +
    # North Arrow
    annotate("segment",
        x = 488000, xend = 488000,
        y = 7659000, yend = 7661000,
        colour = "white",
        linewidth = 1.5,
        arrow = arrow()
    ) +
    annotate("text",
        x = 489500, y = 7660000,
        label = "N",
        size = 22 / .pt,
        colour = "white"
    ) +
    coord_sf(
        xlim = c(478000, 517000),
        ylim = c(7659000, 7678000)
    ) +
    theme_map()

# Generate overview map / globe insert
source("scripts/figures/fig_1_globe_insert.R")

# Combine with the conceptual sketches
plot_grid(
    ggdraw() + 
    draw_plot(site_map) + 
    draw_image("figures/overview_map.png", scale = 0.5, x = 0.35, y = -0.2),
    ggdraw() + 
    draw_image("figures/sketches_site_plot_layout.png", clip = "on", scale = 1.225),
    nrow = 2,
    labels = c("a)", ""),
    label_x = 0.025,
    label_y = 0.975,
    label_size = 18,
    label_colour = "white") %>%
save_plot("figures/figure1_map_sketch.png", ., nrow = 2,
    base_asp = (517000 - 478000) / (7678000 - 7659000), base_height = 5, bg = "white")
