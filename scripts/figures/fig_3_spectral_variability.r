# Figure 3 Spectral variability and species richness
# Jakob J. Assmann jakob.assmann@uzh.ch February 2026

# Dependencies
library(tidyverse)
library(sf)
library(ggplot2)
library(terra)
library(tidyterra)
library(cowplot)

## Load data and model

# Load richness and multispectral cv data
multispec_cv <- read_csv("data/combined_obs_multispec_cv.csv")

# Load plot locations
plot_coords <- read_sf("data/geometries/plot_locations.gpkg")

# Buffer plots by 0 5 m
plots_buffered <- plot_coords %>%
    st_buffer(0.5) %>%
    mutate(site = substr(Name, 1, 2) %>% gsub("S", "Site ", .))

# Load linear model (spp_richness_lm) and generate predictions
load("data/analysis/lm_nspp_multispectral_cv.Rda")
multispec_cv <- bind_cols(
    multispec_cv,
    predict(spp_richness_lm, interval = "confidence")
)

# Load file names of multispec drone rasters
multispec_rasts <- data.frame(
    file_name = list.files(
        "data/drone_data/multispec/2026/",
        full.names = T,
        pattern = ".tif$"
    )
) %>%
    mutate(
        site = gsub(".*(site_[1-3]).*", "\\1", file_name),
        band = gsub(".*_(.*)\\..*", "\\1", file_name)
    ) %>%
    mutate(
        wavelength = case_when(
            band == "blue" ~ 475,
            band == "green" ~ 560,
            band == "red" ~ 668,
            band == "rededge" ~ 717,
            band == "nir" ~ 842,
            TRUE ~ NA
        )
    )


# Plott richness (nspp) vs. spectral variation (multispec cv)
(n_spp_vs_spect_cv_plot <- ggplot(
    multispec_cv,
    aes(x = spectral_cv, y = n_spp)
) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#cdcdcd66") +
    geom_line(aes(y = fit), colour = "#962c7e", linewidth = 2) +
    geom_point(aes(shape = site_id, color = spectral_cv), size = 5) +
    scale_color_viridis_c(
        option = "magma",
        limits = c(0.2, 1),
        na.value = "white",
        guide = "none"
    ) +
    scale_x_continuous(
        limits = c(0.2, 0.9),
        breaks = seq(0.2, 0.9, 0.1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_continuous(limits = c(0, 32), expand = c(0, 0)) +
    labs(
        x = "Spectral Variability (CV)",
        y = "Species Richness (n)",
        colour = "",
        shape = "",
    ) +
    theme_cowplot(18) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.125),
        legend.title = element_blank(),
        legend.box.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        legend.box.background = element_rect(colour = "black")
    ))

# Prepare examples for plots at extremes of CV scale
# Get IDs for extreme plots
plot_id_min <- multispec_cv[
    which(multispec_cv$spectral_cv == min(multispec_cv$spectral_cv)),
]$plot_id
plot_id_max <- multispec_cv[
    which(multispec_cv$spectral_cv == max(multispec_cv$spectral_cv)),
]$plot_id
# Get geometeries
plot_min <- plots_buffered %>% filter(Name == plot_id_min)
plot_max <- plots_buffered %>% filter(Name == plot_id_max)
# Histograms to help with colour scaling
crop(rast(multispec_rasts[c(17, 14, 13), 1]), st_buffer(plot_min, 0.3)) %>%
    hist()
crop(rast(multispec_rasts[c(11, 8, 7), 1]), st_buffer(plot_max, 0.3)) %>%
    hist()

# Example 1: Drone Image (true colour)
example_1_drone <- ggplot() +
    geom_spatraster_rgb(
        data = crop(
            rast(multispec_rasts[c(17, 14, 13), 1]),
            st_buffer(plot_min, 0.3)
        ),
        zlim = c(0.1, 0.3),
        max_col_value = 0.3,
        interpolate = FALSE
    ) +
    geom_sf(
        data = st_geometry(plot_min),
        color = "white",
        fill = NA,
        linewidth = 2
    ) +
    annotate(
        "segment",
        x = 482788.8,
        xend = 482789.3,
        y = 7667984.05,
        yend = 7667984.05,
        linewidth = 3,
        colour = "#595959"
    ) +
    annotate(
        "text",
        x = 482789.05,
        y = 7667984.2,
        label = "0.5 m",
        colour = "black",
        size = 14 / .pt
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_map() +
    theme(axis.title = element_blank())

# Example 2: Drone Image (true colour)
example_2_drone <- ggplot() +
    geom_spatraster_rgb(
        data = crop(
            rast(multispec_rasts[c(11, 8, 7), 1]),
            st_buffer(plot_max, 0.3)
        ),
        zlim = c(0.01, 0.10),
        max_col_value = 0.1,
        interpolate = FALSE
    ) +
    geom_sf(
        data = st_geometry(plot_max),
        color = "white",
        fill = NA,
        linewidth = 2
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_map() +
    theme(axis.title = element_blank())

# Example 1: Statistics
example_1_stats <- plot_grid(
    # Species Richness of example plot 1
    ggplot(filter(multispec_cv, plot_id == plot_id_min)) +
        geom_col(aes(x = " ", y = n_spp)) +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
        theme_cowplot(18) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ),
    # Spectral Variation of example plot 1
    ggplot(filter(multispec_cv, plot_id == plot_id_min)) +
        geom_col(aes(x = "", y = spectral_cv)) +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_cowplot(18) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ),
    align = "v"
)

# Example 2: Statistics
example_2_stats <- plot_grid(
    # Species Richness of example plot 2
    ggplot(filter(multispec_cv, plot_id == plot_id_max)) +
        geom_col(aes(x = " ", y = n_spp)) +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
        theme_cowplot(18) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ),
    # Spectral variation of example plot 2
    ggplot(filter(multispec_cv, plot_id == plot_id_max)) +
        geom_col(aes(x = "", y = spectral_cv)) +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme_cowplot(18) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ),
    align = "v"
)

# Composite of drone imagery and stats for both examples
(plots_example <- plot_grid(
    # Label plot
    plot_grid(
        ggplot() +
            labs(x = "Drone\nImagery") +
            theme_nothing(18) +
            theme(axis.title.x = element_text(vjust = 2, size = 18)),
        ggplot() +
            labs(x = "Species\nRichness (n)") +
            theme_nothing(18) +
            theme(axis.title.x = element_text(vjust = 2, size = 18)),
        ggplot() +
            labs(x = "Spectral\nVariability (CV)") +
            theme_nothing(18) +
            theme(axis.title.x = element_text(vjust = 2, size = 18)),
        nrow = 1,
        rel_widths = c(0.4, 0.3, 0.3),
        align = ""
    ),
    # Example plot 1
    plot_grid(
        example_1_drone,
        example_1_stats,
        nrow = 1,
        rel_widths = c(0.4, 0.6)
    ),
    # Example Plot 2
    plot_grid(
        example_2_drone,
        example_2_stats,
        nrow = 1,
        rel_widths = c(0.4, 0.6)
    ),
    NULL,
    ncol = 1,
    nrow = 4,
    rel_heights = c(0.13, 0.45, 0.45, 0.07)
))
# Combined plot with examples and all
combined_plot_nspp_cv <- plot_grid(
    # n spp vs CV
    n_spp_vs_spect_cv_plot,
    # Example plots
    plots_example,
    ncol = 2,
    nrow = 1,
    rel_widths = c(1.1, 0.9),
    align = "",
    labels = c("a)", "b)"),
    label_size = 18,
    hjust = c(-0.5, 0.75)
)

# Add icons
(combined_plot_nspp_cv <- ggdraw() +
    draw_plot(combined_plot_nspp_cv) +
    # Plant observation icons
    # draw_image("figures/glasses.png", scale = 0.065, x = -0.41, y = +0.425) +
    draw_image(
        "figures/vascular_plant_bw.png",
        scale = 0.12,
        x = -0.425,
        y = +0.425
    ) +
    draw_image("figures/moss_bw.png", scale = 0.12, x = -0.37, y = +0.425) +
    # Drone Icon
    draw_image("figures/drone.png", scale = 0.125, x = 0.02, y = -0.275))

# Write out as file
save_plot(
    "figures/figure_3.png",
    combined_plot_nspp_cv,
    base_asp = 1.2,
    ncol = 2,
    base_height = 6,
    bg = "white"
)
