## Fig 4 Landscape-level species richness
# Jakob Assmann jakob.assmann@gu.se Updated 30 July 2026

library(tidyverse)
library(ggplot2)
library(terra)
library(tidyterra)
library(sf)
library(cowplot)
library(rayshader)
library(magick)

# Load species richness model
load("data/analysis/lm_nspp_multispectral_cv.Rda")

# Load plot locations
plot_coords <- read_sf("data/geometries/plot_locations.gpkg")

# Load site boundaries
site_bounds <- read_sf("data/geometries/site_boundaries.gpkg") %>%
    mutate(
        site = case_when(
            site == "Site 1" ~ "site_1",
            site == "Site 2" ~ "site_2",
            site == "Site 3" ~ "site_3"
        )
    )

# Load drone rasters
# RGB
rgb_rasts <- list.files(
    "data/drone_data/rgb/",
    full.names = T,
    pattern = ".tif$"
) %>%
    .[!grepl("crop", .)] %>%
    map(rast)
rgb_rasts_agg <- list.files(
    "data/drone_data/rgb/",
    full.names = T,
    pattern = ".tif$"
) %>%
    .[!grepl("crop", .)] %>%
    map(function(file_name) {
        cat(file_name, "\n")
        site_id <- gsub(".*(site_[0-3]).*", "\\1", file_name)
        cat(site_id, "\n")
        site_rast <- crop(rast(file_name), filter(site_bounds, site == site_id))
        site_rast <- aggregate(
            site_rast,
            fact = round(1 / res(site_rast)),
            fun = mean,
            na.rm = T
        ) %>%
            crop(filter(site_bounds, site == site_id))
        return(site_rast)
    })
# Multispec file names only
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
# DSMs
dsm_rasts <- list.files(
    "data/drone_data/dsm/",
    full.names = T,
    pattern = ".tif$"
) %>%
    map(function(file_name) {
        cat(file_name, "\n")
        site_id <- gsub(".*(site_[0-3]).*", "\\1", file_name)
        cat(site_id, "\n")
        site_rast <- crop(rast(file_name), filter(site_bounds, site == site_id))
        site_rast
        site_rast <- aggregate(
            site_rast,
            fact = round(1 / res(site_rast)),
            fun = mean
        ) %>%
            setNames("dsm") %>%
            crop(filter(site_bounds, site == site_id))
        return(site_rast)
    })

# Aggregate multispec rasters to CV (1 m rsolution)
spec_cv_rast <- multispec_rasts %>%
    filter(band != "ndvi") %>%
    arrange(site) %>%
    group_by(site) %>%
    group_split() %>%
    map(function(site_rasts) {
        site_rast <- rast(site_rasts$file_name) %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
        site_rast <- aggregate(
            site_rast,
            fact = c(round(1 / res(site_rast)), nlyr(site_rast)),
            fun = function(x) {
                sd(x, na.rm = T) / mean(x, na.rm = T)
            }
        ) %>%
            setNames("spectral_cv") %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
    })

# Load plot-level data
multispec_cv <- read_csv("data/combined_obs_multispec_cv.csv")

# # Override buggy cowplot::get_legend if needed
# get_legend <- function(plot) {
#     legends <- get_plot_component(plot, "guide-box", return_all = T)
#     legends_new <- list()
#     for (i in 1:length(legends)) {
#         if (!("zeroGrob" %in% class(legends[[i]]))) {
#             legends_new <- c(legends_new, list(legends[[i]]))
#         }
#     }
#     if (length(legends_new) > 0) {
#         do.call(cbind, legends_new)
#     } else {
#         return(legends[[1]])
#     }
# }
# (test_plot <- ggplot() +
#     geom_point(
#         aes(
#             x = 1:10,
#             y = 1:10,
#             colour = 1:10,
#             size = 1:10
#         ),
#         shape = 21
#     )) +
#     guides(
#         colour = guide_colourbar(position = "bottom"),
#         size = guide_legend(position = "top")
#     )
# ggdraw() + draw_grob(cowplot::get_legend(test_plot))
# get_legend(ggplot())

# Generate water / snow masks (NDVI < 0)
water_masks <- multispec_rasts %>%
    filter(band == "ndvi") %>%
    arrange(site) %>%
    group_by(site) %>%
    group_split() %>%
    map(function(site_rasts) {
        site_rast <- rast(site_rasts$file_name) %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
        m <- c(-1, 0, NA, 0, 1, 1)
        m <- matrix(m, ncol = 3, byrow = TRUE)
        site_rast <- aggregate(
            site_rast,
            fact = c(round(1 / res(site_rast)), nlyr(site_rast)),
            fun = mean
        ) %>%
            classify(m) %>%
            setNames("water") %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
    })


# Generate plots
spectra_div_maps <- map(1:3, function(x) {
    # Subset data
    plots <- full_join(
        multispec_cv,
        select(plot_coords, Name, geom),
        by = c("plot_id" = "Name")
    ) %>%
        filter(grepl(x, site_id))
    cv_rast <- spec_cv_rast[[x]] %>% mask(water_masks[[x]])

    # Visualise
    div_plot <- ggplot() +
        geom_spatraster(data = cv_rast) +
        # geom_spatraster_contour(data = cv_rast) +
        scale_fill_viridis_c(
            option = "magma",
            limits = c(0.2, 1),
            na.value = "white",
            guide = guide_colorbar(
                frame.colour = "black",
                ticks.colour = "black"
            )
        ) +
        geom_sf(
            data = st_as_sf(plots),
            aes(fill = spectral_cv, size = n_spp),
            shape = 21,
            colour = "dodgerblue",
            stroke = 2
        ) +
        scale_size_continuous(
            limits = c(0, 30),
            breaks = seq(0, 30, 10),
            range = c(0.1, 5)
        ) +
        labs(
            title = paste("Site", x),
            fill = "Spectral Variability (CV)",
            size = "Plot Species Richness (n)"
        ) +
        theme_map(18) +
        theme(
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            plot.title = element_text(hjust = 0.5)
        )
})

spect_div_plot <- plot_grid(
    spectra_div_maps[[1]] + theme(legend.position = "none"),
    spectra_div_maps[[2]] + theme(legend.position = "none"),
    spectra_div_maps[[3]] + theme(legend.position = "none"),
    ggdraw() + draw_grob(get_legend(spectra_div_maps[[1]]), x = 0.125, y = 0),
    ncol = 2
)

save_plot(
    "figures/figure_SX_spect_div_maps.png",
    spect_div_plot,
    base_asp = 1,
    ncol = 2,
    nrow = 2,
    bg = "white"
)

# # Combined with scatter
# save_plot(
#     "figures/figure_3_spect_div_maps_with_scat.png",
#     plot_grid(
#         n_spp_vs_spect_cv_plot +
#             labs(title = " ") +
#             theme_cowplot(18) +
#             theme(
#                 legend.position = "none",
#                 plot.margin = margin(7, 14, 7, 7, "points")
#             ),
#         spectra_div_maps[[1]] + theme(legend.position = "none"),
#         spectra_div_maps[[2]] + theme(legend.position = "none"),
#         ggplot() +
#             annotate(
#                 "segment",
#                 x = 0.275,
#                 xend = 0.725,
#                 y = 0.35,
#                 yend = 0.35,
#                 linewidth = 1.5
#             ) +
#             annotate(
#                 "text",
#                 x = 0.5,
#                 y = 0.4,
#                 label = "150 m",
#                 size = 20 / .pt
#             ) +
#             annotate(
#                 "segment",
#                 x = 0.45,
#                 xend = 0.45,
#                 y = 0.50,
#                 yend = 0.65,
#                 arrow = arrow(),
#                 linewidth = 1.5
#             ) +
#             annotate(
#                 "text",
#                 x = 0.525,
#                 y = 0.575,
#                 label = "N",
#                 size = 32 / .pt
#             ) +
#             coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#             theme_nothing(),
#         spectra_div_maps[[3]] + theme(legend.position = "none"),
#         ggdraw() + draw_grob(get_legend(spectra_div_maps[[1]]), x = 0.2, y = 0),
#         ncol = 3,
#         nrow = 2,
#         labels = c("a)", "b)", "", "", "", ""),
#         label_size = 18
#     ),
#     base_asp = 1,
#     ncol = 3,
#     nrow = 2,
#     base_height = 5,
#     bg = "white"
# )

# Contour plots of predicted species richness
(est_rich_maps <- map(1:3, function(x) {
    # Subset data
    plots <- full_join(
        multispec_cv,
        select(plot_coords, Name, geom),
        by = c("plot_id" = "Name")
    ) %>%
        filter(grepl(x, site_id))
    cv_rast <- spec_cv_rast[[x]] %>% mask(water_masks[[x]])

    # Predict number of species
    n_spp_rast <- predict(cv_rast, spp_richness_lm)
    # round to nearest integer
    n_spp_rast <- round(n_spp_rast)

    # Classify raster
    breaks <- matrix(
        c(
            0,
            5,
            1,
            5,
            10,
            2,
            10,
            15,
            3,
            15,
            20,
            4,
            20,
            25,
            5,
            25,
            30,
            6,
            30,
            35,
            7
        ),
        ncol = 3,
        byrow = T
    )
    n_spp_rast <- classify(n_spp_rast, breaks, right = FALSE)
    levels(n_spp_rast) <- data.frame(
        id = 1:7,
        n_species = c(
            "0-4",
            "5-9",
            "10-14",
            "15-19",
            "20-24",
            "25-29",
            "30-34"
        )
    )

    # Visualise
    est_rich_plot <- ggplot() +
        geom_spatraster(data = n_spp_rast, show.legend = c(fill = TRUE)) +
        # scale_fill_viridis_c(limits = c(0, 35), breaks = seq(0,35,5)) +
        # geom_spatraster_contour_filled(data = (predict(cv_rast, spp_richness_lm)),
        #                         breaks = seq(0,35,5)) +
        scale_fill_viridis_d(
            na.value = "white",
            limits = factor(c(levels(n_spp_rast)[[1]]$n_species[-1], NA)),
            labels = c(
                "5 - 9",
                "10 - 14",
                "15 - 19",
                "20 - 24",
                "25 - 29",
                "30 - 34",
                "Water/Snow"
            )
        ) +
        labs(
            fill = "Estimated\nSpecies\nRichness (n)" #,
        ) +
        theme_map(18) +
        theme(
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            plot.title = element_blank(),
            plot.margin = unit(c(0, 0, 0.1, 0), "inch"),
            legend.key = element_rect(colour = "black")
            #       plot.title = element_text(hjust = 0.5))
        )
}))

# Generate RGB plots
rgb_plots <- map(1:3, function(x) {
    ggplot() +
        geom_spatraster_rgb(data = crop(rgb_rasts[[x]], site_bounds[x, ])) +
        labs(title = paste0("Site ", x)) +
        theme_map(18) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            plot.margin = unit(c(0.1, 0, 0, 0), "inch")
        )
})

est_rich_plot <- plot_grid(
    plot_grid(
        rgb_plots[[1]],
        est_rich_maps[[1]] + theme(legend.position = "none"),
        align = "hv",
        ncol = 1,
        labels = c("a)", "d)"),
        label_size = 18,
        hjust = -0.8
    ),
    plot_grid(
        rgb_plots[[2]],
        est_rich_maps[[2]] + theme(legend.position = "none"),
        align = "hv",
        ncol = 1,
        labels = c("b)", "e)"),
        label_size = 18,
        hjust = -0.8
    ),
    plot_grid(
        rgb_plots[[3]],
        est_rich_maps[[3]] + theme(legend.position = "none"),
        align = "hv",
        ncol = 1,
        labels = c("c)", "f)"),
        label_size = 18,
        hjust = -0.8
    ),
    plot_grid(
        ggplot() +
            annotate(
                "text",
                x = 0.18 + 0.29,
                y = 0.425,
                label = "150 m",
                # angle = 90,
                size = 18 / .pt
            ) +
            annotate(
                "segment",
                x = 0.18,
                xend = 0.765,
                y = 0.5,
                yend = 0.5,
                linewidth = 2
            ) +
            annotate(
                "text",
                x = 0.18,
                y = 0.78,
                label = "True Colour\nImage",
                size = 18 / .pt,
                hjust = 0
            ) +
            annotate(
                "segment",
                x = 0.3,
                xend = 0.3,
                y = 0.1,
                yend = 0.25,
                arrow = arrow(),
                linewidth = 1.5
            ) +
            annotate(
                "text",
                x = 0.45,
                y = 0.175,
                label = "N",
                size = 32 / .pt
            ) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            theme_map(18) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "inch")),
        ggdraw() +
            draw_grob(
                cowplot::get_legend(est_rich_maps[[3]]),
                x = 0.2,
                y = -0.05
            ),
        nrow = 2,
        ncol = 1
    ),
    rel_widths = c(1, 1, 1, 0.8),
    align = "hv",
    ncol = 4
)

# Save
save_plot(
    "figures/figure_4_est_spp_richness_maps.png",
    est_rich_plot,
    base_asp = 0.9,
    ncol = 4,
    nrow = 2,
    bg = "white"
)

## 3D plots
Sys.setenv(RAYON_NUM_THREADS = 20)
rgb_plots_3d <- map(1:3, function(x) {
    # Get elevation raster
    elmat <- raster_to_matrix(dsm_rasts[[x]])
    # Prepare rgb raster
    rgb_rast <- rgb_rasts_agg[[x]]
    rgb_rast[[1]][] <- rgb_rast[[1]][] / 255
    rgb_rast[[2]][] <- rgb_rast[[2]][] / 255
    rgb_rast[[3]][] <- rgb_rast[[3]][] / 255
    #plotRGB(rgb_rast, scale = 1)

    rgb_mat <- as.array(rgb_rast)[,, 1:3]
    rgb_mat <- rgb_mat^1.4
    if (sum(is.na(rgb_mat)) > 0) {
        rgb_mat[is.na(rgb_mat)] <- 1
    }
    #plot_map(rgb_mat)

    # Initate 3D plot
    plot_3d(
        rgb_mat,
        elmat,
        zscale = 0.5,
        zoom = 0.7, #0.85,
        phi = 45,
        theta = 45,
        windowsize = 2000
    )
    # Set outfile and render in high quality
    out_file <- tempfile(fileext = ".png")
    render_highquality(
        out_file,
        samples = 12, # Up to 128 for best quality
        lightdirection = c(60, 150),
        lightintensity = 375,
        clear = TRUE,
        parallel = TRUE,
        interactive = FALSE,
        ground_size = 0
    )
    rgl::close3d()

    # Load image for post-processing
    img <- image_read(out_file)
    # Turn black background to white
    mask <- image_transparent(img, "black")
    img_masked <- image_background(mask, "white", flatten = TRUE)
    # Downscale for antialisaing
    img_scaled <- image_resize(
        img_masked,
        "1000x1000",
        filter = "Lanczos"
    )
    # Write out imagge
    image_write(img_scaled, out_file)

    # Re-load and convert to ggplot object
    (rast_plot_3d <- (ggdraw() + draw_image(out_file)))

    # tidy up
    file.remove(out_file)

    # Return file
    return(rast_plot_3d)
})
n_spp_plots_3d <- map(1:3, function(x) {
    n_spp_rast <- predict(spec_cv_rast[[x]], spp_richness_lm)

    # Classify raster
    breaks <- matrix(
        c(
            0,
            5,
            1,
            5,
            10,
            2,
            10,
            15,
            3,
            15,
            20,
            4,
            20,
            25,
            5,
            25,
            30,
            6,
            30,
            35,
            7
        ),
        ncol = 3,
        byrow = T
    )
    n_spp_rast <- classify(n_spp_rast, breaks)
    levels(n_spp_rast) <- data.frame(
        id = 1:7,
        n_species = c(
            "0-5",
            "5-10",
            "10-15",
            "15-20",
            "20-25",
            "25-30",
            "30-35"
        )
    )
    # Get elevation raster
    elmat <- raster_to_matrix(dsm_rasts[[x]])
    # Prepare water mask
    water_mask <- water_masks[[x]]
    # water_mask_inv_array <- water_mask
    # water_mask_inv_array[] <- NA
    # water_mask_inv_array[] <- is.na(water_mask[])
    # water_mask_inv_array <- as.array(water_mask_inv_array)
    water_mask_inv <- is.na(water_mask[])
    water_col <- col2rgb("white")
    # Prepare color scale
    viridis_cols <- col2rgb(c(
        "#000000",
        "#470052",
        "#494084",
        "#34768d",
        "#22a785",
        "#72d25c",
        "#f6eb46"
    )) /
        255
    # Prepare raster layers
    n_spp_mat_red <- n_spp_rast[[1]]
    n_spp_mat_green <- n_spp_rast[[1]]
    n_spp_mat_blue <- n_spp_rast[[1]]
    n_spp_mat_red[] <- NA
    n_spp_mat_green[] <- NA
    n_spp_mat_blue[] <- NA
    n_spp_mat_red[] <- viridis_cols[1, as.numeric(n_spp_rast$n_species[])]
    n_spp_mat_green[] <- viridis_cols[2, as.numeric(n_spp_rast$n_species[])]
    n_spp_mat_blue[] <- viridis_cols[3, as.numeric(n_spp_rast$n_species[])]
    n_spp_mat_red[water_mask_inv] <- water_col[1] / 255
    n_spp_mat_green[water_mask_inv] <- water_col[2] / 255
    n_spp_mat_blue[water_mask_inv] <- water_col[3] / 255
    n_spp_mat_rgb <- c(n_spp_mat_red, n_spp_mat_green, n_spp_mat_blue)
    plotRGB(n_spp_mat_rgb, scale = 1)
    n_spp_mat_rgb <- as.array(n_spp_mat_rgb)
    # Gamma correction
    n_spp_mat_rgb <- n_spp_mat_rgb^1.4
    plot_map(n_spp_mat_rgb)

    # Initate 3D plot
    plot_3d(
        n_spp_mat_rgb,
        elmat,
        zscale = 0.5,
        zoom = 0.7, #.85
        phi = 45,
        theta = 45,
        windowsize = 2000
    )
    # Set outfile and render in high quality
    out_file <- tempfile(fileext = ".png")
    render_highquality(
        out_file,
        samples = 12, # Up to 128 for best quality
        lightdirection = c(60, 150),
        lightintensity = 375,
        clear = TRUE,
        parallel = TRUE,
        interactive = FALSE,
        ground_size = 0
    )
    rgl::close3d()

    # Load image for post-processing
    img <- image_read(out_file)
    # Turn black background to white
    mask <- image_transparent(img, "black")
    img_masked <- image_background(mask, "white", flatten = TRUE)
    # Downscale for antialisaing
    img_scaled <- image_resize(
        img_masked,
        "1000x1000",
        filter = "Lanczos"
    )
    # Write out imagge
    image_write(img_scaled, out_file)

    # Re-load and convert to ggplot object
    (rast_plot_3d <- (ggdraw() + draw_image(out_file)))

    # tidy up
    file.remove(out_file)

    # Return file
    return(rast_plot_3d)
})

# Annotate plots and save
{
    # Create copies for annotation
    rgb_plots_3d_annotated <- rgb_plots_3d
    n_spp_plots_3d_annotated <- n_spp_plots_3d

    # Scale bar
    rgb_plots_3d_annotated[[3]] <-
        rgb_plots_3d_annotated[[3]] +
        draw_label(
            "150 m",
            x = 0.75,
            y = 0.1925,
            hjust = 0.5,
            size = 42,
            angle = 37,
            fontface = "bold",
            vjust = 1.5
        ) +
        draw_line(
            x = c(0.65, 0.84),
            y = c(0.1, 0.29),
            # arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
            color = "black",
            linewidth = 2.5
        )
    # North Arrow
    n_spp_plots_3d_annotated[[3]] <-
        n_spp_plots_3d_annotated[[3]] +
        draw_label(
            "N",
            x = 0.84,
            y = 0.14,
            hjust = 0.5,
            size = 42,
            #angle = 37,
            fontface = "bold"
        ) +
        draw_line(
            x = c(0.77, 0.85),
            y = c(0.14, 0.23),
            arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
            color = "black",
            linewidth = 2.5
        )
    plot_grid(
        plotlist = c(
            rgb_plots_3d_annotated,
            list(NULL),
            n_spp_plots_3d_annotated,
            list(
                ggdraw() +
                    draw_grob(
                        cowplot::get_legend(
                            est_rich_maps[[3]] +
                                theme_map(32) +
                                theme(
                                    legend.key = element_rect(colour = "black")
                                )
                        ),
                        x = 0,
                        y = 0,
                        scale = 1
                    )
            )
        ),
        labels = c(
            " a)  Site 1",
            " b)  Site 2",
            " c)  Site 3",
            "",
            " d)",
            " e)",
            " f)",
            ""
        ),
        rel_widths = c(1, 1, 1, 0.4),
        hjust = 0,
        nrow = 2,
        ncol = 4,
        label_size = 36
    ) %>%
        save_plot(
            "figures/figure_4_drone_images.png",
            .,
            base_asp = 1.195,
            nrow = 2,
            ncol = 4,
            base_height = 6,
            bg = "white"
        )
}
