# Extract spectral diversity from drone data for plots
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(terra)
library(tidyterra)
library(sf)
library(NbClust)
library(pblapply)
library(cowplot)
library(factoextra)
library(tidyterra)
library(ggnewscale)

# Set up parallel envrionment
cl <- parallel::makeCluster(5)

# Set site colours
site_colours <- data.frame(
        site = c("Site 1", "Site 2", "Site 3"),
        value = c(
            "#F2D64B",
            "#5F1773",
            "#F2507B"
        ))

# Load plot locations
plot_coords <- read_sf("data/geometries/plot_locations.gpkg")

# Load site boundaries
site_bounds <- read_sf("data/geometries/site_boundaries.gpkg") %>%
    mutate(site = case_when(
        site == "Site 1" ~ "site_1",
        site == "Site 2" ~ "site_2",
        site == "Site 3" ~ "site_3"
    ))

# Load drone rasters
# RGB
rgb_rasts <- list.files("data/drone_data/rgb/", full.names = T, pattern = ".tif$") %>%
    .[!grepl("crop", .)] %>%
    map(rast)
# Multispec file names only
multispec_rasts <- data.frame(
    file_name = list.files("data/drone_data/multispec/", full.names = T, pattern = ".tif$")
) %>%
    mutate(
        site = gsub(".*(site_[1-3]).*", "\\1", file_name),
        band = gsub(".*_(.*)\\..*", "\\1", file_name)
    ) %>%
    mutate(wavelength = case_when(
        band == "blue" ~ 475,
        band == "green" ~ 560,
        band == "red" ~ 668,
        band == "rededge" ~ 717,
        band == "nir" ~ 842,
        TRUE ~ NA
    ))

# Buffer plots by 0 5 m
plots_buffered <- plot_coords %>%
    st_buffer(0.5) %>%
    mutate(site = substr(Name, 1, 2) %>% gsub("S", "Site ", .))

# Extract RGB data for a 0.5 m buffer around each plot
plot_rgb <- map(1:3, function(x) {
    site_plots <- filter(plots_buffered, site == paste0("Site ", x))
    plots_rgb <- terra::extract(rgb_rasts[[x]], site_plots) %>%
        select(1:4) %>%
        setNames(c("ID", "R", "G", "B")) %>%
        pivot_longer(2:3, names_to = "band", values_to = "value") %>%
        group_by(ID) %>%
        summarise(spectral_cv = (sd(value) / mean(value)))
    plots_rgb$plot_id <- site_plots$Name
    plots_rgb$site_id <- site_plots$site
    plots_rgb <- relocate(plots_rgb, site_id, plot_id)
    return(plots_rgb)
}) %>% bind_rows()

# Extract multispec data for a 0.5 m buffer around each plot
plot_multispec <- map(1:3, function(x) {
    # Subset geometries to site
    site_plots <- filter(plots_buffered, site == paste0("Site ", x))
    # Extract raster values
    plots_multispec <- filter(multispec_rasts, grepl(x, site)) %>%
        split(., 1:nrow(.)) %>%
        map(function(y) {
            terra::extract(rast(y$file_name), site_plots) %>%
                select(1:2) %>%
                setNames(c("ID", y$band))
        }) %>%
        reduce(bind_cols) %>%
        tibble() %>%
        select(c(1, seq(2, 12, 2)))
    # Match plot names
    plots_multispec <- site_plots %>%
        mutate(., ID = 1:nrow(.)) %>%
        select(ID, Name) %>%
        full_join(plots_multispec, by = c("ID" = "ID...1")) %>%
        rename(plot_id = Name) %>%
        select(-ID) %>%
        st_drop_geometry()
    # Set site_id
    plots_multispec$site_id <- unique(site_plots$site)
    # Tidy up
    plots_multispec <- relocate(plots_multispec, site_id, plot_id)
    return(plots_multispec)
}) %>% bind_rows()

# Load presence / absence vascular observations
vascular_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    ) %>%
    # Add site_id string
    rename(site_id = site) %>%
    select(-date, -full_name) %>%
    mutate(phylum = "Vascular plants")

# Load presence / absence moss observations
bryo_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv") %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # convert al sp. to species
    mutate(species = (gsub("sp\\.", "species", species))) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) & !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    arrange(site_id, plot_id, short_name) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(phylum = "Bryophytes")

# Combine dataframes
combined_obs <- bind_rows(vascular_obs, bryo_obs) %>%
    arrange(site_id, plot_id, short_name, phylum) %>%
    relocate(site_id) %>%
    group_by(plot_id) %>%
    summarise(n_spp = n_distinct(short_name)) %>%
    full_join(plot_rgb)

# Point framing
point_framing <- read_csv("data/botanical_obs/pointframing.csv") %>%
    select(-observer, -scribe, -date) %>%
    group_by(plot_id, full_name) %>%
    summarise(sum_hits = n()) %>%
    split(., .$plot_id) %>%
    map(function(x) {
        x <- mutate(x,
            p = sum_hits / sum(x$sum_hits),
            p_ln_p = p * log(p)
        ) %>%
            summarise(shan_div = -1 * sum(p_ln_p))
    }) %>%
    bind_rows()

# Combine
combined_obs <- full_join(combined_obs, point_framing)

## RGB based anlaysis

# Correlate the the species richness with the spectral diversity
ggplot() +
    geom_point(data = combined_obs, aes(x = n_spp, spectral_cv)) +
    theme_classic()

# Correlate the the cover richness with the spectral diversity
ggplot(data = combined_obs, mapping = aes(x = shan_div, spectral_cv)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic()

## Multispectral analysis
multispec_cv <- plot_multispec %>%
    select(-ndvi) %>%
    pivot_longer(3:7, names_to = "band", values_to = "reflectance") %>%
    group_by(site_id, plot_id) %>%
    summarise(spectral_cv = sd(reflectance) / mean(reflectance)) %>%
    full_join(., select(combined_obs, plot_id, n_spp, shan_div))

# Test for correlation spectral cv and shannon diversity (point frame)
ggplot(multispec_cv, aes(x = spectral_cv, y = shan_div)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic()
hist(multispec_cv$shan_div)
hist(multispec_cv$spectral_cv)
cor(multispec_cv$spectral_cv, multispec_cv$shan_div)
cor.test(multispec_cv$spectral_cv, multispec_cv$shan_div)
summary(lm(shan_div ~ spectral_cv, multispec_cv))

# Test for correlation spectral cv and n_spp (abundance)
(n_spp_vs_spect_cv_plot <- ggplot(multispec_cv, aes(x = spectral_cv, y = n_spp)) +
    geom_smooth(method = "lm", colour = "grey20", fill = "grey60") +
    geom_point(aes(colour = site_id), size = 3) +
    scale_color_manual(
        values = site_colours$value,
        labels = c("Site 1", "Site 2", "Site 3")
    ) +
    scale_x_continuous(limits = c(0.2, 0.9),
                    breaks = seq(0.2,0.9,0.1),
                    expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 32), expand = c(0,0)) +
    labs(
        x = "Spectral Variation (CV)", y = "Species Richness (n)",
        colour = ""
    ) +
    theme_cowplot())
save_plot("figures/figure_3_n_spp_vs_spec_cv.png",
            n_spp_vs_spect_cv_plot,
            base_asp = 1.3)

hist(multispec_cv$n_spp)
hist(multispec_cv$spectral_cv)
cor(multispec_cv$spectral_cv, multispec_cv$n_spp)
cor.test(multispec_cv$spectral_cv, multispec_cv$n_spp)
spp_richness_lm <- lm(n_spp ~ spectral_cv, multispec_cv)
summary(spp_richness_lm)

# Test landscape level abundance
full_join(
    bind_rows(vascular_obs, bryo_obs) %>%
        arrange(site_id, plot_id, short_name, phylum) %>%
        relocate(site_id) %>%
        group_by(site_id, short_name) %>%
        summarise(p = n_distinct(short_name) / 18) %>%
        mutate(p_ln_p = p * log(p)) %>%
        group_by(site_id) %>%
        summarise(shan_div = -1 * sum(p_ln_p)),
    plot_multispec %>%
        group_by(site_id) %>%
        select(-ndvi) %>%
        pivot_longer(3:7, names_to = "band", values_to = "reflectance") %>%
        summarise(spectral_cv = sd(reflectance) / mean(reflectance))
)

# # Test clustering at sites
# set.seed(46)
# multispec_rasts %>%
#     split(., .$site) %>%
#     map(function(site_rasts) {
#         site_rast <- rast(site_rasts$file_name) %>%
#             crop(filter(site_bounds, site == unique(site_rasts$site))) %>%
#             setNames(site_rasts$band)
#         clust_results <- map(1:5, function(x) spatSample(site_rast, 1000)) %>%
#             pblapply(., FUN = function(y) {
#                 NbClust::NbClust(data = y, method = "kmeans", max.nc = 6)
#             }, cl = cl) %>%
#             map(function(x) DescTools::Mode(x$Best.nc[1, ])[1])
#         clust_results
#     })
# plot(site_rast[[1]])

# plot_multispec_sum <- plot_multispec %>%
#     group_by(site_id, plot_id) %>%
#     summarise(
#         blue = mean(blue),
#         green = mean(green),
#         # ndvi = mean(ndvi),
#         nir = mean(nir),
#         red = mean(red),
#         rededge = mean(rededge)
#     ) %>%
#     group_split() %>%
#     map(function(site_vals) {
#         site_vals <- as.data.frame(site_vals)
#         row.names(site_vals) <- site_vals$plot_id
#         pca <- prcomp(site_vals[, 3:ncol(site_vals)])
#         plot_grid(
#             fviz_pca_ind(pca),
#             fviz_pca_var(pca)
#         )
#     }) %>%
#     plot_grid(plotlist = ., nrow = 3)

# Test for correlation with sd of ndv
ndvi_sd <- plot_multispec %>%
    group_by(site_id, plot_id) %>%
    summarise(sd_ndvi = sd(ndvi)) %>%
    full_join(point_framing)
ggplot(ndvi_sd, aes(x = sd_ndvi, y = shan_div)) +
    geom_point() +
    geom_smooth(method = "lm")
cor(ndvi_sd$sd_ndvi, ndvi_sd$shan_div)
cor.test(ndvi_sd$sd_ndvi, ndvi_sd$shan_div)


ndvi_sd <- full_join(
    bind_rows(vascular_obs, bryo_obs) %>%
        group_by(site_id, plot_id) %>%
        summarise(rich = n_distinct(short_name)),
    plot_multispec %>% group_by(site_id, plot_id) %>%
        summarise(sd_ndvi = sd(ndvi))
)

ggplot(ndvi_sd, aes(x = sd_ndvi, y = rich)) +
    geom_point() +
    geom_smooth(method = "lm")
cor(ndvi_sd$sd_ndvi, ndvi_sd$rich)
cor.test(ndvi_sd$sd_ndvi, ndvi_sd$rich)

# Get overall site spectral variability
multispec_rasts %>%
    split(., .$site) %>%
    map(function(site_rasts) {
        site_rasts <- filter(site_rasts, band != "ndvi")
        site_rast <- rast(site_rasts$file_name) %>%
            crop(filter(site_bounds, site == unique(site_rasts$site))) %>%
            setNames(site_rasts$band)
        global_sd <- sd(as.vector(site_rast), na.rm = T)
        return(data.frame(
            site_id = unique(site_rasts$site),
            global_sd = global_sd
        ))
    }) %>%
    bind_rows() %>%
    full_join(bind_rows(vascular_obs, bryo_obs) %>%
        arrange(site_id, plot_id, short_name, phylum) %>%
        relocate(site_id) %>%
        group_by(site_id, short_name) %>%
        summarise(p = n_distinct(short_name) / 18) %>%
        mutate(p_ln_p = p * log(p)) %>%
        group_by(site_id) %>%
        summarise(shan_div = -1 * sum(p_ln_p)))


## Generate landscape level maps based on spectral variation (acrcross all bands)

# Aggregate rasters to CV (1 m rsolution)
spec_cv_rast <- multispec_rasts %>%
    filter(band != "ndvi") %>%
    arrange(site) %>%
    group_by(site) %>%
    group_split() %>%
    map(function(site_rasts) {
        site_rast <- rast(site_rasts$file_name) %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
        site_rast <- aggregate(site_rast,
            fact = c(round(1 / res(site_rast)), nlyr(site_rast)),
            fun = function(x) {
                sd(x, na.rm = T) / mean(x, na.rm = T)
            }
        ) %>%
            setNames("spectral_cv") %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
    })

# Override buggy cowplot::get_legend if needed
get_legend <- function(plot) {
    legends <- get_plot_component(plot, "guide-box", return_all = T)
    legends_new <- list()
    for (i in 1:length(legends)) {
        if (!("zeroGrob" %in% class(legends[[i]]))) legends_new <- c(legends_new, list(legends[[i]]))
    }
    if (length(legends_new) > 0) {
        do.call(cbind, legends_new)
    } else {
        return(legends[[1]])
    }
}
(test_plot <- ggplot() +
    geom_point(
        aes(
            x = 1:10, y = 1:10,
            colour = 1:10,
            size = 1:10
        ),
        shape = 21
    )) +
    guides(
        colour = guide_colourbar(position = "bottom"),
        size = guide_legend(position = "top")
    )
ggdraw() + draw_grob(get_legend(test_plot))
get_legend(ggplot())


# Generate water / snow masks (NDVI < 0)
water_masks <- multispec_rasts %>%
    filter(band == "ndvi") %>%
    arrange(site) %>%
    group_by(site) %>%
    group_split() %>%
    map(function(site_rasts) {
        site_rast <- rast(site_rasts$file_name) %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
        m <- c(-1, 0, NA,
               0, 1, 1)
        m <- matrix(m, ncol=3, byrow=TRUE)
        site_rast <- aggregate(site_rast,
            fact = c(round(1 / res(site_rast)), nlyr(site_rast)),
            fun = mean) %>%
            classify(m) %>%
            setNames("water") %>%
            crop(filter(site_bounds, site == unique(site_rasts$site)))
    })

# Generate plots
spectra_div_maps <- map(1:3, function(x) {
    # Subset data
    plots <- full_join(multispec_cv, select(plot_coords, Name, geom),
        by = c("plot_id" = "Name")
    ) %>%
        filter(grepl(x, site_id))
    cv_rast <- spec_cv_rast[[x]] %>% mask(water_masks[[x]])

    # Determine site colour
    site_colour <- site_colours %>%
        filter(grepl(x, site))

    # Visualise
    div_plot <- ggplot() +
        geom_spatraster(data = cv_rast) +
        # geom_spatraster_contour(data = cv_rast) +
        scale_fill_viridis_c(option = "magma", 
            limits = c(0.2, 1), na.value = "white",
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        geom_sf(
            data = st_as_sf(plots), aes(fill = spectral_cv, size = n_spp),
            shape = 21, colour = "dodgerblue", stroke = 2
        ) +
        scale_size_continuous(limits = c(0, 30), breaks = seq(0,30,10), range = c(0.1,5)) +
        labs(
            title = site_colour$site,
            fill = "Spectral Variation (CV)", 
            size = "Plot Species Richness (n)") +
        theme_map(18) +
        theme(panel.border = element_rect(colour = site_colour$value,
                                          fill = NA,
                                          linewidth = 2),
              plot.title = element_text(hjust = 0.5))
})

spect_div_plot <- plot_grid(spectra_div_maps[[1]] + theme(legend.position = "none"),
    spectra_div_maps[[2]] + theme(legend.position = "none"),
    spectra_div_maps[[3]] + theme(legend.position = "none"),
    ggdraw() + draw_grob(get_legend(spectra_div_maps[[1]]), 
                            x = 0.35 , y = 0),
    ncol = 2
    )

save_plot("figures/figure_3_spect_div_maps.png", spect_div_plot,
    base_asp = 1,
    ncol = 2, nrow = 2,
    bg = "white"
)

# Combined with scatter
save_plot("figures/figure_3_spect_div_maps_with_scat.png",
plot_grid(
    n_spp_vs_spect_cv_plot +
        labs(title = " ") +
        theme_cowplot(18) +
        theme(legend.position = "none", plot.margin = margin(7,14,7,7, "points")),
    spectra_div_maps[[1]] + theme(legend.position = "none"),
    spectra_div_maps[[2]] + theme(legend.position = "none"),
    ggplot() +
        annotate("segment", x = 0.275, xend = 0.725,
                            y = 0.35, yend = 0.35, linewidth = 1.5) +
        annotate("text", x = 0.5, y = 0.4, label = "150 m", size = 20/.pt) +
        annotate("segment", x = 0.45, xend = 0.45,
                            y = 0.50, yend = 0.65,
                arrow = arrow(), linewidth = 1.5) +
        annotate("text", x = 0.525, y = 0.575, label = "N", size = 32/.pt) +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
        theme_nothing(), 
    spectra_div_maps[[3]] + theme(legend.position = "none"),
    ggdraw() + draw_grob(get_legend(spectra_div_maps[[1]]), 
                            x = 0.2 , y = 0), 
    ncol = 3,
    nrow = 2,
    labels = c("a)", "b)", "", "", "", ""),
    label_size = 18),
    base_asp = 1,
    ncol = 3,
    nrow = 2,
    base_height = 5,
    bg = "white"
)

# Contour plots
est_rich_maps <- map(1:3, function(x) {
    # Subset data
    plots <- full_join(multispec_cv, select(plot_coords, Name, geom),
        by = c("plot_id" = "Name")
    ) %>%
        filter(grepl(x, site_id))
    cv_rast <- spec_cv_rast[[x]] %>% mask(water_masks[[x]])

    # Determine site colour
    site_colour <- site_colours %>%
        filter(grepl(x, site))

    # Predict number of species
    n_spp_rast <- predict(cv_rast, spp_richness_lm)

    # Classify raster 
    breaks <- matrix(c(0,5,1,
             5,10,2,
             10,15,3,
             15,20,4,
             20,25,5,
             25,30,6,
             30,35,7), ncol = 3, byrow = T)
    n_spp_rast <- classify(n_spp_rast, breaks) 
    levels(n_spp_rast)  <- data.frame(id = 1:7,
                                            n_species = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35"))
    
    # Visualise
    est_rich_plot <- ggplot() +
    geom_spatraster(data = n_spp_rast) +
    # scale_fill_viridis_c(limits = c(0, 35), breaks = seq(0,35,5)) +
    # geom_spatraster_contour_filled(data = (predict(cv_rast, spp_richness_lm)),
    #                         breaks = seq(0,35,5)) +
    scale_fill_viridis_d(na.value = "white", 
    limits = factor(c(levels(n_spp_rast)[[1]]$n_species[-1], NA)),
    labels = c(
     "5-10",
     "10-15",
     "15-20",
     "20-25", 
     "25-30",
     "30-35",
     "Water/Snow")
    ) +
    #     # geom_sf(
        #     data = st_as_sf(plots), aes(size = n_spp),
        #     fill = NA,
        #     shape = 21, colour = "magenta", stroke = 2
        # ) +
        scale_size_continuous(limits = c(0, 35)) +
        labs(
            #title = site_colour$site,
        fill = "Estimated Species\nRichness (n)" #,
        # size = "Species Richness (n)"
        ) +
    theme_map(18)  +
        theme(
            panel.border = element_rect(colour = site_colour$value,
                                          fill = NA,
                                          linewidth = 2),
            plot.margin = unit(c(0,0,0.1,0), "inch"),
            legend.key = element_rect(colour = "black")
        #       plot.title = element_text(hjust = 0.5))
        )
})

# Generate RGB plots
rgb_plots <- map(1:3, function(x){
    ggplot() + 
    geom_spatraster_rgb(data = crop(rgb_rasts[[x]], site_bounds[x,])) + 
    labs(title = paste0("Site ", x)) +
    theme_map(18) +
    theme(plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = site_colours$value[x],
                                          fill = NA,
                                          linewidth = 2),
        plot.margin = unit(c(0.1,0,0,0), "inch"))
})

est_rich_plot <- plot_grid(
    plot_grid(
        rgb_plots[[1]],
        est_rich_maps[[1]] + theme(legend.position = "none"),
        align = "hv",
        ncol = 1),
    plot_grid(
        rgb_plots[[2]],
    est_rich_maps[[2]] + theme(legend.position = "none"),
    align = "hv",
    ncol = 1),
    plot_grid(
        rgb_plots[[3]],
        est_rich_maps[[3]] + theme(legend.position = "none"),
        align = "hv",
        ncol = 1),
    plot_grid(
        ggplot() +
        annotate("text", x = 0.075, y = 0.25, label = "150 m", angle = 90,
        size = 15/.pt) +
        annotate("segment", x= 0, xend = 0, y = 0.025, yend = 0.885 - 0.43, 
                linewidth = 1) +
        annotate("text", x = 0.5, y = 0.75, label = "True colour imagery",
                 size = 18/.pt) +
        annotate("segment", x = 0.4, xend = 0.4, y = 0.15, yend = 0.3,
                arrow = arrow(), linewidth = 1.5) +
        annotate("text", x = 0.55, y = 0.225, label = "N", size = 32/.pt) +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
        theme_map(18) +
        theme(plot.margin = unit(c(0,0,0,0), "inch")),
    ggdraw() + draw_grob(get_legend(est_rich_maps[[3]]), 
                           x = 0.2 , y = 0
                           ),
            nrow = 2, ncol = 1),
    rel_widths = c(1,1,1,0.8),
    align = "hv",
    ncol = 4)
    

# Save
save_plot("figures/figure_4_est_spp_richness_maps.png",
    est_rich_plot, 
    base_asp = 0.9,
    ncol = 4,
    nrow = 2,
    bg = "white"
)
