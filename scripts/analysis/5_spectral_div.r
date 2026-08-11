# Extract spectral diversity from drone data for plots
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(terra)
library(tidyterra)
library(sf)
library(NbClust)
library(pbapply)
library(cowplot)
library(factoextra)
library(ggnewscale)
library(broom)
library(gt)

# Set up parallel envrionment
cl <- parallel::makeCluster(5)

# Set site colours (not used anymore)
site_colours <- data.frame(
    site = c("Site 1", "Site 2", "Site 3"),
    value = c(
        "#F2D64B",
        "#5F1773",
        "#F2507B"
    )
)

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
}) %>%
    bind_rows()

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
}) %>%
    bind_rows()

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
    mutate(phylum = "Vascular plants") %>%
    # Remove uncertain species-level ids
    filter(species != "species")

# Load presence / absence moss observations
bryo_obs <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb_updated_Aug_2026.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus)) %>%
    filter(!grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    arrange(site_id, plot_id, short_name) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(phylum = "Bryophytes")

# Combine dataframes
combined_obs <- bind_rows(vascular_obs, bryo_obs) %>%
    arrange(site_id, plot_id, short_name, phylum) %>%
    distinct() %>%
    relocate(site_id) %>%
    group_by(plot_id) %>%
    summarise(
        n_spp = sum(phylum %in% c("Vascular plants", "Bryophytes")),
        n_spp_vasc = sum(phylum == "Vascular plants"),
        n_spp_bryo = sum(phylum == "Bryophytes")
    ) %>%
    full_join(plot_rgb)

# Point framing
point_framing <- read_csv("data/botanical_obs/pointframing.csv") %>%
    select(-observer, -scribe, -date) %>%
    group_by(plot_id, full_name) %>%
    summarise(sum_hits = n()) %>%
    split(., .$plot_id) %>%
    map(function(x) {
        x <- mutate(x, p = sum_hits / sum(x$sum_hits), p_ln_p = p * log(p)) %>%
            summarise(shan_div = -1 * sum(p_ln_p))
    }) %>%
    bind_rows()

# Combine
combined_obs <- full_join(combined_obs, point_framing)

## RGB based anlaysis

# # Correlate the the species richness with the spectral diversity
# ggplot() +
#     geom_point(data = combined_obs, aes(x = n_spp, spectral_cv)) +
#     theme_classic()

# # Correlate the the cover richness with the spectral diversity
# ggplot(data = combined_obs, mapping = aes(x = shan_div, spectral_cv)) +
#     geom_point() +
#     geom_smooth(method = "lm") +
#     theme_classic()

## Multispectral analysis
multispec_cv <- plot_multispec %>%
    select(-ndvi) %>%
    pivot_longer(3:7, names_to = "band", values_to = "reflectance") %>%
    group_by(site_id, plot_id) %>%
    summarise(spectral_cv = sd(reflectance) / mean(reflectance)) %>%
    full_join(
        .,
        select(combined_obs, plot_id, n_spp, n_spp_vasc, n_spp_bryo, shan_div)
    )

# Save data
write_csv(multispec_cv, "data/combined_obs_multispec_cv.csv")

# # Test for correlation spectral cv and shannon diversity (point frame)
# ggplot(multispec_cv, aes(x = spectral_cv, y = shan_div)) +
#     geom_point() +
#     geom_smooth(method = "lm") +
#     theme_classic()
# hist(multispec_cv$shan_div)
# hist(multispec_cv$spectral_cv)
# cor(multispec_cv$spectral_cv, multispec_cv$shan_div)
# cor.test(multispec_cv$spectral_cv, multispec_cv$shan_div)
# summary(lm(shan_div ~ spectral_cv, multispec_cv))

# Test for correlation spectral cv and n_spp (abundance)
ggplot(multispec_cv, aes(x = spectral_cv, y = n_spp)) +
    geom_smooth(
        method = "lm",
        colour = "#962c7e",
        fill = "#cdcdcd11",
        linewidth = 2
    ) +
    geom_point(aes(shape = site_id, color = spectral_cv), size = 5) +
    scale_color_viridis_c(
        option = "magma",
        limits = c(0.2, 1),
        na.value = "white",
        guide = "none"
    ) +
    theme_cowplot()

# Assess for normality
hist(multispec_cv$n_spp)
hist(multispec_cv$n_spp_vasc)
hist(multispec_cv$n_spp_bryo)
hist(multispec_cv$spectral_cv)

## Relationship between plot-level spectral cv and species richness

# Test Pearson correlation
cor(multispec_cv$spectral_cv, multispec_cv$n_spp)
cor.test(multispec_cv$spectral_cv, multispec_cv$n_spp)

#         Pearson's product-moment correlation

# data:  multispec_cv$spectral_cv and multispec_cv$n_spp
# t = 5.9817, df = 52, p-value = 2.058e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4471244 0.7739084
# sample estimates:
#       cor
# 0.6384485

spp_richness_lm <- lm(n_spp ~ spectral_cv, multispec_cv)
summary(spp_richness_lm)
# Call:
# lm(formula = n_spp ~ spectral_cv, data = multispec_cv)
# Residuals:
#     Min      1Q  Median      3Q     Max
# -11.812  -2.862   0.332   3.193   8.125

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    4.372      2.288   1.911   0.0616 .
# spectral_cv   25.975      4.342   5.982 2.06e-07 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 4.754 on 52 degrees of freedom
# Multiple R-squared:  0.4076,    Adjusted R-squared:  0.3962
# F-statistic: 35.78 on 1 and 52 DF,  p-value: 2.058e-07

# Save linear model and export table
save(spp_richness_lm, file = "data/analysis/lm_nspp_multispectral_cv.Rda")
spp_richness_lm %>%
    summary() %>%
    tidy() %>%
    # rename() %>%
    mutate(
        across(2:4, function(x) round(x, 2)),
        across(5, function(x) round(x, 3))
    ) %>%
    gt() %>%
    gtsave("tables/6_lm_nspp_multispectral.html")

## Rrepeat analysis with vascular / bryophytes only
ggplot(multispec_cv, aes(x = spectral_cv)) +
    geom_point(aes(y = n_spp), colour = "blue") +
    geom_smooth(aes(y = n_spp), method = "lm", colour = "blue") +
    geom_point(aes(y = n_spp_vasc), colour = "red") +
    geom_smooth(aes(y = n_spp_vasc), method = "lm", colour = "red") +
    geom_point(aes(y = n_spp_bryo), colour = "green") +
    geom_smooth(aes(y = n_spp_bryo), method = "lm", colour = "green") +
    theme_cowplot()

# Export models
spp_richness_vasc_lm <- lm(n_spp_vasc ~ spectral_cv, multispec_cv)
summary(spp_richness_vasc_lm)

# Call:lm(formula = n_spp_vasc ~ spectral_cv, data = multispec_cv)
# Residuals:
#     Min      1Q  Median      3Q     Max
# -3.6915 -1.8285 -0.5537  2.0628  8.3421

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    1.926      1.196   1.611    0.113
# spectral_cv   11.908      2.269   5.248 2.88e-06 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 2.484 on 52 degrees of freedom
# Multiple R-squared:  0.3462,    Adjusted R-squared:  0.3336
# F-statistic: 27.54 on 1 and 52 DF,  p-value: 2.88e-06
save(
    spp_richness_vasc_lm,
    file = "data/analysis/lm_nspp_vasc_multispectral_cv.Rda"
)
spp_richness_vasc_lm %>%
    summary() %>%
    tidy() %>%
    # rename() %>%
    mutate(
        across(2:4, function(x) round(x, 2)),
        across(5, function(x) round(x, 3))
    ) %>%
    gt() %>%
    gtsave("tables/7_lm_nspp_vasc_multispectral.html")
spp_richness_bryo_lm <- lm(n_spp_bryo ~ spectral_cv, multispec_cv)
save(
    spp_richness_bryo_lm,
    file = "data/analysis/lm_nspp_bryo_multispectral_cv.Rda"
)
summary(spp_richness_bryo_lm)
# Residuals:
#    Min     1Q Median     3Q    Max
# -8.890 -2.658 -0.262  3.022  8.954

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    2.446      1.926   1.270 0.209622
# spectral_cv   14.067      3.655   3.849 0.000326 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 4.001 on 52 degrees of freedom
# Multiple R-squared:  0.2217,    Adjusted R-squared:  0.2068
# F-statistic: 14.82 on 1 and 52 DF,  p-value: 0.0003262
spp_richness_bryo_lm %>%
    summary() %>%
    tidy() %>%
    # rename() %>%
    mutate(
        across(2:4, function(x) round(x, 2)),
        across(5, function(x) round(x, 3))
    ) %>%
    gt() %>%
    gtsave("tables/8_lm_nspp_bryo_multispectral.html")

## OLD code

# # Test the reverse relationship to assess the added contributions of bryophytes
# # to the relationship
# spp_richness_sep_rev_lm <- lm(
#     spectral_cv ~ n_spp_bryo + n_spp_vasc,
#     multispec_cv
# )
# spp_richness_vasc_rev_lm <- lm(
#     spectral_cv ~ n_spp_vasc,
#     multispec_cv
# )
# summary(spp_richness_vasc_rev_lm)
# summary(spp_richness_sep_rev_lm)
# anova(spp_richness_vasc_rev_lm, spp_richness_sep_rev_lm)
# summary(spp_richness_sep_rev_lm)$adj.r.squared
# summary(spp_richness_vasc_rev_lm)$adj.r.squared
# summary(spp_richness_sep_rev_lm)$adj.r.squared -
#     summary(spp_richness_vasc_rev_lm)$adj.r.squared

# # Test landscape level abundance
# full_join(
#     bind_rows(vascular_obs, bryo_obs) %>%
#         arrange(site_id, plot_id, short_name, phylum) %>%
#         relocate(site_id) %>%
#         group_by(site_id, short_name) %>%
#         summarise(p = n_distinct(short_name) / 18) %>%
#         mutate(p_ln_p = p * log(p)) %>%
#         group_by(site_id) %>%
#         summarise(shan_div = -1 * sum(p_ln_p)),
#     plot_multispec %>%
#         group_by(site_id) %>%
#         select(-ndvi) %>%
#         pivot_longer(3:7, names_to = "band", values_to = "reflectance") %>%
#         summarise(spectral_cv = sd(reflectance) / mean(reflectance))
# )

# # # Test clustering at sites
# # set.seed(46)
# # multispec_rasts %>%
# #     split(., .$site) %>%
# #     map(function(site_rasts) {
# #         site_rast <- rast(site_rasts$file_name) %>%
# #             crop(filter(site_bounds, site == unique(site_rasts$site))) %>%
# #             setNames(site_rasts$band)
# #         clust_results <- map(1:5, function(x) spatSample(site_rast, 1000)) %>%
# #             pblapply(., FUN = function(y) {
# #                 NbClust::NbClust(data = y, method = "kmeans", max.nc = 6)
# #             }, cl = cl) %>%
# #             map(function(x) DescTools::Mode(x$Best.nc[1, ])[1])
# #         clust_results
# #     })
# # plot(site_rast[[1]])

# # plot_multispec_sum <- plot_multispec %>%
# #     group_by(site_id, plot_id) %>%
# #     summarise(
# #         blue = mean(blue),
# #         green = mean(green),
# #         # ndvi = mean(ndvi),
# #         nir = mean(nir),
# #         red = mean(red),
# #         rededge = mean(rededge)
# #     ) %>%
# #     group_split() %>%
# #     map(function(site_vals) {
# #         site_vals <- as.data.frame(site_vals)
# #         row.names(site_vals) <- site_vals$plot_id
# #         pca <- prcomp(site_vals[, 3:ncol(site_vals)])
# #         plot_grid(
# #             fviz_pca_ind(pca),
# #             fviz_pca_var(pca)
# #         )
# #     }) %>%
# #     plot_grid(plotlist = ., nrow = 3)

# # Test for correlation with sd of ndv
# ndvi_sd <- plot_multispec %>%
#     group_by(site_id, plot_id) %>%
#     summarise(sd_ndvi = sd(ndvi)) %>%
#     full_join(point_framing)
# ggplot(ndvi_sd, aes(x = sd_ndvi, y = shan_div)) +
#     geom_point() +
#     geom_smooth(method = "lm")
# cor(ndvi_sd$sd_ndvi, ndvi_sd$shan_div)
# cor.test(ndvi_sd$sd_ndvi, ndvi_sd$shan_div)

# ndvi_sd <- full_join(
#     bind_rows(vascular_obs, bryo_obs) %>%
#         group_by(site_id, plot_id) %>%
#         summarise(rich = n_distinct(short_name)),
#     plot_multispec %>%
#         group_by(site_id, plot_id) %>%
#         summarise(sd_ndvi = sd(ndvi))
# )

# ggplot(ndvi_sd, aes(x = sd_ndvi, y = rich)) +
#     geom_point() +
#     geom_smooth(method = "lm")
# cor(ndvi_sd$sd_ndvi, ndvi_sd$rich)
# cor.test(ndvi_sd$sd_ndvi, ndvi_sd$rich)

# # Get overall site spectral variability
# multispec_rasts %>%
#     split(., .$site) %>%
#     map(function(site_rasts) {
#         site_rasts <- filter(site_rasts, band != "ndvi")
#         site_rast <- rast(site_rasts$file_name) %>%
#             crop(filter(site_bounds, site == unique(site_rasts$site))) %>%
#             setNames(site_rasts$band)
#         global_sd <- sd(as.vector(site_rast), na.rm = T)
#         return(data.frame(
#             site_id = unique(site_rasts$site),
#             global_sd = global_sd
#         ))
#     }) %>%
#     bind_rows() %>%
#     full_join(
#         bind_rows(vascular_obs, bryo_obs) %>%
#             arrange(site_id, plot_id, short_name, phylum) %>%
#             relocate(site_id) %>%
#             group_by(site_id, short_name) %>%
#             summarise(p = n_distinct(short_name) / 18) %>%
#             mutate(p_ln_p = p * log(p)) %>%
#             group_by(site_id) %>%
#             summarise(shan_div = -1 * sum(p_ln_p))
#     )

# End of file
