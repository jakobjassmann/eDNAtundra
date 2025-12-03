# Figure 1 - Vascular plant richness detection
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

# Load data
obs_combined <- read_csv("data/analysis/combined_observations.csv")

# Plot level detection of species
plot_level_detection <- obs_combined %>%
    # Remove all genus level only ids
    filter(species != "species" | genus == "Salix") %>%
    # eDNA detection of Salix at the genus level
    group_by(plot_id, genus) %>%
    group_split() %>%
    map(function(x) {
        # If it's not a Salix species skip
        if (unique(x$genus) != "Salix") {
            return(x)
        }
        # If it's just one observaiton of a salix species skip
        if (nrow(x) == 1) {
            return(x)
        }
        # If it's just botanical detections keep only species ids
        if (sum(grepl("eDNA", x$source)) == 0) {
            x <- filter(x, species != "species")
            return(x)
        }
        # In all other cases, assign a random species of the sub set
        # to the eDNA detection
        x[x$source == "eDNA", ]$species <- filter(x, source != "eDNA") %>%
            pull(species) %>%
            sample(1)
        x <- mutate(x, short_name = paste(genus, species))
        # Filter any remaining entries wihtout species id
        x <- filter(x, species != "species")
        return(x)
    }) %>%
    bind_rows() %>%
    group_by(plot_id, species) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(plot_id, short_name, source) %>%
    mutate(source = factor(source, levels = rev(c("botanical eDNA", "botanical", "eDNA")))) %>%
    arrange(plot_id, source)

# Plot richness
plot_richness <- plot_level_detection %>%
    group_by(plot_id, source) %>%
    summarise(n = n()) %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id))

# Site level detection of species
site_level_detection <- obs_combined %>%
    # Treating salix as one species removing all other species-level ids
    filter(species != "species" | genus == "Salix") %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    select(-plot_id) %>%
    # eDNA detection of Salix at the genus level
    group_by(site_id, genus) %>%
    group_split() %>%
    map(function(x) {
        # If it's not a Salix species skip
        if (unique(x$genus) != "Salix") {
            return(x)
        }
        # If it's just one observaiton of a salix species skip
        if (nrow(x) == 1) {
            return(x)
        }
        # If it's just botanical detections keep only species ids
        if (sum(grepl("eDNA", x$source)) == 0) {
            x <- filter(x, species != "species")
            return(x)
        }
        # In all other cases, assign a random species of the sub set
        # to the eDNA detection
        x[x$source == "eDNA", ]$species <- filter(x, source != "eDNA") %>%
            pull(species) %>%
            sample(1)
        x <- mutate(x, short_name = paste(genus, species))
        # Filter any remaining entries wihtout species id
        x <- filter(x, species != "species")
        return(x)
    }) %>%
    bind_rows() %>%
    group_by(site_id, species) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(site_id, short_name, source) %>%
    mutate(source = factor(source, levels = rev(c("botanical eDNA", "botanical", "eDNA")))) %>%
    arrange(site_id, short_name, source)

# Quick qc:
print(site_level_detection, n = 110)

# Summarise source groupings and visualise
site_richness <- site_level_detection %>%
    group_by(site_id, source) %>%
    summarise(n = n()) %>%
    mutate(site_id = gsub("S", "Site ", site_id))

# pal <- c("#60BFBF", "#BF638E", "#A6936F")
# pal <- c("#BF612A", "#020202", "#1F1A40")
# pal <- c("#8BA37CCC", "#275FA1CC", "#614C4CCC")
# pal <- c("#A12B93FF", "#4BA12BFF", "#332732FF")
pal <- c("#F28F16", "#A62121", "#223F59")
# pal <- c("#8BA37C", "#0490F2", "#614C4C")
# pal <- c("#8BA37C", "#BF638E", "#0490F2")
# pal_fill <- pal
# pal_colour <- pal
# pal_fill <- c("white", "#BF638E", "#8BA37C")
# pal_colour <- c("#858585", "#BF638E", "#8BA37C")
# pal_fill <-  c("#FFFFFF", "#4BA12B", "#332732")
# pal_colour <-  c("#A12B93", "#4BA12B", "#332732")

# Plot level plot
(richness_detection <- ggplot(plot_richness) +
    geom_col(aes(
        x = plot_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Vascular plant species richness (n)",
        x = "Plot ID"
    ) +
    scale_colour_manual(
        values = pal, name = "Observation",
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observation",
        labels = c("eDNA only", "botanical only",  "botanical and eDNA")
    ) +
    scale_y_continuous(limits = c(0, 20)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)))

# Site level plot
(site_richness_plot <- ggplot(site_richness) +
    geom_col(aes(
        x = site_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Vascular plant species richness (n)"
    ) +
    scale_colour_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only",  "botanical and eDNA")
    )+
    scale_y_continuous(limits = c(0, 40)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank()))


save_plot("figures/figure_2a.png",
    site_richness_plot,
    #plot_grid(richness_detection, site_richness_plot),
    #base_asp = 2 * 1.6,
    bg = "white"
)
