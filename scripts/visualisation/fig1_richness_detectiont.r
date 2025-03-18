# Figure 1 - Vascular plant richness detection
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)

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
    map(function(x){
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>% bind_rows() %>%
    distinct(plot_id, short_name, source) %>%
    # Remove eDNA only
    filter(source != "eDNA")

# Plot richness
plot_richness <- plot_level_detection %>%
    group_by(plot_id, source) %>%
    summarise(n = n()) %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
mutate(site_id = gsub("S", "Site ", site_id))

# SSite level detection of species
site_level_detection <- obs_combined %>%
# Treating salix as one species removing all other species-level ids
    filter(species != "species" | genus == "Salix") %>%
    mutate(site_id = substr(plot_id, 1,2)) %>%
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
    map(function(x){
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>% bind_rows() %>%
    distinct(site_id, short_name, source) %>%
    arrange(site_id, short_name, source)

# Quick qc:
print(site_level_detection, n = 110)
# Summarise source groupings and visualise
site_richness <- site_level_detection %>%
group_by(site_id, source) %>%
summarise(n = n()) %>%
filter(source != "eDNA") %>%
mutate(site_id = gsub("S", "Site ", site_id))


# Plot level plot
(richness_detection <- ggplot(plot_richness) +
  geom_col(aes(
    x = plot_id, y = n,
    fill = paste(site_id, source),
    colour = site_id)) +
  labs(y = "Vascular plant species richness (n)", fill = "Site",
       x = "Plots") +
scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  scale_fill_manual(values = c("white", "#F2D64B", "white", "#5F1773", "white", "#F2507B")) +
scale_y_continuous(limits = c(0,15)) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")) 
# Site level plot
(site_richness_plot <- ggplot(site_richness) +
  geom_col(aes(
    x = site_id, y = n,
    fill = paste(site_id, source),
    colour = site_id)) +
  labs(y = "Vascular plant species richness (n)", fill = "Site",
       x = "Sites") +
scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  scale_fill_manual(values = c("white", "#F2D64B", "white", "#5F1773", "white", "#F2507B")) +
scale_y_continuous(limits = c(0,40)) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")) 
save_plot("figures/figure_1_vasc_sp_detection.png", 
plot_grid(richness_detection, site_richness_plot),
          base_asp = 2* 1.6,
          bg = "white") 
