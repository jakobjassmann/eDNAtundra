# Figure 2 - Vascular plant richness detection
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

## Panel a) Bryophyte vs. Moss Species Richness

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
relocate(site_id)

# Caluclate number of species per site
(combined_obs_nsp_site <- combined_obs %>%
select(-plot_id) %>%
group_by(site_id, phylum) %>%
summarise(nsp = n_distinct(short_name)) %>%
bind_rows() %>%
arrange(site_id, phylum))

# Visualise at the site level
(bryo_vs_vasc_plot <- ggplot(combined_obs_nsp_site) +
    geom_col(aes(x= site_id, y = nsp , fill = phylum)) +
    labs(x = "",
         y = "Number of Species (n)",
         title = "Landscape Richness") +
    scale_fill_manual(values = c("#EB7075", "#a62121"),
                      name = "Botanical observation",
                      labels = c("Bryophytes", "Vascular Plants")) +
    scale_y_continuous(limits = c(0, 120),
                       breaks = seq(0,120,20)) +
    theme_cowplot() +
    theme(legend.margin = margin(0.25, 0.35, 0.25, 0.35, unit = "cm"),
          legend.position = "right",
          legend.box.background = element_rect(colour = "black"),
          legend.justification.right = c(0.5, 0.25),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank()))

# Add icons
# (bryo_vs_vasc_plot_icons <- ggdraw() +
#     draw_plot(bryo_vs_vasc_plot) +
#     draw_image("figures/moss_bw.png", scale = 0.1, x= 0.275, y = 0.175) +
#     draw_image("figures/vascular_plant_bw.png", scale = 0.15, x = 0.375, y = 0.190) )

# save_plot("figures/fig_2b.png",
#     bryo_vs_vasc_plot_icons,
#     bg = "white")

## Vascular Plant DNA detection rates

# Load data
obs_combined <- read_csv("data/analysis/combined_observations.csv")

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

# Set colours
pal <- c("#F28F16", "#A62121", "#223F59")

# Site level plot
(site_richness_plot <- ggplot(site_richness) +
    geom_col(aes(
        x = site_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Number of Species (n)",
        title = "Vascular Plant Detection"
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
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "right", 
        #legend.position.inside = c(0.06, 0.9),
        legend.justification.right = c(0.5, 0.25),
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(0.25,0.50,0.25,0.35, unit = "cm"),
        plot.title = element_text(hjust = 0.5),))

# # Add icons
# (site_richness_plot_icons <- ggdraw() + 
#     draw_plot(site_richness_plot) +
#     # Vascular plants
#     draw_image("figures/vascular_plant_bw.png", scale = 0.15, x = 0.375, y = 0.200) +
#     # DNA 
#     draw_image("figures/dna.png", scale = 0.13, x = 0.275, y = 0.195))

# save_plot("figures/figure_2b.png",
#     site_richness_plot_icons,
#     #plot_grid(richness_detection, site_richness_plot),
#     #base_asp = 2 * 1.6,
#     bg = "white"
# )

# Assemble figure 2 (one colum)
(panel_plot <- (plot_grid(
           bryo_vs_vasc_plot,
           site_richness_plot,
           #bryo_site_richness_plot_icons,
           ncol = 1,
           nrow = 2,
           labels = c("a)", "b)"), # "c)"),
           align = "hv")))

# Add icons
(panel_plot_icons <- ggdraw() + 
    draw_plot(panel_plot) +
    # Bryophytes
    draw_image("figures/moss_bw.png", scale = 0.12, x= 0.25, y = 0.290) +
    # Vascular plants I
    draw_image("figures/vascular_plant_bw.png", scale = 0.1, x = 0.375, y = 0.305) +
    # Vascular plants II
    draw_image("figures/vascular_plant_bw.png", scale = 0.1, x = 0.375, y = -0.18) +
    # DNA 
    draw_image("figures/dna.png", scale = 0.1, x = 0.25, y = -0.18))
    

save_plot("figures/figure_2.png", panel_plot_icons, ncol = 1, nrow = 2, bg = "white")

## Bryophyte detection

# Load botanical observations
botanical_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv") %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # convert al sp. to species
    mutate(species = (gsub("sp\\.", "species", species))) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) & !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species))

# Load eDNA data
eDNA_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_bryo_taxa_in_sample.csv") %>%
    mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    mutate(site_id = substr(plot_id, 1,2))

# Check theoretical species level id
combined_ref_species <- full_join(
    eDNA_obs %>%
        distinct(genus, species, short_name) %>%
        mutate(eDNA = TRUE),
    botanical_obs %>%
        distinct(genus, species, short_name) %>%
        mutate(botanical = TRUE)
) %>%
    mutate(source = case_when(
        botanical & is.na(eDNA) ~ "botanical",
        eDNA & is.na(botanical) ~ "eDNA",
        botanical & eDNA ~ "botanical & eDNA",
        TRUE ~ NA
    )) %>%
    select(-eDNA, -botanical) %>%
    distinct() %>%
    print(n = 100)

# Combine reference list including:
#   botanical observations,
#   eDNA observations if the species was observed (anywhere) in the botanical obs
combined_ref <- # Species in botanical or both
    combined_ref_species %>% filter(source %in% c("botanical & eDNA", "botanical")) %>%
    mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct() %>%
    print(n = 100)

# Filter eDNA obs to subsets containing only taxa observed in botanical obs
eDNA_obs_sub <- eDNA_obs %>%
    filter(short_name %in% combined_ref$short_name |
        genus %in% combined_ref$genus) %>%
    # Remove extraction controls
    filter(!is.na(subsample_id)) %>%
    select(site_id, plot_id, genus, species) %>%
    mutate(source = "eDNA") %>%
    # If detection to species level was not possible set to "species"
    mutate(species = case_when(!(species %in% combined_ref$species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct()

# combine botanical obs and eDNA detection
obs_combined <- bind_rows(
    # Botanical observations species level id
    mutate(botanical_obs, source = "botanical"),
    # Botanical observations genus level id
    botanical_obs %>% 
        distinct(site_id, plot_id, genus) %>%
        mutate(species = "species") %>%
        mutate(short_name = paste0(genus, " ", species),
                source = "botanical"),
    # eDNA detection
    eDNA_obs_sub) %>%
    arrange(site_id, plot_id, short_name, source) %>%
    distinct() %>%
    mutate(source = factor(source, levels = c("botanical", "eDNA", "botanical & eDNA"))) %>%
    mutate(site_id = case_when(site_id == "S1" ~ "Site 1",
    site_id == "S2" ~ "Site 2",
    site_id == "S3" ~ "Site 3",
    TRUE ~ NA))

# Site level detection of species
site_level_detection <- obs_combined %>%
    # Remove all ids not to the species level
    filter(species != "species") %>%
    # Keep only speces detected in the botanical observations
    filter(short_name %in% unique(botanical_obs$short_name)) %>%
    select(-plot_id) %>%
    group_by(site_id, short_name) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(site_id, short_name, source) %>%
    mutate(source = factor(source, levels = rev(c("botanical eDNA", "botanical", "eDNA")))) %>%
    arrange(site_id, short_name, source)

site_level_detection %>% filter(source == "eDNA")

# Summarise source groupings and visualise
(bryo_site_richness <- site_level_detection %>%
    group_by(site_id, source) %>%
    summarise(n = n()))

# Set colors 
pal <- c("#F28F16", "#EB7075", "#223F59") # "#A62121"

(bryo_site_richness_plot <- ggplot(bryo_site_richness) +
    geom_col(aes(
        x = site_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Number of Species (n)",
        title = "Bryophyte Detection"
    ) +
    scale_colour_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only",  "botanical and eDNA")
    )+
    scale_y_continuous(limits = c(0, 60)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "right", 
        #legend.position.inside = c(0.06, 0.9),
        #legend.justification.right = c(0.5, 0.105),
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(0.25,0.75,0.25,0.25, unit = "cm"),
        plot.title = element_text(hjust = 0.5)))

# Add icons
(bryo_site_richness_plot_icons <- ggdraw() + 
    draw_plot(bryo_site_richness_plot) +
    # Vascular plants
    draw_image("figures/moss_bw.png", scale = 0.1, x = 0.31, y = 0.225) +
    # DNA 
    draw_image("figures/dna.png", scale = 0.13, x = 0.21, y = 0.230))

save_plot("figures/figure_S1.png",
    bryo_site_richness_plot_icons,
    #plot_grid(richness_detection, site_richness_plot),
    #base_asp = 2 * 1.6,
    bg = "white"
)

