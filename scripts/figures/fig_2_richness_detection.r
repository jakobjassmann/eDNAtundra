# Figure 2 - Vascular plant richness detection
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025, updated Julyl 2026

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

## a) Vascular Plant Detection

# Load data
obs_combined <- read_csv("data/analysis/detection_species_genus.csv")

# Prepare data frame with site-level detection for plotting
site_richness <- obs_combined %>%
    # add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    relocate(site_id, short_name) %>%
    select(-plot_id, -genus, -species) %>%
    # aggregate detections by site
    group_by(site_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    # Aggregate into detection classes
    group_by(site_id) %>%
    summarise(
        `Botanical observations only` = sum(botanical & !eDNA & !eDNA_harm),
        `Botanical and eDNA (Harmonised)` = sum(botanical & !eDNA & eDNA_harm),
        `Botanical and eDNA` = sum(botanical & eDNA),
        `eDNA only` = sum(!botanical & eDNA)
    ) %>%
    ungroup() %>%
    pivot_longer(
        cols = `Botanical observations only`:`eDNA only`,
        names_to = "source",
        values_to = "n"
    ) %>%
    mutate(
        source = factor(
            source,
            levels = c(
                "eDNA only",
                "Botanical observations only",
                "Botanical and eDNA (Harmonised)",
                "Botanical and eDNA"
            )
        )
    )

# Plot site-level detection
(site_richness_plot <- ggplot(site_richness) +
    geom_col(
        aes(
            x = site_id,
            y = n,
            fill = source
            # colour = source
        )
        # linewidth = 1.5,
        # lineend = "square"
    ) +
    labs(
        y = "Number of Species (n)"
    ) +
    scale_fill_manual(
        values = c("#F28F16", "#A62121", "#016E63", "#223F59"),
        name = "Vascular plants observed with"
    ) +
    # scale_colour_manual(
    #     values = c("#F28F16", "#A62121", "#223F59", "#223F59"),
    #     name = "Vascular plants observed with"
    # ) +
    scale_y_continuous(limits = c(0, 50)) +
    guides(
        fill = guide_legend(nrow = 4)
    ) +
    theme_cowplot() +
    theme(axis.title.x = element_blank()) +
    theme(
        legend.position = "right",
        #legend.position.inside = c(0.06, 0.9),
        legend.justification.right = c(0, 0.25),
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(0.25, 0.50, 0.25, 0.35, unit = "cm"),
        plot.title = element_text(hjust = 0.5),
    ))

## b) Bryophyte vs. Vascular Plant Species Richness

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
bryo_obs <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # convert al sp. to species
    mutate(species = (gsub("sp\\.", "species", species))) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    arrange(site_id, plot_id, short_name) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(phylum = "Bryophytes")

# Combine dataframes
combined_obs <- bind_rows(vascular_obs, bryo_obs) %>%
    arrange(site_id, plot_id, short_name, phylum) %>%
    relocate(site_id) %>%
    # Keep only id's to species level
    filter(species != "species")

# Caluclate number of species per site
(combined_obs_nsp_site <- combined_obs %>%
    select(-plot_id) %>%
    group_by(site_id, phylum) %>%
    summarise(nsp = n_distinct(short_name)) %>%
    bind_rows() %>%
    arrange(site_id, phylum))

# Visualise at the site level
(bryo_vs_vasc_plot <- ggplot(combined_obs_nsp_site) +
    geom_col(aes(x = site_id, y = nsp, fill = phylum)) +
    labs(x = "", y = "Number of Species (n)") +
    scale_fill_manual(
        values = c("#EB7075", "#a62121"),
        name = "Botanical observations of",
        labels = c("Bryophytes", "Vascular plants")
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    theme_cowplot() +
    theme(
        legend.margin = margin(0.25, 0.35, 0.25, 0.35, unit = "cm"),
        legend.position = "right",
        legend.box.background = element_rect(colour = "black"),
        legend.justification.right = c(0, 0.25),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank()
    ))

## Assemble figure and annotate with icons

# Assemble multipanel plot
(panel_plot <- (plot_grid(
    site_richness_plot,
    bryo_vs_vasc_plot,
    ncol = 1,
    nrow = 2,
    labels = c("a)", "b)"), # "c)"),
    align = "hv"
)))

# Add icons
(panel_plot_icons <- ggdraw() +
    draw_plot(panel_plot) +
    # DNA
    draw_image("figures/dna.png", scale = 0.1, x = 0.14, y = 0.355) +
    # Vascular plant I
    draw_image(
        "figures/vascular_plant_bw.png",
        scale = 0.1,
        x = 0.255,
        y = 0.36
    ) +
    # Bryophytes
    draw_image("figures/moss_bw.png", scale = 0.12, x = 0.135, y = -0.19) +
    # Vascular plants II
    draw_image(
        "figures/vascular_plant_bw.png",
        scale = 0.1,
        x = 0.255,
        y = -0.175
    ))

# Save plot
save_plot(
    "figures/figure_2.png",
    panel_plot_icons,
    ncol = 1,
    nrow = 2,
    base_asp = 1.9,
    bg = "white"
)

## End of File
