# Moss vs. vascular diversity
# Jakob J. Assmann jakob.assmann@uzh.ch 19 March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

# Load presence / absence vascular observations
vascular_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    # Keep only species with confirmed species-level ID
    filter(species != "species") %>%
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
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    filter(!grepl("sp.", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    arrange(site_id, plot_id, short_name) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(phylum = "Bryophytes")

# Combine dataframes
combined_obs <- bind_rows(vascular_obs, bryo_obs) %>%
    arrange(site_id, plot_id, short_name, phylum) %>%
    relocate(site_id)

# Caluclate number of species per plot
combined_obs_nsp_plot <- combined_obs %>%
    group_by(site_id, plot_id, phylum) %>%
    summarise(nsp = n_distinct(short_name)) %>%
    # Add 0 for empty observations
    group_split() %>%
    map(function(x) {
        if (nrow(x) == 2) {
            return(x)
        }
        x <- bind_rows(
            x,
            mutate(
                x,
                phylum = case_when(
                    x$phylum == "Bryophytes" ~ "Vascular plants",
                    x$phylum == "Vascular plants" ~ "Bryophytes"
                ),
                nsp = 0
            )
        )
        return(x)
    }) %>%
    bind_rows() %>%
    arrange(site_id, plot_id, phylum)

# Are mossess more diverse than vascular plants at the plot level?
combined_obs_nsp_plot %>%
    group_by(site_id, plot_id) %>%
    arrange(site_id, plot_id, desc(nsp)) %>%
    slice(1) %>%
    ungroup() %>%
    summarise(bryo_lt_vasc = sum(phylum == "Bryophytes"))
n_distinct(pull(combined_obs_nsp_plot, plot_id))
# No! Bryophyte diverstiy only exceeds vascular diversity for 27 / 54
combined_obs_nsp_plot %>%
    group_by(site_id, plot_id) %>%
    arrange(site_id, plot_id, desc(nsp)) %>%
    pivot_wider(names_from = phylum, values_from = nsp) %>%
    summarise(ratio = Bryophytes / (Bryophytes + `Vascular plants`)) %>%
    ungroup() %>%
    summarise(prec = 100 * mean(ratio))
# # A tibble: 1 × 1
#    prec
#   <dbl>
# 1  48.8

# How is it at site level?
# Caluclate number of species per site
(combined_obs_nsp_site <- combined_obs %>%
    select(-plot_id) %>%
    group_by(site_id, phylum) %>%
    summarise(nsp = n_distinct(short_name)) %>%
    bind_rows() %>%
    arrange(site_id, phylum))
# Bryophytes exeed vascular plant diversity dramatically in all cases
# # A tibble: 6 × 3
#   site_id phylum            nsp
#   <chr>   <chr>           <int>
# 1 Site 1  Bryophytes         50
# 2 Site 1  Vascular plants    30
# 3 Site 2  Bryophytes         32
# 4 Site 2  Vascular plants    31
# 5 Site 3  Bryophytes         57
# 6 Site 3  Vascular plants    38
# Calculate percentage
combined_obs_nsp_site %>%
    pivot_wider(names_from = phylum, values_from = nsp) %>%
    mutate(ratio = Bryophytes / (Bryophytes + `Vascular plants`)) %>%
    ungroup() %>%
    summarise(perc = round(100 * mean(ratio)))
# A tibble: 1 × 1
#    perc
#   <dbl>
# 1    58

# What about Shannon diversity?
(combined_obs_div_site <- combined_obs %>%
    group_by(site_id, phylum, short_name) %>%
    summarise(p_ln_p = (n() / 13) * log(n() / 13)) %>%
    group_by(site_id, phylum) %>%
    summarise(shan_div = (-1) * sum(p_ln_p)) %>%
    arrange(site_id, phylum))
# Same here!
#   site_id phylum          shan_div
#  site_id phylum          shan_div
#   <chr>   <chr>              <dbl>
# 1 Site 1  Bryophytes         13.1
# 2 Site 1  Vascular plants     6.75
# 3 Site 2  Bryophytes          8.82
# 4 Site 2  Vascular plants     7.77
# 5 Site 3  Bryophytes         14.9
# 6 Site 3  Vascular plants     7.77

# Visualise at the site level
(bryo_vs_vasc_plot <- ggplot(combined_obs_nsp_site) +
    geom_col(aes(x = site_id, y = nsp, fill = phylum)) +
    labs(x = "", y = "Number of Species (n)") +
    scale_fill_manual(
        values = c("#a6212166", "#a62121"),
        labels = c("Bryophytes", "Vascular Plants")
    ) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 20)) +
    theme_cowplot() +
    theme(
        legend.title = element_blank(),
        legend.margin = margin(0.25, 0.25, 0.25, 0.25, unit = "cm"),
        legend.position = "right",
        legend.box.background = element_rect(colour = "black"),
        legend.justification.right = c(0.5, 0.05),
        axis.title.x = element_blank()
    ))

(bryo_vs_vasc_plot_icons <- ggdraw() +
    draw_plot(bryo_vs_vasc_plot) +
    draw_image("figures/moss_bw.png", scale = 0.275, x = 0.35, y = 0.325) +
    draw_image("figures/vascular_plant_bw.png", scale = 0.3, x = 0.34, y = 0))

save_plot("figures/fig_2b.png", bryo_vs_vasc_plot_icons, bg = "white")
