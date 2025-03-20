# Figure 4: Byophyte vs. vascular diversity
# Jakob J. Assmann jakob.asmann@uzh.ch 19 March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

# Source analysis script to get the data object
source("scripts/analysis/4_bryophyte_vs_vascular_diversity.r")

(richness_plot <- plot_grid(
    #Plot level
    ggplot() +
  geom_col(data = filter(combined_obs_nsp_plot, phylum == "Bryophytes"),
  aes(
    x = plot_id, y = nsp,
    fill = site_id,
    colour = site_id)) +
    geom_col(data = filter(combined_obs_nsp_plot, phylum != "Bryophytes"),
  aes(
    x = plot_id, y = nsp), alpha = 0, color = "black") +
  labs(y = "Brophyte species richness (n)", fill = "Site",
       x = "Plots") +
scale_colour_manual(values = c("#F2D64B", 
                               "#5F1773",
                               "#F2507B")) +
scale_fill_manual(values = c("#F2D64B",
                               "#5F1773", 
                               "#F2507B")) +
scale_y_continuous(limits = c(0,20)) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"),
# Site level
ggplot() +
  geom_col(data = filter(combined_obs_nsp_site, phylum == "Bryophytes"),
  aes(
    x = site_id, y = nsp,
    fill = site_id,
    colour = site_id)) +
    geom_col(data = filter(combined_obs_nsp_site, phylum != "Bryophytes"),
  aes(
    x = site_id, y = nsp), alpha = 0, color = "black", linewidth = 1) +
  labs(y = "Bryophyte species richness (n)", fill = "Site",
       x = "Site") +
scale_colour_manual(values = c("#F2D64B", 
                               "#5F1773",
                               "#F2507B")) +
scale_fill_manual(values = c("#F2D64B",
                               "#5F1773", 
                               "#F2507B")) +
scale_y_continuous(limits = c(0,70), breaks = seq(0,70, 10)) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")))
save_plot("figures/fig_4_bryo_vs_vasc_richness.png", 
    richness_plot,
    bg = "white",
    base_asp = 1.6*2
)
