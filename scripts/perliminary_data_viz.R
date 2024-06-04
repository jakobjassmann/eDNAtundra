# Preliminary data analysis and visualisation for the botanical
# observations data 
# Jakob J. Assmann 2023-08-15 j.assmann@uzh.ch
# Updated 4 June 2024

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(gt)

# Load data
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") 

# generate separate df containing only species where ID was certain
botanical_obs_certain <- botanical_obs %>%
    filter(!grepl("cf", species))

# Calculate species richness for each site
plot_richness <- botanical_obs %>%
    group_by(site, plot_id) %>%
    summarise(species_richness = n_distinct(species)) 
site_richness <- botanical_obs_certain %>%
    group_by(site) %>%
    summarise(species_richness = n_distinct(species)) 

# Visualise

# Plot richness
plot_richness_plot <- ggplot(plot_richness) +
    geom_col(aes(
        x = plot_id, y = species_richness,
        fill = site
    )) +
    labs(y = "Vascular plant species richness (n)", fill = "Site",
    x = "Plot") +
scale_y_continuous(limits = c(0,20)) +
theme_cowplot() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))

# Site richness
site_richness_plot <- ggplot(site_richness) +
    geom_col(aes(
        x = site, y = species_richness,
        fill = site
    )) +
    labs(
        x = "",
        y = "Vascular plant species richness (n)"
    ) +
scale_y_continuous(limits = c(0, 50)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90), legend.position = "none")

# write out plots
plot_grid(site_richness_plot, plot_richness_plot) %>%
    save_plot("figures/vascular_richness_plots.png", ., ncol = 2, bg = "white")

# Identify common species
species_frequencies <- data.frame(
    species = names(table(botanical_obs$species)),
    n_plots = as.vector(table(botanical_obs$species))
) %>%
    mutate(spec_freq = n_plots / 54) %>%
        arrange(desc(n_plots))
common_species <- species_frequencies %>%
    filter(n_plots > 10) %>%
    pull(species)

# Generate a plot specific contingency table
plot_species <- table(
    botanical_obs_certain$plot_id,
    botanical_obs_certain$species
)
site_species <- table(
    botanical_obs_certain$site,
    botanical_obs_certain$species
)

# Correspondence analysis plot level
CA(plot_species)
# -> no significant association between plot and speices

# Correspondence analysis site level
(ca_site_species <- CA(site_species))
# -> significant association between site and species
# chi square = 136.8995 (p-value =  0.00579884)

# Visualise site level correspondence analysis
get_eigenvalue(ca_site_species)
fviz_screeplot(ca_site_species)
fviz_ca_biplot(ca_site_species)

# Simplify objects
ca_site_sp_site_dims <- ca_site_species$row$coord %>%
    as.data.frame() %>%
    mutate(., site = row.names(.))
ca_site_sp_spp_dims <- ca_site_species$col$coord %>%
    as.data.frame() %>%
    mutate(., species = row.names(.))

ca_plot_site <- ggplot() +
    geom_point(
        data = ca_site_sp_spp_dims,
        aes(x = `Dim 1`, y = `Dim 2`)
    ) +
    geom_text_repel(
        data = ca_site_sp_spp_dims
        # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
        ,
        aes(
            x = `Dim 1`, y = `Dim 2`,
            label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
        ),
        max.overlaps = 20
    ) +    
    geom_point(
        data = ca_site_sp_site_dims,
        aes(x = `Dim 1`, y = `Dim 2`,
        color = site),
        size = 2,
    ) + 
    geom_text(
        data = ca_site_sp_site_dims,
        aes(
            x = `Dim 1`, y = `Dim 2`, label = site,
            color = site
        ),
        nudge_x = 0.1,
        size = 5
    ) +
    scale_x_continuous(
        limits = c(-1.5, 1.5),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
        geom_hline(yintercept = 0, alpha = 0.5) +
        scale_y_continuous(
            limits = c(-1.5, 1.5),
            breaks = seq(-1.5, 1.5, 0.5)
        ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.8, 
           label = "Vascular plant species", size = 18 / .pt,
           fontface = "bold", colour = "darkgreen") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -1.8, label = "p-value = 0.006") +
    labs(x = paste0("Dim 1 (", round(ca_site_species$eig[1, 2], 2), "%)"), y = paste0("Dim 2 (", round(ca_site_species$eig[2, 2], 2), "%)")) +
    theme_cowplot() +
        theme(
            legend.position = "none",
            panel.border = element_rect(color = "black")
        )
save_plot("figures/ca_vascular_plot_site.png", ca_plot_site,
    bg = "white",
    base_height = 7,
    base_asp = 1.6
)


## Mosses!

rm(list = ls())

# Load moss presence absence
moss_obs <- read_csv("data/botanical_obs/Moss_data_EC_Summer2023_UZH.csv") %>%
  # Remove redundant plot_id column
  select(-Plot_ID) %>%
  # Remove some empty columns
  filter(!is.na(Site))

# Make columnames lowercase for easier handling
colnames(moss_obs) <- c("site", "plot_id", "genus", "sp", "species", "presence", "date_entered") 

# Calculate species richness for each site
plot_richness <- moss_obs %>%
  group_by(site, plot_id) %>%
  summarise(species_richness = n_distinct(species)) 
site_richness <- moss_obs %>%
  group_by(site) %>%
  summarise(species_richness = n_distinct(species)) 

# Plot richness
plot_richness_plot <- ggplot(plot_richness) +
  geom_col(aes(
    x = plot_id, y = species_richness,
    fill = site
  )) +
  labs(y = "Moss species richness (n)", fill = "Site",
       x = "Plot") +
  scale_y_continuous(limits = c(0,20)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5))

# Site richness
site_richness_plot <- ggplot(site_richness) +
  geom_col(aes(
    x = site, y = species_richness,
    fill = site
  )) +
  labs(
    x = "",
    y = "Moss species richness (n)"
  ) +
  scale_y_continuous(limits = c(0, 70)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")

# write out plots
plot_grid(site_richness_plot, plot_richness_plot) %>%
  save_plot("figures/moss_richness_plots.png", ., ncol = 2, bg = "white")

# Identify common species
species_frequencies <- data.frame(
  species = names(table(moss_obs$species)),
  n_plots = as.vector(table(moss_obs$species))
) %>%
  mutate(spec_freq = n_plots / 54) %>%
  arrange(desc(n_plots))
common_species <- species_frequencies %>%
  filter(n_plots > 10) %>%
  pull(species)
species_frequencies %>%
  mutate(spec_freq = round(spec_freq, 2)) %>%
  select(Species = species,
         `Number of Plots` = n_plots,
         `Proportion of Plots` = spec_freq) %>%
  gt() %>%
    gtsave("tables/moss_species_frequency.html")

# Generate a plot specific contingency table
plot_species <- table(
  moss_obs$plot_id,
  moss_obs$species
)
site_species <- table(
  moss_obs$site,
  moss_obs$species
)

# Correspondence analysis plot level
CA(plot_species)
# -> Again significant association between plot and species!

# Correspondence analysis site level
(ca_site_species <- CA(site_species))
# -> no significant association between site and species
# chi square = 148.9782 (p-value =  0.3709957)

# Visualise site level correspondence analysis
get_eigenvalue(ca_site_species)
fviz_screeplot(ca_site_species)
fviz_ca_biplot(ca_site_species)

# Simplify objects
ca_site_sp_site_dims <- ca_site_species$row$coord %>%
  as.data.frame() %>%
  mutate(., site = row.names(.))
ca_site_sp_spp_dims <- ca_site_species$col$coord %>%
  as.data.frame() %>%
  mutate(., species = row.names(.))

ca_plot_site <- ggplot() +
  geom_point(
    data = ca_site_sp_spp_dims,
    aes(x = `Dim 1`, y = `Dim 2`)
  ) +
  geom_text_repel(
    data = ca_site_sp_spp_dims
    # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
    ,
    aes(
      x = `Dim 1`, y = `Dim 2`,
      label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
    ),
    max.overlaps = 20
  ) +    
  geom_point(
    data = ca_site_sp_site_dims,
    aes(x = `Dim 1`, y = `Dim 2`,
        color = site),
    size = 2,
  ) + 
  geom_text(
    data = ca_site_sp_site_dims,
    aes(
      x = `Dim 1`, y = `Dim 2`, label = site,
      color = site
    ),
    nudge_x = 0.1,
    size = 5,
    fontface = "bold"
  ) +
  scale_x_continuous(
    limits = c(-2, 1.5),
    breaks = seq(-2, 1.5, 0.5)
  ) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  scale_y_continuous(
    limits = c(-2, 1),
    breaks = seq(-2, 1, 0.5)
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.8, 
           label = "Moss species", size = 18 / .pt,
           fontface = "bold", colour = "darkgreen") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -1.8, label = "p-value = 0.371") +
  labs(x = paste0("Dim 1 (", round(ca_site_species$eig[1, 2], 2), "%)"), y = paste0("Dim 2 (", round(ca_site_species$eig[2, 2], 2), "%)")) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black")
  )
save_plot("figures/ca_moss_plot_site.png", ca_plot_site,
          bg = "white",
          base_height = 7,
          base_asp = 1.6
)

