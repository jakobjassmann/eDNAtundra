# Preliminary eDNA data viz

# Dependencies
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(gt)
library(vegan)

# Load out list per plot 
load("C:/Users/jakob/repositories/eDNAtundra/data/sequences/ITS2/its2_custom_run2/2_Filtered.project/2_Filtered/Amplicon_analysis/Downstream/R_objects/data.RData")
sample_sums(data)

# Convert OTU table to tibble annotated with species
otu_presence <- otu_table(data) %>% 
  as.data.frame() %>%
  mutate(OTU = rownames(.)) %>%
  relocate(OTU) %>%
  remove_rownames() %>%
  pivot_longer(cols = 2:last_col(), names_to = "plot_id", values_to = "count") %>%
  filter(count > 0) %>%
  select(-count) %>%
  mutate(plot_id = substr(plot_id, 1, 6)) %>%
  full_join(
    tax_table(data) %>%
      as.data.frame()  %>%
      mutate(OTU = rownames(.)) %>%
      relocate(OTU) %>%
      remove_rownames()) %>%
  relocate(plot_id) %>%
  mutate(plot_sub = plot_id, 
         plot_id = substr(plot_id, 1, 4),
         is_pooled = grepl(".*P", plot_sub),
         Species = gsub(".*_(.*)", "\\1", Species))

# Load species lists
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") 

# Calculate species richness and detection rates
detection_df <- distinct(botanical_obs,plot_id) %>%
  pull(plot_id) %>%
  map(function(current_plot){
    data.frame(
      plot_id = current_plot,
      n_species = filter(botanical_obs, plot_id == current_plot) %>% nrow(),
      n_genera = filter(botanical_obs, plot_id == current_plot) %>% 
        mutate(genus = gsub("(.*) .*", "\\1", species)) %>% distinct() %>% nrow(),
      its2_n_otus = filter(otu_presence, plot_id == current_plot) %>% nrow(),
      its2_n_sp_detect = filter(otu_presence, plot_id == current_plot) %>%
        select(Genus, Species) %>%
        na.omit() %>%
        distinct() %>%
        split(1:nrow(.)) %>%
        sapply(function(x){
          x <- as.data.frame(x)
          species_in_plot <- filter(botanical_obs, plot_id == current_plot) %>%
            pull(species)
          any(grepl(x$Genus, species_in_plot) & grepl(x$Species, species_in_plot))
        }) %>% sum(),
      its2_n_gen_detect = filter(otu_presence, plot_id == current_plot) %>%
        select(Genus) %>%
        na.omit() %>%
        distinct() %>%
        split(1:nrow(.)) %>%
        sapply(function(x){
          x <- as.data.frame(x)
          species_in_plot <- filter(botanical_obs, plot_id == current_plot) %>%
            pull(species)
          any(grepl(x$Genus, species_in_plot))
        }) %>% sum()
      )
  }) %>%
  bind_rows() %>%
  mutate(its2_sp_detec_r = round(its2_n_sp_detect / n_species, 2),
         its2_gen_detec_r = round(its2_n_gen_detect / n_genera, 2))

# its2 mean detection rate for species and plot
(mean_rates <- summarise(detection_df, 
                         its2_sp = mean(its2_sp_detec_r) %>% round(2),
                         its2_gen = mean(its2_gen_detec_r) %>% round(2)))

# Current detection rates at species level and genus level are 47% and 66% respectively
# That is not great!

# Do pooled samples tell us the same story as the indiviudal samples combined?
otu_presence <- otu_presence %>% group_by(plot_id) %>% summarise(n_subs = n_distinct(plot_sub)) %>%
  full_join(otu_presence, .)
pooling_stats <- otu_presence %>%
  filter(n_subs == 4) %>%
  split(., .$plot_id) %>%
  map(function(x){
    not_pooled <- filter(x, is_pooled == FALSE)
    pooled <- filter(x, is_pooled == TRUE)
    data.frame(plot_id = pooled$plot_id,
               n_in_sample_1 = filter(not_pooled, grepl("1", plot_sub)) %>% 
                 pull(OTU) %>% n_distinct(),
               n_in_sample_2 = filter(not_pooled, grepl("2", plot_sub)) %>% 
                 pull(OTU) %>% n_distinct(),
               n_in_sample_3 = filter(not_pooled, grepl("3", plot_sub)) %>% 
                 pull(OTU) %>% n_distinct(),
               n_overlap_samples = (filter(not_pooled, grepl("1", plot_sub)) %>% 
                                      pull(OTU) %>%
                                      unique())[
                                        (filter(not_pooled, grepl("1", plot_sub)) %>% 
                                      pull(OTU) %>%
                                      unique()) %in% 
                 (filter(not_pooled, grepl("2", plot_sub)) %>% 
                    pull(OTU) %>% 
                    unique())] %in% 
                 (filter(not_pooled, grepl("3", plot_sub)) %>% 
                    pull(OTU) %>% 
                    unique()) %>% sum(), 
               n_in_not_pooled = length(unique(not_pooled$OTU)),
               n_in_pooled = length(unique(pooled$OTU)),
               n_detected = sum((unique(pooled$OTU) %in% unique(not_pooled$OTU))),
               per_detected = round(sum(unique(pooled$OTU) %in% unique(not_pooled$OTU))/ length(unique(not_pooled$OTU)), 2),
               n_pooled_only = sum(!(unique(pooled$OTU) %in% unique(not_pooled$OTU)))
               ) %>% distinct()
               
  }) %>% bind_rows() 
mean(pooling_stats$per_detected) %>% round(2)
# 65% of all OTUs of the three subsamples in the pooled sample - that is pretty poor
# Working with sum of sub-samples better
pooling_stats <- otu_presence %>%
  filter(n_subs == 4) %>%
  filter(!is.na(Species)) %>%
  mutate(species = paste0(Genus, Species)) %>%
  split(., .$plot_id) %>%
  map(function(x){
    not_pooled <- filter(x, is_pooled == FALSE)
    pooled <- filter(x, is_pooled == TRUE)
    data.frame(plot_id = pooled$plot_id,
               n_in_sample_1 = filter(not_pooled, grepl("1", plot_sub)) %>% 
                 pull(species) %>% n_distinct(),
               n_in_sample_2 = filter(not_pooled, grepl("2", plot_sub)) %>% 
                 pull(species) %>% n_distinct(),
               n_in_sample_3 = filter(not_pooled, grepl("3", plot_sub)) %>% 
                 pull(species) %>% n_distinct(),
               n_overlap_samples = (filter(not_pooled, grepl("1", plot_sub)) %>% 
                                      pull(species) %>%
                                      unique())[
                                        (filter(not_pooled, grepl("1", plot_sub)) %>% 
                                           pull(species) %>%
                                           unique()) %in% 
                                          (filter(not_pooled, grepl("2", plot_sub)) %>% 
                                             pull(species) %>% 
                                             unique())] %in% 
                 (filter(not_pooled, grepl("3", plot_sub)) %>% 
                    pull(species) %>% 
                    unique()) %>% sum(), 
               n_in_not_pooled = length(unique(not_pooled$species)),
               n_in_pooled = length(unique(pooled$species)),
               n_detected = sum((unique(pooled$species) %in% unique(not_pooled$species))),
               per_detected = round(sum(unique(pooled$species) %in% unique(not_pooled$species))/ length(unique(not_pooled$species)), 2),
               n_pooled_only = sum(!(unique(pooled$species) %in% unique(not_pooled$species)))
    ) %>% distinct()
    
  }) %>% bind_rows() 
# Deviation of species in replicast
mean(pooling_stats$per_detected) %>% round(2) # Pooling is really poor for species capture wiht ITS
pooling_stats %>% 
  mutate(average_sample_detection = ((n_in_sample_1 +  
                                              n_in_sample_2 +
                                              n_in_sample_3) / n_in_not_pooled) / 3) %>%
  summarise(mean(average_sample_detection))

# Correspondence analysis of OTUs at site level
ca_otu_site <- otu_presence %>%
  # Remove pooled samples
  filter(is_pooled == FALSE) %>%
  # Keep only one occurance per plot
  split(., .$plot_id) %>%
  map(function(x) distinct(x, plot_id, OTU)) %>%
  bind_rows() %>%
  # Create site identifier
  mutate(site = substr(plot_id, 1, 2)) %>%
  list() %>%
  # Convert to table and do CA analysis
  map(function(x) table(pull(x, site), pull(x, OTU))) %>%
  .[[1]] %>%
  CA()
ca_otu_site

# Correspondence analysis of species ids at site level
ca_otu_site <- otu_presence %>%
  # Remove pooled samples
  filter(is_pooled == FALSE) %>%
  # Keep only one occurance per plot
  split(., .$plot_id) %>%
  map(function(x) distinct(x, plot_id, OTU)) %>%
  bind_rows() %>%
  # Create site identifier
  mutate(site = substr(plot_id, 1, 2)) %>%
  list() %>%
  # Convert to table and do CA analysis
  map(function(x) table(pull(x, site), pull(x, OTU))) %>%
  .[[1]] %>%
  CA()

# Bray curtis dissimelarity
otu_presence %>%
  # Remove pooled samples
  filter(is_pooled == FALSE) %>%
  # Keep only one occurance per plot
  split(., .$plot_id) %>%
  map(function(x) distinct(x, plot_id, OTU)) %>%
  bind_rows() %>%
  # Create site identifier
  mutate(site = substr(plot_id, 1, 2)) %>%
  list() %>%
  # Convert to table and do CA analysis
  map(function(x) table(pull(x, site), pull(x, OTU))) %>%
  .[[1]] %>%
  vegdist()
#           S1        S2
# S2 0.2041643          
# S3 0.1962848 0.2520908



# Repeat for molecularly identified species 
ca_species_site <- otu_presence %>%
  filter(!is.na(Species)) %>%
  mutate(species = paste0(Genus, " ", Species)) %>%
  # Remove pooled samples
  filter(is_pooled == FALSE) %>%
  # Keep only one occurance per plot
  split(., .$plot_id) %>%
  map(function(x) distinct(x, plot_id, species)) %>%
  bind_rows() %>%
  # Create site identifier
  mutate(site = substr(plot_id, 1, 2)) %>%
  list() %>%
  # Convert to table and do CA analysis
  map(function(x) table(pull(x, site), pull(x, species))) %>%
  .[[1]] %>%
  CA()

# Bray Curtis distance
otu_presence %>%
  filter(!is.na(Species)) %>%
  mutate(species = paste0(Genus, " ", Species)) %>%
  # Remove pooled samples
  filter(is_pooled == FALSE) %>%
  # Keep only one occurance per plot
  split(., .$plot_id) %>%
  map(function(x) distinct(x, plot_id, species)) %>%
  bind_rows() %>%
  # Create site identifier
  mutate(site = substr(plot_id, 1, 2)) %>%
  list() %>%
  # Convert to table and do CA analysis
  map(function(x) table(pull(x, site), pull(x, species))) %>%
  .[[1]] %>%
  vegdist()
#           S1        S2
# S2 0.2217573          
# S3 0.2058824 0.2432432

# Plot CA results by species in a nice way

# Simplify objects
ca_site_sp_site_dims <- ca_species_site$row$coord %>%
  as.data.frame() %>%
  mutate(., site = row.names(.))
ca_site_sp_spp_dims <- ca_species_site$col$coord %>%
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
        fill = site),
    colour = "white",
    shape = 21,
    size = 6,
  ) + 
  geom_text(
    data = ca_site_sp_site_dims,
    aes(
      x = `Dim 1`, y = `Dim 2`, label = site,
      color = site
    ),
    nudge_x = 0.2,
    size = 24 /.pt
  ) +
  scale_x_continuous(
    limits = c(-1.5, 1.5),
    breaks = seq(-1.5, 1.5, 0.5)
  ) +
  scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  scale_y_continuous(
    limits = c(-1.5, 1.5),
    breaks = seq(-1.5, 1.5, 0.5)
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.8, 
           label = "eDNA ITS2 vascular species", size = 20 / .pt,
           fontface = "bold", colour = "black") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -1.8, label = "p-value = 0.006",
           size = 20/.pt) +
  labs(x = paste0("Dim 1 (", round(ca_species_site$eig[1, 2], 2), "%)"), y = paste0("Dim 2 (", round(ca_species_site$eig[2, 2], 2), "%)")) +
  theme_cowplot(20) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black")
  )
save_plot("figures/ca_its2_plot_site.png", ca_plot_site,
          bg = "white",
          base_height = 7,
          base_asp = 1.6
)
