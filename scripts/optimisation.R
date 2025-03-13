# Sampling scheme opimisation

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(gt)
library(pbapply)

# Load data
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv")

# generate separate df containing only species where ID was certain
botanical_obs_certain <- botanical_obs %>%
  filter(!grepl("cf", species))

# Species per site
n_spp_site <- group_by(botanical_obs_certain, site) %>%
    summarise(n = n_distinct(species))
# Species richness per site depending on random vs. systematic plots
n_spp_per_plot_type <- botanical_obs_certain %>%
  mutate(plot_type = case_when(grepl("S.S.", plot_id) ~ "systematic",
                               grepl("S.R.", plot_id) ~ "stratified random")) %>%
  group_by(site, plot_type) %>%
  group_split() %>%
  map(function(x) data.frame(site = unique(x$site),
                             plot_type = unique(x$plot_type),
                             n_spp = n_distinct(x$species))) %>%
  bind_rows()
n_spp_per_plot_type %>%
  mutate(perc_total = n_spp /   (bind_rows(n_spp_site, n_spp_site) %>%
           arrange(site) %>% pull(n)) ) %>%
  pivot_wider(names_from = plot_type,  values_from = perc_total) %>%
  summarise_at(c("system", "strat_random"), ~ mean(.x, na.rm = T))
# On average sysematic sampling detected 63% of the total species
# Stratified ranodm sampling detected 91% of the total species at the site
t.test(n_spp ~ plot_type, n_spp_per_plot_type)

# Not significnat

# How much would we on average detect using 3 sampling plots?
prec_detected <- botanical_obs_certain %>%
  mutate(plot_type = case_when(grepl("S.S.", plot_id) ~ "systematic",
                               grepl("S.R.", plot_id) ~ "stratified random")) %>%
  group_by(site, plot_type) %>%
  group_split() %>%
  map(function(x){
    cat(unique(x$site), " ", unique(x$plot_type), "\n")
    return_df <- data.frame()
    max_plots <- n_distinct(x$plot_id)
    for(i in 1:max_plots){
      cat("sampling ", i, "plots\n")
      combinations <- combn(unique(x$plot_id), i) %>%
        t() %>%
        split(., 1:nrow(.))
      return_df <- bind_rows(return_df,
                data.frame(
                  site = unique(x$site),
                  plot_type = unique(x$plot_type),
                  n_sampled = i,
                  perc_detected = pblapply(combinations, x = x, function(y, x){
                    return((filter(x, plot_id %in% y) %>%
                       summarise(n = n_distinct(species)) %>%
                         pull(n)) /
                      (filter(n_spp_site, site %in% x$site) %>% pull(n)))
                  }) %>% unlist() ))
    }
    return(return_df)
  }) %>%
  bind_rows()

# Visualise
rarefraction <- ggplot() +
  geom_line(aes(x = n_sampled , y = perc_detected * 100, colour = site,
                 linetype = plot_type, group = paste0(site, "_", plot_type)),
            data = prec_detected %>%
              group_by(site, plot_type, n_sampled ) %>%
              summarise(perc_detected = mean(perc_detected)),
            linewidth = 1) +
  geom_ribbon(aes(x = n_sampled , ymin = min_perc * 100, ymax = max_perc * 100,
                  fill = site,
                linetype = plot_type, group = paste0(site, "_", plot_type)),
            data = prec_detected %>%
              group_by(site, plot_type, n_sampled ) %>%
              summarise(min_perc = min(perc_detected),
                        max_perc = max(perc_detected)),
            alpha = 0.1) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(limits = c(0,100), breaks = 100* seq(0,1,0.1)) +
  scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  labs(x = "Number of plot sampled", y = "Vascular species detected (%)",
       fill = "Site", colour = "Site", 
       linetype = "Sampling Strategy") +
  theme_cowplot()
save_plot("figures/rarefraction_vascular.png", rarefraction,
          base_asp = 2,
          bg = "white")    
