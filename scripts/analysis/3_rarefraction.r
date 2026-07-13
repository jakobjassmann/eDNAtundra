# Rarefraction analysis to optimise sampling
# Jakob J. Assmann jakob.assmann@uzh.zh 19 March 2025

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)
library(pbapply)
library(parallel)

# Load combined observations
obs_combined <- read_csv("data/analysis/detection_species_genus.csv")

# Add polot type and prepare data for analsysis
obs_combined <- obs_combined %>%
    # Keep only botanical and (strict species-level) eDNA observations
    select(plot_id, short_name, botanical, eDNA) %>%
    mutate(
        plot_type = case_when(
            grepl("S.S.", plot_id) ~ "systematic",
            grepl("S.R.", plot_id) ~ "stratified random"
        )
    ) %>%
    # Add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    relocate(site_id, plot_id, plot_type) %>%
    # Add labeld combinding plot_id and observation type
    pivot_longer(botanical:eDNA, names_to = "source", values_to = "present") %>%
    filter(present) %>%
    mutate(plot_source = paste(plot_id, source))


# Generate all possible plot combinations splitting startificaiton systems
all_plot_combos <- obs_combined %>%
    distinct(site_id, plot_type, plot_id) %>%
    arrange(site_id, plot_type, plot_id) %>%
    group_by(site_id, plot_type) %>%
    group_split() %>%
    map(function(x) {
        # A combinaiton can have between 1-9 plots (9 = all)
        map(1:9, function(i) {
            # Get all possible plot combinations for i plots
            plot_combos <- combn(x$plot_id, i)
            # Get all possible source combinations for i plots
            source_type_combos <- sapply(1:i, function(y) {
                c(rep("eDNA", y), rep("botanical", i - y))
            })
            # Expand all source combinaitons for all plot combinations
            combos_dfs <- map(1:ncol(plot_combos), function(j) {
                combos_df <- data.frame(
                    # Expand combinations
                    plot_source = paste(plot_combos[, j], source_type_combos),
                    # Carry over site and stratification ids
                    site_id = unique(x$site_id),
                    plot_type = unique(x$plot_type),
                    # number of plots samples
                    n_sampled = i,
                    # Assging unique id to ech set of i combinations
                    id = paste0(
                        # Site
                        gsub(" ", "_", unique(x$site_id)),
                        "_",
                        # Stratification type
                        gsub(" ", "_", unique(x$plot_type)),
                        "_",
                        # Number of plots in combination
                        i,
                        "_",
                        # ID for each combo based
                        j,
                        "_",
                        # ID for each combo based on source types
                        sort(rep(1:i, i))
                    )
                )
            }) %>%
                bind_rows()
        }) %>%
            bind_rows()
    }) %>%
    bind_rows()

# Check total number of combinations
nrow(all_plot_combos)
# This requires parallel processing!

# prepare parallel envrionment
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, library(tidyverse))
n_max <- obs_combined %>%
    filter(source == "botanical" & present) %>%
    distinct(site_id, short_name) %>%
    group_by(site_id) %>%
    summarise(n_max = n())
clusterExport(cl, varlist = c("obs_combined", "n_max"))

# Calculate proportion of correctly identfied species per plot
fract_detected <- all_plot_combos %>%
    split(.$id) %>%
    pblapply(cl = cl, FUN = function(combo) {
        # Look up maximum number of species to detect
        n_max <- filter(n_max, site_id == unique(combo$site_id)) %>%
            pull(n_max)
        n_obs <- filter(obs_combined, plot_source %in% combo$plot_source) %>%
            distinct(short_name) %>%
            nrow()
        # Return proportion as a data frame
        return(data.frame(
            site_id = unique(combo$site_id),
            plot_type = unique(combo$plot_type),
            n_sampled = unique(combo$n_sampled),
            fract_detected = n_obs / n_max
        )) %>%
            bind_rows()
    }) %>%
    bind_rows()

# Save calculatoins
write_csv(fract_detected, "data/analysis/fract_detected.csv")

# Visualise
(rarefraction <- ggplot() +
    geom_line(
        aes(
            x = n_sampled,
            y = fract_detected * 100,
            colour = plot_type,
            linetype = plot_type
            # group = paste0(site_id, "_", plot_type)
        ),
        data = fract_detected %>%
            group_by(site_id, plot_type, n_sampled) %>%
            summarise(fract_detected = mean(fract_detected)),
        linewidth = 1
    ) +
    geom_ribbon(
        aes(
            x = n_sampled,
            ymin = min_perc * 100,
            ymax = max_perc * 100,
            fill = plot_type,
            linetype = plot_type,
            group = paste0(site_id, "_", plot_type)
        ),
        data = fract_detected %>%
            group_by(site_id, plot_type, n_sampled) %>%
            summarise(
                min_perc = min(fract_detected),
                max_perc = max(fract_detected)
            ),
        alpha = 0.1
    ) +
    scale_x_continuous(breaks = 1:9) +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 110, 10)) +
    scale_colour_manual(values = c("#2A0081", "#BA1650")) + # yellow #F5C700
    scale_fill_manual(values = c("#2A0081", "#BA1650")) +
    facet_wrap(~site_id) +
    labs(
        x = "Number of plots sampled",
        y = "Vascular species detected (%)",
        fill = "Sampling Strategy",
        colour = "Sampling Strategy",
        linetype = "Sampling Strategy"
    ) +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "white")))
save_plot(
    "figures/figure_3_optimisation.png",
    rarefraction,
    base_asp = 4,
    bg = "white"
)
