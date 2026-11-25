# Rarefraction analysis to optimise sampling
# Jakob J. Assmann jakob.assmann@uzh.zh 19 March 2025

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)
library(pbapply)
library(parallel)
library(RcppAlgos)

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

# Get all combinations of plots, sampling schemes and data sources
all_plots <- expand.grid(
    plot_id = unique(obs_combined$plot_id),
    source = unique(obs_combined$source)
) %>%
    mutate(
        plot_type = case_when(
            grepl("S.S.", plot_id) ~ "systematic",
            grepl("S.R.", plot_id) ~ "stratified random"
        )
    ) %>%
    mutate(
        site_id = gsub("S(.).*", "Site \\1", plot_id),
        combo_id = paste(site_id, plot_id, plot_type, source),
    ) %>%
    relocate(site_id)

# For each site, sample plots with the following crietria
# - Nine sets of samples
#   - Stratified set: eDNA, botanical and both
#   - Systematic set: eDNA, botanical and both
#   - Mixed set: eDNA, botanical and both
# - Samples contain 1-18 plots (where possible)
# - Samples contain 50 samples for each number of plots (where possible)

# Helper function
get_sample <- function(plots_subset, max_plots, max_sample = max_plots) {
    pblapply(1:max_plots, x = plots_subset, function(y, x) {
        {
            # Check if the whole existing set needs to be included
            # and set max counter
            z_max <- if_else(
                max(c(choose(max_plots, y), choose(max_sample, y))) >= 100,
                100,
                max(c(choose(max_plots, y), choose(max_sample, y)))
            )
            # Empty list
            z <- list()
            cat(z_max, "\n")
            while (length(z) < z_max) {
                sample <- data.frame(
                    combo_id = sample(x$combo_id, y, replace = FALSE)
                )
                if (!(sample %in% z)) {
                    z <- c(z, sample)
                }
            }
            return(z)
        } %>%
            do(c) %>%
            unlist(recursive = FALSE)
    })
}
all_combos <- bind_rows(
    # Stratified set eDNA
    all_plots %>%
        filter(plot_type == "stratified random", source == "eDNA") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("stratified eDNA", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # Stratified set botanical
    all_plots %>%
        filter(plot_type == "stratified random", source == "botanical") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("stratified botanical", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # Stratified set eDNA and botancial
    all_plots %>%
        filter(
            plot_type == "stratified random",
            source == "botanical" | source == "eDNA"
        ) %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9, max_sample = 18) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("stratified both", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # Systematic set eDNA
    all_plots %>%
        filter(plot_type == "systematic", source == "eDNA") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("systematic eDNA", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # Systematic set botanical
    all_plots %>%
        filter(plot_type == "systematic", source == "botanical") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("systematic botanical", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # Systematic set eDNA and botancial
    all_plots %>%
        filter(
            plot_type == "systematic",
            source == "botanical" | source == "eDNA"
        ) %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 9, max_sample = 18) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("systematic both", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # All eDNA
    all_plots %>%
        filter(source == "eDNA") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 18) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("all eDNA", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # All botanical
    all_plots %>%
        filter(source == "botanical") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 18) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("all botanical", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample)),
    # All eDNA and botancial
    all_plots %>%
        filter(source == "botanical" | source == "eDNA") %>%
        split(.$site_id) %>%
        map(get_sample, max_plots = 18, max_sample = 36) %>%
        unlist(recursive = FALSE) %>%
        unlist(recursive = FALSE) %>%
        set_names(paste0("all both", 1:length(.))) %>%
        tibble(sample = .) %>%
        mutate(sample_id = names(sample))
)

# prepare parallel envrionment
stopCluster(cl)
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, library(tidyverse))
n_max <- obs_combined %>%
    filter(source == "botanical" & present) %>%
    distinct(site_id, short_name) %>%
    group_by(site_id) %>%
    summarise(n_max = n())
clusterExport(cl, varlist = c("all_plots", "obs_combined", "n_max"))

# Calculate proportion of correctly identfied species per plot
n_detected <- all_combos %>%
    split(.$sample_id) %>%
    pblapply(
        cl = cl,
        FUN = function(combo) {
            plots_df <- filter(all_plots, combo_id %in% unlist(combo$sample))
            current_site <- unique(plots_df$site_id)
            # Look up maximum number of species to detect
            n_max_site <- filter(n_max, site_id == current_site) %>%
                pull(n_max)
            n_obs <- filter(
                obs_combined,
                (site_id == current_site) &
                    (plot_id %in% unique(plots_df$plot_id)) &
                    (source %in% unique(plots_df$source))
            ) %>%
                distinct(short_name) %>%
                nrow()
            # Return proportion as a data frame
            return(tibble(
                sample = sort(unique(combo$sample_id)),
                site_id = current_site,
                plot_type = paste0(
                    sort(unique(plots_df$plot_type)),
                    collapse = " and "
                ),
                method = paste0(
                    sort(unique(plots_df$source)),
                    collapse = " and "
                ),
                n_sampled = nrow(plots_df),
                n_detected = n_obs,
                n_max_site = n_max_site,
                fract_detected = n_obs / n_max_site
            )) %>%
                bind_rows()
        }
    ) %>%
    bind_rows()

# Save calculations
write_csv(n_detected, "data/analysis/n_detected.csv")
# n_detected <- read_csv("data/analysis/n_detected")

# Prepare data frame for plotting
n_detected_plot <- n_detected %>%
    group_by(site_id, method, plot_type, n_sampled) %>%
    summarise(
        mean_detected = mean(n_detected),
        min_detected = min(n_detected),
        max_detected = max(n_detected),
        n_max_site = unique(n_max_site),
        eighty_perc = unique(n_max_site) * 0.8
    ) %>%
    mutate(
        panel_id = factor(
            paste(site_id, plot_type),
            levels = c(
                "Site 1 systematic",
                "Site 2 systematic",
                "Site 3 systematic",
                "Site 1 stratified random",
                "Site 2 stratified random",
                "Site 3 stratified random",
                "Site 1 stratified random and systematic",
                "Site 2 stratified random and systematic",
                "Site 3 stratified random and systematic"
            )
        )
    )

# Visualise
(rarefraction <- ggplot(
    data = n_detected_plot
) +
    geom_hline(
        aes(yintercept = eighty_perc),
        linetype = "dashed",
        colour = "grey70"
    ) +
    geom_hline(
        aes(yintercept = n_max_site),
        linetype = "dashed",
        colour = "grey70"
    ) +
    geom_text(
        aes(y = eighty_perc),
        x = 19.5,
        label = "80%",
        size = 10 / .pt,
        vjust = 1.5,
        hjust = 1,
        colour = "grey70"
    ) +
    geom_text(
        aes(y = n_max_site),
        x = 19.5,
        label = "100%",
        size = 10 / .pt,
        vjust = 1.5,
        hjust = 1,
        colour = "grey70"
    ) +
    geom_line(
        aes(
            x = n_sampled,
            y = mean_detected,
            colour = method
            #linetype = plot_type
            # group = paste0(site_id, "_", plot_type)
        ),
        linewidth = 1
    ) +
    geom_ribbon(
        aes(
            x = n_sampled,
            ymin = min_detected,
            ymax = max_detected,
            fill = method
            #linetype = plot_type
        ),
        alpha = 0.1
    ) +
    scale_x_continuous(limits = c(1, 19), breaks = 1:18) +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_colour_manual(values = c("#F28F16", "#A62121", "#223F59")) + # yellow #F5C700
    scale_fill_manual(values = c("#F28F16", "#A62121", "#223F59")) +
    facet_wrap(~panel_id, scales = "free") +
    labs(
        x = "Number of plots sampled",
        y = "Vascular species detected (n)",
        fill = "Method",
        colour = "Method"
        #linetype = "Sampling Strategy"
    ) +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "white")))
save_plot(
    "figures/figure_SX_rarefaction.png",
    rarefraction,
    base_asp = 2.5,
    base_height = 7,
    bg = "white"
)

obs_combined %>%
    filter(site_id == "Site 1") %>%
    distinct(plot_id, source, plot_type)
