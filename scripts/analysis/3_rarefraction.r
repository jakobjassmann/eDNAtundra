# Rarefraction analysis to optimise sampling
# Jakob J. Assmann jakob.assmann@uzh.zh 19 March 2025

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)
library(pbapply)
library(parallel)

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
    mutate(source = "botanical")


# Load presence / absence from eDNA
trnl_obs <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa_in_sample.csv") %>%
    mutate(short_name = paste(genus, species), sq = "trnl")
its2_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa_in_sample.csv") %>%
    mutate(short_name = paste(genus, species), sq = "its2")
eDNA_obs <- bind_rows(trnl_obs, its2_obs) %>%
    group_by(sample_id, short_name) %>%
    group_split() %>%
    map(function(x) {
        # Arrange by count
        x <- arrange(x, count)
        # Get unique markers
        sqs <- sort(unique(x$sq)) %>% paste(collapse = " ")
        # Slice tor row with most counts
        x <- slice(x, 1)
        # Update markers
        x$sq <- sqs
        return(x)
    }) %>%
    bind_rows() %>%
    filter(!is.na(plot_id)) %>%
    select(-subsample_id, -sq, -count, -sample_id) %>%
    mutate(site_id = substr(plot_id, 1,2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(source = "eDNA") 

# combine
obs_combined <- bind_rows(vascular_obs, 
    eDNA_obs) %>%
    mutate(plot_type = case_when(
        grepl("S.S.", plot_id) ~ "systematic",
        grepl("S.R.", plot_id) ~ "stratified random"
    )) %>%
    # Collapse plot_id and source
    mutate(plot_source = paste(plot_id, source)) %>%
    # Keep only ids to species level except for Salix
    filter(species != "species" | genus == "Salix")


# Generate all possible plot combinations splitting startificaiton systems
all_plot_combos <- obs_combined  %>%
    distinct(site_id, plot_type, plot_id) %>%
    arrange(site_id, plot_type, plot_id) %>%
    group_by(site_id, plot_type) %>%
    group_split() %>%
    map(function(x){
        # A combinaiton can have between 1-9 plots (9 = all)
        map(1:9, function(i){
            # Get all possible plot combinatoins for i plots
        plot_combos <- combn(x$plot_id, i)
        # Get all possible source combinations for i plots
        source_type_combos <- sapply(1:i, function(y) c(rep("eDNA", y), rep("botanical",i-y)))
        # Expand all source combinaitons for all plot combinations
        combos_dfs <- map(1:ncol(plot_combos), function(j){
            combos_df <- data.frame(
                # Expand combinations
                plot_source = paste(plot_combos[,j], source_type_combos),
                # Carry over site and stratification ids
                site_id = unique(x$site_id),
                plot_type = unique(x$plot_type),
                # number of plots samples
                n_sampled = i,
                # Assging unique id to ech set of i combinations
                id = paste0(
                    # Site
                    gsub(" ", "_", unique(x$site_id)), "_",
                    # Stratification type
                    gsub(" ", "_", unique(x$plot_type)), "_",
                    # Number of plots in combination
                    i, "_",
                    # ID for each combo based
                    j, "_",
                    # ID for each combo based on source types
                    sort(rep(1:i, i))))}) %>%
            bind_rows()
        }) %>% bind_rows()
    }) %>% bind_rows()

# Check total number of combinations
nrow(all_plot_combos)
# This requires parallel processing!

# prepare parallel envrionment
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, library(tidyverse))
n_max <- obs_combined %>%
         group_by(site_id, plot_type) %>% 
         summarise(n_max == n())
    nrow()
clusterExport(cl, varlist = c("obs_combined", "n_max"))

# Calculate proportion of correctly identfied species per plot
fract_detected <- all_plot_combos %>%
split(.$id) %>%
pblapply(cl = cl, FUN = function(combo){
    # Calculate maximum species to be detected
    n_max <- filter(obs_combined, 
        site_id == unique(combo$site_id),
        source == "botanical") %>%
        distinct(short_name) %>% 
        # Remove the only Salix sp in the botanical obs (Site 1)
        filter(short_name != "Salix species") %>% 
        nrow()
    n_obs <- filter(obs_combined, plot_source %in% combo$plot_source) %>%
    distinct(genus, species, short_name) %>%
    split(.$genus) %>%
    map(function(x) {
        # If it's not a Salix species skip
        if (unique(x$genus) != "Salix") {
            return(x)
        }
        # If it's just one observaiton of a Salix species skip
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
    }) %>% bind_rows() %>%
    distinct() %>%
    nrow()
    # Return proportion as a data frame
    return(data.frame(
        site_id = unique(combo$site_id),
    plot_type = unique(combo$plot_type),
    n_sampled = unique(combo$n_sampled),
    fract_detected = n_obs / n_max
    )) %>% bind_rows()
}) %>% bind_rows()

# Save calculatoins
write_csv(fract_detected, "data/analysis/fract_detected.csv")

# Visualise
rarefraction <- ggplot() +
    geom_line(
        aes(
            x = n_sampled, y = fract_detected * 100, colour = site_id,
            linetype = plot_type, group = paste0(site_id, "_", plot_type)
        ),
        data = fract_detected %>%
            group_by(site_id, plot_type, n_sampled) %>%
            summarise(fract_detected = mean(fract_detected)),
        linewidth = 1
    ) +
    geom_ribbon(
        aes(
            x = n_sampled, ymin = min_perc * 100, ymax = max_perc * 100,
            fill = site_id,
            linetype = plot_type, group = paste0(site_id, "_", plot_type)
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
    scale_y_continuous(limits = c(0, 104), breaks = 100 * seq(0, 1, 0.1)) +
    scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    labs(
        x = "Number of plot sampled", y = "Vascular species detected (%)",
        fill = "Site", colour = "Site",
        linetype = "Sampling Strategy"
    ) +
    theme_cowplot()
save_plot("figures/figure_3_optimisation.png", rarefraction,
    base_asp = 2,
    bg = "white"
)

