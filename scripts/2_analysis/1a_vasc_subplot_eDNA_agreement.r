# eDNA detection performance and quality control
# Jakob J. Assmann 18 March 2025 jakob.assmann@uh.ch

# Dependencies
library(tidyverse)

# Load presence / absence from eDNA
trnl_obs <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa.csv") %>%
    mutate(short_name = paste(genus, species), sq = "trnl")
its2_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa.csv") %>%
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
    bind_rows()

# Consistency of detection in sub-samples ----
eDNA_obs %>%
    filter(subsample_id != "P" & !is.na(subsample_id)) %>%
    group_by(plot_id, short_name) %>%
    summarise(n_detected = n(), prop_detected = n() / 3) %>%
    ungroup() %>%
    summarise(
        mean_n = mean(n_detected),
        sd_n = sd(n_detected),
        mean_prop = mean(prop_detected),
        sd_prop = sd(prop_detected)
    )
# Samples are consistently detected across subplots
# Mean detection in 2.1 plots (sd = 0.880 plots)
# Mean detection rate = 70.0% sd = 29.3%

# Consistency of detection in pooled samples ----
eDNA_obs %>%
    filter(!is.na(subsample_id)) %>%
    group_by(plot_id) %>%
    group_split() %>%
    map(function(x) {
        pooled <- filter(x, subsample_id == "P")
        individual <- filter(x, subsample_id != "P") %>%
            distinct(short_name)
        n_detected <- sum(individual$short_name %in% pooled$short_name)
        return(data.frame(
            plot_id = unique(x$plot_id),
            prop_detected = n_detected / nrow(individual)
        ))
    }) %>%
    bind_rows() %>%
    summarise(mean_prop = mean(prop_detected), sd_prop = sd(prop_detected))
# On average 71% of the species are detected in the pooled samples (sd = 13%)
# Useful for streamlining the protocol.
