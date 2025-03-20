# eDNA detection rates for vascular plants

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)

# Set colours
col_edna <- "#156082"
col_bot <- "#00CC66"

# Load presence/ absence botanical observations and prepare for comaprison
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    )

# Load species list from botanical observaitons
botanical_ref <- read_csv("data/botanical_obs/species_list_cleaned.csv") %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    )

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
    bind_rows()

# Load reference databases
its2_ref <- read_csv("data/sequences/reference/its2_reference_Assmann_et_al_2025.csv") %>% mutate(sq = "its2")
trnl_ref <- read_csv("data/sequences/reference/trnl_reference_Assmann_et_al_2025.csv") %>% mutate(sq = "trnl")
ref_combined <- bind_rows(its2_ref, trnl_ref) %>%
    distinct(realised_taxon, sq) %>%
    dplyr::rename(short_name = realised_taxon) %>%
    mutate(
        genus = gsub("([A-Za-z]*) ([A-Za-z]*)", "\\1", short_name),
        species = gsub("([A-Za-z]*) ([A-Za-z]*)", "\\2", short_name)
    ) %>%
    group_by(short_name) %>%
    group_split() %>%
    map(function(x) {
        # Get unique markers
        sqs <- sort(unique(x$sq)) %>% paste(collapse = " ")
        # Slice tor row with most counts
        x <- slice(x, 1)
        # Update markers
        x$sq <- sqs
        return(x)
    }) %>%
    bind_rows() %>%
    relocate(short_name, genus, species)

## Theoretical detection rates ----

# Species-level detection
theo_detect <- botanical_ref %>%
    filter(species != "species") %>%
    inner_join(ref_combined, by = join_by(short_name)) %>%
    select(short_name, genus = genus.x, species = species.x, sq) %>%
    filter(!is.na(short_name)) %>%
    mutate(detect_level = "species")

# Check n species
nrow(theo_detect)
# In theory we can detect 48 species to the species level

# Genus level detection
theo_detect <- botanical_ref %>%
    filter(!(short_name %in% theo_detect$short_name)) %>%
    filter(genus %in% ref_combined$genus) %>%
    select(-full_name) %>%
    mutate(detect_level = "genus", sq = "its2 trnl") %>%
    bind_rows(., theo_detect)

theo_detect %>% filter(detect_level == "genus")
botanical_ref %>% filter(!(short_name %in% theo_detect$short_name))
# All Salix speices except for Salix arctica can only be detected at the
# genus level. Woodsia cannot be detected at all.
# For the puprose of the analyis it seems like no Salix arcitca
# Are specifically detected with eDNA so we just treat all
# willow species as one species

## Joined analysis of botanical and edna data ----

# Combined datasets
obs_combined <- bind_rows(
    eDNA_obs %>%
        filter(!is.na(subsample_id)) %>%
        select(plot_id, genus, species, short_name) %>%
        mutate(source = "eDNA"),
    botanical_obs %>%
        select(plot_id, genus, species, short_name) %>%
        mutate(source = "botanical")
) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct()

# Save file
write_csv(
    obs_combined,
    "data/analysis/combined_observations.csv"
)

# Frequencies of species detected with each method
species_frequencies <- obs_combined %>%
    filter(species != "species" | genus == "Salix") %>%
    group_by(short_name, source) %>%
    summarise(n_detected = n()) %>%
    arrange(desc(short_name))

# Add 0s for no detection
species_frequencies <- species_frequencies %>%
    group_by(short_name) %>%
    filter(n() == 1) %>%
    mutate(
        source = case_when(
            source == "botanical" ~ "eDNA",
            source == "eDNA" ~ "botanical"
        ),
        n_detected = 0
    ) %>%
    bind_rows(species_frequencies)

# re order factor for plotting
species_frequencies$short_name <- factor(species_frequencies$short_name,
    levels = sort(unique(species_frequencies$short_name), decreasing = TRUE)
)
species_frequencies$source <- gsub("botanical", "botanical survey", species_frequencies$source)

# Visualise
sp_detection_plot <- ggplot(
    species_frequencies,
    aes(
        x = short_name,
        y = n_detected,
        fill = source,
        group = source
    )
) +
    geom_col(position = "dodge") +
    scale_x_discrete() +
    scale_fill_manual(values = c(col_bot, col_edna)) +
    labs(
        x = "", y = "species detected (n plots)",
        fill = "data source:"
    ) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "top")
save_plot("figures/fig_s1_sp_detection.png",
    sp_detection_plot,
    base_height = 10,
    base_asp = 3 / 5,
    bg = "white"
)

# Plot level detection rate ----

# Species level
plot_level_detection <- obs_combined %>%
    # Remove all ids without species id except for salix
    filter(species != "species" | genus == "Salix") %>%
    # eDNA detection of Salix at the genus level
    group_by(plot_id, genus) %>%
    group_split() %>%
    map(function(x) {
        # If it's not a Salix species skip
        if (unique(x$genus) != "Salix") {
            return(x)
        }
        # If it's just one observaiton of a salix species skip
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
    }) %>%
    bind_rows() %>%
    group_by(plot_id, species) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(plot_id, short_name, source)

# Calculate mean detection rates
plot_level_detection %>%
    group_by(plot_id) %>%
    group_split() %>%
    map(function(x) {
        n_matches <- filter(x, source == "botanical eDNA") %>% nrow()
        return(data.frame(
            plot_id = unique(x$plot_id),
            bot_in_edna = n_matches / nrow(filter(x, grepl("botanical", source))),
            edna_in_bot = n_matches / nrow(filter(x, grepl("eDNA", source)))
        ))
    }) %>%
    bind_rows() %>%
    summarise(
        bot_in_edna_mean = mean(bot_in_edna),
        bot_in_edna_sd = sd(bot_in_edna),
        edna_in_bot_mean = mean(edna_in_bot),
        edna_in_bot_sd = sd(edna_in_bot)
    )
#   bot_in_edna_mean bot_in_edna_sd edna_in_bot_mean edna_in_bot_sd
# 1        0.7595775      0.1821415         0.544813      0.1726605
# Mean detection botany in eDNA: 76.0% (sd = 18.2%)
# Mean detection of eDNA in botany: 54.5% (sd = 17.3%)

# Summarise source groupings and visualise
plot_level_detection <- plot_level_detection %>%
    group_by(plot_id, source) %>%
    summarise(n = n()) %>%
    mutate(site = substr(plot_id, 1, 2))

# Update source levels for plotting
plot_level_detection$source <- gsub(" ", " & ", plot_level_detection$source)
plot_level_detection$source <- case_when(
    plot_level_detection$source == "botanical" ~ "botanical only",
    plot_level_detection$source == "eDNA" ~ "eDNA only",
    TRUE ~ plot_level_detection$source
)
plot_level_detection$source <- factor(plot_level_detection$source,
    levels = sort(unique(plot_level_detection$source), decreasing = TRUE)
)

plot_level_detect_plot <- ggplot(
    plot_level_detection,
    aes(x = plot_id, y = n, fill = source, group = source)
) +
    geom_col() +
    scale_y_continuous(limits = c(0, 20)) +
    scale_fill_manual(values = c(col_edna, col_bot, "gray25")) +
    labs(
        x = "", y = "species detected (n)",
        fill = "detected in:"
    ) +
    theme_cowplot() +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
save_plot("figures/fig_s2_sp_detect_per_plot.png",
    plot_level_detect_plot,
    base_height = 6,
    bg = "white"
)


# Genus level
plot_level_detection <- obs_combined %>%
    group_by(plot_id, genus) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(plot_id, genus, source)

# Calculate mean detection rates
plot_level_detection %>%
    group_by(plot_id) %>%
    group_split() %>%
    map(function(x) {
        n_matches <- filter(x, source == "botanical eDNA") %>% nrow()
        return(data.frame(
            plot_id = unique(x$plot_id),
            bot_in_edna = n_matches / nrow(filter(x, grepl("botanical", source))),
            edna_in_bot = n_matches / nrow(filter(x, grepl("eDNA", source)))
        ))
    }) %>%
    bind_rows() %>%
    summarise(
        bot_in_edna_mean = mean(bot_in_edna),
        bot_in_edna_sd = sd(bot_in_edna),
        edna_in_bot_mean = mean(edna_in_bot),
        edna_in_bot_sd = sd(edna_in_bot)
    )
#   bot_in_edna_mean bot_in_edna_sd edna_in_bot_mean edna_in_bot_sd
# 1        0.9226257      0.1083756        0.6073142      0.1894256
# Mean detection botany in eDNA: 92.2% (sd = 10.8%)
# Mean detection of eDNA in botany: 60.7% (sd = 18.9%)

# Summarise source groupings and visualise
plot_level_detection <- plot_level_detection %>%
    group_by(plot_id, source) %>%
    summarise(n = n()) %>%
    mutate(site = substr(plot_id, 1, 2))

# Update source levels for plotting
plot_level_detection$source <- gsub(" ", " & ", plot_level_detection$source)
plot_level_detection$source <- case_when(
    plot_level_detection$source == "botanical" ~ "botanical only",
    plot_level_detection$source == "eDNA" ~ "eDNA only",
    TRUE ~ plot_level_detection$source
)
plot_level_detection$source <- factor(plot_level_detection$source,
    levels = sort(unique(plot_level_detection$source), decreasing = TRUE)
)

plot_level_detect_plot <- ggplot(
    plot_level_detection,
    aes(x = plot_id, y = n, fill = source, group = source)
) +
    geom_col() +
    scale_y_continuous(limits = c(0, 15)) +
    scale_fill_manual(values = c(col_edna, col_bot, "gray25")) +
    labs(
        x = "", y = "genera detected (n)",
        fill = "detected in:"
    ) +
    theme_cowplot() +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
save_plot("figures/fig_s2_genus_detect_per_plot.png",
    plot_level_detect_plot,
    base_height = 6,
    bg = "white"
)

## Site level detection ----
# Species level
site_level_detection <- obs_combined %>%
    # Treating salix as one species removing all other species-level ids
    filter(species != "species" | genus == "Salix") %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    select(-plot_id) %>%
    # eDNA detection of Salix at the genus level
    group_by(site_id, genus) %>%
    group_split() %>%
    map(function(x) {
        # If it's not a Salix species skip
        if (unique(x$genus) != "Salix") {
            return(x)
        }
        # If it's just one observaiton of a salix species skip
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
    }) %>%
    bind_rows() %>%
    group_by(site_id, species) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(site_id, short_name, source) %>%
    arrange(site_id, short_name) 

# Print list for qc
site_level_detection %>%
group_by(site_id) %>%
group_split() %>% map(print,n = 100)

# Calculate mean detection rates
site_level_detection %>%
    group_by(site_id) %>%
    group_split() %>%
    map(function(x) {
        n_matches <- filter(x, source == "botanical eDNA") %>% nrow()
        return(data.frame(
            plot_id = unique(x$site_id),
            bot_in_edna = n_matches / nrow(filter(x, grepl("botanical", source))),
            edna_in_bot = n_matches / nrow(filter(x, grepl("eDNA", source)))
        ))
    }) %>%
    bind_rows() %>%
    summarise(
        bot_in_edna_mean = mean(bot_in_edna),
        bot_in_edna_sd = sd(bot_in_edna),
        edna_in_bot_mean = mean(edna_in_bot),
        edna_in_bot_sd = sd(edna_in_bot)
    )
# site levle detection
#   bot_in_edna_mean bot_in_edna_sd edna_in_bot_mean edna_in_bot_sd
# 1        0.7896563      0.0594771        0.8380561     0.06004709
# Mean detection botany in eDNA: 79.0% (sd = 5.9%)
# Mean detection of eDNA in botany: 83.8% (sd = 6.0%)

# Summarise source groupings and visualise
site_level_detection <- site_level_detection %>%
    group_by(site_id, source) %>%
    summarise(n = n())

# Update source levels for plotting
site_level_detection$source <- gsub(" ", " & ", site_level_detection$source)
site_level_detection$source <- case_when(
    site_level_detection$source == "botanical" ~ "botanical only",
    site_level_detection$source == "eDNA" ~ "eDNA only",
    TRUE ~ site_level_detection$source
)
site_level_detection$source <- factor(site_level_detection$source,
    levels = sort(unique(site_level_detection$source), decreasing = TRUE)
)

site_level_detect_plot <- ggplot(
    site_level_detection,
    aes(x = site_id, y = n, fill = source, group = source)
) +
    geom_col() +
    scale_y_continuous(limits = c(0, 40)) +
    scale_fill_manual(values = c(col_edna, col_bot, "gray25")) +
    labs(
        x = "", y = "species detected (n)",
        fill = "detected in:"
    ) +
    theme_cowplot() +
    theme(legend.position = "top")
save_plot("figures/fig_s3_sp_detect_per_site.png",
    site_level_detect_plot,
    base_height = 6,
    bg = "white"
)
