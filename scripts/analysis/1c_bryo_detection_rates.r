# eDNA detection rates for vascular plants

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)

# Set colours
col_edna <- "#156082"
col_bot <- "#00CC66"

# Load botanical observations
botanical_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv") %>%
mutate(plot_id = paste0(Site, Plot_ID)) %>%
select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
# convert al sp. to species
mutate(species = (gsub("sp\\.","species", species))) %>%
# Remove all uncertain ids
filter(!grepl("cf", genus) & !grepl("cf", species)) %>%
mutate(short_name = paste(genus, species)) 

# Load eDNA data
eDNA_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_bryo_taxa_in_sample.csv") %>%
mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species)) %>%
mutate(short_name = paste(genus, species))

# Check theoretical genus level id
combined_ref_genus <- full_join(eDNA_obs %>% 
distinct(genus) %>%
mutate(botanical = TRUE), 
botanical_obs %>% 
distinct(genus) %>%
mutate(eDNA = TRUE)) %>%
mutate(source = case_when(botanical & is.na(eDNA) ~ "botanical",
eDNA & is.na(botanical) ~ "eDNA",
botanical & eDNA ~ "botanical & eDNA",
TRUE ~ NA)) %>%
select(-eDNA, -botanical) %>%
distinct() %>% print(n = 100)
combined_ref_genus %>% filter(source == "botanical & eDNA")

# Check theoretical species level id
combined_ref_species <- full_join(eDNA_obs %>% 
distinct(genus, species, short_name) %>%
mutate(botanical = TRUE), 
botanical_obs %>% 
distinct(genus, species, short_name) %>%
mutate(eDNA = TRUE)) %>%
mutate(source = case_when(botanical & is.na(eDNA) ~ "botanical",
eDNA & is.na(botanical) ~ "eDNA",
botanical & eDNA ~ "botanical & eDNA",
TRUE ~ NA)) %>%
select(-eDNA, -botanical) %>%
distinct() %>% print(n = 100)

# this looks  a bit messy, what about the matches
combined_ref_species %>% filter(source == "botanical & eDNA")
combined_ref_species %>% filter(source == "botanical & eDNA") %>% distinct(genus)

# With a bit of luck we might combine the data from the matching species and 
# genuses to get a somehwat meaningful analysis
combined_ref <- bind_rows(combined_ref_genus %>% 
filter(source == "botanical & eDNA"),
combined_ref_species %>% filter(source == "botanical & eDNA")) %>%
mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species))  %>%
mutate(short_name = paste(genus, species))

# Resolution of the combined reference
distinct(combined_ref, genus) %>% tally()
distinct(botanical_obs, genus) %>% tally()
# 18 / 47 genusses
combined_ref %>% filter(species != "species") %>% distinct(species) %>% tally()
distinct(botanical_obs, genus, species) %>% tally()
# 16 species / 73 species

# Filter datasets to subsets containing only the detectable taxa
botanical_obs_sub <- botanical_obs %>% 
    filter(short_name %in% combined_ref$short_name |
    genus %in% combined_ref$genus) %>%
    select(plot_id, genus, species) %>%
    mutate(source = "botanical") %>%
    mutate(species = case_when(!(species %in% combined_ref$species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct()
eDNA_obs_sub <- eDNA_obs %>% 
    filter(short_name %in% combined_ref$short_name |
    genus %in% combined_ref$genus) %>%
    filter(!is.na(subsample_id)) %>%
    select(plot_id, genus, species) %>%
    mutate(source = "eDNA") %>%
    mutate(species = case_when(!(species %in% combined_ref$species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct()

# combine and collapse where detected
obs_combined <- bind_rows(botanical_obs_sub, eDNA_obs_sub) 

# Frequencies of species detected with each method
species_frequencies <- obs_combined %>%
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
save_plot("figures/fig_s5_bryo_sp_detection.png",
    sp_detection_plot,
    base_height = 10,
    base_asp = 3 / 5,
    bg = "white"
)

# Collapse species list where detected
obs_combined <- obs_combined %>%
    group_by(plot_id, species) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(plot_id, short_name, source)

# Detection per site
obs_combined %>%
filter(!grepl("species", short_name)) %>%
mutate(site_id = substr(plot_id, 1,2),
site_id = gsub("S", "Site ", site_id)) %>%
group_by(site_id) %>%
summarise(n = n_distinct(short_name))
