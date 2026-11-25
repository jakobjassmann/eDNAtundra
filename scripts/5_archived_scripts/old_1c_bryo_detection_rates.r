# eDNA detection rates for vascular plants

# Dependencies
library(tidyverse)
library(cowplot)
library(ggplot2)

# Set colours
col_edna <- "#F28F16"
col_bot <- "#A62121"

# Load botanical observations
botanical_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv") %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # convert al sp. to species
    mutate(species = (gsub("sp\\.", "species", species))) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) & !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species))

# Load eDNA data
eDNA_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_bryo_taxa_in_sample.csv") %>%
    mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    mutate(site_id = substr(plot_id, 1,2))

# Check theoretical genus level id
combined_ref_genus <- full_join(
    eDNA_obs %>%
        distinct(genus) %>%
        mutate(eDNA = TRUE),
    botanical_obs %>%
        distinct(genus) %>%
        mutate(botanical = TRUE)
) %>%
    mutate(source = case_when(
        botanical & is.na(eDNA) ~ "botanical",
        eDNA & is.na(botanical) ~ "eDNA",
        botanical & eDNA ~ "botanical & eDNA",
        TRUE ~ NA
    )) %>%
    select(-eDNA, -botanical) %>%
    distinct() %>%
    print(n = 100)
combined_ref_genus %>% group_by(source) %>% tally()
combined_ref_genus %>% filter(source == "botanical & eDNA")
# 18 genus were detected in both
# A tibble: 18 × 2
#    genus              source
#    <chr>              <chr>
#  1 Orthotrichum       botanical & eDNA
#  2 Encalypta          botanical & eDNA
#  3 Plagiomnium        botanical & eDNA
#  4 Distichium         botanical & eDNA
#  5 Bryoerythrophyllum botanical & eDNA
#  6 Ptychostomum       botanical & eDNA
#  7 Timmia             botanical & eDNA
#  8 Isopterygiopsis    botanical & eDNA
#  9 Didymodon          botanical & eDNA
# 10 Drepanocladus      botanical & eDNA
# 11 Orthothecium       botanical & eDNA
# 12 Tomentypnum        botanical & eDNA
# 13 Myurella           botanical & eDNA
# 14 Ditrichum          botanical & eDNA
# 15 Platydictya        botanical & eDNA
# 16 Hypnum             botanical & eDNA
# 17 Tortella           botanical & eDNA
# 18 Meesia             botanical & eDNA

# Check theoretical species level id
combined_ref_species <- full_join(
    eDNA_obs %>%
        distinct(genus, species, short_name) %>%
        mutate(eDNA = TRUE),
    botanical_obs %>%
        distinct(genus, species, short_name) %>%
        mutate(botanical = TRUE)
) %>%
    mutate(source = case_when(
        botanical & is.na(eDNA) ~ "botanical",
        eDNA & is.na(botanical) ~ "eDNA",
        botanical & eDNA ~ "botanical & eDNA",
        TRUE ~ NA
    )) %>%
    select(-eDNA, -botanical) %>%
    distinct() %>%
    print(n = 100)

# this looks  a bit messy, what about the matches
combined_ref_species %>% filter(source == "botanical & eDNA")
# 17 species detecteable at species level
#    <chr>              <chr>            <chr>                             <chr>
#  1 Distichium         capillaceum      Distichium capillaceum            botani…
#  2 Bryoerythrophyllum recurvirostrum   Bryoerythrophyllum recurvirostrum botani…
#  3 Encalypta          procera          Encalypta procera                 botani…
#  4 Ptychostomum       pseudotriquetrum Ptychostomum pseudotriquetrum     botani…
#  5 Distichium         inclinatum       Distichium inclinatum             botani…
#  6 Timmia             norvegica        Timmia norvegica                  botani…
#  7 Isopterygiopsis    pulchella        Isopterygiopsis pulchella         botani…
#  8 Orthothecium       chryseum         Orthothecium chryseum             botani…
#  9 Tomentypnum        nitens           Tomentypnum nitens                botani…
# 10 Myurella           julacea          Myurella julacea                  botani…
# 11 Ditrichum          flexicaule       Ditrichum flexicaule              botani…
# 12 Encalypta          rhaptocarpa      Encalypta rhaptocarpa             botani…
# 13 Ptychostomum       species          Ptychostomum species              botani…
# 14 Tortella           fragilis         Tortella fragilis                 botani…
# 15 Meesia             uliginosa        Meesia uliginosa                  botani…
# 16 Hypnum             bambergeri       Hypnum bambergeri                 botani…
# 17 Meesia             triquetra        Meesia triquetra                  botani…

# Combine reference list including:
#   botanical observations,
#   eDNA observations if the genus was observed (anywhere) in the botanical obs
#   eDNA observations if the species was observed (anywhere) in the botanical obs
combined_ref <- bind_rows(
    # Genuses in botanical or both
    combined_ref_genus %>%
        filter(source %in% c("botanical & eDNA", "botanical")),
    # Species in botanical or both
    combined_ref_species %>% filter(source %in% c("botanical & eDNA", "botanical"))
) %>%
    mutate(species = case_when(is.na(species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct() %>%
    print(n = 100)

# Confirm resolution for genus
distinct(combined_ref, genus) %>% tally()
distinct(botanical_obs, genus) %>% tally()
combined_ref %>%
    filter(source == "botanical & eDNA") %>%
    distinct(genus) %>%
    tally()
# 18 / 47 genusses
combined_ref %>%
    filter(species != "species") %>%
    distinct(genus, species) %>%
    tally()
botanical_obs %>%
    filter(species != "species") %>%
    distinct(genus, species) %>%
    tally()
combined_ref %>%
    filter(species != "species") %>%
    filter(source == "botanical & eDNA") %>%
    distinct(genus, species) %>%
    tally()
# 16 species / 68 species id to species level vs botanical obs

# Filter datasets to subsets containing only the detectable taxa
# botanical_obs_sub <- botanical_obs %>%
#     filter(short_name %in% combined_ref$short_name |
#         genus %in% combined_ref$genus) %>%
#     select(plot_id, genus, species) %>%
#     mutate(source = "botanical") %>%
#     mutate(species = case_when(!(species %in% combined_ref$species) ~ "species", TRUE ~ species)) %>%
#     mutate(short_name = paste(genus, species)) %>%
#     distinct()

# Filter eDNA obs to subsets containing only taxa observed in botanical obs
eDNA_obs_sub <- eDNA_obs %>%
    filter(short_name %in% combined_ref$short_name |
        genus %in% combined_ref$genus) %>%
    # Remove extraction controls
    filter(!is.na(subsample_id)) %>%
    select(site_id, plot_id, genus, species) %>%
    mutate(source = "eDNA") %>%
    # If detection to species level was not possible set to "species"
    mutate(species = case_when(!(species %in% combined_ref$species) ~ "species", TRUE ~ species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct()

# combine botanical obs and eDNA detection
obs_combined <- bind_rows(
    # Botanical observations species level id
    mutate(botanical_obs, source = "botanical"),
    # Botanical observations genus level id
    botanical_obs %>% 
        distinct(site_id, plot_id, genus) %>%
        mutate(species = "species") %>%
        mutate(short_name = paste0(genus, " ", species),
                source = "botanical"),
    # eDNA detection
    eDNA_obs_sub) %>%
    arrange(site_id, plot_id, short_name, source) %>%
    distinct() %>%
    mutate(source = factor(source, levels = c("botanical", "eDNA", "botanical & eDNA"))) %>%
    mutate(site_id = case_when(site_id == "S1" ~ "Site 1",
    site_id == "S2" ~ "Site 2",
    site_id == "S3" ~ "Site 3",
    TRUE ~ NA))

# Frequencies of species / genus detected with each method (genus level)
species_frequencies <- obs_combined %>%
    # filter(species == "species") %>%
    group_by(short_name, source) %>%
    summarise(n_detected = n()) %>%
    arrange(desc(short_name)) %>%
    print(n = 100)

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
# species_frequencies$source <- gsub("botanical", "botanical survey", species_frequencies$source)

# Visualise Genus detection
(genus_detection_plot <- ggplot(
    filter(species_frequencies, grepl("species", short_name)),
    aes(
        x = gsub(" species", "", short_name),
        y = n_detected,
        fill = source,
        group = source
    ),
    colour = "#black"
) +
    geom_col(position = "dodge") +
    scale_x_discrete() +
    scale_y_continuous(limits = c(0, 60)) +
    scale_fill_manual(values = c(col_bot, col_edna)) +
    labs(
        x = "", y = "n plots genus detected",
        fill = "Data source"
    ) +
    coord_flip() +
    theme_cowplot() +
       theme(legend.position = "inside", 
    legend.position.inside = c(0.7, 0.06)))
save_plot("figures/fig_s_bryo_genus_detection.png",
    genus_detection_plot,
    base_height = 10,
    base_asp = 3 / 5,
    bg = "white"
)

# Visualise species detection
(sp_detection_plot <- ggplot(
    filter(species_frequencies, !grepl("species", short_name)),
    aes(
        x = short_name,
        y = n_detected,
        fill = source,
        group = source
    )
) +
    geom_col(position = "dodge") +
    scale_x_discrete() +
    scale_y_continuous(limits = c(0, 40)) +
    scale_fill_manual(values = c(col_bot, col_edna)) +
    labs(
        x = "", y = "n plots species detected",
        fill = "Data source"
    ) +
    coord_flip() +
    theme_cowplot() +
    theme(legend.position = "inside", 
    legend.position.inside = c(0.55, 0.06)))
save_plot("figures/fig_s_bryo_sp_detection.png",
    sp_detection_plot,
    base_height = 11,
    base_asp = 3 / 6,
    bg = "white"
)

# # Collapse species list where detected
# obs_combined <- obs_combined %>%
#     group_by(plot_id, species) %>%
#     group_split() %>%
#     map(function(x) {
#         srcs <- unique(x$source)
#         slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
#     }) %>%
#     bind_rows() %>%
#     distinct(plot_id, short_name, source)

# Set colour palette (botanical, both)
pal <- c("#A62121", "#223F59")

# Number of genera per plot
bryo_detect_plot <- obs_combined %>%
    filter(species == "species") %>%
    mutate(source = as.character(source)) %>%
    group_by(plot_id, short_name) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows()  %>%
    mutate(source = factor(source, levels = c("botanical", "eDNA", "botanical eDNA")))

filter(obs_combined, plot_id == "S1R1", species == "species")
filter(bryo_detect_plot, plot_id == "S1R1")


# Genus per plot
(genus_plot <- ggplot(obs_combined ) +
    geom_col(aes(
        x = plot_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Bryophyte genus richness (n)",
        x = "Plot ID"
    ) +
    scale_colour_manual(
        values = pal, name = "Observation",
        labels = c("botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observation",
        labels = c("botanical only",  "botanical and eDNA")
    ) +
    scale_y_continuous(limits = c(0, 20)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)))

# Site level detection of species
site_level_detection <- obs_combined %>%
    # Remove all ids not to the species level
    filter(species != "species") %>%
    # Keep only speces detected in the botanical observations
    filter(short_name %in% unique(botanical_obs$short_name)) %>%
    select(-plot_id) %>%
    group_by(site_id, short_name) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(site_id, short_name, source) %>%
    mutate(source = factor(source, levels = rev(c("botanical eDNA", "botanical", "eDNA")))) %>%
    arrange(site_id, short_name, source)

site_level_detection %>% filter(source == "eDNA")

# Summarise source groupings and visualise
(bryo_site_richness <- site_level_detection %>%
    group_by(site_id, source) %>%
    summarise(n = n()))

# Set colors 
pal <- c("#F28F16", "#dba6a6", "#223F59") # "#A62121"

(bryo_site_richness_plot <- ggplot(bryo_site_richness) +
    geom_col(aes(
        x = site_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Number of Species (n)"
    ) +
    scale_colour_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only",  "botanical and eDNA")
    )+
    scale_y_continuous(limits = c(0, 60)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "right", 
        #legend.position.inside = c(0.06, 0.9),
        legend.justification.right = c(0.5, 0.105),
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(0.25,0.25,0.25,0.25, unit = "cm")))

# Add icons
bryo_site_richness_plot <- ggdraw() + 
    draw_plot(bryo_site_richness_plot) +
    # Vascular plants
    draw_image("figures/moss_bw.png", scale = 0.2, x = 0.225, y = 0.2) +
    # DNA 
    draw_image("figures/dna.png", scale = 0.35, x = 0.4, y = 0.19) 

save_plot("figures/figure_2b.png",
    bryo_site_richness_plot,
    #plot_grid(richness_detection, site_richness_plot),
    #base_asp = 2 * 1.6,
    bg = "white"
)


# Calculate mean detection
bryo_site_richness %>%
split(.$site_id) %>%
map(function(x){
    x <- filter(x, source != "eDNA") 
    filter(x, source == "botanical eDNA") %>% pull(n) %>% {./sum(x$n)}
}) %>% unlist() %>% mean()

# Site level detection of genus
genus_site_level_detection <- obs_combined %>%
    select(-plot_id) %>%
    distinct(site_id, genus, source) %>%
    # Keep only genera detected in the botanical observations
    filter(genus %in% unique(botanical_obs$genus)) %>%
    group_by(site_id, genus) %>%
    group_split() %>%
    map(function(x) {
        srcs <- unique(x$source)
        slice(x, 1) %>% mutate(source = paste(sort(srcs), collapse = " "))
    }) %>%
    bind_rows() %>%
    distinct(site_id, genus, source) %>%
    mutate(source = factor(source, levels = rev(c("botanical eDNA", "botanical", "eDNA")))) %>%
    arrange(site_id, genus, source)
# Summarise source groupings and visualise
(bryo_site_genus_richness <- genus_site_level_detection %>%
    group_by(site_id, source) %>%
    summarise(n = n()) )
 
 genus_site_level_detection %>% arrange(site_id, source, genus) %>%
 print(n = 120)

(bryo_site_richness_plot <- ggplot(bryo_site_genus_richness) +
    geom_col(aes(
        x = site_id,
        y = n,
        fill = source, 
        colour = source
    )) +
    labs(
        y = "Number of Genera (n)"
    ) +
    scale_colour_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    scale_fill_manual(
        values = pal, name = "Observed in",
        labels = c("eDNA only", "botanical only",  "botanical and eDNA")
    )+
    scale_y_continuous(limits = c(0, 50)) +
    theme_cowplot() +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "right", 
        #legend.position.inside = c(0.06, 0.9),
        legend.justification.right = c(0.5, 0.105),
        legend.box.background = element_rect(colour = "black"),
        legend.box.margin = margin(0.25,0.25,0.25,0.25, unit = "cm")))

# Calculate mean detection
bryo_site_genus_richness %>%
split(.$site_id) %>%
map(function(x){
    x <- filter(x, source != "eDNA") 
    filter(x, source == "botanical eDNA") %>% pull(n) %>% {./sum(x$n)}
}) %>% unlist() %>% mean()
