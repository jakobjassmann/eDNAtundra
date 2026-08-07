# eDNA detection for bryophytes
# Jakob J. Assmann jakob.assmann@gu.se 5 August 2026

# Dependencies
library(tidyverse)
library(gt)

#### Load and prepare botanical and eDNA datasets ----

# Load presence/ absence botanical observations and prepare for comaprison
botanical_obs_orignal <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb_updated_Aug_2026.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(plot_id, genus = Genus, species = Species, synonym) %>%
    # convert al sp. to species
    mutate(
        species = (gsub("sp\\.", "species", species)),
        synonym = (gsub("sp\\.", "species", synonym))
    ) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    # Remove any duplicates (n = 2) from subspecies-ids
    distinct()

# Save all name and synonym combinations for later
botanical_obs_synonyms <- botanical_obs_orignal %>%
    distinct(genus, species, short_name, synonym) %>%
    arrange(short_name) %>%
    filter(species != "species")

# Override original species names with synonyms present in DS-BITSR
botanical_obs <- botanical_obs_orignal %>%
    mutate(
        short_name = synonym,
        genus = gsub("([A-Za-z]*) ([A-Za-z]*)", "\\1", synonym),
        species = gsub("([A-Za-z]*) ([A-Za-z]*)", "\\2", synonym)
    ) %>%
    select(-synonym)

# Load presence / absence from eDNA
eDNA_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_bryo_taxa.csv") %>%
    mutate(species = if_else(is.na(species), "species", species)) %>%
    mutate(short_name = paste(genus, species), eDNA = TRUE) %>%
    select(-asv_id, -count, -sq) %>%
    distinct()

# Tally unqiue detections
nrow(botanical_obs)
nrow(filter(botanical_obs, species != "species"))
nrow(eDNA_obs)
# botanical 516 vs edna 120 (23%)
# These are very few detections, how many are at species level?
nrow(filter(eDNA_obs, species != "species"))
# even less 77 (15%)

# How many species were successfuly detected at least once?
eDNA_obs %>%
    filter(species != "species") %>%
    distinct(genus, species) %>%
    arrange(genus, species)
# 17
# How many of these are in our botanical observations?
eDNA_obs %>%
    filter(species != "species") %>%
    distinct(short_name) %>%
    arrange(short_name) %>%
    filter(short_name %in% botanical_obs_synonyms$synonym)
# 12

# Which are not?
eDNA_obs %>%
    filter(species != "species") %>%
    distinct(short_name) %>%
    arrange(short_name) %>%
    filter(!(short_name %in% botanical_obs_synonyms$synonym))
# # A tibble: 7 × 1
#   short_name
#   <chr>
# 1 Campylium bambergeri
# 2 Distichium hagenii
# 3 Encalypta longicollis
# 4 Flexitrichum flexicaule
# 5 Hymenostylium recurvirostrum
# 6 Orthothecium intricatum
# 7 Timmiella anomala

#### Join botanical and edna data ----

# Combine datasets and add species level detection
obs_combined <- full_join(
    eDNA_obs %>%
        select(plot_id, genus, species, short_name, eDNA) %>%
        # Summarise at plot level
        group_by(plot_id, genus, species, short_name) %>%
        summarise(eDNA = any(eDNA)) %>%
        ungroup(),
    botanical_obs %>%
        select(plot_id, genus, species, short_name) %>%
        filter(species != "species") %>%
        mutate(botanical = TRUE)
) %>%
    mutate(short_name = paste(genus, species))

# Replace detection NAs with FALSE
obs_combined <- obs_combined %>%
    mutate(
        eDNA = if_else(is.na(eDNA), FALSE, eDNA),
        botanical = if_else(is.na(botanical), FALSE, botanical)
    )

# Add genus detection
obs_combined <- obs_combined %>%
    group_by(plot_id) %>%
    group_split() %>%
    map(function(x) {
        # Get genera with eDNA level detection
        genus_detected <- select(x, genus, eDNA) %>%
            filter(eDNA) %>%
            distinct() %>%
            group_by(genus) %>%
            summarise(
                eDNA_genus = any(eDNA)
            )

        # Update detection level
        x <- full_join(x, genus_detected) %>%
            filter(species != "species")

        return(x)
    }) %>%
    bind_rows() %>%
    mutate(
        eDNA_genus = if_else(is.na(eDNA_genus), FALSE, eDNA_genus)
    )

# revert synonyms and fill columns for species not in botanical observations
obs_combined <- obs_combined %>%
    mutate(synonym = short_name) %>%
    select(-genus, -species, -short_name) %>%
    full_join(botanical_obs_synonyms) %>%
    mutate(
        short_name = if_else(is.na(short_name), synonym, short_name),
        genus = if_else(
            is.na(genus),
            gsub("([A-Za-z]*) ([A-Za-z]*)", "\\1", synonym),
            genus
        ),
        species = if_else(
            is.na(species),
            gsub("([A-Za-z]*) ([A-Za-z]*)", "\\2", synonym),
            species
        )
    )


# Arrange columns
obs_combined <- obs_combined %>%
    relocate(
        plot_id,
        genus,
        species,
        short_name,
        synonym,
        botanical,
        eDNA,
        eDNA_genus
    )

# Save file
write_csv(
    obs_combined,
    "data/analysis/bryo_detection_species_genus.csv"
)
# obs_combined <- read_csv("data/analysis/bryo_detection_species_genus.csv")

# For QC with Erin, pull out those plots where there is at least one species
# that was detected only through DS-BITSR
obs_combined_new <-
    filter(
        obs_combined,
        plot_id %in%
            (obs_combined %>%
                group_by(plot_id) %>%
                summarise(
                    new = any(
                        short_name %in%
                            (eDNA_obs %>%
                                filter(species != "species") %>%
                                distinct(short_name) %>%
                                arrange(short_name) %>%
                                filter(
                                    !(short_name %in%
                                        botanical_obs_synonyms$short_name)
                                ) %>%
                                pull(short_name))
                    )
                ) %>%
                filter(new) %>%
                pull(plot_id))
    ) %>%
    mutate(
        new_species = short_name %in%
            (eDNA_obs %>%
                filter(species != "species") %>%
                distinct(short_name) %>%
                arrange(short_name) %>%
                filter(
                    !(short_name %in%
                        botanical_obs_synonyms$short_name)
                ) %>%
                pull(short_name))
    )
write_csv(
    obs_combined_new,
    "data/analysis/bryo_detection_plots_with_new_species.csv"
)

#### Summary statistics ----

# Taxa detection across the different methods
obs_combined %>%
    group_by(short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    relocate(short_name, botanical) %>%
    bind_rows(
        .,
        mutate(
            summarise(., across(botanical:eDNA_genus, sum)),
            short_name = "Total"
        )
    ) %>%
    mutate(across(botanical:eDNA_genus, function(x) {
        case_when(x == 1 ~ "X", x == 0 ~ "O", TRUE ~ as.character(x))
    })) %>%
    gt() %>%
    cols_label(
        short_name = "species",
        botanical = html("botanical<br>survey"),
        eDNA = html("eDNA at<br>species-level"),
        eDNA_genus = html("eDNA at<br>genus-level"),
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            eDNA,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/bryo_1_detection_across_methods.html")

## Detections rates at the site level
obs_combined %>%
    # add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    relocate(site_id, short_name) %>%
    select(-plot_id, -genus, -species) %>%
    # aggregate detections by site
    group_by(site_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    # keep only species detected in botanical observations (reference)
    filter(botanical) %>%
    summarise(
        across(botanical:eDNA_genus, sum)
    ) %>%
    ungroup() %>%
    bind_rows(
        .,
        mutate(
            summarize(., round(across(botanical:eDNA_genus, mean), 1)),
            site_id = "Mean"
        )
    ) %>%
    mutate(
        across(eDNA:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, drop0trailing = TRUE),
                " (",
                format(
                    round(100 * x / botanical),
                    width = 2
                ),
                "%)"
            )
        }),
        botanical = format(botanical, width = 2, drop0trailing = TRUE),
    ) %>%
    gt() %>%
    cols_label(
        site_id = "site",
        botanical = html("botanical<br>survey"),
        eDNA = html("eDNA at<br>species-level"),
        eDNA_genus = html("eDNA at<br>genus-level"),
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            eDNA,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/bryo_2_detections_per_site.html")

## Detections rates at the plot level
obs_combined %>%
    select(-genus, -species) %>%
    # aggregate detections by site
    group_by(plot_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    # keep only species detected in botanical observations (reference)
    filter(botanical) %>%
    summarise(
        across(botanical:eDNA_genus, sum)
    ) %>%
    ungroup() %>%
    bind_rows(
        .,
        mutate(
            summarize(., round(across(botanical:eDNA_genus, mean), 1)),
            plot_id = "Mean"
        )
    ) %>%
    mutate(
        across(eDNA:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, drop0trailing = TRUE),
                " (",
                format(
                    round(100 * x / botanical),
                    width = 3
                ),
                "%)"
            )
        }),
        botanical = format(botanical, width = 2, drop0trailing = TRUE),
    ) %>%
    gt() %>%
    cols_label(
        plot_id = "Plot",
        botanical = html("botanical<br>survey"),
        eDNA = html("eDNA at<br>species-level"),
        eDNA_genus = html("eDNA at<br>genus-level"),
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            eDNA,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/bryo_3_detections_per_plot.html")

## Sites per species by method
obs_combined %>%
    # add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    relocate(site_id, short_name) %>%
    select(-plot_id, -genus, -species) %>%
    # aggregate detections by site
    group_by(site_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    group_by(short_name) %>%
    arrange(short_name) %>%
    summarise(across(botanical:eDNA_genus, sum)) %>%
    mutate(
        across(eDNA:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0)
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    cols_label(
        short_name = "species",
        botanical = html("botanical<br>survey"),
        eDNA = html("eDNA at<br>species-level"),
        eDNA_genus = html("eDNA at<br>genus-level"),
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            eDNA,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/bryo_4_sites_per_species.html")

## Plot per species by method
obs_combined %>%
    select(-genus, -species) %>%
    group_by(short_name) %>%
    arrange(short_name) %>%
    summarise(across(botanical:eDNA_genus, sum)) %>%
    mutate(
        across(eDNA:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0)
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    cols_label(
        short_name = "species",
        botanical = html("botanical<br>survey"),
        eDNA = html("eDNA at<br>species-level"),
        eDNA_genus = html("eDNA at<br>genus-level"),
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            eDNA,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/bryo_5_plots_per_species.html")

## End of file
