# eDNA detection rates for vascular plants

# Dependencies
library(tidyverse)
library(gt)

#### Load and prepare botanical and eDNA datasets ----

# Load presence/ absence botanical observations and prepare for comaprison
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    ) %>%
    # Remove any duplicated rows
    distinct()

# Load species list from botanical observaitons
botanical_ref <- read_csv("data/botanical_obs/species_list_cleaned.csv") %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    )

# Load presence / absence from eDNA
trnl_obs <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa.csv") %>%
    mutate(species = if_else(is.na(species), "species", species)) %>%
    mutate(short_name = paste(genus, species), trnl = TRUE) %>%
    select(-asv_id, -count, -sq) %>%
    distinct()
its2_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa.csv") %>%
    mutate(species = if_else(is.na(species), "species", species)) %>%
    mutate(short_name = paste(genus, species), its2 = TRUE) %>%
    select(-asv_id, -count, -sq) %>%
    distinct()
eDNA_obs <- full_join(trnl_obs, its2_obs) %>%
    mutate(eDNA = TRUE)

# Tally unqiue detection by each marker
eDNA_obs %>%
    summarise(
        trnl_only = sum(is.na(its2)),
        its2_only = sum(is.na(trnl)),
        both = sum(!is.na(trnl) & !is.na(its2))
    )
# # A tibble: 1 × 3
#   trnl_only its2_only  both
#       <int>     <int> <int>
# 1       704       670   968

#### Join botanical and edna data ----

# Combine datasets and add species level detection
obs_combined <- full_join(
    eDNA_obs %>%
        select(plot_id, genus, species, short_name, trnl, its2, eDNA) %>%
        # Summarise at plot level
        group_by(plot_id, genus, species, short_name) %>%
        summarise(trnl = any(trnl), its2 = any(its2), eDNA = any(eDNA)) %>%
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
        trnl = if_else(is.na(trnl), FALSE, trnl),
        its2 = if_else(is.na(its2), FALSE, its2),
        eDNA = if_else(is.na(eDNA), FALSE, eDNA),
        botanical = if_else(is.na(botanical), FALSE, botanical)
    )

# Add genus detection
obs_combined <- obs_combined %>%
    group_by(plot_id) %>%
    group_split() %>%
    map(function(x) {
        # Get genera with eDNA level detection
        genus_detected <- select(x, genus, trnl, its2, eDNA) %>%
            filter(eDNA) %>%
            distinct() %>%
            group_by(genus) %>%
            summarise(
                eDNA_genus_trnl = any(trnl),
                eDNA_genus_its2 = any(its2),
                eDNA_genus = any(its2 | trnl)
            )

        # Update detection level
        x <- full_join(x, genus_detected) %>%
            filter(species != "species")

        return(x)
    }) %>%
    bind_rows() %>%
    mutate(
        eDNA_genus_trnl = if_else(
            is.na(eDNA_genus_trnl),
            FALSE,
            eDNA_genus_trnl
        ),
        eDNA_genus_its2 = if_else(
            is.na(eDNA_genus_its2),
            FALSE,
            eDNA_genus_its2
        ),
        eDNA_genus = if_else(is.na(eDNA_genus), FALSE, eDNA_genus)
    )

# Add harmonised detection for confused taxa
obs_combined <- obs_combined %>%
    mutate(
        eDNA_harm_trnl = case_when(
            !trnl &
                eDNA_genus_trnl &
                (short_name %in%
                    c(
                        "Carex borealipolaris",
                        "Carex myosuroides",
                        "Carex capillaris",
                        "Carex scirpoidea",
                        "Carex membranacea",
                        "Carex fuliginosa",
                        "Oxytropis arctobia",
                        "Oxytropis maydelliana",
                        "Salix arctica",
                        "Salix polaris",
                        "Salix reticulata",
                        "Salix richardsonii"
                    )) ~ TRUE,
            TRUE ~ trnl
        ),
        eDNA_harm_its2 = case_when(
            !its2 &
                eDNA_genus_its2 &
                (short_name %in%
                    c(
                        "Carex aquatilis",
                        "Carex bigelowii",
                        "Eriophorum angustifolium",
                        "Eriophorum triste",
                        "Salix reticulata",
                        "Salix richardsonii"
                    )) ~ TRUE,
            TRUE ~ its2
        )
    ) %>%
    mutate(eDNA_harm = eDNA_harm_trnl | eDNA_harm_its2)

# Arrange columns
obs_combined <- obs_combined %>%
    relocate(
        plot_id,
        genus,
        species,
        short_name,
        botanical,
        trnl,
        its2,
        eDNA,
        eDNA_harm_trnl,
        eDNA_harm_its2,
        eDNA_harm,
        eDNA_genus_trnl,
        eDNA_genus_its2,
        eDNA_genus
    )

# Save file
write_csv(
    obs_combined,
    "data/analysis/detection_species_genus.csv"
)
# obs_combined <- read_csv("data/analysis/detection_species_genus.csv")

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
    tab_spanner(label = "species-level", columns = c(trnl, its2, eDNA)) %>%
    tab_spanner(
        label = "harmonised",
        columns = c(eDNA_harm_trnl, eDNA_harm_its2, eDNA_harm)
    ) %>%
    tab_spanner(
        label = "genus-level",
        columns = c(eDNA_genus_trnl, eDNA_genus_its2, eDNA_genus)
    ) %>%
    cols_label(
        short_name = "Species",
        botanical = html("botanical<br>survery"),
        trnl = "trnl",
        its2 = "ITS2",
        eDNA = "both",
        eDNA_harm_trnl = "trnl",
        eDNA_harm_its2 = "ITS2",
        eDNA_harm = "both",
        eDNA_genus_trnl = "trnl",
        eDNA_genus_its2 = "ITS2",
        eDNA_genus = "both",
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            trnl,
            its2,
            eDNA,
            eDNA_harm_trnl,
            eDNA_harm_its2,
            eDNA_harm,
            eDNA_genus_trnl,
            eDNA_genus_its2,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/1_Detection_across_methods.html")

## Detections rates at the site level
obs_combined %>%
    # add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    relocate(site_id, short_name) %>%
    select(-plot_id, -genus, -species) %>%
    # keep only species detected in botanical observations (reference)
    filter(botanical) %>%
    # aggregate detections by site
    group_by(site_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    summarise(
        across(botanical:eDNA_genus, sum)
    ) %>%
    ungroup() %>%
    bind_rows(
        .,
        mutate(
            summarize(., across(botanical:eDNA_genus, mean)),
            site_id = "Mean"
        )
    ) %>%
    mutate(
        across(trnl:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0),
                " (",
                format(
                    round(100 * x / botanical),
                    width = 2
                ),
                "%)"
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    tab_spanner(label = "species-level", columns = c(trnl, its2, eDNA)) %>%
    tab_spanner(
        label = "harmonised",
        columns = c(eDNA_harm_trnl, eDNA_harm_its2, eDNA_harm)
    ) %>%
    tab_spanner(
        label = "genus-level",
        columns = c(eDNA_genus_trnl, eDNA_genus_its2, eDNA_genus)
    ) %>%
    cols_label(
        site_id = "Site",
        botanical = html("botanical<br>survery"),
        trnl = "trnl",
        its2 = "ITS2",
        eDNA = "both",
        eDNA_harm_trnl = "trnl",
        eDNA_harm_its2 = "ITS2",
        eDNA_harm = "both",
        eDNA_genus_trnl = "trnl",
        eDNA_genus_its2 = "ITS2",
        eDNA_genus = "both",
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            trnl,
            its2,
            eDNA,
            eDNA_harm_trnl,
            eDNA_harm_its2,
            eDNA_harm,
            eDNA_genus_trnl,
            eDNA_genus_its2,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/2_detections_per_site.html")

## Detections rates at the plot level
obs_combined %>%
    select(-genus, -species) %>%
    # keep only species detected in botanical observations (reference)
    filter(botanical) %>%
    # aggregate detections by site
    group_by(plot_id, short_name) %>%
    summarise(across(botanical:eDNA_genus, any)) %>%
    summarise(
        across(botanical:eDNA_genus, sum)
    ) %>%
    ungroup() %>%
    bind_rows(
        .,
        mutate(
            summarize(., across(botanical:eDNA_genus, mean)),
            plot_id = "Mean"
        )
    ) %>%
    mutate(
        across(trnl:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0),
                " (",
                format(
                    round(100 * x / botanical),
                    width = 2
                ),
                "%)"
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    tab_spanner(label = "species-level", columns = c(trnl, its2, eDNA)) %>%
    tab_spanner(
        label = "harmonised",
        columns = c(eDNA_harm_trnl, eDNA_harm_its2, eDNA_harm)
    ) %>%
    tab_spanner(
        label = "genus-level",
        columns = c(eDNA_genus_trnl, eDNA_genus_its2, eDNA_genus)
    ) %>%
    cols_label(
        plot_id = "Plot",
        botanical = html("botanical<br>survery"),
        trnl = "trnl",
        its2 = "ITS2",
        eDNA = "both",
        eDNA_harm_trnl = "trnl",
        eDNA_harm_its2 = "ITS2",
        eDNA_harm = "both",
        eDNA_genus_trnl = "trnl",
        eDNA_genus_its2 = "ITS2",
        eDNA_genus = "both",
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            trnl,
            its2,
            eDNA,
            eDNA_harm_trnl,
            eDNA_harm_its2,
            eDNA_harm,
            eDNA_genus_trnl,
            eDNA_genus_its2,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/3_detections_per_plot.html")

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
        across(trnl:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0)
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    tab_spanner(label = "species-level", columns = c(trnl, its2, eDNA)) %>%
    tab_spanner(
        label = "harmonised",
        columns = c(eDNA_harm_trnl, eDNA_harm_its2, eDNA_harm)
    ) %>%
    tab_spanner(
        label = "genus-level",
        columns = c(eDNA_genus_trnl, eDNA_genus_its2, eDNA_genus)
    ) %>%
    cols_label(
        botanical = html("botanical<br>survery"),
        trnl = "trnl",
        its2 = "ITS2",
        eDNA = "both",
        eDNA_harm_trnl = "trnl",
        eDNA_harm_its2 = "ITS2",
        eDNA_harm = "both",
        eDNA_genus_trnl = "trnl",
        eDNA_genus_its2 = "ITS2",
        eDNA_genus = "both",
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            trnl,
            its2,
            eDNA,
            eDNA_harm_trnl,
            eDNA_harm_its2,
            eDNA_harm,
            eDNA_genus_trnl,
            eDNA_genus_its2,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/4_sites_per_species.html")

## Plot per species by method
obs_combined %>%
    select(-genus, -species) %>%
    group_by(short_name) %>%
    arrange(short_name) %>%
    summarise(across(botanical:eDNA_genus, sum)) %>%
    mutate(
        across(trnl:eDNA_genus, function(x) {
            paste0(
                format(x, width = 2, digits = 2, nsmall = 0)
            )
        }),
        botanical = format(botanical, width = 2),
    ) %>%
    gt() %>%
    tab_spanner(label = "species-level", columns = c(trnl, its2, eDNA)) %>%
    tab_spanner(
        label = "harmonised",
        columns = c(eDNA_harm_trnl, eDNA_harm_its2, eDNA_harm)
    ) %>%
    tab_spanner(
        label = "genus-level",
        columns = c(eDNA_genus_trnl, eDNA_genus_its2, eDNA_genus)
    ) %>%
    cols_label(
        botanical = html("botanical<br>survery"),
        trnl = "trnl",
        its2 = "ITS2",
        eDNA = "both",
        eDNA_harm_trnl = "trnl",
        eDNA_harm_its2 = "ITS2",
        eDNA_harm = "both",
        eDNA_genus_trnl = "trnl",
        eDNA_genus_its2 = "ITS2",
        eDNA_genus = "both",
    ) %>%
    cols_align(
        align = "center",
        columns = c(
            botanical,
            trnl,
            its2,
            eDNA,
            eDNA_harm_trnl,
            eDNA_harm_its2,
            eDNA_harm,
            eDNA_genus_trnl,
            eDNA_genus_its2,
            eDNA_genus
        )
    ) %>%
    tab_options(
        table_body.hlines.color = "white",
        table_body.vlines.color = "white"
    ) %>%
    gtsave("tables/5_plots_per_species.html")

## End of file
