# Quick script to check the completness of the DS-BITSR for our moss species list
# 4 August 2026 Jakob J. Assmann jakob.assmann@gu.se

# Dependencies
library(tidyverse)
library(Biostrings)

# Load bryophyte species list and filter unique species
bryo_obs <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb_updated_Aug_2026.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(
        site_id = Site,
        plot_id,
        genus = Genus,
        species = Species,
        short_name = Full_Name,
        synonym
    ) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    filter(!grepl("sp\\.", species)) %>%
    distinct()

# Save distinct species and their synonyms for later
bryo_obs_distinct <- distinct(bryo_obs, short_name, synonym) %>%
    arrange(short_name)

# Use synynoms for DS-BITSR compariosn
bryo_obs <- bryo_obs %>%
    mutate(short_name = synonym) %>%
    distinct(short_name) %>%
    arrange(short_name)

# Load species names in DS-BITSR
ds_bitsr <- readDNAStringSet(
    "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta"
) %>%
    names() %>%
    gsub("ID_", "", .) %>%
    gsub("_", " ", .) %>%
    unique() %>%
    sort()

# Remove any sub-species level IDs for matching
ds_bitsr[grepl("subsp\\.", ds_bitsr)] <- gsub(
    "(.) subsp\\..*",
    "\\1",
    ds_bitsr[grepl("subsp\\.", ds_bitsr)]
)
ds_bitsr[grepl("var\\.", ds_bitsr)] <- gsub(
    "(.) subsp\\..*",
    "\\1",
    ds_bitsr[grepl("var\\.", ds_bitsr)]
)

# Check which species are present in reference database
ds_bitsr_complete_check <- tibble(
    bryo_obs,
    ds_bitsr = bryo_obs$short_name %in% ds_bitsr
)

# Quick check of direct matches
sum(ds_bitsr_complete_check$ds_bitsr)
# [1] 52

# No match
sum(!ds_bitsr_complete_check$ds_bitsr)
# 15

# revert to names used in our taxonomy
ds_bitsr_complete_check <- ds_bitsr_complete_check %>%
    dplyr::rename(synonym = short_name) %>%
    full_join(bryo_obs_distinct) %>%
    mutate(synonym = if_else(synonym == short_name, "", synonym)) %>%
    arrange(short_name)

# Export CSV of matches
write_excel_csv(
    ds_bitsr_complete_check %>%
        select(short_name, synonym, ds_bitsr) %>%
        setNames(c("Name", "Synonym", "In_DS-BITSR")),
    "data/external/reference_libraries/DS-BITSR/completness_check.csv"
)

# EOF
