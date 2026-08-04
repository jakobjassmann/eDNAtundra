# Quick script to check the completness of the DS-BITSR for our moss species list
# 4 August 2026 Jakob J. Assmann jakob.assmann@gu.se

# Dependencies
library(tidyverse)
library(Biostrings)
library(readxl)

# Load bryophyte species list and filter unique species
bryo_obs <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    filter(!grepl("sp.", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    distinct(short_name) %>%
    sort()

# Load species names in DS-BITSR
ds_bitsr <- readDNAStringSet(
    "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta"
) %>%
    names() %>%
    gsub("ID_", "", .) %>%
    gsub("_", " ", .) %>%
    unique() %>%
    sort()

# Check which species are present in reference database
ds_bitsr_complete_check <- tibble(
    bryo_obs,
    ds_bitsr = bryo_obs$short_name %in% ds_bitsr
)

# Quick check of direct matches
sum(ds_bitsr_complete_check$ds_bitsr)
# [1] 41

# Export CSV of matches
write_excel_csv(
    ds_bitsr_complete_check %>%
        setNames(c("Name", "In_DS-BITSR")),
    "data/external/reference_libraries/DS-BITSR/completness_check.csv"
)

# EOF
