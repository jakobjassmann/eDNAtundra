# Script to perpare species list -- 2025-03-09
# Jakob J. Assmann jakob.assmann@uzh.ch

# Dependencies
library(tidyverse)
library(readxl)

# Load file from Saarela et al. 2020
species_list <- read_xlsx("data/external/species_list/phytokeys-141-001-s003.xlsx")
species_list <- species_list %>%
    select(2:13) %>%
    set_names(., slice(., 1) %>% as.character()) %>%
    slice(c(-1, -2)) %>%
    filter(., !is.na(`Cambridge Bay`) | !is.na(`Ovayok Territorial Park`)) %>%
    select(1:3) %>%
    mutate(
        genus = gsub("([A-Za-z]*) (.*)", "\\1", Taxon),
        species = gsub("([A-Za-z]*) (.*)", "\\2", Taxon),
        full_name = Taxon
    ) %>%
    select(-Taxon) 

# Write to file
write_csv(
    species_list,
    "data/external/species_list/species_list_Saarela_et_al_2020.csv"
)
