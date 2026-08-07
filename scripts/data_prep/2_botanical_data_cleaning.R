# Botanical data prep / cleaning
# Jakob J. Assmann 14 August 2023

# Dependencies
library(tidyverse)
library(readxl)

## Presence/Absence Observations -----

# Load XLS spread sheet
botanical_obs <- read_xlsx(
    "data/botanical_obs/2025-03-11_botanical_ops_final.xlsx"
) %>%
    filter(Cleaned == "Yes") %>%
    filter(Group == "Vascular") %>%
    select(
        date = "Date",
        observer = "Observer",
        scribe = "Scribe",
        site = "Site",
        plot_id = "Plot_ID",
        genus = "Genus",
        species = "Species",
        full_name = `Full name`
    ) %>%
    mutate(plot_id = paste0("S", gsub("Site ", "", site), plot_id)) %>%
    filter(!grepl("NA", full_name) & !is.na(full_name))

# Check unique species
distinct(botanical_obs, full_name) %>%
    arrange(full_name) %>%
    print(n = 101)

# Remove all entries without an ID to genus and set unknown species to sp.
n_before <- nrow(botanical_obs)
botanical_obs <- botanical_obs %>%
    filter(!grepl("Eudicot", genus)) %>%
    filter(!grepl("Monocot", genus)) %>%
    filter(!grepl("Fabaceae", genus)) %>%
    filter(!grepl("Poaceae", genus)) %>%
    filter(!grepl("Cyperaceae", genus)) %>%
    filter(!grepl("Brassicaceae", genus)) %>%
    filter(!grepl("cf", genus)) %>%
    mutate(
        species = case_when(
            grepl("cf", species) ~ "species",
            grepl("^sp\\.", species) ~ "species",
            TRUE ~ species
        ),
        full_name = case_when(
            is.na(species) ~ genus,
            TRUE ~ paste(genus, species)
        )
    )
n_after <- nrow(botanical_obs)
n_before - n_after
# Removed 29 rows

# Check unique species again
distinct(botanical_obs, full_name) %>%
    arrange(full_name) %>%
    print(n = 100)

# Remove subspecies names in species column
botanical_obs <- botanical_obs %>%
    mutate(
        species = sapply(str_extract_all(full_name, "\\S+"), `[`, 2)
    )

# write out cleaned data file
write_csv(botanical_obs, "data/botanical_obs/presence_absence_cleaned.csv")

# Write out a species list
botanical_obs %>%
    distinct(genus, species, full_name) %>%
    mutate(
        short_name = gsub("([A-Za-z]*) ([A-Za-z]*).*", "\\1 \\2", full_name)
    ) %>%
    filter(!grepl("species", species)) %>%
    arrange(short_name) %>%
    write_csv(., "data/botanical_obs/species_list_cleaned.csv")

## Point fraiming ----

# Load point frame data
point_frame_obs <- read_xlsx(
    "data/botanical_obs/2025-03-11_botanical_ops_final.xlsx",
    sheet = "Pointframing"
)

# Rows before filtering
nrow(point_frame_obs)

# Reformat column names
point_frame_obs <- point_frame_obs %>%
    filter(`Cleaned?` == "Yes") %>%
    select(
        date = "Date",
        observer = "Observer",
        scribe = "Scribe",
        site = "Site",
        plot_id = "Plot_ID",
        point = "Point_No.",
        hit = "Hit_No.",
        genus = "Genus",
        species = "Species"
    ) %>%
    # Add site id to plot id
    mutate(plot_id = paste0("S", gsub("Site ", "", site), plot_id)) %>%
    # Remove any NAs (there are none)
    filter(!is.na(genus)) %>%
    # Set "species" for all uncertain species ids to NA
    mutate(
        species = case_when(
            grepl("cf", species) ~ "species",
            grepl("^sp\\.", species) ~ "species",
            TRUE ~ species
        ),
        genus = case_when(
            grepl("^cf", genus) ~ "Unknown",
            grepl("NA", genus) ~ "Unknown",
            grepl("Cyperaceae", genus) ~ "Unknown",
            grepl("Unknown", genus) ~ "Unknown",
            TRUE ~ genus
        ),
        species = case_when(
            grepl("Unknown", genus) ~ "species",
            grepl("Carex", genus) & is.na(species) ~ "species",
            TRUE ~ species
        ),
        full_name = case_when(
            is.na(species) ~ genus,
            TRUE ~ paste(genus, species)
        )
    )
nrow(point_frame_obs)

# Check unique species
distinct(point_frame_obs, full_name) %>%
    arrange(full_name) %>%
    print(n = 250)

# Remove subspecies names in species column
point_frame_obs <- point_frame_obs %>%
    mutate(
        species = sapply(str_extract_all(species, "\\S+"), `[`, 1)
    )

write_csv(point_frame_obs, "data/botanical_obs/pointframing.csv")

## QC: Compare with species list from Saarela et al. 2020 ----
species_list <- read_csv(
    "data/external/species_list/species_list_Saarela_et_al_2020.csv"
)

# Presence / absence
botanical_obs %>%
    # Remove all rows with uncertain species ids
    filter(species != "species") %>%
    distinct(genus, species) %>%
    mutate(obs = "botanical obs") %>%
    full_join(
        .,
        select(species_list, genus, species) %>% mutate(list = "species list"),
        by = join_by(genus == genus, species == species)
    ) %>%
    arrange(genus, species) %>%
    print(n = 200)
# Check observations not in species list...
botanical_obs %>%
    filter(species != "species") %>%
    distinct(genus, species) %>%
    mutate(obs = "botanical obs") %>%
    full_join(
        .,
        select(species_list, genus, species) %>% mutate(list = "species list"),
        by = join_by(genus == genus, species == species)
    ) %>%
    arrange(genus, species) %>%
    filter(!is.na(obs), is.na(list)) %>%
    print(n = 200)
# Not matching due to ID differences:
# Carex bigelowii Torr. ex Schwein. not id'ed to sub-species, but only listed as
# subspecies in Saarela et al. 2020
# Not matching because previously not recorded
# Tofieldia coccinea
# Note: Previously spelled Tofielda coccinea Richardson (genus spelled wrong), should be : Tofieldia -> corrected in file version 09-03-2025

# Point framing
point_frame_obs %>%
    filter(species != "species") %>%
    mutate(species = gsub("([A-Za-z]*) (.*)", "\\1", species)) %>%
    distinct(genus, species) %>%
    mutate(obs = "botanical obs") %>%
    full_join(
        .,
        select(species_list, genus, species) %>%
            mutate(
                list = "species list",
                species = gsub("([A-Za-z]*) (.*)", "\\1", species)
            ),
        by = join_by(genus == genus, species == species)
    ) %>%
    arrange(genus, species) %>%
    filter(!is.na(obs), is.na(list)) %>%
    print(n = 200)
# Only mosses that are not listed.
# Check genera with no species ID
point_frame_obs %>%
    filter(species == "species") %>%
    distinct(genus) %>%
    print(n = 100)
# Only Carex is a vascular plant genus (and this is correct)
# Check other surfaces
point_frame_obs %>%
    filter(is.na(species)) %>%
    distinct(genus)

## QC: Check agreement with pan-Arctic species list ----

# Load pan-Arctic species list

pasl <- read_xlsx("data/external/species_list/pasl_breen_uaf_original.xlsx")

# Prepare observations to PASL format
(not_in_pasl <- botanical_obs %>%
    mutate(
        species = gsub("([A-Za-z]*) (.*)", "\\1", species),
        full_name = paste(genus, species),
        obs = "botanical obs"
    ) %>%
    distinct(full_name, obs) %>%
    left_join(
        .,
        select(pasl, ABBREVIAT) %>% mutate(list = "PASL"),
        by = join_by(full_name == ABBREVIAT)
    ) %>%
    arrange(full_name) %>%
    distinct() %>%
    # Only show species not listed in PASL
    filter(is.na(list)))

# Species with no match in PASL:
# 1 Carex borealipolaris -> listed as Kobresia sibirica (synonym)
# 2 Cassiope tetragona -> listed as Cassiope tetragona [s. tetragona]
# 3 Hedysarum boreale  -> listed as Hedysarum boreale subsp. mackenzii
# 4 Sabulina rubella   -> listed as Minuartia rubella (synonym)

## QC: Check Mosses ----

mosses <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv")
(not_in_pasl <- mosses %>%
    distinct(Genus, Species) %>%
    mutate(
        Species = case_when(
            Species == "sp." ~ "species",
            TRUE ~ Species
        ),
        full_name = paste(Genus, Species),
        obs = "obs"
    ) %>%
    select(full_name, obs) %>%
    print(n = 100) %>%
    left_join(
        .,
        select(pasl, ABBREVIAT) %>% mutate(list = "PASL"),
        by = join_by(full_name == ABBREVIAT)
    ) %>%
    arrange(full_name) %>%
    distinct() %>%
    filter(is.na(list)))
# Species not in PASL
#  1 Drepanocladus arcticus              -> listed as Campylium arcticum (synonym)
#  2 Encalypta vittiana                  -> not in PASL incl. no known synonyms - listed in Flora of North America
#  3 Eurynchiastrum pulchellum           -> listed as Hypnum pulchellum (basionym)
#  4 Neoorthocaulis floerkei             -> listed as Orthocaulis floerkei (synonym)
#  5 Ptychostomum pallens                -> not in PASL incl. no know synonyms - listed in Flora of North America
#  6 Ptychostomum species                -> genus not listed in PASL incl. no know synonyms
#  7 Scapania gymnostophila              -> listed as Scapania gymnostomophila (synonym)
#  8 Stegonia latifolia var. pilifera    -> sp., but not var. listed in PASL
#  9 Tortella tortuosa var. fragilifolia -> sp., but not var. listed in PASL
# 10 Tortella tortuosa var. tortuosa     -> sp., but not var. listed in PASL

# Update 5 August 2028
# Erin spotted some spelling mistakes and synonyms that might have caused
# mismatch with DS-BITSR even though species were detected.
# Here I correct for those and add a synonym colum
moss_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv")

# Spelling errors (2 species affected)
moss_obs <- moss_obs %>%
    mutate(
        Full_Name = case_when(
            # Missing "h" in genus for Eurhynchiastrum pulchellum
            Full_Name ==
                "Eurynchiastrum pulchellum" ~ "Eurhynchiastrum pulchellum",
            # Missing "mo" in genus for Scapania gymnostomophila
            Full_Name == "Scapania gymnostophila" ~ "Scapania gymnostomophila",
            TRUE ~ Full_Name
        )
    ) %>%
    mutate(
        Genus = if_else(
            Full_Name == "Eurhynchiastrum pulchellum",
            "Eurhynchiastrum",
            Genus
        ),
        Species = if_else(
            Full_Name == "Scapania gymnostomophila",
            "gymnostomophila",
            Species
        )
    ) %>%
    # Remove any duplicates
    distinct()

# Remove identifications at the sub-species level
moss_obs <- moss_obs %>%
    mutate(
        Species = case_when(
            Full_Name == "Stegonia latifolia var. pilifera" ~ "latifolia",
            Full_Name == "Tortella tortuosa var. tortuosa" ~ "tortuosa",
            Full_Name == "Tortella tortuosa var. fragilifolia" ~ "tortuosa",
            TRUE ~ Species
        ),
        Full_Name = case_when(
            Full_Name ==
                "Stegonia latifolia var. pilifera" ~ "Stegonia latifolia",
            Full_Name ==
                "Tortella tortuosa var. tortuosa" ~ "Tortella tortuosa",
            Full_Name ==
                "Tortella tortuosa var. fragilifolia" ~ "Tortella tortuosa",
            TRUE ~ Full_Name
        )
    ) %>%
    # Remove any newly created duplicates
    distinct()

# Add synonyms for mathching with DS-BITSR
moss_obs <- moss_obs %>%
    mutate(
        synonym = case_when(
            Full_Name == "Ditrichum flexicaule" ~ "Flexitrichum flexicaule",
            Full_Name == "Ditrichum gracile" ~ "Flexitrichum gracile",
            Full_Name == "Hypnum bambergeri" ~ "Campylium bambergeri",
            Full_Name == "Hypnum procerrimum" ~ "Pseudostereodon procerrimus",
            Full_Name == "Hypnum revolutum" ~ "Roaldia revoluta",
            Full_Name == "Hypnum vaucheri" ~ "Buckia vaucheri",
            Full_Name == "Orthotrichum speciosum" ~ "Lewinskya speciosa",
            Full_Name ==
                "Pseudocalliergon brevifolium" ~ "Drepanocladus brevifolius",
            # Add the original name for those that don't have synonyms
            TRUE ~ Full_Name
        )
    )

# Write out update
write_csv(
    moss_obs,
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb_updated_Aug_2026.csv"
)
#
# End of file ----
