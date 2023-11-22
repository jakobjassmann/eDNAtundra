# Get ITS2 and rbcL data for vascular plant species list from BOLD
# Jakob J. Assmann jakob.assmann@uzh.ch 22 November 2023

# Dependencies
# install.packages("bold")
library(bold)
library(tidyverse)

species <- "Salix arctica"
marker <- "*rbcL*"

# Define helper funciton to get sequence for a species and marker
get_bold <- function(species, marker) {
    # Status
    cat("\nSpecies:", species, "\n")
    cat("Marker:", marker, "\n\n")

    # Grab database entries
    meta_data <- bold_seqspec(taxon = species)
    meta_data <- filter(meta_data, grepl(marker, markercode))
    
    # Summary stats
    cat("Total entries: ", nrow(meta_data), "\n")
    meta_data_can <- filter(meta_data, country == "Canada")
    cat("Canada entries: ", nrow(meta_data_can), "\n")
    meta_data_nun <- filter(meta_data_can, province_state == "Nunavut")
    cat("Nunavut entries: ", nrow(meta_data_nun), "\n")
    meta_data_vic <- filter(meta_data_can, region == "Victoria Island")
    cat("Victoria Island entries: ", nrow(meta_data_vic), "\n\n")
    
    # Placeholders
    selection <- NULL
    return_seq <- NULL
    
    # Select from Victoria Island data
    if (nrow(meta_data_vic) >= 1) {
        cat("Victoria Island\n")
        meta_data_vic %>%
            select(exactsite, institution_storing, country, province_state, region) %>%
            tibble() %>%
            print(n = nrow(meta_data_vic))
        selection <- readline("\nSelect line number [break = n]:")
        while (!(as.numeric(selection) <= nrow(meta_data_vic))) {
            if (selection == "n") {
                selection <- NULL
                break()
            }
            selection <- readline("\nSelect line number [break = n]:")
        }
        return_seq <- meta_data_vic[selection, ]
    }
    
    
    # Select from Nunavut data
    if (nrow(meta_data_nun) >= 1 && is.null(selection)) {
        cat("Nunavut\n")
        meta_data_nun %>%
            select(exactsite, institution_storing, country, province_state, region) %>%
            tibble() %>%
            print(n = nrow(meta_data_nun))
        selection <- readline("\nSelect line number [break = n]:")
        while (!(as.numeric(selection) <= nrow(meta_data_nun))) {
            if (selection == "n") {
                selection <- NULL
                break()
            }
            selection <- readline("\nSelect line number [break = n]:")
        }
        return_seq <- meta_data_nun[selection, ]
    }

    # Select from Canada data    
    if (nrow(meta_data_can) >= 1 && is.null(selection)) {
        cat("Canada\n")
        meta_data_can %>%
            select(exactsite, institution_storing, country, province_state, region) %>%
            tibble() %>%
            print(n = nrow(meta_data_can))
        selection <- readline("\nSelect line number [break = n]:")
        while(!(as.numeric(selection) <= nrow(meta_data_can))) {
            if (selection == "n") {
                selection <- NULL
                break()
            }
            selection <- readline("\nSelect line number [break = n]:")
        }
        return_seq <- meta_data_can[selection, ]
    }

    # Select from global data  
    if (nrow(meta_data) >= 1 && is.null(selection)) {
        cat("Global\n")
        meta_data %>%
            select(exactsite, institution_storing, country, province_state, region) %>%
            tibble() %>%
            print(n = nrow(meta_data))
        selection <- readline("\nSelect line number [break = n]:")
        while (!(as.numeric(selection)<= nrow(meta_data))) {
            if (selection == "n") {
                selection <- NULL
                break()
            }
            selection <- readline("\nSelect line number [break = n]:")
        }
        return_seq <- meta_data[selection, ]
    }
    cat("\nDone.\n")
    return(return_seq)
}

# Get species list and trnl sequences
species_list <- read_csv("data/willerslev_db_subset_final.csv")


## rbcl 
# Placeholders
rbcl_seqs <- NULL
for (i in 1:nrow(species_list)) {
    species <- species_list$species_name[i]
    marker <- "*rbcL*"
    rbcl_seqs <- get_bold(species, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, marker = marker) %>%
        bind_rows(rbcl_seqs, .)
}
write_csv(rbcl_seqs, "data/bold_rbcl_seqs_first.csv")
# rbcl_seqs <- read_csv("data/bold_rbcl_seqs_first.csv")
# Synonym search required for missing species: 
species_list$species_name[!(species_list$species_name %in% rbcl_seqs$species)]
# "Carex borealipolaris"    "Carex myosuroides"  "Chamaenerion latifolium"
# and "Carex simpliciuscula" and "Sabulina rubella" (only global entries)
species_with_synonyms <- filter(species_list, !is.na(synonym)) %>%
    filter(species_name %in% c(
        "Carex borealipolaris",
        "Carex myosuroides",
        "Chamaenerion latifolium",
        "Carex simpliciuscula",
        "Sabulina rubella"
    ))
# Remove any old entries
rbcl_seqs <- filter(rbcl_seqs, !(species %in%
    c(
        "Carex borealipolaris",
        "Carex myosuroides",
        "Chamaenerion latifolium",
        "Carex simpliciuscula",
        "Sabulina rubella"
    )))
for (i in 1:nrow(species_with_synonyms)) {
    species <- species_with_synonyms$species_name[i]
    synonym <- species_with_synonyms$synonym[i]
    marker <- "*rbcL*"
    rbcl_seqs <- get_bold(synonym, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, synonym = synonym, marker = marker) %>%
        bind_rows(rbcl_seqs, .)
}
# Check final outputs
rbcl_seqs
distinct(rbcl_seqs, species)
write_csv(rbcl_seqs, "data/bold_rbcl_seqs_final.csv")


## Repeat for ITS2
# Placeholders
its_seqs <- NULL
for(i in 45:nrow(species_list)) {
    species <- species_list$species_name[i]
    marker <- "*ITS2*"
    its_seqs <- get_bold(species, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, marker = marker) %>%
        bind_rows(its_seqs, .)
}
write_csv(its_seqs, "data/bold_its_seqs_first.csv")
#its_seqs <- read_csv("data/bold_its_seqs_first.csv")

# species with issues
its_seqs$species[is.na(its_seqs$processid)]
# and "Carex simpliciuscula", "Sabulina rubella"
species_to_repeat <- c(
    "Carex bigelowii",
    "Stellaria longipes",
    "Carex myosuroides",
    "Carex simpliciuscula",
    "Sabulina rubella"
)
# remove any already existing ones
its_seqs <- filter(its_seqs, !(species %in% species_to_repeat))
# Get synonyms
species_to_repeat <- species_list %>% filter(species_name %in% species_to_repeat)

for(i in 1:nrow(species_to_repeat[1:2,])) {
    species <- species_to_repeat[1:2,]$species_name[i]
    marker <- "*ITS2*"
    its_seqs <- get_bold(species, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, marker = marker) %>%
        bind_rows(its_seqs, .)
}

for(i in 3:nrow(species_to_repeat)) {
    species <- species_to_repeat$species_name[i]
    synonym <- species_to_repeat$synonym[i]
    marker <- "*ITS2*"
    its_seqs <- get_bold(synonym, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, synonym = synonym, marker = marker) %>%
        bind_rows(its_seqs, .)
}

# Checks
tally(its_seqs)
distinct(its_seqs, species) %>% tally()

# 51... which ones are missing?
(species_to_repeat <- species_list[!(species_list$species_name %in% its_seqs$species), ])
# Repeat:
for(i in 1:nrow(species_to_repeat)) {
    species <- species_to_repeat$species_name[i]
    synonym <- species_to_repeat$synonym[i]
    marker <- "*ITS2*"
    its_seqs <- get_bold(synonym, marker) %>%
        tibble() %>%
        select(
            -contains("image"),
            -contains("copyright"),
            -contains("catalognum")
        ) %>%
        mutate(species = species, synonym = synonym, marker = marker) %>%
        bind_rows(its_seqs, .)
}

# Final Checks
tally(its_seqs)
distinct(its_seqs, species) %>% tally()
species_list[!(species_list$species_name %in% its_seqs$species), ]

# Write out
write_csv(its_seqs, "data/bold_its_seqs_final.csv")
