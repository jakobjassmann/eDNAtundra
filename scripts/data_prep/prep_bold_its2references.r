# Quick script to help with updated reference library 12 March 2025
# Jakob J. Assmann jakob.assmann@uzh.ch

# Dependencies
library(tidyverse)
library(msa)

# get species list
species_list <- read_csv("data/botanical_obs/species_list_cleaned.csv")

# load DS_VASCAN database from BOLD
ds_vascan <- read_tsv("data/external/reference_libraries/BOLD/DS-VASCAN.tsv")
names(ds_vascan)

# Define helper function to filter DS_VASCAN for a species
filter_bold <- function(short_name) {
        # Status
        cat(short_name, "\n")
        # Filter DS_VASCAN dataset for species and ITS2 
        bold_sub <- filter(
            ds_vascan,
            grepl(short_name, species),
            grepl("ITS2", marker_code),
            grepl("ITS_S2F", primers_forward),
            grepl("ITS4", primers_reverse)
        ) %>%
            relocate(processid, species, subspecies, inst, site, marker_code, primers_forward, primers_reverse, nuc) %>%
            mutate(nuc = gsub("-", "", nuc)) %>%
            mutate(short_name == short_name)
        # Check agreement between obtained reference sequences 
        if (nrow(bold_sub) > 1) {
            seq_agreement <- msa(bold_sub$nuc, type = "dna") %>%
                consensusMatrix() %>%
                (function(x) colSums(x == nrow(bold_sub))) %>%
                (function(x) round(100 * sum(x == 1) / length(x), 2)) %>%
                 paste0(" | agreement: ", ., "%")
        } else {
        # If there was only one or none don't specify agreement
            seq_agreement <- ""
        }
        # Status
        cat("   Records: ", nrow(bold_sub), seq_agreement, "\n")
        return(bold_sub)
    }
# Filter DS_VASCAN for species in species list
its2_bold <- species_list %>%
    pull(short_name) %>%
    map(filter_bold) %>%
    bind_rows()
its2_bold %>% print(n = 100)
ds_vascan %>%
    filter(processid == "MKPCH424-09", marker_code == "ITS2") %>%
    as.data.frame()


(test <- filter_bold("Cardamine digitata"))
test$country_iso
test$region
test$site
test$inst
test %>% pull(nuc)
test$species %>% sort()
msa(test$nuc, type = "dna")
print(msa(test$nuc , type = "dna"), show="complete")

print(msa(c(test$nuc, 
     filter(ref_seqs, grepl("Cardamine digitata", species)) %>% pull(its2_seq)), 
     type = "dna"), show="complete")
