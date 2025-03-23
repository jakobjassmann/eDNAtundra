# Bryophytes: Retrieve presence absence for each plot based on its2
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(dada2)
library(Biostrings)
library(ShortRead)

# Load sequences and taxonomic assignments
load("data/sequences/its2/edna/cleaned/its2_seqtab_nochim.Rda")

##  Assign species ----
its2_taxa <- assignTaxonomy(its2_seqtab_nochim,
                           "data/sequences/reference/its2_reference_bryo_BOLD.fasta",
                           verbose = T,
                           multithread = TRUE)
its2_taxa %>% data.frame() %>% remove_rownames() %>% filter(!is.na(Genus)) %>% 
  distinct() %>% 
  arrange(Genus, Species) %>% View()

its2_bryo_taxa <- its2_taxa

# Save object
save(its2_bryo_taxa, file = "data/sequences/its2/edna/cleaned/its2_bryo_taxa.Rda")

# Collapse data into a dataframe of species presence / absence per plot
tax_counts <- right_join(
  # Prepare sequence ocunts
  t(its2_seqtab_nochim) %>%
    data.frame() %>%
    mutate(., sq = row.names(.)) %>%
    remove_rownames() %>%
    relocate(sq) %>%
    pivot_longer(.,
                 cols = 2:ncol(.),
                 names_to = "sample_id",
                 values_to = "count"
    ) %>%
    mutate(
      sample_id = gsub("\\.", "-", sample_id),
      plot_id = case_when(
        grepl("S[0-9]{1}[SR]{1}[0-9]{1}", sample_id) ~ substr(sample_id, 1, 4),
        TRUE ~ NA
      ),
      subsample_id = case_when(
        grepl("S[0-9]{1}[SR]{1}[0-9]{1}", sample_id) ~ substr(sample_id, 6, 6),
        TRUE ~ NA
      )
    ),
  # Preprare taxonomic assignments
  data.frame(its2_taxa) %>%
    mutate(., sq = row.names(.)) %>%
    remove_rownames() %>%
    filter(!is.na(Genus))
) %>%
  relocate(sample_id, plot_id, subsample_id, Genus, Species) %>%
  dplyr::rename(genus = Genus, species = Species) %>%
  arrange(sample_id) %>%
  # Retain only presences
  filter(count > 0)

# Write out data table
write_csv(tax_counts, "data/sequences/its2/edna/cleaned/its2_bryo_taxa_in_sample.csv")
