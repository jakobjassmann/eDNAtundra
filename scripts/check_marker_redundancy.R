# A quick script to check for redundancy in our marker sequences for the 
# species observed in the eDNA project
# Jakob J. Assmann jakob.assmann@uzh.ch 22 November 2023
# Inspired by phind https://www.phind.com/search?cache=bob7dh3msjjta8lgz70mu84w

# Dependencies
library(tidyverse)
# BiocManager::install("msa")
# BiocManager::install("ggmsa")
library(msa)
library(seqinr)
library(ggmsa)
library(tidysq)

# Read sequences
trnl <- read_csv("data/sequences/willerslev_db_subset_final.csv")
rbcl <- read_csv("data/sequences/bold_rbcl_seqs_final.csv")
its <- read_csv("data/sequences/bold_its_seqs_final.csv")

# Convert to DNA string sets
trnl_seqs <- DNAStringSet(trnl$trnL_seq[c(-1,-2)])
rbcl_seqs <- DNAStringSet(rbcl$nucleotides)
its_seqs <- DNAStringSet(its$nucleotides)

# Alignment
trnl_aligned <- msa(trnl_seqs)
rbcl_aligned <- msa(rbcl_seqs)
its_aligned <- msa(its_seqs)

# Calculate distance matrices
trnl_dist_matrix <- dist.alignment(msaConvert(trnl_aligned, "seqinr::alignment"), matrix = "identity")
rbcl_dist_matrix <- dist.alignment(msaConvert(rbcl_aligned, "seqinr::alignment"), matrix = "identity")
its_dist_matrix <- dist.alignment(msaConvert(its_aligned, "seqinr::alignment"), matrix = "identity")

# Find indices with perfect alignment
trnl_matches <- data.frame(i = NULL, matches = NULL, n_matches = NULL)
for(i in 1:nrow(trnl)){
    mat_row <- as.matrix(trnl_dist_matrix)[i, ]
    match_ids <- as.numeric(which(mat_row == 0))
    trnl_matches <- trnl_matches %>%
    bind_rows(tibble(i = i+2, matches = list(match_ids+2), n_matches = length(match_ids)))
}

rbcl_matches <- data.frame(i = NULL, matches = NULL, n_matches = NULL)
for (i in 1:nrow(rbcl)) {
    mat_row <- as.matrix(rbcl_dist_matrix)[i, ]
    match_ids <- as.numeric(which(mat_row == 0))
    rbcl_matches <- rbcl_matches %>%
        bind_rows(tibble(i = i, matches = list(match_ids), n_matches = length(match_ids)))
}

its_matches <- data.frame(i = NULL, matches = NULL, n_matches = NULL)
for(i in 1:nrow(its)){
    mat_row <- as.matrix(its_dist_matrix)[i, ]
    match_ids <- as.numeric(which(mat_row == 0))
    its_matches <- its_matches %>%
    bind_rows(tibble(i = i, matches = list(match_ids), n_matches = length(match_ids)))
}

# combine all into one data frame to visualise
combined_data <- trnl %>%
    select(species = species_name, synonym, trnL_seq, trnL_seq_lenght, trnL_unique) %>%
    mutate(., index = 1:nrow(.)) %>%
    relocate(index) %>%
    full_join(rbcl %>% select(species, rbcl_seq = nucleotides)) %>%
    full_join(rbcl_matches %>% select(index = i, rbcl_n_matches = n_matches)) %>%
    mutate(rbcl_n_matches = rbcl_n_matches - 1) %>%
    mutate(rbcl_seq_length = gsub("-", "", rbcl_seq) %>% sq() %>% get_sq_lengths()) %>%
    full_join(its %>% select(species, its2_seq = nucleotides)) %>%
    full_join(its_matches %>% select(index = i, its2_n_matches = n_matches)) %>%
    mutate(its2_n_matches = its2_n_matches - 1) %>%
    mutate(its2_seq_length =  gsub("-", "", its2_seq) %>% sq() %>% get_sq_lengths())
write_csv(combined_data, "data/sequences/species_list_with_all_markers.csv")

# Quick stats summary
combined_data %>%
    summarise(
        trnl_resolution = sum(trnL_unique == "Unique", na.rm = T) / 53,
        rbcl_resolution = sum(rbcl_n_matches == 0) / 53,
        its2_resolution = sum(its2_n_matches == 0) / 53
    )

# What proportion of species is poorly resolved in all?
combined_data %>%
    filter(trnL_unique == "Not unique", rbcl_n_matches > 0, its2_n_matches > 0)
# rbcl and its
combined_data %>%
    filter(rbcl_n_matches > 0, its2_n_matches > 0)
# trnl and its
combined_data %>%
    filter(trnL_unique == "Not unique", its2_n_matches > 0)
# trnl and rbcl
combined_data %>%
    filter(trnL_unique == "Not unique", rbcl_n_matches > 0)
