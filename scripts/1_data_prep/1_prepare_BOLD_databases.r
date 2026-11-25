# Load libraries
library(tidyverse)
library(Biostrings)

## DS-ZACKPLAN

# Load data
lookup_table <- read_tsv(
    "data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN.tsv"
)

# Filter to species and its2 sequences
its2_seqs <- lookup_table %>%
    filter(marker_code == "ITS2") %>%
    select(species, genus, nuc) %>%
    mutate(
        species = case_when(
            is.na(species) ~ paste(genus, "sp."),
            TRUE ~ species
        )
    ) %>%
    mutate(species = paste("ID", species)) %>%
    mutate(species = gsub(" ", "_", species)) %>%
    mutate(nuc = gsub("-", "", nuc)) %>%
    na.omit()

# Create DNAStringSet
sequences <- DNAStringSet(its2_seqs$nuc)
names(sequences) <- its2_seqs$species

# Write FASTA
writeXStringSet(
    sequences,
    filepath = "data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN_its2.fasta"
)

# Create BLAST database in terminal (external command)
# makeblastdb -in data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN_its2.fasta -dbtype nucl -out data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN_its2_db

## DS-BITSR

# Load data
lookup_table <- read_tsv(
    "data/external/reference_libraries/DS-BITSR/DS-BITSR.tsv"
)

# Filter to species and its2 sequences
its2_seqs <- lookup_table %>%
    filter(marker_code == "ITS2") %>%
    select(species, genus, nuc) %>%
    mutate(
        species = case_when(
            is.na(species) ~ paste(genus, "sp."),
            TRUE ~ species
        )
    ) %>%
    mutate(species = paste("ID", species)) %>%
    mutate(species = gsub(" ", "_", species)) %>%
    mutate(nuc = gsub("-", "", nuc)) %>%
    na.omit()

# Create DNAStringSet
sequences <- DNAStringSet(its2_seqs$nuc)
names(sequences) <- its2_seqs$species

# Write FASTA
writeXStringSet(
    sequences,
    filepath = "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta"
)


# Create BLAST database in terminal (external command)
# makeblastdb -in data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta -dbtype nucl -out data/external/reference_libraries/DS-BITSR/DS-BITSR_its2_db

## DS-VASCAN

# Load data
lookup_table <- read_tsv(
    "data/external/reference_libraries/DS-VASCAN/DS-VASCAN.tsv"
)

# Filter to species and its2 sequences
its2_seqs <- lookup_table %>%
    filter(marker_code == "ITS2") %>%
    select(species, genus, nuc) %>%
    mutate(
        species = case_when(
            is.na(species) ~ paste(genus, "sp."),
            TRUE ~ species
        )
    ) %>%
    mutate(species = paste("ID", species)) %>%
    mutate(species = gsub(" ", "_", species)) %>%
    mutate(nuc = gsub("-", "", nuc)) %>%
    na.omit()

# Create DNAStringSet
sequences <- DNAStringSet(its2_seqs$nuc)
names(sequences) <- its2_seqs$species

# Write FASTA
writeXStringSet(
    sequences,
    filepath = "data/external/reference_libraries/DS-VASCAN/DS-VASCAN_its2.fasta"
)

# Create BLAST database in terminal (external command)
# makeblastdb -in data/external/reference_libraries/DS-VASCAN/DS-VASCAN_its2.fasta -dbtype nucl -out data/external/reference_libraries/DS-VASCAN/DS-VASCAN_its2_db
