# Quick script to check trimmed reads

# Dependencies
library(tidyverse)
library(readxl)
library(dada2)
library(Biostrings)
library(ShortRead)
library(msa)

# Look up table
sample_look_up <- read_xlsx("data/sequences/Assmann_et_al_-_metadata.xlsx") %>%
    filter(species != "NA") %>%
    select(sample = SampleID, species)


## Obtain majority reads from next gen sequencing of single species samples ----
# Get single species files
trnl_reads <- list.files("data/sequences/trnl/trnL_single_species_run1/2_Filtered.project/2_Filtered/Merged_reads",
    pattern = "fastq",
    full.names = T
)
its2_reads <- list.files("data/sequences/ITS2/its2_single_species_run1/2_Filtered.project/2_Filtered/Merged_reads",
    pattern = "fastq",
    full.names = T
)

# Get sample names
trnl_sample_names <- substr(basename(trnl_reads), 1, 2)
its2_sample_names <- substr(basename(its2_reads), 1, 2)

# De-replicate reads
trnl_derep <- derepFastq(trnl_reads)
its2_derep <- derepFastq(its2_reads)

# Get most common sequence using an anonymous helper function
# trnl
trnl_seqs <- lapply(1:length(trnl_sample_names), function(i) {
    data.frame(
        sample = trnl_sample_names[i],
        trnl_seq = names(trnl_derep[[i]]@.Data[[1]][1]),
        trnl_count = trnl_derep[[i]]@.Data[[1]][1],
        trnl_perc = 100 * round(trnl_derep[[i]]@.Data[[1]][1] / sum(trnl_derep[[1]]@.Data[[1]]), 2)
    )
}) %>%
    bind_rows() %>%
    remove_rownames()
# its2
its2_seqs <- lapply(1:length(its2_sample_names), function(i) {
    data.frame(
        sample = its2_sample_names[i],
        its2_seq = names(its2_derep[[i]]@.Data[[1]][1]),
        its2_count = its2_derep[[i]]@.Data[[1]][1],
        its2_perc = 100 * round(its2_derep[[i]]@.Data[[1]][1] / sum(its2_derep[[1]]@.Data[[1]]), 2)
    )
}) %>%
    bind_rows() %>%
    remove_rownames()

# Merge datasets
seqs <- full_join(trnl_seqs, its2_seqs, by = join_by(sample == sample))

# Add species id
ref_seqs <- full_join(sample_look_up, seqs, by = join_by(sample == sample))

# Samples without trnl / its2 reference
ref_seqs %>% filter(is.na(trnl_seq)) # 17 samples without a trnl sequence
ref_seqs %>% filter(is.na(its2_seq)) # 10 samples without a its2 sequence
ref_seqs %>% filter(is.na(trnl_seq), is.na(its2_seq)) # 8 samples without either

# Check for duplictes / uncertain ids
ref_seqs %>% filter(grepl("cf", species))

# Only Sabulina rubella is dublicated
(sabulina_sp <- ref_seqs %>% filter(grepl("Sabulina", species)))
# Compare sequences between both samples
msa(sabulina_sp$trnl_seq, type = "dna") %>%
    consensusMatrix() %>%
    (function(x) colSums(x == 2)) %>%
    (function(x) !any(x != 1)) # trnl match!
msa(sabulina_sp$its2_seq, type = "dna") %>%
    consensusMatrix() %>%
    (function(x) colSums(x == 2)) %>%
    (function(x) !any(x != 1)) # ITS match!
# Samples are identical ... remove cf sample
ref_seqs <- ref_seqs %>%
    filter(sample != sabulina_sp$sample[2])

# Sequences with low coverage below the majortity (< 60 %)
# trnl
ref_seqs %>% filter(trnl_perc < 60) %>% pull(species)
# Carex maritima
# Cassiope tetragona
# Pedicularis lanata
# We check below whether these sequences agree with Willerslev et al. 2014 (they do!)
# its2

# Print all as a check 
ref_seqs %>% select(-sample, -trnl_count, -trnl_perc, -its2_count, -its2_perc) %>% 
arrange(species)  %>% print(n = 100)

# Write oute reference
write_csv(ref_seqs, "data/sequences/reference/Assmann_et_al_CB_2023.csv")

## Compare sequences to previously obtained references
online_ref_seqs <- read_csv("data/sequences/reference/Assmann_et_al_eDNA-tundra_species_list_with_markers_updated_27_Nov_23.csv") %>% select(species, synonym, trnl_seq_online = trnL_seq, its2_seq_online = its2_seq)

# Look up and add sample ID
online_ref_seqs <- online_ref_seqs %>%
    mutate(sample = sapply(species, function(x) {
        online_genus <- gsub("([A-Za-z]*) ([A-Za-z]*)", "\\1", x)
        online_species <- gsub("([A-Za-z]*) ([A-Za-z]*)", "\\2", x)
        sample_id <- ref_seqs %>%
            filter(grepl(online_genus, ref_seqs$species) &
                grepl(online_species, ref_seqs$species)) %>%
            pull(sample)
        if (length(sample_id) == 0) {
            return("None")
        }
        return(sample_id)
    }, USE.NAMES = FALSE))

# Merge online references with sequenced
ref_seqs_combined <- full_join(ref_seqs %>% dplyr::rename(full_name = species),
    online_ref_seqs,
    by = join_by(sample == sample)
) %>%
    select(
        sample, species, full_name, trnl_seq, trnl_seq_online,
        its2_seq, its2_seq_online
    )

# Check sequence availability
ref_seqs_combined %>% select(full_name, trnl_seq, trnl_seq_online, its2_seq, its2_seq_online) %>%
print(n = 100)

# Compare sample sequences and online sequences
ref_seqs_combined <- ref_seqs_combined %>%
    mutate(
        trnl_match = mapply(trnl_seq, trnl_seq_online,
            FUN = function(trnl_seq, trnl_online) {
                if (is.na(trnl_seq) | is.na(trnl_online)) {
                    return(NA)
                }
                msa(c(trnl_seq, trnl_online), type = "dna") %>%
                    consensusMatrix() %>%
                    (function(x) colSums(x == 2)) %>%
                    (function(x) round(100 * sum(x == 1) / length(x), 2))
            }
        ),
        its2_match = mapply(its2_seq, its2_seq_online,
            FUN = function(its2_seq, its2_online) {
                if (is.na(its2_seq) | is.na(its2_online)) {
                    return(NA)
                }
                msa(c(its2_seq, its2_online), type = "dna") %>%
                    consensusMatrix() %>%
                    (function(x) colSums(x == 2)) %>%
                    (function(x) round(100 * sum(x == 1) / length(x), 2))
            }
        )
    )

# Filter out species with less than 97% match for trnl
ref_seqs_combined %>%
    select(species, full_name, trnl_match, trnl_seq, trnl_seq_online) %>%
    filter(trnl_match < 97) %>%
    print() %>%
    split(., .$species) %>%
    map(function(x) {
        print(paste0(x$species, " - match: ", x$trnl_match))
        print(x$trnl_seq)
        print(x$trnl_seq_online)
        return(NULL)
    })
# - Three species have differences in local vs Willerslev et al. 2014
# Draba nivalis
# Sabulina rubella
# Bistorta vivipara

# Check how the species above with low percentage of reads for the trnl seq (above)
# match the references from Willerslev et al. 2014
ref_seqs_combined %>% filter(species %in% c("Carex marina", 
"Cassiope tetragona", "Pedicularis lanata")) %>% select(species, trnl_match)
# All agree 100%!

# Filter out species with less than 97% match for its2
ref_seqs_combined %>%
    select(species, full_name, its2_match, its2_seq, its2_seq_online) %>%
    filter(its2_match < 97) %>%
    # arrange(full_name) %>%
    print(n = 100) %>%
    split(., .$species) %>%
    map(function(x) {
        cat(paste0("\n", x$species, " - match: ", x$its2_match))
        cat("\n", x$its2_seq)
        cat("\n", x$its2_seq_online)
        return(NULL)
    })

# There is a lot of difference - perhaps we should not trus the sequences
# derived from the online databases, but just use the barcodes we derived.

## Compare with reference database generated by Microsynth ----
# Load microsynth reference
microsynth_seqs <- readDNAStringSet("C:/Users/jakob/repositories/eDNAtundra/data/sequences/trnl/trnL_eDNA_run1_site_vs_single_species_and_arctic_gh/customdb2_single_species_and_arctic_gh.fasta")

# Convert to same format aus our database
microsynth_seqs <- data.frame(
    species = microsynth_seqs@ranges@NAMES %>%
        gsub(".*s:([A-Za-z]*)_([A-Za-z]*);", "\\1 \\2", .),
    trnl_micro = as.character(microsynth_seqs, use.names = FALSE)
)

# Merge with our sequences
ref_seqs_combined <- left_join(ref_seqs_combined, microsynth_seqs,
by = join_by(species == species))

ref_seqs_combined %>% select(species, trnl_micro) %>% distinct() %>% print(n = 100)

## Compare availability of sequences with the botanical data
species_list <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>% 
    distinct(full_name) %>%
    filter(!grepl("species", full_name))

(!(species_list$full_name %in% ref_seqs_combined$full_name)) %>%
species_list$full_name[.]
ref_seqs_combined %>% select(species, full_name) %>% distinct() %>% print(n = 100)


## Compare with updated BOLD data base 
bold_data <- read_csv("data/external/reference_libraries/BOLD/bold_sequences_2025_03_12.csv")

ref_seqs_bold <- full_join(ref_seqs, bold_data, join_by(species == full_name)) %>%
select(-contains("trnl"), -species.y, its2_bold = its2) %>%
mutate(its2_match = mapply(its2_seq, its2_bold,
            FUN = function(its2_seq, its2_bold) {
                if (is.na(its2_seq) | is.na(its2_bold)) {
                    return(NA)
                }
                msa(c(its2_seq, its2_bold), type = "dna") %>%
                    consensusMatrix() %>%
                    (function(x) colSums(x == 2)) %>%
                    (function(x) round(100 * sum(x == 1) / length(x), 2))
            }
        ))
ref_seqs_bold %>% select(species, its2_count, its2_perc, its2_match) %>% 
arrange(species) %>% print(n = 100)
# Cardamine digitata sequences actually match, except BOLD barcode ist too long (exceeding base 322)
# Carxie marina strong match (8 base pair differences), except BOLD barcode is too short.
# Barcodes are quite different, bold barcode is too short. 

msa(c(ref_seqs_bold %>% 
    filter(grepl("Carex myosuroides", species)) %>% select(its2_seq, its2_bold) %>%
    unlist()),
    type = "dna") %>% print(., show="complete")


# Check updated Willerslev et al. from March 2025
willerslev <- read_csv("data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2024_subset.csv")


ref_seqs_willer <- full_join(ref_seqs, willerslev, join_by(species == full_name))

ref_seqs_willer <-ref_seqs_willer %>%
    mutate(
        trnl_match = mapply(trnl_seq, trnl,
            FUN = function(trnl_seq, trnl) {
                if (is.na(trnl_seq) | is.na(trnl)) {
                    return(NA)
                }
                msa(c(trnl_seq, trnl), type = "dna") %>%
                    consensusMatrix() %>%
                    (function(x) colSums(x == 2)) %>%
                    (function(x) round(100 * sum(x == 1) / length(x), 2))
            }
        ))
ref_seqs_willer %>% select(species, trnl_count, trnl_perc, trnl_match) %>% 
arrange(species) %>% print(n = 100)
# Exact match as previously. 
# - Three species have differences in local vs Willerslev et al. 2014
# Draba nivalis
# Sabulina rubella
# Bistorta vivipara

