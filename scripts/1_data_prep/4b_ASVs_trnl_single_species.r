# Prepare trnl ASVs for the single species samples
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(readxl)
library(dada2)
library(Biostrings)
library(ShortRead)
library(msa)
library(pwalign)

# Look up table for species in sample IDs
sample_look_up <- read_xlsx("data/sequences/Assmann_et_al_-_metadata.xlsx") %>%
    filter(species != "NA") %>%
    select(sample = SampleID, species)

# Get file names for trnl reads
trnl_fwd <- list.files(
    "data/sequences/trnl/single_species/trimmed",
    pattern = ".*_R1.*fastq",
    full.names = T
)
trnl_rev <- list.files(
    "data/sequences/trnl/single_species/trimmed",
    pattern = ".*_R2.*fastq",
    full.names = T
)

# Set out file paths
out_path <- "data/sequences/trnl/single_species"
trnl_fwd_filtered <- file.path(out_path, "filtered", basename(trnl_fwd))
trnl_rev_filtered <- file.path(out_path, "filtered", basename(trnl_rev))
dir.create(file.path(out_path, "filtered"), showWarnings = FALSE)

# Create out dir for data QC
dir.create("data/sequences/trnl/single_species/dada2_qc/")

## Visual QC of a couple of samples
(fwd_quality_sample_plot <- plotQualityProfile(trnl_fwd[17:18]))
ggsave(
    "data/sequences/trnl/single_species/dada2_qc/fwd_quality_sample_plot.png",
    fwd_quality_sample_plot
)
(rev_quality_sample_plot <- plotQualityProfile(trnl_rev[17:18]))
ggsave(
    "data/sequences/trnl/single_species/dada2_qc/rev_quality_sample_plot.png",
    rev_quality_sample_plot
)
# Short reads as expected, but looks not too bad.

# Filter and trim based on quality
out <- filterAndTrim(
    fwd = trnl_fwd,
    filt = trnl_fwd_filtered,
    rev = trnl_rev,
    filt.rev = trnl_rev_filtered,
    maxN = 0, # Discard sequences with 1 N or more (missing read)
    maxEE = c(2, 2), # After truncation, reads with higher than 2 "expected errors" will be discarded
    truncQ = 2, # Truncate  at the first instance of a quality score less than or equal to 2
    minLen = 10,
    rm.phix = TRUE, # RUE Discard reads that match against the phiX genome (ILLUMINA control)
    compress = TRUE,
    multithread = TRUE
)
# Check results
out
write_csv(
    as.data.frame(out),
    "data/sequences/trnl/single_species/dada2_qc/filterAndTrim_out.csv"
)

## Identify ASVs ----

# Learn errors
err_fwd <- learnErrors(trnl_fwd_filtered, multithread = TRUE)
# 40168026 total bases in 657299 reads from 47 samples will be used for learning the error rates.
err_rev <- learnErrors(trnl_rev_filtered, multithread = TRUE)
# 40140036 total bases in 657299 reads from 47 samples will be used for learning the error rates.

# Visualise the estimated error rates
fwd_errors_plot <- plotErrors(err_fwd, nominalQ = TRUE)
ggsave(
    "data/sequences/trnl/single_species/dada2_qc/fwd_errors_plot.png",
    fwd_errors_plot
)
rew_errors_plot <- plotErrors(err_rev, nominalQ = TRUE)
ggsave(
    "data/sequences/trnl/single_species/dada2_qc/rew_errors_plot.png",
    rew_errors_plot
)
# It's perhaps not super ideal, but looks good enough?!

# Dereplicate
exists <- file.exists(trnl_fwd_filtered)
trnl_few_derep <- derepFastq(trnl_fwd_filtered[exists], verbose = TRUE)
trnl_rev_derep <- derepFastq(trnl_rev_filtered[exists], verbose = TRUE)

# Run the data core algorithm to infer ASV
trnl_fwd_dada <- dada(trnl_few_derep, err = err_fwd, multithread = TRUE)
trnl_rev_dada <- dada(trnl_rev_derep, err = err_rev, multithread = TRUE)

# Merge paired-end reads
trnl_merged <- mergePairs(
    trnl_fwd_dada,
    trnl_few_derep,
    trnl_rev_dada,
    trnl_rev_derep,
    verbose = TRUE
)

# Convert to sequence table of ASVs
trnl_seqtab <- makeSequenceTable(trnl_merged)
dim(trnl_seqtab) # n samples n sequences
table(nchar(getSequences(trnl_seqtab))) # distribution of sequence lengths

# Remove chimeras
trnl_seqtab_nochim <- removeBimeraDenovo(
    trnl_seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
)
dim(trnl_seqtab_nochim) # n samples n non-chimeric sequences
sum(trnl_seqtab_nochim) / sum(trnl_seqtab) # prop. non-chimeric sequences
table(nchar(getSequences(trnl_seqtab_nochim))) # distribution of non-chimeric sequence lengths

# Track sequences through pipeline (from tutorial)
getN <- function(x) sum(getUniques(x))
track <- cbind(
    out,
    sapply(trnl_fwd_dada, getN),
    sapply(trnl_rev_dada, getN),
    sapply(trnl_merged, getN),
    rowSums(trnl_seqtab_nochim)
)
colnames(track) <- c(
    "input",
    "filtered",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim"
)
head(track)
View(track)
write_csv(
    as.data.frame(track),
    "data/sequences/trnl/single_species/dada2_qc/trnl_dada2_track.csv"
)

# Save cleaned sequences
dir.create("data/sequences/trnl/single_species/cleaned")
save(
    trnl_seqtab_nochim,
    file = "data/sequences/trnl/single_species/cleaned/trnl_seqtab_nochim.Rda"
)
# Load filtered table (if restarting)
# load("data/sequences/trnl/single_species/cleaned/trnl_seqtab_nochim.Rda")

# Clean up ASV table and filter out low abundance and low length ASVs ----
trnl_seqtab_nochim_df <- trnl_seqtab_nochim %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "sequence", values_to = "count") %>%
    group_by(sample) %>%
    # calculate percentage of reads for each sequence
    mutate(perc = 100 * count / sum(count)) %>%
    # Add sequence length
    mutate(sequence_length = nchar(sequence)) %>%
    ungroup() %>%
    # Filter ASVs with zero counts and add sequence length
    filter(count > 0) %>%
    # Filter ASVs with less than 10% of reads (arbitrary threshold)
    filter(perc >= 10) %>%
    # filter ASVs with less than 10 BP (to short for a complete trnl sequence)
    filter(sequence_length >= 10) %>%
    arrange(sample, desc(count)) %>%
    # Shorten sample names to match with look up table
    mutate(sample = substr(sample, 1, 2))
print(trnl_seqtab_nochim_df, n = 200)

# Helper function to calculate pairwise identities using pid from Biostrings
pairwise_pid <- function(sequence_df) {
    # Check whether data frame has only one row, if so return itself with NA for pid value
    if (nrow(sequence_df) == 1) {
        sequence_df$pid_value <- 100
        sequence_df$nmis_value <- 0
        sequence_df$length_value <- 0
        return(sequence_df)
    }

    # Calculate pairwise identities for all sequences relative to the most common sequence
    pid_value <- numeric()
    nmis_value <- numeric()
    length_value <- numeric()

    # Get sequences
    sequences <- sequence_df$sequence
    for (i in 1:length(sequences)) {
        seq1 <- sequences[1]
        seq2 <- sequences[i]
        alignment <- pairwiseAlignment(seq1, seq2, type = "global")
        # Percentage difference
        pid_value[i] <- pid(alignment)
        # N BP difference
        nmis_value[i] <- nmismatch(alignment)
        # difference in length
        length_value[i] <- abs(nchar(seq1) - nchar(seq2))
    }

    # update data frame
    sequence_df$pid_value <- pid_value
    sequence_df$nmis_value <- nmis_value
    sequence_df$length_value <- length_value

    return(sequence_df)
}

# Apply helper function to calculate pairwise identity with most common sequene in sample
trnl_seqtab_nochim_df <- trnl_seqtab_nochim_df %>%
    group_by(sample) %>%
    group_split() %>%
    map(pairwise_pid) %>%
    bind_rows() %>%
    ungroup() %>%
    # Add unique id
    mutate(., ID = 1:nrow(.))
trnl_seqtab_nochim_df %>% select(-sequence) %>% print(n = 200)

# Export sequences as FASTA to query with BOLD Plant Database
trnl_to_check_seqs <- DNAStringSet(trnl_seqtab_nochim_df$sequence)
names(trnl_to_check_seqs) <- trnl_seqtab_nochim_df$ID
writeXStringSet(
    trnl_to_check_seqs,
    "data/sequences/trnl/single_species/trnl_to_check.fasta"
)

# Prepare Willerselv db for matching
willerslev <- readDNAStringSet(
    "data/external/reference_libraries/Willerslev_et_al_2014/references.arctic-gh.fasta"
)

# phrase names
names(willerslev) <- names(willerslev) %>%
    gsub(".*scientific_name=(.*); taxonomic_rank.*", "\\1", .) %>%
    paste("ID", ., sep = "_") %>%
    gsub("GROUP_[0-9]", "species", .) %>%
    gsub("species_species", "species", .)
# Write FASTA
writeXStringSet(
    willerslev,
    filepath = "data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2014.fasta"
)
# Create BLAST database in terminal (external command)
# makeblastdb -in data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2014.fasta -dbtype nucl -out data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2014_db

# Query in terminal (external command)
# blastn -query data/sequences/trnl/single_species/trnl_to_check.fasta -db data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2014_db -outfmt '6 qseqid sseqid pident length' -out data/sequences/trnl/single_species/Willerslev_trnl_to_check_matches.txt

# Read and format the BLAST results
trnl_Willerslev_results <- read_tsv(
    "data/sequences/trnl/single_species/Willerslev_trnl_to_check_matches.txt",
    col_names = c("query_id", "species", "identity", "length")
)
print(trnl_Willerslev_results, n = 250)

# Filter for the best match (highest identity)
trnl_Willerslev_results <- trnl_Willerslev_results %>%
    group_by(query_id) %>%
    arrange(desc(identity)) %>%
    dplyr::slice(1) %>%
    mutate(species = gsub("_", " ", species)) %>% # Replace underscores with spaces in species names
    mutate(species = gsub("ID ", "", species)) %>%
    mutate(
        Willerslev_full_name = species,
        Willerslev_genus = sapply(str_extract_all(species, "\\S+"), `[`, 1),
        Willerslev_species = sapply(str_extract_all(species, "\\S+"), `[`, 2)
    ) %>%
    dplyr::rename(
        ID = query_id,
        Willerslev_identity = identity,
        Willerslev_length = length
    ) %>%
    select(-species)

# Load ITS2 confirmed IDs
its2_ids <- read_csv(
    "data/sequences/its2/single_species/its2_annotated_manual_qc.csv"
) %>%
    filter(qc_include == "yes") %>%
    select(sample, its2_genus = genus, its2_species = species) %>%
    distinct()

# merge with trnl sequences
trnl_annotated <- trnl_seqtab_nochim_df %>%
    left_join(., sample_look_up) %>%
    mutate(full_name = species) %>%
    select(-species) %>%
    mutate(
        genus = sapply(str_extract_all(full_name, "\\S+"), `[`, 1),
        species = sapply(str_extract_all(full_name, "\\S+"), `[`, 2)
    ) %>%
    # add BOLD matches
    full_join(trnl_Willerslev_results, by = c("ID" = "ID")) %>%
    left_join(its2_ids, by = c("sample" = "sample")) %>%
    # Annotate matches
    mutate(
        Willerslev_match_species = case_when(
            genus == Willerslev_genus ~ "yes",
            TRUE ~ "no"
        ),
        its2_match_species = case_when(
            genus == its2_genus ~ "yes",
            TRUE ~ "no"
        )
    )

# Add
# Export to CSV and manually annotate
write_excel_csv(
    trnl_annotated,
    "data/sequences/trnl/single_species/trnl_annotated.csv"
)

# Load manualy quality controlled annotated data and clean up
trnl_annotated <- read_csv(
    "data/sequences/trnl/single_species/trnl_annotated_manual_qc.csv"
) %>%
    filter(qc_include == "yes") %>%
    mutate(source = paste("Assmann et al. sample", sample)) %>%
    select(genus, species, full_name, sequence, sequence_length, source)
View(trnl_annotated)

# Write out
write_excel_csv(
    trnl_annotated,
    "data/sequences/trnl/single_species/trnl_single_species_cleaned.csv"
)

### EOF
###############################################################################
