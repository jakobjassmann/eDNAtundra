# Prepare ITS2 ASVs for the single species samples
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

# Get file names for its2 reads
its2_fwd <- list.files(
    "data/sequences/its2/single_species/trimmed",
    pattern = ".*_R1.*fastq",
    full.names = T
)
its2_rev <- list.files(
    "data/sequences/its2/single_species/trimmed",
    pattern = ".*_R2.*fastq",
    full.names = T
)

# Set out file paths
out_path <- "data/sequences/its2/single_species"
its2_fwd_filtered <- file.path(out_path, "filtered", basename(its2_fwd))
its2_rev_filtered <- file.path(out_path, "filtered", basename(its2_rev))
dir.create(file.path(out_path, "filtered"), showWarnings = FALSE)

# Create out dir for data QC
dir.create("data/sequences/its2/single_species/dada2_qc/")

## Visual QC of a couple of samples
(fwd_quality_sample_plot <- plotQualityProfile(its2_fwd[17:18]))
ggsave(
    "data/sequences/its2/single_species/dada2_qc/fwd_quality_sample_plot.png",
    fwd_quality_sample_plot
)
(rev_quality_sample_plot <- plotQualityProfile(its2_rev[17:18]))
ggsave(
    "data/sequences/its2/single_species/dada2_qc/rev_quality_sample_plot.png",
    rev_quality_sample_plot
)
# Some samples have failed to amplify (C3) others look good (C4)

# Filter and trim based on quality
out <- filterAndTrim(
    fwd = its2_fwd,
    filt = its2_fwd_filtered,
    rev = its2_rev,
    filt.rev = its2_rev_filtered,
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
    "data/sequences/its2/single_species/dada2_qc/filterAndTrim_out.csv"
)

## Identify ASVs ----

# Learn errors
err_fwd <- learnErrors(its2_fwd_filtered, multithread = TRUE)
# 104739891 total bases in 465587 reads from 20 samples will be used for learning the error rates.
err_rev <- learnErrors(its2_rev_filtered, multithread = TRUE)
# 104742748 total bases in 465587 reads from 20 samples will be used for learning the error rates.

# Visualise the estimated error rates
fwd_errors_plot <- plotErrors(err_fwd, nominalQ = TRUE)
ggsave(
    "data/sequences/its2/single_species/dada2_qc/fwd_errors_plot.png",
    fwd_errors_plot
)
rew_errors_plot <- plotErrors(err_rev, nominalQ = TRUE)
ggsave(
    "data/sequences/its2/single_species/dada2_qc/rew_errors_plot.png",
    rew_errors_plot
)
# It's perhaps not super ideal, but looks good enough?!

# Dereplicate
exists <- file.exists(its2_fwd_filtered)
its2_few_derep <- derepFastq(its2_fwd_filtered[exists], verbose = TRUE)
its2_rev_derep <- derepFastq(its2_rev_filtered[exists], verbose = TRUE)

# Run the data core algorithm to infer ASV
its2_fwd_dada <- dada(its2_few_derep, err = err_fwd, multithread = TRUE)
its2_rev_dada <- dada(its2_rev_derep, err = err_rev, multithread = TRUE)

# Merge paired-end reads
its2_merged <- mergePairs(
    its2_fwd_dada,
    its2_few_derep,
    its2_rev_dada,
    its2_rev_derep,
    verbose = TRUE
)

# Convert to sequence table of ASVs
its2_seqtab <- makeSequenceTable(its2_merged)
dim(its2_seqtab) # n samples n sequences
table(nchar(getSequences(its2_seqtab))) # distribution of sequence lengths

# Remove chimeras
its2_seqtab_nochim <- removeBimeraDenovo(
    its2_seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
)
dim(its2_seqtab_nochim) # n samples n non-chimeric sequences
sum(its2_seqtab_nochim) / sum(its2_seqtab) # prop. non-chimeric sequences
table(nchar(getSequences(its2_seqtab_nochim))) # distribution of non-chimeric sequence lengths

# Track sequences through pipeline (from tutorial)
getN <- function(x) sum(getUniques(x))
track <- cbind(
    out,
    sapply(its2_fwd_dada, getN),
    sapply(its2_rev_dada, getN),
    sapply(its2_merged, getN),
    rowSums(its2_seqtab_nochim)
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
    "data/sequences/its2/single_species/dada2_qc/its2_dada2_track.csv"
)

# Save cleaned sequences
dir.create("data/sequences/its2/single_species/cleaned")
save(
    its2_seqtab_nochim,
    file = "data/sequences/its2/single_species/cleaned/its2_seqtab_nochim.Rda"
)
# Load filtered table (if restarting)
# load("data/sequences/its2/single_species/cleaned/its2_seqtab_nochim.Rda")

# Clean up ASV table and filter out low abundance and low length ASVs ----
its2_seqtab_nochim_df <- its2_seqtab_nochim %>%
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
    # filter ASVs with less than 160 BP (to short for a complete ITS2 sequence)
    filter(sequence_length >= 160) %>%
    arrange(sample, desc(count)) %>%
    # Shorten sample names to match with look up table
    mutate(sample = substr(sample, 1, 2))
print(its2_seqtab_nochim_df, n = 200)

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
its2_seqtab_nochim_df <- its2_seqtab_nochim_df %>%
    group_by(sample) %>%
    group_split() %>%
    map(pairwise_pid) %>%
    bind_rows() %>%
    ungroup() %>%
    # Add unique id
    mutate(., ID = 1:nrow(.))
its2_seqtab_nochim_df %>% select(-sequence) %>% print(n = 200)

# Export sequences as FASTA to query with BOLD Plant Database
its2_to_check_seqs <- DNAStringSet(its2_seqtab_nochim_df$sequence)
names(its2_to_check_seqs) <- its2_seqtab_nochim_df$ID
writeXStringSet(
    its2_to_check_seqs,
    "data/sequences/its2/single_species/its2_to_check.fasta"
)

# Load FASTA queried on BOLD ("plant public", "genus and species") on 2026-06-10 at 09:46
its2_BOLD_id <- read_csv(
    "data/sequences/its2/single_species/Barcode ID - BOLD - db_public.plantsmi_0.9mo_100maxh_50order_3-e620f0f9bdca403fab732742ba3f108b.csv"
) %>%
    group_by(`Query ID`) %>%
    filter(!is.na(Species)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(
        ID = `Query ID`,
        PID = `PID [BIN]`,
        species = Species,
        bold_perc_match = `ID%`
    ) %>%
    mutate(
        bold_full_name = species,
        bold_genus = sapply(str_extract_all(species, "\\S+"), `[`, 1),
        bold_species = sapply(str_extract_all(species, "\\S+"), `[`, 2)
    ) %>%
    select(ID, PID, bold_genus, bold_species, bold_full_name, bold_perc_match)

# Query in terminal (external command)
# blastn -query data/sequences/its2/single_species/its2_to_check.fasta -db data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN_its2_db -outfmt '6 qseqid sseqid pident length' -out data/sequences/its2/single_species/DS-ZACKPLAN_its2_to_check_matches.txt

# Read and format the BLAST results
its2_DS_ZACKPLAN_results <- read_tsv(
    "data/sequences/its2/single_species/DS-ZACKPLAN_its2_to_check_matches.txt",
    col_names = c("query_id", "species", "identity", "length")
)

# Filter for the best match (highest identity)
its2_DS_ZACKPLAN_best_match <- its2_DS_ZACKPLAN_results %>%
    group_by(query_id) %>%
    arrange(desc(identity)) %>%
    dplyr::slice(1) %>%
    mutate(species = gsub("_", " ", species)) %>% # Replace underscores with spaces in species names
    mutate(species = gsub("ID ", "", species)) %>%
    mutate(
        DS_ZACKPLAN_full_name = species,
        DS_ZACKPLAN_genus = sapply(str_extract_all(species, "\\S+"), `[`, 1),
        DS_ZACKPLAN_species = sapply(str_extract_all(species, "\\S+"), `[`, 2)
    ) %>%
    dplyr::rename(
        ID = query_id,
        DS_ZACKPLAN_identity = identity,
        DS_ZACKPLAN_length = length
    ) %>%
    select(-species) %>%
    # remove matches that had less than 140 BP for the aligned sequence
    filter(DS_ZACKPLAN_length >= 140)

# merge with ITS2 sequences
its2_annotated <- its2_seqtab_nochim_df %>%
    left_join(., sample_look_up) %>%
    mutate(
        full_name = species,
        genus = sapply(str_extract_all(species, "\\S+"), `[`, 1),
        species = sapply(str_extract_all(species, "\\S+"), `[`, 2)
    ) %>%
    # add BOLD matches
    full_join(its2_BOLD_id, by = c("ID" = "ID")) %>%
    # add DS-ZACKPLAN matches
    left_join(its2_DS_ZACKPLAN_best_match, by = c("ID" = "ID")) %>%
    # Annotate matches
    mutate(
        bold_match_species = case_when(
            genus == bold_genus & species == bold_species ~ "yes",
            TRUE ~ "no"
        ),
        DS_ZACKPLAN_match_species = case_when(
            genus == DS_ZACKPLAN_genus & species == DS_ZACKPLAN_species ~ "yes",
            TRUE ~ "no"
        )
    )

# Export to CSV and manually annotate
write_excel_csv(
    its2_annotated,
    "data/sequences/its2/single_species/its2_annotated.csv"
)

# Load manualy quality controlled annotated data and clean up
its2_annotated <- read_csv(
    "data/sequences/its2/single_species/its2_annotated_manual_qc.csv"
) %>%
    filter(qc_include == "yes") %>%
    mutate(source = paste("Assmann et al. sample", sample)) %>%
    select(genus, species, full_name, sequence, sequence_length, source)
View(its2_annotated)

# Write out
write_csv(
    its2_annotated,
    "data/sequences/its2/single_species/its2_single_species_cleaned.csv"
)

### EOF
###############################################################################

### Old code (delete once ready!)

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
        trnl_perc = 100 *
            round(
                trnl_derep[[i]]@.Data[[1]][1] / sum(trnl_derep[[1]]@.Data[[1]]),
                2
            )
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
        its2_perc = 100 *
            round(
                its2_derep[[i]]@.Data[[1]][1] / sum(its2_derep[[1]]@.Data[[1]]),
                2
            )
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

# Check for duplicates / uncertain ids
ref_seqs %>% filter(grepl("cf", species))

# Only Sabulina rubella is duplicated
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
ref_seqs %>%
    select(-sample, -trnl_count, -trnl_perc, -its2_count, -its2_perc) %>%
    arrange(species) %>%
    print(n = 100)

# Write oute reference
write_csv(ref_seqs, "data/sequences/reference/Assmann_et_al_CB_2023.csv")
#ref_seqs <- read_csv("data/sequences/reference/Assmann_et_al_CB_2023.csv")

## Compare sequences to previously obtained references
online_ref_seqs <- read_csv(
    "data/sequences/reference/Assmann_et_al_eDNA-tundra_species_list_with_markers_updated_27_Nov_23.csv"
) %>%
    select(
        species,
        synonym,
        trnl_seq_online = trnL_seq,
        its2_seq_online = its2_seq
    )

# Look up and add sample ID
online_ref_seqs <- online_ref_seqs %>%
    mutate(
        sample = sapply(
            species,
            function(x) {
                online_genus <- gsub("([A-Za-z]*) ([A-Za-z]*)", "\\1", x)
                online_species <- gsub("([A-Za-z]*) ([A-Za-z]*)", "\\2", x)
                sample_id <- ref_seqs %>%
                    filter(
                        grepl(online_genus, ref_seqs$species) &
                            grepl(online_species, ref_seqs$species)
                    ) %>%
                    pull(sample)
                if (length(sample_id) == 0) {
                    return("None")
                }
                return(sample_id)
            },
            USE.NAMES = FALSE
        )
    )

# Merge online references with sequenced
ref_seqs_combined <- full_join(
    ref_seqs %>% dplyr::rename(full_name = species),
    online_ref_seqs,
    by = join_by(sample == sample)
) %>%
    select(
        sample,
        species,
        full_name,
        trnl_seq,
        trnl_seq_online,
        its2_seq,
        its2_seq_online
    )

# Check sequence availability
ref_seqs_combined %>%
    select(full_name, trnl_seq, trnl_seq_online, its2_seq, its2_seq_online) %>%
    print(n = 100)

# Compare sample sequences and online sequences
ref_seqs_combined <- ref_seqs_combined %>%
    mutate(
        trnl_match = mapply(
            trnl_seq,
            trnl_seq_online,
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
        its2_match = mapply(
            its2_seq,
            its2_seq_online,
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
ref_seqs_combined %>%
    filter(
        species %in%
            c("Carex marina", "Cassiope tetragona", "Pedicularis lanata")
    ) %>%
    select(species, trnl_match)
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
microsynth_seqs <- readDNAStringSet(
    "C:/Users/jakob/repositories/eDNAtundra/data/sequences/trnl/trnL_eDNA_run1_site_vs_single_species_and_arctic_gh/customdb2_single_species_and_arctic_gh.fasta"
)

# Convert to same format aus our database
microsynth_seqs <- data.frame(
    species = microsynth_seqs@ranges@NAMES %>%
        gsub(".*s:([A-Za-z]*)_([A-Za-z]*);", "\\1 \\2", .),
    trnl_micro = as.character(microsynth_seqs, use.names = FALSE)
)

# Merge with our sequences
ref_seqs_combined <- left_join(
    ref_seqs_combined,
    microsynth_seqs,
    by = join_by(species == species)
)

ref_seqs_combined %>%
    select(species, trnl_micro) %>%
    distinct() %>%
    print(n = 100)

## Compare availability of sequences with the botanical data
species_list <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    distinct(full_name) %>%
    filter(!grepl("species", full_name))

(!(species_list$full_name %in% ref_seqs_combined$full_name)) %>%
    species_list$full_name[.]
ref_seqs_combined %>%
    select(species, full_name) %>%
    distinct() %>%
    print(n = 100)


## Compare with updated BOLD data base
bold_data <- read_csv(
    "data/external/reference_libraries/BOLD/bold_sequences_2025_03_12.csv"
)

ref_seqs_bold <- full_join(
    ref_seqs,
    bold_data,
    join_by(species == full_name)
) %>%
    select(-contains("trnl"), -species.y, its2_bold = its2) %>%
    mutate(
        its2_match = mapply(
            its2_seq,
            its2_bold,
            FUN = function(its2_seq, its2_bold) {
                if (is.na(its2_seq) | is.na(its2_bold)) {
                    return(NA)
                }
                msa(c(its2_seq, its2_bold), type = "dna") %>%
                    consensusMatrix() %>%
                    (function(x) colSums(x == 2)) %>%
                    (function(x) round(100 * sum(x == 1) / length(x), 2))
            }
        )
    )
ref_seqs_bold %>%
    select(species, its2_count, its2_perc, its2_match) %>%
    arrange(species) %>%
    print(n = 100)
# Cardamine digitata sequences actually match, except BOLD barcode ist too long (exceeding base 322)
# Carxie marina strong match (8 base pair differences), except BOLD barcode is too short.
# Barcodes are quite different, bold barcode is too short.

msa(
    c(
        ref_seqs_bold %>%
            filter(grepl("Carex myosuroides", species)) %>%
            select(its2_seq, its2_bold) %>%
            unlist()
    ),
    type = "dna"
) %>%
    print(., show = "complete")


# Check updated Willerslev et al. from March 2025
willerslev <- read_csv(
    "data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2024_subset.csv"
)


ref_seqs_willer <- full_join(
    ref_seqs,
    willerslev,
    join_by(species == full_name)
)

ref_seqs_willer <- ref_seqs_willer %>%
    mutate(
        trnl_match = mapply(trnl_seq, trnl, FUN = function(trnl_seq, trnl) {
            if (is.na(trnl_seq) | is.na(trnl)) {
                return(NA)
            }
            msa(c(trnl_seq, trnl), type = "dna") %>%
                consensusMatrix() %>%
                (function(x) colSums(x == 2)) %>%
                (function(x) round(100 * sum(x == 1) / length(x), 2))
        })
    )
ref_seqs_willer %>%
    select(species, trnl_count, trnl_perc, trnl_match) %>%
    arrange(species) %>%
    print(n = 100)
# Exact match as previously.
# - Three species have differences in local vs Willerslev et al. 2014
# Draba nivalis
# Sabulina rubella
# Bistorta vivipara
