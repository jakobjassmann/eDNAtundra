# Retrieve presence absence for each plot based on its2
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(dada2)
library(Biostrings)
library(ShortRead)
# library(pbapply)
# library(parallel)

## Data prep ---

# Get file names for all its2 reads
its2_fwd <- list.files(
    "data/sequences/its2/edna/trimmed",
    pattern = ".*_R1.*fastq",
    full.names = T
)
its2_rev <- list.files(
    "data/sequences/its2/edna/trimmed",
    pattern = ".*_R2.*fastq",
    full.names = T
)

# Set out path
out_path <- "data/sequences/its2/edna"

# Get sample ids
sample_ids <- gsub("([A-Za-z0-9]*-?[123P]?)-.*", "\\1", basename(its2_fwd))

# Create out dir for data QC
dir.create("data/sequences/its2/edna/dada2_qc/")

## Visual QC of a couple of samples
(fwd_quality_sample_plot <- plotQualityProfile(its2_fwd[17:18]))
ggsave(
    "data/sequences/its2/edna/dada2_qc/fwd_quality_sample_plot.png",
    fwd_quality_sample_plot
)
(rev_quality_sample_plot <- plotQualityProfile(its2_rev[17:18]))
ggsave(
    "data/sequences/its2/edna/dada2_qc/rev_quality_sample_plot.png",
    rev_quality_sample_plot
)
# These reads look a lot better then the trnl ones (overall)

## Filtering ---

# Set out file paths
its2_fwd_filtered <- file.path(out_path, "filtered", basename(its2_fwd))
its2_rev_filtered <- file.path(out_path, "filtered", basename(its2_rev))

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
    "data/sequences/its2/edna/dada2_qc/filterAndTrim_out.csv"
)

## Identify ASVs ----

# Learn errors
err_fwd <- learnErrors(its2_fwd_filtered, multithread = TRUE)
err_rev <- learnErrors(its2_rev_filtered, multithread = TRUE)

# Visulaise estimated error rates
fwd_errors_plot <- plotErrors(err_fwd, nominalQ = TRUE)
ggsave(
    "data/sequences/its2/edna/dada2_qc/fwd_errors_plot.png",
    fwd_errors_plot
)
rew_errors_plot <- plotErrors(err_rev, nominalQ = TRUE)
ggsave(
    "data/sequences/its2/edna/dada2_qc/rew_errors_plot.png",
    rew_errors_plot
)
# There are some  dips in the error rate functions at high quality scores.
# Looks good enough_

# Dereplicate
exists <- file.exists(its2_fwd_filtered)
its2_few_derep <- derepFastq(its2_fwd_filtered[exists], verbose = TRUE)
its2_rev_derep <- derepFastq(its2_rev_filtered[exists], verbose = TRUE)

# Add sample names
names(its2_few_derep) <- sample_ids[exists]
names(its2_rev_derep) <- sample_ids[exists]

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

# Reove chimeras
its2_seqtab_nochim <- removeBimeraDenovo(
    its2_seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
)
dim(its2_seqtab_nochim) # n samples n non-chimeric sequences
sum(its2_seqtab_nochim) / sum(its2_seqtab) # prop. non-chimeric sequences
table(nchar(getSequences(its2_seqtab_nochim))) # distribution of non-chimeric sequence lengths

# Save
dir.create("data/sequences/its2/edna/cleaned")
save(
    its2_seqtab_nochim,
    file = "data/sequences/its2/edna/cleaned/its2_seqtab_nochim.Rda"
)

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
rownames(track) <- sample_ids[exists]
head(track)
View(track)
write_csv(
    as.data.frame(track),
    "data/sequences/its2/edna/dada2_qc/its2_dada2_track.csv"
)

##  Assign species ----

its2_taxa <- assignSpecies(
    its2_seqtab_nochim,
    "data/sequences/reference/its2_reference_Assmann_et_al_2025.fasta",
    verbose = T
)
its2_taxa %>%
    data.frame() %>%
    remove_rownames() %>%
    filter(!is.na(Genus)) %>%
    distinct() %>%
    arrange(Genus, Species) %>%
    View()

# Save objects

save(its2_taxa, file = "data/sequences/its2/edna/cleaned/its2_taxa.Rda")
write_csv(
    track %>% data.frame(),
    "data/sequences/its2/edna/cleaned/its2_track.csv"
)

## Post processing to generate presence absence data frames ----

# Load sequneces and taxonomic assignments
load("data/sequences/its2/edna/cleaned/its2_seqtab_nochim.rda")
load("data/sequences/its2/edna/cleaned/its2_taxa.rda")

# Collapse data into a datafram of species presence / absence per plot
tax_counts <- right_join(
    # Prepare sequence ocunts
    t(its2_seqtab_nochim) %>%
        data.frame() %>%
        mutate(., sq = row.names(.)) %>%
        remove_rownames() %>%
        relocate(sq) %>%
        pivot_longer(
            .,
            cols = 2:ncol(.),
            names_to = "sample_id",
            values_to = "count"
        ) %>%
        mutate(
            sample_id = gsub("\\.", "-", sample_id),
            plot_id = case_when(
                grepl("S[0-9]{1}[SR]{1}[0-9]{1}", sample_id) ~ substr(
                    sample_id,
                    1,
                    4
                ),
                TRUE ~ NA
            ),
            subsample_id = case_when(
                grepl("S[0-9]{1}[SR]{1}[0-9]{1}", sample_id) ~ substr(
                    sample_id,
                    6,
                    6
                ),
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
write_csv(
    tax_counts,
    "data/sequences/its2/edna/cleaned/its2_taxa_in_sample.csv"
)

## Quick quality checks ----

# Check for NAs
tax_counts %>%
    filter(is.na(species))
# No NAs -> good!

# Check for n taxa per sample
tax_counts %>% group_by(sample_id) %>% tally() %>% pull(n) %>% hist()
read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    group_by(plot_id) %>%
    distinct(full_name) %>%
    tally() %>%
    pull(n) %>%
    hist()
# Distribution look more similar than what we saw for trnl

# Check counts
tax_counts %>% pull(count) %>% hist()

# How many counts are below 16?
tax_counts %>% filter(count < 16) %>% nrow()
# A bit more than for trnl, but still...
# This is fewer than expected -> it seems like the dada2 pipline filter worked
# rather well. Check later whether count is associated with presence.

# Overall looking good.

# End of file
