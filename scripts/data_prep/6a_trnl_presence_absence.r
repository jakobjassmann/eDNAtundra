# Retrieve presence absence for each plot based on trnl
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(dada2)
library(Biostrings)
library(ShortRead)
# library(pbapply)
# library(parallel)

## Data prep ---

# Get file names for all trnl reads
trnl_fwd <- list.files("data/sequences/trnl/edna/trimmed",
    pattern = ".*_R1.*fastq",
    full.names = T
)
trnl_rev <- list.files("data/sequences/trnl/edna/trimmed",
    pattern = ".*_R2.*fastq",
    full.names = T
)

# Set out path
out_path <- "data/sequences/trnl/edna"
dir.create(out_path)

# Get sample ids
sample_ids <- gsub("([A-Za-z0-9]*-?[123P]?)-.*", "\\1", basename(trnl_fwd))

## Visual QC of a couple of samples
plotQualityProfile(trnl_fwd[17:18])
plotQualityProfile(trnl_rev[17:18])

## Filtering ---

# Set out file paths
trnl_fwd_filtered <- file.path(out_path, "filtered", basename(trnl_fwd))
trnl_rev_filtered <- file.path(out_path, "filtered", basename(trnl_rev))

# Filter and trim based on quality
out <- filterAndTrim(
    fwd = trnl_fwd, filt = trnl_fwd_filtered,
    rev = trnl_rev, filt.rev = trnl_rev_filtered,
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

## Identify ASVs ----

# Learn errors
err_fwd <- learnErrors(trnl_fwd_filtered, multithread = TRUE)
err_rev <- learnErrors(trnl_rev_filtered, multithread = TRUE)

# Visulaise estimated error rates
plotErrors(err_fwd, nominalQ = TRUE)
plotErrors(err_rev, nominalQ = TRUE)

# There are some strange dips in the error rate functions at high quality scores.
# Not sure what to do with it for now, proceeding and then reading up on it later.

# Dereplicate
exists <- file.exists(trnl_fwd_filtered)
trnl_few_derep <- derepFastq(trnl_fwd_filtered[exists], verbose = TRUE)
trnl_rev_derep <- derepFastq(trnl_rev_filtered[exists], verbose = TRUE)

# Add sample names
names(trnl_few_derep) <- sample_ids[exists]
names(trnl_rev_derep) <- sample_ids[exists]

# Run the data core algorithm to infer ASV
trnl_fwd_dada <- dada(trnl_few_derep, err = err_fwd, multithread = TRUE)
trnl_rev_dada <- dada(trnl_rev_derep, err = err_rev, multithread = TRUE)

# Merge paired-end reads
trnl_merged <- mergePairs(trnl_fwd_dada, trnl_few_derep,
    trnl_rev_dada, trnl_rev_derep,
    verbose = TRUE
)

# Convert to sequence table of ASVs
trnl_seqtab <- makeSequenceTable(trnl_merged)
dim(trnl_seqtab) # n samples n sequences
table(nchar(getSequences(trnl_seqtab))) # distribution of sequence lengths

# Reove chimeras
trnl_seqtab_nochim <- removeBimeraDenovo(trnl_seqtab,
    method = "consensus",
    multithread = TRUE, verbose = TRUE
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
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_ids[exists]
head(track)
View(track)

##  Assign species ----

trnl_taxa <- assignSpecies(trnl_seqtab_nochim,
    "data/sequences/reference/trnl_reference_Assmann_et_al_2025.fasta",
    verbose = T
)
trnl_taxa %>%
    data.frame() %>%
    remove_rownames() %>%
    filter(!is.na(Genus)) %>%
    View()

# Save objects
dir.create("data/sequences/trnl/edna/cleaned")
save(trnl_seqtab_nochim, file = "data/sequences/trnl/edna/cleaned/trnl_seqtab_nochim.Rda")
save(trnl_taxa, file = "data/sequences/trnl/edna/cleaned/trnl_taxa.Rda")
write_csv(track %>% data.frame(), "data/sequences/trnl/edna/cleaned/trnl_track.csv")

## Post processing to generate presence absence data frames

# Load sequneces and taxonomic assignments
load("data/sequences/trnl/edna/cleaned/trnl_seqtab_nochim.Rda")
load("data/sequences/trnl/edna/cleaned/trnl_taxa.Rda")

# Collapse data into a datafram of species presence / absence per plot
tax_counts <- right_join(
    # Prepare sequence ocunts
    t(trnl_seqtab_nochim) %>%
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
    data.frame(trnl_taxa) %>%
        mutate(., sq = row.names(.)) %>%
        remove_rownames() %>%
        filter(!is.na(Genus))
) %>%
    relocate(sample_id, plot_id, subsample_id, Genus, Species) %>%
    dplyr::rename(genus = Genus, species = Species) %>%
    arrange(sample_id)

# Write out data table
write_csv(tax_counts, "data/sequences/trnl/edna/cleaned/trnl_taxa_in_sample.csv")

## Quick quality checks ----

# Why are some species NA?
tax_counts %>%
filter(is.na(species)) %>%
dplyr::slice(2) %>%
pull(sq)

# Compare to botanical obs
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv")

botanical_obs %>%
    filter(plot_id %in% unique(substr(species_detected_df_filtered$plot_id, 1, 4))) %>%
    group_by(plot_id) %>%
    summarise(n_spp = n())

species_detected_df_filtered %>%
    mutate(plot_id = substr(plot_id, 1, 4)) %>%
    distinct(plot_id, Genus, Species) %>%
    group_by(plot_id) %>%
    summarise(n_spp = n())

test_plot <- unique(substr(species_detected_df_filtered$plot_id, 1, 4))[3]
botanical_obs %>%
    filter(plot_id %in% test_plot) %>%
    arrange(plot_id, full_name) %>%
    select(plot_id, full_name)
species_detected_df_filtered %>%
    mutate(plot_id = substr(plot_id, 1, 4)) %>%
    filter(plot_id %in% test_plot) %>%
    group_by(plot_id, Genus, Species) %>%
    summarise(total_count = sum(Count))
arrange(Genus, Species)

# End of file