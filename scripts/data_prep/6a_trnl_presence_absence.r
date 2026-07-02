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
trnl_fwd <- list.files(
    "data/sequences/trnl/edna/trimmed",
    pattern = ".*_R1.*fastq",
    full.names = T
)
trnl_rev <- list.files(
    "data/sequences/trnl/edna/trimmed",
    pattern = ".*_R2.*fastq",
    full.names = T
)

# Set out path
out_path <- "data/sequences/trnl/edna"
dir.create(out_path)

# Get sample ids
sample_ids <- gsub("([A-Za-z0-9]*-?[123P]?)-.*", "\\1", basename(trnl_fwd))

# Create out dir for data QC
dir.create("data/sequences/trnl/edna/dada2_qc/")

## Visual QC of a couple of samples
(fwd_quality_sample_plot <- plotQualityProfile(trnl_fwd[17:18]))
ggsave(
    "data/sequences/trnl/edna/dada2_qc/fwd_quality_sample_plot.png",
    fwd_quality_sample_plot
)
(rev_quality_sample_plot <- plotQualityProfile(trnl_rev[17:18]))
ggsave(
    "data/sequences/trnl/edna/dada2_qc/rev_quality_sample_plot.png",
    rev_quality_sample_plot
)
# Looking good.

## Filtering ---

# Set out file paths
trnl_fwd_filtered <- file.path(out_path, "filtered", basename(trnl_fwd))
trnl_rev_filtered <- file.path(out_path, "filtered", basename(trnl_rev))

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
    "data/sequences/trnl/edna/dada2_qc/filterAndTrim_out.csv"
)


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

# Save
dir.create("data/sequences/trnl/edna/cleaned")
save(
    trnl_seqtab_nochim,
    file = "data/sequences/trnl/edna/cleaned/trnl_seqtab_nochim.Rda"
)

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
rownames(track) <- sample_ids[exists]
head(track)
View(track)
write_csv(
    as.data.frame(track),
    "data/sequences/trnl/edna/dada2_qc/trnl_dada2_track.csv"
)

##  Assign species ----

# Build data frame of ASVs and filter our short sequences
trnl_asvs <- tibble(
    id = 1:length(dimnames(trnl_seqtab_nochim)[[2]]),
    sequence = dimnames(trnl_seqtab_nochim)[[2]],
    width = nchar(sequence)
) %>%
    ### Denoising Step 1 !!!
    # Set all sequences shorter than 10 to NA (min widht in ref db is 11)
    mutate(
        sequence = case_when(width < 10 ~ NA, TRUE ~ sequence)
    )

# Export to a fasta
sequences <- trnl_asvs %>%
    na.omit() %>%
    pull(sequence) %>%
    DNAStringSet()
names(sequences) <- trnl_asvs %>%
    na.omit() %>%
    pull(id)
writeXStringSet(
    sequences,
    filepath = "data/sequences/trnl/edna/cleaned/trnl_asvs.fasta"
)

# Blast against reference databse (on linux command line)
blastn -query data/sequences/trnl/edna/cleaned/trnl_asvs.fasta \
       -db data/sequences/trnl/reference/Assmann_et_al_trnl_ref_db \
       -out data/sequences/trnl/edna/cleaned/trnl_blast_res.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95 \
       -word_size 7 \
       -dust no \
       -evalue 1000 \
       -qcov_hsp_perc 100
# Load results
blast_res <- read_delim("data/sequences/trnl/edna/cleaned/trnl_blast_res.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# Filter top hips
top_hits <- blast_res %>%
  # Get all hits tied for the maximum score
  group_by(qseqid) %>%
  # And keep only thos with 100% match
  filter(pident == 100) %>%
  filter(bitscore == max(bitscore)) %>% 
  # Get taxonomy from match id
  separate(sseqid, c("clade", "family", "genus", "species"), sep = ";") 
print(top_hits, n = 100)

# Determine least common ancestor
lca_taxonomy <- top_hits %>%
  group_by(qseqid) %>%
  summarize(
    # Count how many unique taxa exist among the tied top hits
    n_species = n_distinct(species),
    n_genus   = n_distinct(genus),
    
    # Store the highest percent identity for reporting
    max_ident    = max(pident),
    
    # DYNAMIC LCA LOGIC:
    # Level 1: Species Assignment (If there's exactly 1 unique species tied at the top)
    species = if_else(n_species == 1, dplyr::first(species), NA_character_),
    
    # Level 2: Genus Fallback (If multiple species tied, but they all belong to 1 genus)
    genus   = if_else(n_genus == 1, dplyr::first(genus), NA_character_),
    
    # Keep family
    family  = dplyr::first(family)
  ) %>%
  relocate(qseqid, family, genus, species) 
tibble::view(lca_taxonomy)

# Subset ASV table and assign taxonomy 
trnl_taxa <- trnl_seqtab_nochim[, lca_taxonomy$qseqid] %>%
    # Reformat into tibble
    t() %>%
        data.frame() %>%
        mutate(., sq = row.names(.)) %>%
        mutate(asv_id = lca_taxonomy$qseqid) %>%
        remove_rownames() %>%
        relocate(sq, asv_id) %>%
        pivot_longer(
            .,
            cols = 3:ncol(.),
            names_to = "sample_id",
            values_to = "count"
        ) %>%
    # Format sample ids
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
        ) %>% 
        # Merge with taxonomy
        full_join(lca_taxonomy, by = c("asv_id" = "qseqid")) %>%
        # Remove non-zero counts
        filter(count > 0)
tibble::view(trnl_taxa)

# Denoising step 2 !!!
# Remove count of ASVs found in blanks
trnl_taxa %>% filter(is.na(plot_id), grepl("C", sample_id)) %>%
print(n = 100)
# # A tibble: 2 × 12
#   sq            asv_id sample_id count plot_id subsample_id family genus species
#   <chr>          <dbl> <chr>     <int> <chr>   <chr>        <chr>  <chr> <chr>  
# 1 ATCCTGTTTTAT…      1 CD230723      5 NA      NA           Rosac… Dryas integr…
# 2 ATCCTATTTTTC…      2 CJ230729      2 NA      NA           Salic… Salix NA  
# Two blanks had contamination (likely tagjumping?)
# --> Deduct these values from the counts
trnl_taxa <- trnl_taxa %>%
    mutate(count = case_when(asv_id == 1 ~ count - 5,
                             asv_id == 2 ~ count - 2,
                             TRUE ~ count))

# Denosing step 3 !!
# Remove any sample ASV-combo with a count of less or equalt to 5
filter(trnl_taxa, count <= 5)
# Remove two sample asv combinations
trnl_taxa <- filter(trnl_taxa, count > 5)

# Final clean up !!! - Keep only assigned samples
filter(trnl_taxa, is.na(plot_id)) %>% print(n = 30)
# This includes pippeting error 1-3 which were accidentially processed
trnl_taxa <- filter(trnl_taxa, !is.na(plot_id))
# Relocate colums
trnl_taxa <- select(trnl_taxa,
    sample_id, plot_id, subsample_id, family, genus, species, count, asv_id, sq)
# Save
write_csv(trnl_taxa, "data/sequences/trnl/edna/cleaned/trnl_taxa.csv")


## EOF File

### OLD CODE:

## Post processing to generate presence absence data frames ----

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
    data.frame(trnl_taxa) %>%
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
    "data/sequences/trnl/edna/cleaned/trnl_taxa_in_sample.csv"
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
# Botanical obs are slightly less normally distributed.

# Check counts
tax_counts %>% pull(count) %>% hist()

# How many counts are below 16?
tax_counts %>% filter(count < 16) %>% nrow()
# This is fewer than expected -> it seems like the dada2 pipline filter worked
# rather well. Check later whether count is associated with presence.

# Overall looking good.

# End of file
