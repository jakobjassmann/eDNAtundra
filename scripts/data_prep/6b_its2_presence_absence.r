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

# Get file names for all its2 reads
its2_fwd <- list.files("data/sequences/its2/edna/trimmed",
                       pattern = ".*_R1.*fastq",
                       full.names = T
)
its2_rev <- list.files("data/sequences/its2/edna/trimmed",
                       pattern = ".*_R2.*fastq",
                       full.names = T
)

# Set out path
out_path <- "data/sequences/its2/edna"

# Get sample ids
sample_ids <- gsub("([A-Za-z0-9]*-?[123P]?)-.*", "\\1", basename(its2_fwd))

## Visual QC of a couple of samples
plotQualityProfile(its2_fwd[17:18])
plotQualityProfile(its2_rev[17:18])

## Filtering ---

# Set out file paths
its2_fwd_filtered <- file.path(out_path, "filtered", basename(its2_fwd))
its2_rev_filtered <- file.path(out_path, "filtered", basename(its2_rev))

# Filter and trim based on quality
out <- filterAndTrim(fwd = its2_fwd, filt = its2_fwd_filtered,
                     rev = its2_rev, filt.rev = its2_rev_filtered, 
                     maxN = 0, # Discard sequences with 1 N or more (missing read) 
                     maxEE = c(2, 2), # After truncation, reads with higher than 2 "expected errors" will be discarded
                     truncQ = 2, # Truncate  at the first instance of a quality score less than or equal to 2
                     minLen = 10,
                     rm.phix = TRUE, # RUE Discard reads that match against the phiX genome (ILLUMINA control)
                     compress = TRUE,
                     multithread = TRUE)
# Check results
out 

## Identify ASVs ----

# Learn errors
err_fwd <- learnErrors(its2_fwd_filtered, multithread = TRUE)
err_rev <- learnErrors(its2_rev_filtered, multithread = TRUE)

# Visulaise estimated error rates
plotErrors(err_fwd, nominalQ=TRUE)
plotErrors(err_rev, nominalQ=TRUE)

# There are some strange dips in the error rate functions at high quality scores. 
# Not sure what to do with it for now, proceeding and then reading up on it later.

# Dereplicate
exists <- file.exists(its2_fwd_filtered)
its2_few_derep <- derepFastq(its2_fwd_filtered[exists], verbose=TRUE)
its2_rev_derep <- derepFastq(its2_rev_filtered[exists], verbose=TRUE)

# Add sample names
names(its2_few_derep) <- sample_ids[exists]
names(its2_rev_derep) <- sample_ids[exists]

# Run the data core algorithm to infer ASV
its2_fwd_dada <- dada(its2_few_derep, err = err_fwd, multithread = TRUE)
its2_rev_dada <- dada(its2_rev_derep, err = err_rev, multithread = TRUE)

# Merge paired-end reads
its2_merged <- mergePairs(its2_fwd_dada, its2_few_derep,
                      its2_rev_dada, its2_rev_derep, verbose=TRUE)

# Convert to sequence table of ASVs
its2_seqtab <- makeSequenceTable(its2_merged)
dim(its2_seqtab) # n samples n sequences
table(nchar(getSequences(its2_seqtab))) # distribution of sequence lengths

# Reove chimeras
its2_seqtab_nochim <- removeBimeraDenovo(its2_seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(its2_seqtab_nochim) # n samples n non-chimeric sequences
sum(its2_seqtab_nochim)/sum(its2_seqtab) # prop. non-chimeric sequences
table(nchar(getSequences(its2_seqtab_nochim))) # distribution of non-chimeric sequence lengths

# Track sequences through pipeline (from tutorial)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(its2_fwd_dada, getN),
               sapply(its2_rev_dada, getN),
               sapply(its2_merged, getN),
               rowSums(its2_seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_ids[exists]
head(track)
View(track)

##  Assign species ----

its2_taxa <- assignSpecies(its2_seqtab_nochim,
                      "data/sequences/reference/its2_reference_Assmann_et_al_2025.fasta",
                      verbose = T)
its2_taxa %>% data.frame() %>% remove_rownames() %>% filter(!is.na(Genus)) %>% View()

# Save objects
dir.create("data/sequences/its2/edna/cleaned")
write(its2_seqtab_nochim, "data/sequences/its2/edna/cleaned/its2_seqtab_nochim.rda")
write(its2_taxa, "data/sequences/its2/edna/cleaned/its2_taxa.rda")
write_csv(track %>% data.frame(), "data/sequences/its2/edna/cleaned/track.csv")

# End of file
    

