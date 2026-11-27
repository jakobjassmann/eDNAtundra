# Retrieve presence for each plot based on trnl
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(Biostrings)
library(ShortRead)

## Load ASVs ----
load("data/sequences/trnl/edna/cleaned/trnl_seqtab_nochim.Rda")

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
################################################################################
