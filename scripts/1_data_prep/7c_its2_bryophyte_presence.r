# Bryophytes: Retrieve presence absence for each plot based on its2
# Jakob Assmann jakob.assmann@uzh.ch March 2025 updated July 2026

# Note: requirement is that script 6b_its2_presence_absence.r has been run

# Dependencies
library(tidyverse)
library(dada2)
library(Biostrings)
library(ShortRead)

# Blast the its2 asv against the BITSR database (run on linux command line)
blastn -query data/sequences/its2/edna/cleaned/its2_asvs.fasta \
       -db data/external/reference_libraries/DS-BITSR/DS-BITSR_its2_db \
       -out data/sequences/its2/edna/cleaned/its2_bryo_blast_res.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95

# Read test result
blast_test <- read_delim("data/sequences/its2/edna/cleaned/its2_bryo_blast_res.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast_test %>% tibble::view()

# Quality filter test result
top_hits <- blast_test %>%
   # For 140 BP match minimum length for a plant
  filter(length >= 140) %>%
  # Get all hits tied for the maximum score
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>% 
  # Get taxonomy from match id
  mutate(sseqid = gsub("ID_", "", sseqid)) %>%
  separate(sseqid, c("genus", "species"), sep = "_") 
print(top_hits, n = 100)

# Assigns taxa based on a dynamic least common ancestor logic
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
    species = if_else(n_species == 1 & max(pident) >= 98, dplyr::first(species), NA_character_),
    
    # Level 2: Genus Fallback (If multiple species tied, but they all belong to 1 genus)
    genus   = if_else(n_genus == 1, dplyr::first(genus), NA_character_),
  ) %>%
  relocate(qseqid, genus, species) 
tibble::view(lca_taxonomy)
hist(lca_taxonomy$max_ident)

# Subset ASV table and assign taxonomy 
load("data/sequences/its2/edna/cleaned/its2_seqtab_nochim.rda")
its2_taxa <- its2_seqtab_nochim[, lca_taxonomy$qseqid] %>%
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
tibble::view(its2_taxa)

# Denoising step 2 !!!
# Remove count of ASVs found in blanks
its2_taxa %>% filter(is.na(plot_id), grepl("C", sample_id)) %>%
print(n = 100)
# A tibble: 0 × 11
# ℹ 11 variables: sq <chr>, asv_id <dbl>, sample_id <chr>, count <int>,
#   plot_id <chr>, subsample_id <chr>, genus <chr>, species <chr>,
#   n_species <int>, n_genus <int>, max_ident <dbl>
# None this is great!

# Denosing step 3 !!
# Remove any sample ASV-combo with a count of less or equalt to 5
filter(its2_taxa, count <= 5) %>% print(n = Inf)
# Remove 25 sample asv combinations
its2_taxa <- filter(its2_taxa, count > 5)

# Final clean up !!! - Keep only assigned samples
filter(its2_taxa, is.na(plot_id)) %>% print(n = 30)
# This includes pippeting error 2-3 which were accidentially processed
its2_taxa <- filter(its2_taxa, !is.na(plot_id))
# Relocate colums
its2_taxa <- select(its2_taxa,
    sample_id, plot_id, subsample_id, genus, species, count, asv_id, sq)
# Save
write_csv(its2_taxa, "data/sequences/its2/edna/cleaned/its2_bryo_taxa.csv")

# Quick check of plot level distribution
its2_taxa %>% distinct(genus, species, plot_id) %>%
group_by(genus, species) %>% summarise(n = n()) %>% tibble:: view()

# None.

### End of file
################################################################################
