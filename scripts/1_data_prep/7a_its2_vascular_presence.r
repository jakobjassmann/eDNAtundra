# Retrieve presence for each plot based on its2
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(Biostrings)
library(ShortRead)

## Load ASVs
load("data/sequences/its2/edna/cleaned/its2_seqtab_nochim.rda")

##  Assign species ----

# Build data frame of ASVs and filter our short sequences
its2_asvs <- tibble(
    id = 1:length(dimnames(its2_seqtab_nochim)[[2]]),
    sequence = dimnames(its2_seqtab_nochim)[[2]],
    width = nchar(sequence)
) %>%
    ### Denoising Step 1 !!!
    # Set all sequences shorter than 140 BP to NA (min widht in ref db is 310)
    # This removes any remaining fungi, chimeras etc
    mutate(
        sequence = case_when(width < 140 ~ NA, TRUE ~ sequence),
        width = case_when(width < 140 ~ NA, TRUE ~ width)
    )
min(its2_asvs$width, na.rm = T)

# Export to a fasta
sequences <- its2_asvs %>%
    na.omit() %>%
    pull(sequence) %>%
    DNAStringSet()
names(sequences) <- its2_asvs %>%
    na.omit() %>%
    pull(id)
writeXStringSet(
    sequences,
    filepath = "data/sequences/its2/edna/cleaned/its2_asvs.fasta"
)

# Blast the data against the reference database (in shell)
blastn -query data/sequences/its2/edna/cleaned/its2_asvs.fasta \
       -db data/sequences/its2/reference/Assmann_et_al_ITS2_ref_db \
       -out data/sequences/its2/edna/cleaned/its2_blast_res.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95

# Read test result
blast_test <- read_delim("data/sequences/its2/edna/cleaned/its2_blast_res.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast_test %>% tibble::view()

# Quality filter test result
top_hits <- blast_test %>%
   # For 140 BP match minimum length for plant
  filter(length >= 140) %>%
  # Get all hits tied for the maximum score
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>% 
  # Get taxonomy from match id
  separate(sseqid, c("clade", "family", "genus", "species"), sep = ";") 
print(top_hits, n = 100)

# Assigns taxa based on a dynamic least common ancestor logic
lca_taxonomy <- top_hits %>%
  group_by(qseqid) %>%
  summarize(
    # Count how many unique taxa exist among the tied top hits
    n_species = n_distinct(species),
    n_genus   = n_distinct(genus),
    n_family  = n_distinct(family),
    
    # Store the highest percent identity for reporting
    max_ident    = max(pident),
    
    # DYNAMIC LCA LOGIC:
    # Level 1: Species Assignment (If there's exactly 1 unique species tied at the top)
    species = if_else(n_species == 1 & max(pident) >= 98, dplyr::first(species), NA_character_),
    
    # Level 2: Genus Fallback (If multiple species tied, but they all belong to 1 genus)
    genus   = if_else(n_genus == 1, dplyr::first(genus), NA_character_),
    
    # Level 3: Family Fallback (If hits span multiple genera, but they all belong to 1 family)
    family  = if_else(n_family == 1, dplyr::first(family), NA_character_)
  ) %>%
  relocate(qseqid, family, genus, species) 
tibble::view(lca_taxonomy)
hist(lca_taxonomy$max_ident)

# Subset ASV table and assign taxonomy 
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
# # A tibble: 1 × 13
#   sq            asv_id sample_id count plot_id subsample_id family genus species
#   <chr>          <dbl> <chr>     <int> <chr>   <chr>        <chr>  <chr> <chr>  
# 1 TGCAGAATCCCG…      2 CJ230729     73 NA      NA           Rosac… Dryas integr…
# One sample had contamination (likely true conamination or tagjumping?)
# --> Deduct these values from the counts
its2_taxa <- its2_taxa %>%
    mutate(count = case_when(asv_id == 2 ~ count - 73,
                             TRUE ~ count))

# Denosing step 3 !!
# Remove any sample ASV-combo with a count of less or equalt to 5
filter(its2_taxa, count <= 5) %>% print(n = Inf)
# Remove 84 sample asv combinations
its2_taxa <- filter(its2_taxa, count > 5)

# Final clean up !!! - Keep only assigned samples
filter(its2_taxa, is.na(plot_id)) %>% print(n = 30)
# This includes pippeting error 1-3 which were accidentially processed
its2_taxa <- filter(its2_taxa, !is.na(plot_id))
# Relocate colums
its2_taxa <- select(its2_taxa,
    sample_id, plot_id, subsample_id, family, genus, species, count, asv_id, sq)
# Save
write_csv(its2_taxa, "data/sequences/its2/edna/cleaned/its2_taxa.csv")

# Quick check of plot level distribution
its2_taxa %>% distinct(family, genus, species, plot_id) %>%
group_by(family, genus, species) %>% summarise(n = n()) %>% tibble:: view()
# Check for any family onlin IDs
its2_taxa %>% filter(is.na(genus), is.na(species)) %>% print(n= 100)
# None.

### EOF
################################################################################
