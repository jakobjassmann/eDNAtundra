# Quick script to check the completness of the DS-BITSR for our moss species list
# 4 August 2026 Jakob J. Assmann jakob.assmann@gu.se

# Dependencies
library(tidyverse)
library(Biostrings)

# Load bryophyte species list and filter unique species
bryo_obs <- read_csv(
    "data/botanical_obs/moss_data_Erin_cleaned_final_Deb_updated_Aug_2026.csv"
) %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(
        site_id = Site,
        plot_id,
        genus = Genus,
        species = Species,
        short_name = Full_Name,
        synonym
    ) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) | !grepl("cf", species)) %>%
    filter(!grepl("sp\\.", species)) %>%
    distinct()

# Save distinct species and their synonyms for later
bryo_obs_distinct <- distinct(bryo_obs, short_name, synonym) %>%
    arrange(short_name)

# Use synynoms for DS-BITSR compariosn
bryo_obs <- bryo_obs %>%
    mutate(short_name = synonym) %>%
    distinct(short_name) %>%
    arrange(short_name)

# Load species names in DS-BITSR
ds_bitsr <- readDNAStringSet(
    "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta"
) %>%
    names() %>%
    gsub("ID_", "", .) %>%
    gsub("_", " ", .) %>%
    unique() %>%
    sort()

# Remove any sub-species level IDs for matching
ds_bitsr[grepl("subsp\\.", ds_bitsr)] <- gsub(
    "(.) subsp\\..*",
    "\\1",
    ds_bitsr[grepl("subsp\\.", ds_bitsr)]
)
ds_bitsr[grepl("var\\.", ds_bitsr)] <- gsub(
    "(.) subsp\\..*",
    "\\1",
    ds_bitsr[grepl("var\\.", ds_bitsr)]
)

# Check which species are present in reference database
ds_bitsr_complete_check <- tibble(
    bryo_obs,
    ds_bitsr = bryo_obs$short_name %in% ds_bitsr
)

# Quick check of direct matches
sum(ds_bitsr_complete_check$ds_bitsr)
# [1] 53

# No match
sum(!ds_bitsr_complete_check$ds_bitsr)
# 14

# Save DS-BITSR taxonomy
ds_bitrs_spp <- ds_bitsr_complete_check

# revert to names used in our taxonomy
ds_bitsr_complete_check <- ds_bitsr_complete_check %>%
    dplyr::rename(synonym = short_name) %>%
    full_join(bryo_obs_distinct) %>%
    mutate(synonym = if_else(synonym == short_name, "", synonym)) %>%
    arrange(short_name)

# Export CSV of matches
write_excel_csv(
    ds_bitsr_complete_check %>%
        select(short_name, synonym, ds_bitsr) %>%
        setNames(c("Name", "Synonym", "In_DS-BITSR")),
    "data/external/reference_libraries/DS-BITSR/completness_check.csv"
)

#### Test reference database against itself (for relevant species)

# Subset species present in DS-BITSR and prepare for ds_bitsr matching
ds_bitrs_spp <- filter(ds_bitrs_spp, ds_bitsr) %>% pull(short_name)
ds_bitrs_spp <- gsub(" ", "_", ds_bitrs_spp)
ds_bitrs_spp <- paste0("ID_", ds_bitrs_spp)

# Load DS-BITR and subset FASTA
ds_bitsr <- readDNAStringSet(
    "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2.fasta"
)

# Quick sanity check that all are present
ds_bitrs_spp %in% names(ds_bitsr)

# Subset DS-BITSR and eexport
ds_bitsr_sub <- ds_bitsr[names(ds_bitsr) %in% ds_bitrs_spp]

# Generate an indexed fasta for testing
names(ds_bitsr_sub) <- gsub("ID_", "", names(ds_bitsr_sub))
names(ds_bitsr_sub) <- gsub("_", ";", names(ds_bitsr_sub))
names(ds_bitsr_sub) <- paste0(1:length(names(ds_bitsr_sub)), ";", names(ds_bitsr_sub))
writeXStringSet(
    ds_bitsr_sub,
    filepath = "data/external/reference_libraries/DS-BITSR/DS-BITSR_its2_test.fasta"
)
# Test blast the database against itself
blastn -query data/external/reference_libraries/DS-BITSR/DS-BITSR_its2_test.fasta \
       -db data/external/reference_libraries/DS-BITSR/DS-BITSR_its2_db \
       -out data/external/reference_libraries/DS-BITSR/test_blast.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95

# Read test result
blast_test <- read_delim("data/external/reference_libraries/DS-BITSR/test_blast.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast_test %>% separate(qseqid, c("id", "genus_id", "species_id"), sep = ";") %>% view()

# Qulaity filter test result
top_hits <- blast_test %>%
   # For 180 BP match minimum length for plant
  filter(length >= 180) %>%
  # Get all hits tied for the maximum score
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>% 
  # Get taxonomy from match id
  separate(sseqid, c("filler", "genus", "species"), sep = "_") 
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
  relocate(qseqid, genus, species) %>%
  separate(qseqid, c("id", "genus_id", "species_id"), sep = ";") %>%
  arrange(id) %>%
  mutate(
        genus_match = genus_id == genus,
        species_match = species_id == species)
view(lca_taxonomy)

# Get confusion matrix
confused_taxa <- top_hits %>%
  separate(qseqid, c("id", "genus_id", "species_id"), sep = ";") %>%
 filter(id %in% (filter(lca_taxonomy, is.na(species_match)) %>% pull(id))) %>%
 distinct(genus_id, species_id, genus, species)
confused_taxa %>% print(n = 50)
 # A tibble: 37 × 4
#    genus_id          species_id       genus             species         
#    <chr>             <chr>            <chr>             <chr>           
#  1 Hygroamblystegium varium           Hygroamblystegium varium          
#  2 Hygroamblystegium varium           Hygroamblystegium noterophilum    
#  3 Brachythecium     acutum           Brachythecium     albicans        
#  4 Brachythecium     acutum           Brachythecium     acutum          
#  5 Lewinskya         speciosa         Lewinskya         speciosa        
#  6 Lewinskya         speciosa         Lewinskya         affinis         
#  7 Lewinskya         speciosa         Lewinskya         pylaisii        
#  8 Myurella          tenerrima        Myurella          tenerrima       
#  9 Myurella          tenerrima        Myurella          julacea         
# 10 Campylium         stellatum        Campylium         bambergeri      
# 11 Campylium         stellatum        Campylium         stellatum       
# 12 Campylium         stellatum        Campylophyllum    hispidulum      
# 13 Drepanocladus     brevifolius      Drepanocladus     aduncus         
# 14 Drepanocladus     brevifolius      Drepanocladus     polygamus       
# 15 Drepanocladus     brevifolius      Drepanocladus     brevifolius     
# 16 Drepanocladus     brevifolius      Drepanocladus     longifolius     
# 17 Campylium         stellatum        Pseudocampylium   radicale        
# 18 Brachythecium     acutum           Brachythecium     acuminatum      
# 19 Brachythecium     acutum           Brachythecium     campestre       
# 20 Brachythecium     turgidum         Brachythecium     turgidum        
# 21 Brachythecium     turgidum         Brachythecium     asperrimum      
# 22 Brachythecium     coruscum         Brachythecium     turgidum        
# 23 Brachythecium     coruscum         Brachythecium     coruscum        
# 24 Brachythecium     coruscum         Brachythecium     erythrorrhizon  
# 25 Brachythecium     coruscum         Brachythecium     salebrosum      
# 26 Brachythecium     turgidum         Brachythecium     coruscum        
# 27 Brachythecium     turgidum         Brachythecium     erythrorrhizon  
# 28 Brachythecium     turgidum         Brachythecium     salebrosum      
# 29 Campylium         bambergeri       Campylium         bambergeri      
# 30 Campylium         bambergeri       Campylium         stellatum       
# 31 Campylium         bambergeri       Campylophyllum    hispidulum      
# 32 Tetraplodon       mnioides         Tetraplodon       mnioides        
# 33 Tetraplodon       mnioides         Tetraplodon       angustatus      
# 34 Ptychostomum      pseudotriquetrum Ptychostomum      neodamense      
# 35 Ptychostomum      pseudotriquetrum Ptychostomum      pseudotriquetrum
# 36 Plagiomnium       ellipticum       Plagiomnium       ciliare         
# 37 Plagiomnium       ellipticum       Plagiomnium       ellipticum   

# Print only those taxa inluced in the confusion matrix
distinct(confused_taxa, genus_id, species_id) %>%
arrange(genus_id, species_id)
# # A tibble: 12 × 2
#    genus_id          species_id      
#    <chr>             <chr>           
#  1 Brachythecium     acutum          
#  2 Brachythecium     coruscum        
#  3 Brachythecium     turgidum        
#  4 Campylium         bambergeri      
#  5 Campylium         stellatum       
#  6 Drepanocladus     brevifolius     
#  7 Hygroamblystegium varium          
#  8 Lewinskya         speciosa        
#  9 Myurella          tenerrima       
# 10 Plagiomnium       ellipticum      
# 11 Ptychostomum      pseudotriquetrum
# 12 Tetraplodon       mnioides  

# -> in the end only four species are truly confused amongst our subset.
# EOF
