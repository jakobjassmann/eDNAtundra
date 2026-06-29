# Build trl reference libraries combinding the Cambridge Bay samples with Willerslev et al.
# Jakob J. Assmann 2026-06-29

# Dependencies
library(tidyverse)
library(readxl)
library(Biostrings)

# Load species list from vascular plant observatins
species_list <- read_csv("data/botanical_obs/species_list_cleaned.csv")

# Load reference our own trnl sequences
cb_trnl <- read_csv(
    "data/sequences/trnl/single_species/trnl_single_species_cleaned.csv"
)

# Join datastes
assmann_et_al <- full_join(species_list, cb_trnl) %>%
    mutate(source = "Assmann et al")

# Loads fasta with sequences from Willerslev et al. 2014
willerslev <- readDNAStringSet("data/external/reference_libraries/Willerslev_et_al_2014/references.arctic-gh.fasta")

# Check species with missing sequences
filter(assmann_et_al, is.na(sequence))
#    <chr>        <chr>       <chr>     <chr>      <chr>              <dbl> <chr> 
#  1 Arctagrostis latifolia   Arctagro… Arctagros… NA                    NA Assma…
#  2 Astragalus   alpinus     Astragal… Astragalu… NA                    NA Assma…
#  3 Carex        atrofusca   Carex at… Carex atr… NA                    NA Assma…
#  4 Carex        rupestris   Carex ru… Carex rup… NA                    NA Assma…
#  5 Cassiope     tetragona   Cassiope… Cassiope … NA                    NA Assma…
#  6 Equisetum    arvense     Equisetu… Equisetum… NA                    NA Assma…
#  7 Equisetum    variegatum  Equisetu… Equisetum… NA                    NA Assma…
#  8 Oxytropis    arctobia    Oxytropi… Oxytropis… NA                    NA Assma…
#  9 Pedicularis  langsdorff… Pedicula… Pedicular… NA                    NA Assma…
# 10 Salix        arctica     Salix ar… Salix arc… NA                    NA Assma…
# 11 Salix        reticulata  Salix re… Salix ret… NA                    NA Assma…
# 12 Salix        richardson… Salix ri… Salix ric… NA                    NA Assma…
# 13 Saxifraga    cespitosa   Saxifrag… Saxifraga… NA                    NA Assma…

# Fill missing values

# Arctagrostis latifolia
sum(grepl(".*Arctagrostis latifolia.*", names(willerslev)))
# It is included in "Pooideae_GROUP_1"
sum(grepl(".*Pooideae_GROUP_1.*", names(willerslev))) 
assmann_et_al$sequence[grepl(".*Arctagrostis latifolia.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Pooideae_GROUP_1.*", names(willerslev))]
assmann_et_al$source[grepl(".*Arctagrostis latifolia.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Astralagus alpinus
sum(grepl(".*Astragalus alpinus.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Astragalus alpinus.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Astragalus alpinus.*", names(willerslev))]
assmann_et_al$source[grepl(".*Astragalus alpinus.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Carex atrofusca
sum(grepl(".*Carex atrofusca.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Carex atrofusca.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Carex atrofusca.*", names(willerslev))]
assmann_et_al$source[grepl(".*Carex atrofusca.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Carex rupestris
sum(grepl(".*Carex rupestris.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Carex rupestris.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Carex rupestris.*", names(willerslev))]
assmann_et_al$source[grepl(".*Carex rupestris.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Cassiope tetragona 
sum(grepl(".*Cassiope tetragona.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Cassiope tetragona.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Cassiope tetragona.*", names(willerslev))]
assmann_et_al$source[grepl(".*Cassiope tetragona.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Equisetum arvense
sum(grepl(".*Equisetum arvense.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Equisetum arvense.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Equisetum arvense.*", names(willerslev))]
assmann_et_al$source[grepl(".*Equisetum arvense.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Equisetum variegatum
sum(grepl(".*Equisetum variegatum.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Equisetum variegatum.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Equisetum variegatum.*", names(willerslev))]
assmann_et_al$source[grepl(".*Equisetum variegatum.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Oxytropis arctobia
sum(grepl(".*Oxytropis arctobia.*", names(willerslev)))
# Included as Oxytropis arctica (Oxytropis_GROUP_1) in Willerslev et al. 2014 
sum(grepl(".*Oxytropis_GROUP_1.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Oxytropis arctobia.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Oxytropis_GROUP_1.*", names(willerslev))]
assmann_et_al$source[grepl(".*Oxytropis arctobia.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Pedicularis langsdorffii
sum(grepl(".*Pedicularis langsdorffii.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Pedicularis langsdorffii.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Pedicularis langsdorffii.*", names(willerslev))]
assmann_et_al$source[grepl(".*Pedicularis langsdorffii.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Salix arctica
sum(grepl(".*Salix arctica.*", names(willerslev)))
# Included as Salix arctica (Saliceae_GROUP_1) in Willerslev et al. 2014
sum(grepl(".*Saliceae_GROUP_1.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Salix arctica.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Saliceae_GROUP_1.*", names(willerslev))]
assmann_et_al$source[grepl(".*Salix arctica.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Salix reticulata
sum(grepl(".*Salix reticulata.*", names(willerslev)))
# Included as Salix reticulata. (Saliceae_GROUP_1) in Willerslev et al. 2014
sum(grepl(".*Saliceae_GROUP_1.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Salix reticulata.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Saliceae_GROUP_1.*", names(willerslev))]
assmann_et_al$source[grepl(".*Salix reticulata.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Salix richardsonii
sum(grepl(".*Salix richardsonii.*", names(willerslev)))
# Included as Salix richardsonii. (Saliceae_GROUP_1) in Willerslev et al. 2014
sum(grepl(".*Saliceae_GROUP_1.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Salix richardsonii.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Saliceae_GROUP_1.*", names(willerslev))]
assmann_et_al$source[grepl(".*Salix richardsonii.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"
   
# Saxifraga cespitosa 
sum(grepl(".*Saxifraga cespitosa.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Saxifraga cespitosa.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Saxifraga cespitosa.*", names(willerslev))]
assmann_et_al$source[grepl(".*Saxifraga cespitosa.*", assmann_et_al$short_name)] <- 
    "Willerslev et al. 2014"

# Check species with missing sequences
filter(assmann_et_al, is.na(sequence))

### Merge taxnonym with database

# Load file from Saarela et al. 2020 to get clade and family
saarela <- read_xlsx(
    "data/external/species_list/phytokeys-141-001-s003.xlsx"
)
saarela <- saarela %>%
    select(2:13) %>%
    set_names(., dplyr::slice(., 1) %>% as.character()) %>%
    dplyr::slice(c(-1, -2)) %>%
    filter(., !is.na(`Cambridge Bay`) | !is.na(`Ovayok Territorial Park`)) %>%
    select(1:3) %>%
    mutate(
        genus = sapply(str_extract_all(Taxon, "\\S+"), `[`, 1),
        species = sapply(str_extract_all(Taxon, "\\S+"), `[`, 2),
        full_name = Taxon
    ) %>%
    select(-Taxon)

# Join clade and family to database
assmann_et_al <- assmann_et_al %>%
    left_join(saarela) %>%
    select(
        clade = Clade,
        family = Family,
        genus,
        species,
        short_name,
        full_name,
        source,
        sequence_length,
        sequence
    )

# Check for nas
filter(assmann_et_al, is.na(clade))

# Fill manually
saarela %>% filter(genus == "Carex") %>% distinct(Clade, Family)
saarela %>% filter(genus == "Tofieldia") %>% distinct(Clade, Family)

assmann_et_al$clade[which(is.na(assmann_et_al$clade))] <- c(
    saarela %>%
            filter(genus == "Carex") %>%
            distinct(Clade, Family) %>%
            pull(Clade),
    saarela %>%
        filter(genus == "Tofieldia") %>%
        distinct(Clade, Family) %>%
        pull(Clade)
)
assmann_et_al$family[which(is.na(assmann_et_al$family))] <- c(
    saarela %>%
            filter(genus == "Carex") %>%
            distinct(Clade, Family) %>%
            pull(Family),
    saarela %>%
        filter(genus == "Tofieldia") %>%
        distinct(Clade, Family) %>%
        pull(Family)
)

# Check for nas
filter(assmann_et_al, is.na(clade))

### Write out ast csv, fasta and create blast database

# Write out
write_excel_csv(
    assmann_et_al,
    "data/sequences/trnl/reference/Assmann_et_al_trnl_ref.csv"
)

# Create DNAStringSet
sequences <- DNAStringSet(assmann_et_al$sequence)
names(sequences) <- paste(
    assmann_et_al$clade,
    assmann_et_al$family,
    assmann_et_al$genus,
    assmann_et_al$species,
    sep = ";"
)

# Write FASTA
writeXStringSet(
    sequences,
    filepath = "data/sequences/trnl/reference/Assmann_et_al_trnl_ref.fasta"
)

# Create BLAST database in terminal (external command)
# makeblastdb -in data/sequences/trnl/reference/Assmann_et_al_trnl_ref.fasta -dbtype nucl -out data/sequences/trnl/reference/Assmann_et_al_trnl_ref_db

#### Test reference database against itself

# Generate an indexed fasta for testing
names(sequences) <- paste0(1:length(names(sequences)), ";", names(sequences))
writeXStringSet(
    sequences,
    filepath = "data/sequences/trnl/reference/Assmann_et_al_trnl_test.fasta"
)
# Test blast the database against itself
blastn -query data/sequences/trnl/reference/Assmann_et_al_trnl_test.fasta \
       -db data/sequences/trnl/reference/Assmann_et_al_trnl_ref_db \
       -out data/sequences/trnl/reference/test_blast.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95

# Read test result
blast_test <- read_delim("data/sequences/trnl/reference/test_blast.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast_test %>% separate(qseqid, c("id", "clade_id", "family_id", "genus_id", "species_id"), sep = ";") %>% view()
blast_test %>% separate(qseqid, c("id", "clade_id", "family_id", "genus_id", "species_id"), sep = ";") %>% distinct(id) %>% nrow()

# Quality filter test result
top_hits <- blast_test %>%
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
  relocate(qseqid, family, genus, species) %>%
  separate(qseqid, c("id", "clade_id", "family_id", "genus_id", "species_id"), sep = ";") %>%
  arrange(id) %>%
  mutate(family_match = family_id == family,
        genus_match = genus_id == genus,
        species_match = species_id == species)
view(lca_taxonomy)

# Get confusion matrix
confused_taxa <- top_hits %>%
  separate(qseqid, c("id", "clade_id", "family_id", "genus_id", "species_id"), sep = ";") %>%
 filter(id %in% (filter(lca_taxonomy, is.na(species_match)) %>% pull(id))) %>%
 distinct(genus_id, species_id, genus, species)
 # 6 confusion pairs
# # A tibble: 12 × 4
#    genus_id   species_id    genus      species      
#    <chr>      <chr>         <chr>      <chr>        
#  1 Carex      aquatilis     Carex      bigelowii    
#  2 Carex      aquatilis     Carex      aquatilis    
#  3 Carex      bigelowii     Carex      bigelowii    
#  4 Carex      bigelowii     Carex      aquatilis    
#  5 Eriophorum angustifolium Eriophorum triste       
#  6 Eriophorum angustifolium Eriophorum angustifolium
#  7 Eriophorum triste        Eriophorum triste       
#  8 Eriophorum triste        Eriophorum angustifolium
#  9 Salix      reticulata    Salix      richardsonii 
# 10 Salix      reticulata    Salix      reticulata   
# 11 Salix      richardsonii  Salix      richardsonii 
# 12 Salix      richardsonii  Salix      reticulata  

# Print only those taxa inluced in the confusion matrix
distinct(confused_taxa, genus_id, species_id)

# 6 species with confusion
# # A tibble: 6 × 2
#   genus_id   species_id   
#   <chr>      <chr>        
# 1 Carex      aquatilis    
# 2 Carex      bigelowii    
# 3 Eriophorum angustifolium
# 4 Eriophorum triste       
# 5 Salix      reticulata   
# 6 Salix      richardsonii 