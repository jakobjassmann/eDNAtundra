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
) %>% select(-full_name)

# Join datastes
assmann_et_al <- full_join(species_list, cb_trnl) %>%
    mutate(source = "Assmann et al")
view(assmann_et_al)
# Loads fasta with sequences from Willerslev et al. 2014
willerslev <- readDNAStringSet("data/external/reference_libraries/Willerslev_et_al_2014/references.arctic-gh.fasta")

# Check species with missing sequences
filter(assmann_et_al, is.na(sequence))
#  1 Arctagrostis latifolia   Arctagro… Arctagros… NA                    NA Assma…
#  2 Astragalus   alpinus     Astragal… Astragalu… NA                    NA Assma…
#  3 Bistorta     vivipara    Bistorta… Bistorta … NA                    NA Assma…
#  4 Carex        atrofusca   Carex at… Carex atr… NA                    NA Assma…
#  5 Carex        rupestris   Carex ru… Carex rup… NA                    NA Assma…
#  6 Cassiope     tetragona   Cassiope… Cassiope … NA                    NA Assma…
#  7 Equisetum    arvense     Equisetu… Equisetum… NA                    NA Assma…
#  8 Equisetum    variegatum  Equisetu… Equisetum… NA                    NA Assma…
#  9 Oxytropis    arctobia    Oxytropi… Oxytropis… NA                    NA Assma…
# 10 Pedicularis  langsdorff… Pedicula… Pedicular… NA                    NA Assma…
# 11 Salix        arctica     Salix ar… Salix arc… NA                    NA Assma…
# 12 Salix        reticulata  Salix re… Salix ret… NA                    NA Assma…
# 13 Salix        richardson… Salix ri… Salix ric… NA                    NA Assma…
# 14 Saxifraga    cespitosa   Saxifrag… Saxifraga… NA                    NA Assma…

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

# Bistorta vivipara
sum(grepl(".*Bistorta vivipara.*", names(willerslev)))
assmann_et_al$sequence[grepl(".*Bistorta vivipara.*", assmann_et_al$short_name)] <- 
    willerslev[grepl(".*Bistorta vivipara.*", names(willerslev))]
assmann_et_al$source[grepl(".*Bistorta vivipara.*", assmann_et_al$short_name)] <- 
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
       -perc_identity 95 \
       -word_size 7 \
       -dust no \
       -evalue 1000

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
    species = if_else(n_species == 1 & max(pident) == 100, dplyr::first(species), NA_character_),
    
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
print(confused_taxa, n = Inf)
# # A tibble: 48 × 4
#    genus_id   species_id     genus      species       
#    <chr>      <chr>          <chr>      <chr>         
#  1 Carex      aquatilis      Carex      bigelowii     
#  2 Carex      aquatilis      Carex      aquatilis     
#  3 Carex      bigelowii      Carex      bigelowii     
#  4 Carex      bigelowii      Carex      aquatilis     
#  5 Carex      borealipolaris Carex      myosuroides   
#  6 Carex      borealipolaris Carex      borealipolaris
#  7 Carex      capillaris     Carex      scirpoidea    
#  8 Carex      capillaris     Carex      membranacea   
#  9 Carex      capillaris     Carex      fuliginosa    
# 10 Carex      capillaris     Carex      capillaris    
# 11 Carex      fuliginosa     Carex      scirpoidea    
# 12 Carex      fuliginosa     Carex      membranacea   
# 13 Carex      fuliginosa     Carex      fuliginosa    
# 14 Carex      fuliginosa     Carex      capillaris    
# 15 Carex      membranacea    Carex      scirpoidea    
# 16 Carex      membranacea    Carex      membranacea   
# 17 Carex      membranacea    Carex      fuliginosa    
# 18 Carex      membranacea    Carex      capillaris    
# 19 Carex      myosuroides    Carex      myosuroides   
# 20 Carex      myosuroides    Carex      borealipolaris
# 21 Carex      scirpoidea     Carex      scirpoidea    
# 22 Carex      scirpoidea     Carex      membranacea   
# 23 Carex      scirpoidea     Carex      fuliginosa    
# 24 Carex      scirpoidea     Carex      capillaris    
# 25 Eriophorum angustifolium  Eriophorum triste        
# 26 Eriophorum angustifolium  Eriophorum angustifolium 
# 27 Eriophorum triste         Eriophorum triste        
# 28 Eriophorum triste         Eriophorum angustifolium 
# 29 Oxytropis  arctobia       Oxytropis  maydelliana   
# 30 Oxytropis  arctobia       Oxytropis  arctobia      
# 31 Oxytropis  maydelliana    Oxytropis  maydelliana   
# 32 Oxytropis  maydelliana    Oxytropis  arctobia      
# 33 Salix      arctica        Salix      richardsonii  
# 34 Salix      arctica        Salix      reticulata    
# 35 Salix      arctica        Salix      polaris       
# 36 Salix      arctica        Salix      arctica       
# 37 Salix      polaris        Salix      richardsonii  
# 38 Salix      polaris        Salix      reticulata    
# 39 Salix      polaris        Salix      polaris       
# 40 Salix      polaris        Salix      arctica       
# 41 Salix      reticulata     Salix      richardsonii  
# 42 Salix      reticulata     Salix      reticulata    
# 43 Salix      reticulata     Salix      polaris       
# 44 Salix      reticulata     Salix      arctica       
# 45 Salix      richardsonii   Salix      richardsonii  
# 46 Salix      richardsonii   Salix      reticulata    
# 47 Salix      richardsonii   Salix      polaris       
# 48 Salix      richardsonii   Salix      arctica   
print(confused_taxa, n = Inf)

# Print only those taxa inluced in the confusion matrix
distinct(confused_taxa, genus_id, species_id)

# A tibble: 16 × 2
#    genus_id   species_id    
#    <chr>      <chr>         
#  1 Carex      aquatilis     
#  2 Carex      bigelowii     
#  3 Carex      borealipolaris
#  4 Carex      capillaris    
#  5 Carex      fuliginosa    
#  6 Carex      membranacea   
#  7 Carex      myosuroides   
#  8 Carex      scirpoidea    
#  9 Eriophorum angustifolium 
# 10 Eriophorum triste        
# 11 Oxytropis  arctobia      
# 12 Oxytropis  maydelliana   
# 13 Salix      arctica       
# 14 Salix      polaris       
# 15 Salix      reticulata    
# 16 Salix      richardsonii 

# These collapse into the follwing groups
# 1. Carex sp. (bigelowii / aquatilis)
# 2. Carex sp. (scirpoidea / membranacea / myosuroides / scirpoidea)
# 3. Eriophorum sp.
# 4. Oxitropis sp. (arctobia / meydelliana)
# 5. Salix sp. 

