# Build ITS2 reference libraries combinding the Cambridge Bay samples with the published databeses
# Jakob J. Assmann 2026-06-18

#### Preparations

# Dependencies
library(tidyverse)
library(readxl)
library(Biostrings)

# Load species list from vascular plant observations
species_list <- read_csv("data/botanical_obs/species_list_cleaned.csv")

# Load reference our own ITS2 sequences
cb_its2 <- read_csv(
    "data/sequences/its2/single_species/its2_single_species_cleaned.csv"
)

# Join datastes
assmann_et_al <- full_join(species_list, cb_its2) %>%
    mutate(source = "Assmann et al")

# Load DS-ZACKPLAN and prepare for merging
dszackplan <- read_tsv(
    "data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN.tsv"
) %>%
    filter(
        grepl("ITS_S2F", primers_forward) &
            grepl("ITS4", primers_reverse) &
            marker_code == "ITS2"
    ) %>%
    select(short_name = species, processid, inst, sequence = nuc) %>%
    # remove any dashes in sequence and convert to all caps
    mutate(sequence = gsub("-", "", sequence)) %>%
    mutate(sequence = str_to_upper(sequence)) %>%
    # Keep only non-ambigous sequences
    filter(grepl("\\<[ATGC]*\\>", sequence)) %>%
    right_join(species_list) %>%
    arrange(genus, species) %>%
    # Add source and calculate length
    mutate(
        source = paste("DS-ZACKPLAN", processid),
        sequence_length = nchar(sequence)
    )

# Load DS-VASCAN and prepare for merging
dsvascan <- read_tsv(
    "data/external/reference_libraries/DS-VASCAN/DS-VASCAN.tsv"
) %>%
    filter(
        grepl("ITS_S2F", primers_forward) &
            grepl("ITS4", primers_reverse) &
            marker_code == "ITS2"
    ) %>%
    select(short_name = species, processid, inst, sequence = nuc) %>%
    # remove any dashes in sequence and convert to all caps
    mutate(sequence = gsub("-", "", sequence)) %>%
    mutate(sequence = str_to_upper(sequence)) %>%
    # Keep only non-ambigous sequences
    filter(grepl("\\<[ATGC]*\\>", sequence)) %>%
    right_join(species_list) %>%
    arrange(genus, species) %>%
    # Add source and calculate length
    mutate(
        source = paste("DS-VASCAN", processid),
        sequence_length = nchar(sequence)
    )


#### Helper datasets for QC

# Export full DS-VASCAN to help with QC
read_tsv(
    "data/external/reference_libraries/DS-VASCAN/DS-VASCAN.tsv"
) %>%
    filter(
        grepl("ITS_S2F", primers_forward) &
            grepl("ITS4", primers_reverse) &
            marker_code == "ITS2"
    ) %>%
    select(short_name = species, processid, inst, sequence = nuc) %>%
    arrange(short_name) %>%
    write_excel_csv("data/sequences/its2/reference/DS-VASCAN_for_qc.csv")

# Export full DS-ZACKPLAN for helping with qc
read_tsv(
    "data/external/reference_libraries/DS-ZACKPLAN/DS-ZACKPLAN.tsv"
) %>%
    filter(
        grepl("ITS_S2F", primers_forward) &
            grepl("ITS4", primers_reverse) &
            marker_code == "ITS2"
    ) %>%
    select(short_name = species, processid, inst, sequence = nuc) %>%
    arrange(short_name) %>%
    write_excel_csv("data/sequences/its2/reference/DS-ZACKPLAN_for_qc.csv")

# To help with gap filling and qc merge all datasets
all_joined <- assmann_et_al %>%
    full_join(filter(
        dszackplan,
        any(short_name %in% species_list$short_name)
    )) %>%
    full_join(filter(
        dsvascan,
        any(short_name %in% species_list$short_name)
    )) %>%
    arrange(genus, species, source) %>%
    select(
        short_name,
        genus,
        species,
        full_name,
        source,
        processid,
        inst,
        sequence_length,
        sequence
    )

# Export for manual QC
write_excel_csv(
    all_joined,
    "data/sequences/its2/reference/Assmann_et_al_ITS2_ref_for_qc.csv"
)

#### Fill gaps in Assmann et al reference database

# Check species with missing sequences
filter(assmann_et_al, is.na(sequence))
# # A tibble: 12 × 7
#    genus        species     full_name short_name sequence sequence_length source
#    <chr>        <chr>       <chr>     <chr>      <chr>              <dbl> <chr>
#  1 Bistorta     vivipara    Bistorta… Bistorta … NA                    NA Assma…
#  2 Chamaenerion latifolium  Chamaene… Chamaener… NA                    NA Assma…
#  3 Equisetum    arvense     Equisetu… Equisetum… NA                    NA Assma…
#  4 Equisetum    variegatum  Equisetu… Equisetum… NA                    NA Assma…
#  5 Salix        arctica     Salix ar… Salix arc… NA                    NA Assma…
#  6 Salix        polaris     Salix po… Salix pol… NA                    NA Assma…
#  7 Salix        reticulata  Salix re… Salix ret… NA                    NA Assma…
#  8 Salix        richardson… Salix ri… Salix ric… NA                    NA Assma…
#  9 Saxifraga    aizoides    Saxifrag… Saxifraga… NA                    NA Assma…
# 10 Saxifraga    cespitosa   Saxifrag… Saxifraga… NA                    NA Assma…
# 11 Saxifraga    oppositifo… Saxifrag… Saxifraga… NA                    NA Assma…
# 12 Woodsia      glabella    Woodsia … Woodsia g… NA                    NA Assma…

# 1 Bistorta vivipara MKPCH419-09 in DS-VASCAN
# Entry has a photograph that matches B vivipara
# Barcode is complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH419-09")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 2 Chamaenerion latifolium MKPCH428-09 in DS-VASCAN
# Entry has a photograph that matches C. latifolium
# Barcode is complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH428-09")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 3 Equisetum variegatum is neither found in DS-VASCAN or DS-ZACKPLAN for our primer pairs
# (only with ambious sequences in DS-ZACKPLAN)

# 4 Equiesetum arvense is neither found in DS-VASCAN or DS-ZACKPLAN for our primer pairs

# 5 Salix arctica GRAFW1877-13 in DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW1877-13")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 6 Salix polaris is neither found in DS-VASCAN or DS-ZACKPLAN for our primer pairs

# 7 Salix reticulata MKPCH108-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH108-09")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 8 Salix richardsonii MKPCH606-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH606-09")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 9 Saxifraga aizoides MKPCH404-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH404-09")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 10 Saxifraga cespitosa GRAFW2116-13 form DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW2116-13")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 11 Saxifraga oppositifolia GRAFW1873-13 from DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode is complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW1873-13")
) %>%
    select(-processid, -inst) %>%
    bind_rows(assmann_et_al)

# 12 Woodsia glabella is found in neither dataset

# Arrange by name
assmann_et_al <- arrange(assmann_et_al, short_name)

# Check all gaps were filled where possible
print(assmann_et_al, n = 200)
# 4 species without barcode
# Equisetum arvense
# Equisetum variegatum
# Salix polaris
# Woodsia glabella

# Remove NAs
assmann_et_al <- na.omit(assmann_et_al)

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
    rep(
        saarela %>%
            filter(genus == "Carex") %>%
            distinct(Clade, Family) %>%
            pull(Clade),
        4
    ),
    saarela %>%
        filter(genus == "Tofieldia") %>%
        distinct(Clade, Family) %>%
        pull(Clade)
)
assmann_et_al$family[which(is.na(assmann_et_al$family))] <- c(
    rep(
        saarela %>%
            filter(genus == "Carex") %>%
            distinct(Clade, Family) %>%
            pull(Family),
        4
    ),
    saarela %>%
        filter(genus == "Tofieldia") %>%
        distinct(Clade, Family) %>%
        pull(Family)
)

# Check for nas
filter(assmann_et_al, is.na(clade))

# Write out
write_excel_csv(
    assmann_et_al,
    "data/sequences/its2/reference/Assmann_et_al_ITS2_ref.csv"
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
    filepath = "data/sequences/its2/reference/Assmann_et_al_ITS2_ref.fasta"
)
# Create BLAST database in terminal (external command)
# makeblastdb -in data/sequences/its2/reference/Assmann_et_al_ITS2_ref.fasta -dbtype nucl -out data/sequences/its2/reference/Assmann_et_al_ITS2_ref_db

#### Test reference database against itself

# Generate an indexed fasta for testing
names(sequences) <- paste0(1:length(names(sequences)), ";", names(sequences))
writeXStringSet(
    sequences,
    filepath = "data/sequences/its2/reference/Assmann_et_al_ITS2_test.fasta"
)
# Test blast the database against itself
blastn -query data/sequences/its2/reference/Assmann_et_al_ITS2_test.fasta \
       -db data/sequences/its2/reference/Assmann_et_al_ITS2_ref_db \
       -out data/sequences/its2/reference/test_blast.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -max_target_seqs 5 \
       -perc_identity 95

# Read test result
blast_test <- read_delim("data/sequences/its2/reference/test_blast.txt", delim="\t", col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
blast_test %>% separate(qseqid, c("id", "clade_id", "family_id", "genus_id", "species_id"), sep = ";") %>% view()

# Qulaity filter test result
top_hits <- blast_test %>%
   # For 180 BP match minimum length for plant
  filter(length >= 180) %>%
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

## EOF