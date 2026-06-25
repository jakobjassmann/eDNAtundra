# Build ITS2 reference libraries combinding the Cambridge Bay samples with the published databeses
# Jakob J. Assmann 2026-06-18

# Dependencies
library(tidyverse)

# Load species list from vascular plant observatins
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
    select(short_name = species, processid, sequence = nuc) %>%
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
    ) %>%
    select(-processid)

# Load DS-VASCAN and prepare for merging
dsvascan <- read_tsv(
    "data/external/reference_libraries/DS-VASCAN/DS-VASCAN.tsv"
) %>%
    filter(
        grepl("ITS_S2F", primers_forward) &
            grepl("ITS4", primers_reverse) &
            marker_code == "ITS2"
    ) %>%
    select(short_name = species, processid, sequence = nuc) %>%
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
    ) %>%
    select(-processid)

# Check species with missing sequences
filter(assmann_et_al, is.na(sequence))

# Bistorta vivipara MKPCH419-09 in DS-VASCAN
# Entry has a photograph that matches B vivipara
# Barcode appears complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH419-09")
) %>%
    bind_rows(assmann_et_al)

# Chamaenerion latifolium MKPCH428-09 in DS-VASCAN
# Entry has a photograph that matches C. latifolium
# Barcode appears complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH428-09")
) %>%
    bind_rows(assmann_et_al)

# Equisetum variegatum GRAFW1958-13 in DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW1958-13")
) %>%
    bind_rows(assmann_et_al)

# Equiesetum arvense is neither found in DS-VASCAN or DS-ZACKPLAN for our primer pairs

# Salix arctica GRAFW1877-13 in DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW1877-13")
) %>%
    bind_rows(assmann_et_al)

# Salix polaris is neither found in DS-VASCAN or DS-ZACKPLAN for our primer pairs

# Salix reticulata MKPCH108-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH108-09")
) %>%
    bind_rows(assmann_et_al)

# Salix richardsonii MKPCH606-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH606-09")
) %>%
    bind_rows(assmann_et_al)

# Saxifraga aizoides MKPCH404-09 from DS-VASCAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dsvascan,
    source == paste("DS-VASCAN", "MKPCH404-09")
) %>%
    bind_rows(assmann_et_al)

# Saxifraga cespitosa GRAFW2116-13 form DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW2116-13")
) %>%
    bind_rows(assmann_et_al)

# Saxifraga oppositifolia GRAFW1873-13 from DS-ZACKPLAN
# Entry has a photograph that matches species
# Barcode appears complete
assmann_et_al <- filter(
    dszackplan,
    source == paste("DS-ZACKPLAN", "GRAFW1873-13")
) %>%
    bind_rows(assmann_et_al)

# Woodsia glabella is found in neither dataset

# Merge all datasets
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
    "data/sequences/its2/reference/Assmann_et_al_ITS2_ref.csv"
)

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
# Check for unique sequnces, if identical select just one

unique(dszackplan$primers_forward)
any(dszackplan$genus %in% species_list$genus)
