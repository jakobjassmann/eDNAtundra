# Compile reference library files (FASTA)
# Jakob J. Assmann jakob.assmann@uzh.ch 13 March 2025

# Dependencies
library(tidyverse)
library(tidysq)
library(msa)

# trnl ----

# Load trnl sequences from cambridge bay
cb_sqs <- read_csv("data/sequences/reference/Assmann_et_al_CB_2023.csv") %>%
    select(sample, species, trnl_seq)

# Helper function to check uniquness of sample sequences in set of samples
get_set_of_unique_spp <- function(current_sample) {
    # Check if sequence if present, if not return NA
    if (is.na(current_sample$trnl_seq)) {
        return(NA)
    }
    # Compare current sample to all other samples
    matching_samples <- cb_sqs %>%
        # Get all samples that are not NA
        filter(!is.na(trnl_seq)) %>%
        split(.$sample) %>%
        map(function(comp_sample) {
            # Align sequences
            agreement <- msa(c(
                current_sample$trnl_seq,
                comp_sample$trnl_seq
            ), type = "dna", substitutionMatrix = "default") %>%
                # Determine fractional match from consensus matrix
                consensusMatrix() %>%
                (function(x) colSums(x == 2)) %>%
                (function(x) sum(x == 1) / length(x))
            # If match return sample ID
            if (agreement == 1) {
                return(comp_sample$sample)
            }
            # If not return NULL
            return(NULL)
        }) %>%
        unlist()
    # Clean up names
    names(matching_samples) <- NULL
    # Arrange alphabeticall for comparison and turn into a single vector
    matching_samples <- paste(sort(matching_samples), collapse = " ")
    # Return
    return(matching_samples)
}

# Apply to all sequences
cb_sqs$matches <- cb_sqs %>%
    split(.$sample) %>%
    map(get_set_of_unique_spp, .progress = TRUE) %>%
    unlist()

# Get overview on all speices
cb_sqs %>% print(n = 100)

# Get overview on groupings
cb_sqs %>%
    distinct(matches) %>%
    print(n = 20)
# 3 groupings of multiple species -> let's look at these
cb_sqs %>%
    filter(!is.na(trnl_seq)) %>%
    arrange(species) %>%
    group_by(matches) %>%
    filter(matches %in% c("A5 A7", "A8 B5 H5", "B1 B2 B4 B7", "C4 I5", "D3 H7")) %>%
    group_split()
# The Bistorta sequence must be contaminated with Carex -> remove this.
# The other two groupings are all Carex species, we will group these accordingly.
# The second last grouping are the two Drabas and the final grouping the
# two Eriophorum species.

# Generate genus-species binomials for the fasta files
cb_sqs_cleaned <- cb_sqs %>%
    # Remove contaminated Bistorta vivipara sample
    filter(!grepl("Bistorta", species)) %>%
    mutate(
        realised_taxon =
            case_when(
                # Assign Carex species to all clusters
                matches %in% c("A5 A7", "A8 B5 H5", "B1 B2 B4 B7") ~ "Carex species",
                matches %in% c("C4 I5") ~ "Draba species",
                matches %in% c("D3 H7") ~ "Eriophorum species",
                grepl("cf", species) ~ gsub("([A-Za-z]*) ([A-Za-z]*).*", "\\1 species", species),
                TRUE ~ gsub("([A-Za-z]*) ([A-Za-z]*).*", "\\1 \\2", species)
            )
    ) %>%
    distinct(realised_taxon, trnl_seq) %>%
    arrange(realised_taxon) %>%
    filter(!is.na(trnl_seq))
# Print result
cb_sqs_cleaned %>% print(n = 100)

# Load the curated subset of Willserslev et al. 2014
willerslev <- read_csv("data/external/reference_libraries/Willerslev_et_al_2014/Willerslev_et_al_2024_subset.csv") %>%
    distinct(realised_taxon, trnl) %>%
    select(realised_taxon, trnl_seq = trnl)

# Combine and keep only distinct sequences -> choosind Assmann et al. if both are identical
trnl_combined <- bind_rows(
    mutate(cb_sqs_cleaned, source = "Assmann et al"),
    mutate(willerslev, source = "Willerslev et al. 2014")
) %>%
    filter(!is.na(realised_taxon)) %>%
    group_by(realised_taxon, trnl_seq) %>%
    group_split() %>%
    map(function(x) {
        if (nrow(x) > 1) x <- filter(x, source == "Assmann et al")
        return(x)
    }) %>%
    bind_rows()

# Quick check:
length(unique(trnl_combined$trnl_seq))
nrow(distinct(trnl_combined, realised_taxon, source))
trnl_combined %>% group_by(trnl_seq) %>% filter(n()>1)

# Two more duplicates remain: 
# Carex aquatilis (Willerslev) groups with Carex aquatilis and Carex bigelowii (Assmann)
# -> we will merge these into Carex species (Assmann)
# Oxytropis maydelliana (Assmann) groups with Oxytropis species (Willerslev)
# -> we will merge thes into Oxytropis species (Willerslev)
trnl_combined %>% filter((realised_taxon == "Carex aquatilis"))
trnl_combined <- trnl_combined %>% filter(!(realised_taxon == "Carex aquatilis"))
trnl_combined %>% filter((realised_taxon == "Oxytropis maydelliana"))
trnl_combined <- trnl_combined %>% filter(!(realised_taxon == "Oxytropis maydelliana"))

# Show final list
print(trnl_combined, n = 100)

# Quick check:
length(unique(trnl_combined$trnl_seq))
nrow(distinct(trnl_combined, realised_taxon, source))
trnl_combined %>% group_by(trnl_seq) %>% filter(n()>1)
trnl_combined %>% group_by(realised_taxon) %>% filter(n()>1)

# Write to CSV
trnl_combined %>%
    write_csv("data/sequences/reference/trnl_reference_Assmann_et_al_2025.csv")
# trnl_combined <- read_csv("data/sequences/reference/trnl_reference_Assmann_et_al_2025.csv")

# Write to fasta file
trnl_combined %>%
    select(name = realised_taxon, sq = trnl_seq) %>%
    mutate(name = paste0("ID ", name), sq = sq(sq)) %>%
    write_fasta("data/sequences/reference/trnl_reference_Assmann_et_al_2025.fasta")


## ITS2 ---

# Load ITS2 sequences from cambridge bay
cb_sqs <- read_csv("data/sequences/reference/Assmann_et_al_CB_2023.csv") %>%
    select(sample, species, its2_seq)

# Helper function to check uniquness of sample sequences in set of samples
get_set_of_unique_spp <- function(current_sample) {
    # Check if sequence if present, if not return NA
    if (is.na(current_sample$its2_seq)) {
        return(NA)
    }
    # Compare current sample to all other samples
    matching_samples <- cb_sqs %>%
        # Get all samples that are not NA
        filter(!is.na(its2_seq)) %>%
        split(.$sample) %>%
        map(function(comp_sample) {
            # Align sequences
            agreement <- msa(c(
                current_sample$its2_seq,
                comp_sample$its2_seq
            ), type = "dna", substitutionMatrix = "default") %>%
                # Determine fractional match from consensus matrix
                consensusMatrix() %>%
                (function(x) colSums(x == 2)) %>%
                (function(x) sum(x == 1) / length(x))
            # If match return sample ID
            if (agreement == 1) {
                return(comp_sample$sample)
            }
            # If not return NULL
            return(NULL)
        }) %>%
        unlist()
    # Clean up names
    names(matching_samples) <- NULL
    # Arrange alphabeticall for comparison and turn into a single vector
    matching_samples <- paste(sort(matching_samples), collapse = " ")
    # Return
    return(matching_samples)
}

# Apply to all sequences
cb_sqs$matches <- cb_sqs %>%
    split(.$sample) %>%
    map(get_set_of_unique_spp, .progress = TRUE) %>%
    unlist()

# Get overview on all speices
cb_sqs %>% print(n = 100)

# Get overview on groupings
cb_sqs %>%
    distinct(matches) %>%
    print(n = 50)
# 3 groupings of multiple species -> let's look at these
cb_sqs %>%
    filter(!is.na(its2_seq)) %>%
    arrange(species) %>%
    group_by(matches) %>%
    filter(matches %in% c("A5 A7", "B2 I3")) %>%
    group_split()
# First grouping:
# It looks like there might have been some contamination between
# Carex aquatilis and Carex bigelowii (or that these two species are super close)
# They match 99.7% and 97.4% to the respective BOLD references.
# Likely indicating that we sampled C. aquatilis twice (or that there was)
# contamination. In this case, I'll just denote both two "Carex speices"
# Second grouping:
# It looks like C. cf. glaerosa is actually C. fuliginosa. -> we can remove it.

# Generate genus-species binomials for the fasta files
cb_sqs_cleaned <- cb_sqs %>%
    # Remove contaminated Bistorta vivipara sample
    filter(!grepl("cf. glareosa", species)) %>%
    mutate(
        realised_taxon =
            case_when(
                # Assign Carex species to all clusters
                matches %in% c("A5 A7") ~ "Carex species",
                grepl("cf", species) ~ gsub("([A-Za-z]*) ([A-Za-z]*).*", "\\1 species", species),
                TRUE ~ gsub("([A-Za-z]*) ([A-Za-z]*).*", "\\1 \\2", species)
            )
    ) %>%
    distinct(realised_taxon, its2_seq) %>%
    arrange(realised_taxon) %>%
    filter(!is.na(its2_seq))
# Print result
cb_sqs_cleaned %>% print(n = 100)

# Load the curated subset of Willserslev et al. 2014
bold <- read_csv("data/external/reference_libraries/BOLD/bold_sequences_2025_03_12.csv") %>%
    distinct(short_name, its2, bold_id) %>%
    select(realised_taxon = short_name, its2_seq = its2, source = bold_id) %>%
    filter(!is.na(its2_seq))

# Combine then for duplicate species-sequences entries keep only "Assmann et al.", then add a temporary id for redundancy checks
its2_combined <- bind_rows(mutate(cb_sqs_cleaned, source = "Assmann et al"), bold) %>%
    filter(!is.na(realised_taxon)) %>%
    group_by(realised_taxon, its2_seq) %>%
    group_split() %>%
    map(function(x) {
        if (nrow(x) > 1) x <- filter(x, source == "Assmann et al")
        return(x)
    }) %>%
    bind_rows() %>%
    arrange(realised_taxon, source) %>%
    mutate(id = 1:n())

# Check for redundancies
its2_combined$matches <- its2_combined %>%
    split(.$id) %>%
    map(function(current_row) {
        matching_ids <- its2_combined$id[grepl(current_row$its2_seq, its2_combined$its2_seq, ignore.case = TRUE)]
    })
# Manually check all combinations of rows with at least one other match and run msa
its2_combined %>%
    filter(sapply(matches, length) >= 2) %>%
    pull(matches) %>%
    unique() %>%
    map(function(x) {
        matches <- slice(its2_combined, x)
        print(matches)
        cat("1 matches 2", grepl(matches$its2_seq[1], matches$its2_seq[2], ignore.case = TRUE), "\n")
        cat("2 matches 1", grepl(matches$its2_seq[2], matches$its2_seq[1], ignore.case = TRUE), "\n")
        msa(matches$its2_seq, type = "DNA") %>% print(show = "complete")
        readline("Waiting for key input...")
        return(NULL)
    })
# BOLD reference from Cadamine digitata is longer but otherwise identical -> KEEP ONLY Assmann et al.
# Carex marina and Carex membranacea sequences are identical from BOLD:
its2_combined %>% filter(realised_taxon == "Carex marina")
its2_combined %>% filter(realised_taxon == "Carex membranacea")
#    -> both sequences from BOLD are incomplete and not super trustworthy (visual ID is okay)
#    -> there are sequences for both species from Assmann et al.
#    -> remove bold sequences from reference database
# Draba species from Assmann et al. is identical with Draba nivalis from BOLD
cb_sqs_cleaned %>% filter(grepl("Draba", realised_taxon))
its2_combined %>%
    filter(realised_taxon == "Draba nivalis") %>%
    pull(its2_seq) %>%
    msa(type = "DNA") %>%
    print(show = "complete")
#     -> the unidentified Draba species could have been a Draba nivalis
#     -> the Draba nivalis sequence from bold and our sequence differ only by one bas
#     -> this could be true in-species variation or it could be a miss identification
#     -> on either end. As the BOLD entry does not have a photo ID I will
#     -> Remove BOLD sequence to be save (i.e, stick to Draba species)
# Salix arctic and Salix polaris are identical from BOLD -> ASSIGN Salix species
# Salix reticulata and Salix richardsonii are identical from BOLD -> ASSIGN Salix species

# Clean up 1: remove entries
its2_combined <- its2_combined %>%
    filter(!(realised_taxon == "Cadamine digitata" & source != "Assmann et al")) %>%
    filter(!(realised_taxon == "Carex marina" & source != "Assmann et al")) %>%
    filter(!(realised_taxon == "Carex membranacea" & source != "Assmann et al")) %>%
    filter(!(realised_taxon == "Draba nivalis" & source != "Assmann et al"))

# Clean up: Genus level demotions
its2_combined <- its2_combined %>%
    mutate(realised_taxon = case_when(
        grepl("Salix", realised_taxon) & source != "Assmann et al" ~ "Salix species",
        grepl("Draba nivalis", realised_taxon) & source != "Assmann et al" ~ "Draba species",
        TRUE ~ realised_taxon
    ))

# Remove temporary ID for matching
its2_combined <- select(its2_combined, -id)

# De-dublicate once more keeping only the first enty of each group
its2_combined <- its2_combined %>%
    group_by(realised_taxon, its2_seq) %>%
    group_split() %>%
    map(function(x) {
        if (nrow(x) > 1) x <- filter(x, source == "Assmann et al")
        return(x)
    }) %>%
    bind_rows() %>%
    arrange(realised_taxon, source)

# Show final list
print(its2_combined, n = 100)

# Quick sanity check:
length(unique(its2_combined$its2_seq))
nrow(distinct(its2_combined, realised_taxon, source))


# Write as CSV
its2_combined %>%
    write_csv("data/sequences/reference/its2_reference_Assmann_et_al_2025.csv")
# its2_combined <- read_csv("data/sequences/reference/its2_reference_Assmann_et_al_2025.csv")

# Write to fasta file
its2_combined %>%
    select(name = realised_taxon, sq = its2_seq) %>%
    mutate(name = paste0("ID ", name), sq = sq(sq)) %>%
    arrange(name) %>%
    write_fasta("data/sequences/reference/its2_reference_Assmann_et_al_2025.fasta")
