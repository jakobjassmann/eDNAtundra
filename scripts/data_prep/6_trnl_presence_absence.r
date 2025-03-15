# Retrieve presence absence for each plot based on trnl
# Jakob Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(dada2)
library(pbapply)
library(parallel)

# Set up parallel environment
cl <- makeCluster(15)
clusterEvalQ(cl, {
    library(dplyr)
})

# Get file names for all trnl reads
trnl_reads <- list.files("data/sequences/trnl/edna/trimmed",
    pattern = "fastq",
    full.names = T
)

# # Remove any error samples
# trnl_reads <- trnl_reads[!grepl("error", trnl_reads)]

# Get sample names
sample_ids[1:10] <- gsub("([A-Za-z0-9]*)-.*", "\\1", basename(trnl_reads))

# De-replicate reads
trnl_derep <- derepFastq(trnl_reads[1:10])

# Convert to sequence table
trnl_seqtab <- makeSequenceTable(trnl_derep)
str(trnl_seqtab)

taxa <- assignSpecies(trnl_seqtab,
    "data/sequences/reference/trnl_reference_Assmann_et_al_2025.fasta",
    verbose = T
)
sample_derep <- trnl_derep[[2]]

# For each sample -> retrieve taxonomy:
species_detected <- trnl_derep %>%
    pblapply(cl = cl, function(sample_derep) {
        # Convert to sequence table
        sample_seqtab <- dada2::makeSequenceTable(sample_derep)

        # Convert sequence table to data frame for processing later
        sample_seqtab_df <- data.frame(t(sample_seqtab)) %>%
            mutate(., sq = row.names(.))
        names(sample_seqtab_df) <- c("Count", "sq")
        row.names(sample_seqtab_df) <- NULL

        # Assign taxonomy
        taxa <- dada2::assignSpecies(
            sample_seqtab,
            "data/sequences/reference/trnl_reference_Assmann_et_al_2025.fasta"
        )
        # Convert taxnonmy into a data frame of unique species
        taxa_df <- taxa %>%
            data.frame() %>%
            filter(!is.na(Genus)) %>%
            mutate(., sq = row.names(.))

        # Remove rownames
        row.names(taxa_df) <- NULL

        # add counts
        taxa_df <- left_join(taxa_df, sample_seqtab_df, join_by(sq == sq)) %>%
            select(-sq)

        # and return
        return(taxa_df)
    })


# Collapse list into a  dataframe
species_detected_df <- map(1:length(trnl_derep), function(i) {
    species_detected[[i]] %>%
        mutate(plot_id = plot_ids[i])
}) %>%
    bind_rows() %>%
    relocate(plot_id)

# Filter out ids with low read numbers
species_detected_df_filtered <- species_detected_df %>%
filter(Count >= 10)

# Stop cluster
stopCluster(cl)

# Compare to botanical obs
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv")

botanical_obs %>%
    filter(plot_id %in% unique(substr(species_detected_df_filtered$plot_id, 1, 4))) %>%
    group_by(plot_id) %>%
    summarise(n_spp = n())

species_detected_df_filtered %>%
    mutate(plot_id = substr(plot_id, 1, 4)) %>%
    distinct(plot_id, Genus, Species) %>%
    group_by(plot_id) %>%
    summarise(n_spp = n())

test_plot <- unique(substr(species_detected_df_filtered$plot_id, 1, 4))[3]
botanical_obs %>%
    filter(plot_id %in% test_plot) %>%
    arrange(plot_id, full_name) %>%
    select(plot_id, full_name)
species_detected_df_filtered %>%
    mutate(plot_id = substr(plot_id, 1, 4)) %>%
    filter(plot_id %in% test_plot) %>%
    group_by(plot_id, Genus, Species) %>%
    summarise(total_count = sum(Count))
    arrange(Genus, Species)
