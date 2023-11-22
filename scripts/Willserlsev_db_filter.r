# Quick script to check out sequence similarity for the
# species observed Jakob J. Assmann j.assmann@uzh.ch 10 Oct 2023

# Depencencies
library(tidyverse)
library(tidysq)


# Download FASTA data base from Willerslev et al. 2014
temp <- tempfile()
download.file("https://datadryad.org/stash/downloads/file_stream/86383", temp)
willerslev_db <- read_fasta(temp, ignore_case = T) 

# # Load species list and preformat
# botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
#     distinct(species)
# write_csv(botanical_obs, "data/species_list.csv")
# # manual data cleaning #1
# botanical_obs <- read_csv("data/species_list.csv") 
# botanical_obs$error_species <- NA
# botanical_obs$error_genus <- NA
# for (i in 1:nrow(botanical_obs)) {
#     botanical_obs$error_species[i] =
#         grepl(botanical_obs$species[i], botanical_obs$full_name[i])
#     botanical_obs$error_genus[i] =
#         grepl(botanical_obs$genus[i], botanical_obs$full_name[i])
# }
# # Control for spelling mistakes - check#1
# botanical_obs %>% filter(!error_species)
# botanical_obs %>% filter(!error_genus)
# # Filter distinct
# botanical_obs <- botanical_obs %>%
#     distinct(genus, species, subspecies) %>%
#     arrange(genus, species, subspecies)
# # Write out
#write_csv(botanical_obs, "data/species_list_cleaned.csv")

# Load botanical obs
botanical_obs <- read_csv("data/species_list_cleaned.csv") 
tally(botanical_obs)

# Add combined species genus colum
botanical_obs <- mutate(botanical_obs, ge_sp = paste(genus, species))

# Filter fasta data by direct matches for species and genus
willerslev_db_first_pass <- botanical_obs %>%
    split(1:nrow(.)) %>%
    map(function(obs) {
        index <- grepl(obs$ge_sp, willerslev_db$name)
        willerslev_db[index, ] %>%
            mutate(full_name = obs$ge_sp, pass = 1) %>%
            slice(1)
    }) %>%
    bind_rows()
tally(willerslev_db_first_pass)
53-37
# 16 sequences are not matched directly!

# # Identify missing
# botanical_obs_second_pass <- botanical_obs %>%
#     filter(!(ge_sp %in% willerslev_db_first_pass$full_name))

# # Go thrugh the missing sequences by genus
# willerslev_db_second_pass <- botanical_obs_second_pass %>%
#     split(., 1:nrow(.)) %>%
#     map(function(obs) {
#         cat(paste0(rep("#", 60), collapse = ""), "\n")
#         cat("Full name:", obs$ge_sp, "\n\n")
#         cat("Genus:", obs$genus, "| Species:", obs$species, "\n\n")
#         index <- grepl(obs$genus, willerslev_db$name)
#         cat(sum(index), "matching DB entries matching genus.\n\n")
#         cat(paste0(1:sum(index), ") ", paste0(willerslev_db$name[index]), "\n"), "\n")
#         cat("\nChoose option [", 1:sum(index), "or n for none]:\n" )
#         response <- readline()
#         while (!response %in% c("n", 1:sum(index))) {
#             cat("\nChoose option [", 1:sum(index), "or n for none]:\n")
#             response <- readline()
#         }
#         if (response == "n") {
#             return(NULL)
#         } else {
#             return(willerslev_db[index[as.numeric(response)], ] %>%
#                 mutate(full_name = obs$ge_sp, pass = 2))
#         }
#         cat("\n\n")
#     }) %>% bind_rows() %>%
#     bind_rows(willerslev_db_first_pass, .)
# tally(willerslev_db_second_pass)
# 53-42 
# # 11 species for manual assessment
# # save species with known sequences
# save(willerslev_db_second_pass, file = "data/willerslev_db_second_pass.Rda")
# write_csv(willerslev_db_second_pass, file = "data/willerslev_db_second_pass.csv")

# Get not included and set sequence to NA
botanical_obs_third_pass <- botanical_obs %>%
    filter(!(ge_sp %in% willerslev_db_first_pass$full_name))

# Merge with second pass
willerslev_db_third_pass <- botanical_obs_third_pass %>%
    mutate(full_name = ge_sp, pass = 3, sq = sq("AA"), name = NA) %>%
    select(sq, name, full_name, pass) %>%
    bind_rows(willerslev_db_second_pass)

# look up species manually 
willerslev_db_third_pass %>% filter(is.na(name)) %>% pull(full_name)

# Write quick helper function to replace sequence and name description for a species
replace_sq <- function(sp_name, seq, name_desc){
    row_content <- willerslev_db_third_pass %>% filter(full_name == sp_name)
    willerslev_db_third_pass <- willerslev_db_third_pass %>% filter(!(full_name == sp_name)) # nolint: object_usage_linter.
    new_row <- tibble(
        sq = seq,
        name = name_desc,
        full_name = sp_name,
        pass = 3
    )
    willerslev_db_third_pass <- bind_rows(willerslev_db_third_pass, new_row)
    rm(new_row)
    rm(row_content)
    return(willerslev_db_third_pass)
}

# Add colum for synonyms
willerslev_db_third_pass$synoynm <- NA

# Species:      Arctagrostis latifolia
# Family:       Poaceae
# Tribe         Poeae
# Synonyms:     Colpodium latifolium
willerslev_db$name[grepl("Arctagrostis", willerslev_db$name)]
willerslev_db$name[grepl("Colpodium", willerslev_db$name)]
willerslev_db$name[grepl("Poaceae", willerslev_db$name)]
willerslev_db$name[grepl("Po", willerslev_db$name)]
willerslev_db$name[grepl("group", willerslev_db$name)]
willerslev_db$name[grepl("Poeae", willerslev_db$name)]
# "Pooideae_GROUP_1" contains the Poeae tribe select this:
willerslev_db_third_pass <- replace_sq("Arctagrostis latifolia",
    willerslev_db$sq[grepl("Pooideae_GROUP_1", willerslev_db$name)],
    willerslev_db$name[grepl("Pooideae_GROUP_1", willerslev_db$name)])

# Species:      Carex borealipolaris
# Synonym:      Kobresia sibirica 
willerslev_db$name[grepl("borealipolaris", willerslev_db$name)]
willerslev_db$name[grepl("Kobresia sibirica", willerslev_db$name)]
# Match!
willerslev_db_third_pass <- replace_sq("Carex borealipolaris",
    willerslev_db$sq[grepl("Kobresia sibirica", willerslev_db$name)],
    willerslev_db$name[grepl("Kobresia sibirica", willerslev_db$name)])
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Carex borealipolaris"] <- "Kobresia sibirica"

# Species:      Carex simpliciuscula
# Synonym:      Kobresia simpliciuscula
willerslev_db$name[grepl("Kobresia simpliciuscula", willerslev_db$name)]
# Match!
willerslev_db_third_pass <- replace_sq(
    "Carex simpliciuscula",
    willerslev_db$sq[grepl("Kobresia simpliciuscula", willerslev_db$name)],
    willerslev_db$name[grepl("Kobresia simpliciuscula", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Carex simpliciuscula"] <- "Kobresia simpliciuscula"

# Species:      Carex myosuroides
# Synonym:      Kobresia myosuroides
willerslev_db$name[grepl("Kobresia myosuroides", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Carex myosuroides",
    willerslev_db$sq[grepl("Kobresia myosuroides", willerslev_db$name)],
    willerslev_db$name[grepl("Kobresia myosuroides", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Carex myosuroides"] <- "Kobresia myosuroides"

# Species:      Chamaenerion latifolium
# Synonyms:     Epilobium latifolium
#               Chamerion latifolium
willerslev_db$name[grepl("Epilobium latifolium", willerslev_db$name)]
willerslev_db$name[grepl("Epilobium", willerslev_db$name)]
willerslev_db$name[grepl("Chamerion latifolium", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Chamaenerion latifolium",
    willerslev_db$sq[grepl("Chamerion latifolium", willerslev_db$name)],
    willerslev_db$name[grepl("Chamerion latifolium", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Chamaenerion latifolium"] <- "Chamerion latifolium"

# Species:      Dupontia fisheri
willerslev_db$name[grepl("Dupontia", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Dupontia fisheri",
    willerslev_db$sq[grepl("Dupontia", willerslev_db$name)],
    willerslev_db$name[grepl("Dupontia", willerslev_db$name)]
)

# Species:      Eriophorum triste
# Synoynm:      Eriophorum angustifolium subsp.
willerslev_db$name[grepl("Eriophorum angustifolium", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Eriophorum triste",
    willerslev_db$sq[grepl("Eriophorum angustifolium", willerslev_db$name)],
    willerslev_db$name[grepl("Eriophorum angustifolium", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Eriophorum triste"] <- "Eriophorum angustifolium"

# Species:      Oxytropis arctobia
# Synonyms:     Oxytropis arctica
#               Oxytropis nigrescens subsp. arctobia / uniflora
willerslev_db$name[grepl("Oxytropis nigrescens", willerslev_db$name)]
willerslev_db$name[grepl("Oxytropis arctica", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Oxytropis arctobia",
    willerslev_db$sq[grepl("Oxytropis arctica", willerslev_db$name)],
    willerslev_db$name[grepl("Oxytropis arctica", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Oxytropis arctobia"] <- "Oxytropis arctica"

# Species:      Sabulina rubella
# Synonyms:     Minuartia rubella (and many others)
willerslev_db$name[grepl("Minuartia rubella", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Sabulina rubella",
    willerslev_db$sq[grepl("Minuartia rubella", willerslev_db$name)],
    willerslev_db$name[grepl("Minuartia rubella", willerslev_db$name)]
)
willerslev_db_third_pass$synoynm[willerslev_db_third_pass$full_name == "Sabulina rubella"] <- "Minuartia rubella"

# Salix species:
#                   Salix arctica
#                   Salix reticulata
#                   Salix richardsonii
#                   Salix polaris
willerslev_db$name[grepl("Salix", willerslev_db$name)]
willerslev_db_third_pass <- replace_sq(
    "Salix arctica",
    willerslev_db$sq[grepl("Salix", willerslev_db$name)],
    willerslev_db$name[grepl("Salix", willerslev_db$name)]
)
willerslev_db_third_pass <- replace_sq(
    "Salix reticulata",
    willerslev_db$sq[grepl("Salix", willerslev_db$name)],
    willerslev_db$name[grepl("Salix", willerslev_db$name)]
)
willerslev_db_third_pass <- replace_sq(
    "Salix richardsonii",
    willerslev_db$sq[grepl("Salix", willerslev_db$name)],
    willerslev_db$name[grepl("Salix", willerslev_db$name)]
)
willerslev_db_third_pass <- replace_sq(
    "Salix polaris",
    willerslev_db$sq[grepl("Salix", willerslev_db$name)],
    willerslev_db$name[grepl("Salix", willerslev_db$name)]
)
# Species:      Saxifraga tricuspidata
# Synonyms:     Leptasea tricuspidata
willerslev_db$name[grepl("Leptasea tricuspidata", willerslev_db$name)]
willerslev_db$name[grepl("tricuspidata", willerslev_db$name)]
willerslev_db$name[grepl("Saxifraga", willerslev_db$name)]
willerslev_db$name[grepl("Trachyphyllum", willerslev_db$name)]
# Not in database (longer trnL sequences on GenBank)

# Species:      Tofielda coccinea
# !!! Typo in genus !!! Tofieldia coccinea
willerslev_db$name[grepl("Tofieldia coccinea", willerslev_db$name)]
# Check!
willerslev_db_third_pass$full_name[willerslev_db_third_pass$full_name == "Tofielda coccinea"] <- "Tofieldia coccinea"
willerslev_db_third_pass <- replace_sq(
    "Tofieldia coccinea",
    willerslev_db$sq[grepl("Tofieldia coccinea", willerslev_db$name)],
    willerslev_db$name[grepl("Tofieldia coccinea", willerslev_db$name)]
)

# Species:      Woodsia glabella
# Synonyms:     Woodsia alpina var. glabella
#               Woodsia hyperborea var. glabella
#               Woodsia lapponica
#               Woodsia yazawai
willerslev_db$name[grepl("Woodsia", willerslev_db$name)]
# Go up taxonomy from Wikipeida
willerslev_db$name[grepl("Woodsiaceae", willerslev_db$name)]
willerslev_db$name[grepl("Polypodiales", willerslev_db$name)]
willerslev_db$name[grepl("Polypodiopsida", willerslev_db$name)]
willerslev_db$name[grepl("Tracheophyta", willerslev_db$name)]
# Not in data base, but trnl L-F available on GenBank

# Check final sequences (should be 51 known)
willerslev_db_third_pass %>%
    filter(!is.na(name)) %>%
    tally()

# Check unique sequences
willerslev_db_third_pass$is_unique <- willerslev_db_third_pass %>%
    split(., 1:nrow(.)) %>%
    sapply(function(tibble_sub){
        seq <- pull(tibble_sub, sq)
        if (sum(willerslev_db_third_pass$sq == seq) > 1) {
            return("Not unique")
        } else {
            return("Unique")
        }
    })

# add seq_lenght
willerslev_db_third_pass$seq_lenght <- get_sq_lengths(willerslev_db_third_pass$sq)

# Rename colums
willerslev_db_third_pass <- willerslev_db_third_pass %>%
    select(
        species_name = full_name,
        synonym = synoynm,
        willerslev_meta = name,
        trnL_seq = sq,
        trnL_seq_lenght = seq_lenght,
        trnL_unique = is_unique
    )

# Check resolution:
(sum(willerslev_db_third_pass$trnL_unique == "Not unique") - 2) / 51
# 37% of species not unqiue!!!

# Safe Rda
save(willerslev_db_third_pass, file = "data/willerslev_db_subset_final.Rda")

# Write out CSV
write_csv(willerslev_db_third_pass, file = "data/willerslev_db_subset_final.csv")
