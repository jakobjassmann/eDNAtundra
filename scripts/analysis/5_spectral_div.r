# Extract spectral diversity from drone data for plots
# Jakob J. Assmann jakob.assmann@uzh.ch March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(terra)
library(tidyterra)
library(sf)

# Load plot locations
plot_coords <- read_sf("data/geometries/plot_locations.gpkg")

# Load drone rasters
rgb_rasts <- list.files("data/drone_data/rgb/", full.names = T, pattern = ".tif$") %>%
.[!grepl("crop", .)] %>%
map(rast)

# Buffer plots by 0 5 m 
plots_buffered <- plot_coords %>%
st_buffer(0.5) %>%
mutate(site = substr(Name, 1,2) %>% gsub("S", "Site ", .))

# Extract RGB data for a 0.5 m buffer around each plot
plot_rgb <- map(1:3, function(x){
    site_plots <- filter(plots_buffered, site == paste0("Site ", x))
     plots_rgb <- terra::extract(rgb_rasts[[x]], site_plots) %>%
     select(1:4) %>% setNames(c("ID", "R", "G", "B")) %>%
     pivot_longer(2:3, names_to = "band", values_to = "value") %>%
     group_by(ID) %>% summarise(spectral_cv =  (sd(value) / mean(value)))
     plots_rgb$plot_id <- site_plots$Name
     plots_rgb$site_id <- site_plots$site
     plots_rgb <- relocate(plots_rgb, site_id, plot_id)
     return(plots_rgb)
}) %>% bind_rows()

# Load presence / absence vascular observations
vascular_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    ) %>%
    # Add site_id string
    rename(site_id = site) %>%
    select(-date, -full_name) %>%
    mutate(phylum = "Vascular plants") 

# Load presence / absence moss observations
bryo_obs <- read_csv("data/botanical_obs/moss_data_Erin_cleaned_final_Deb.csv") %>%
    mutate(plot_id = paste0(Site, Plot_ID)) %>%
    select(site_id = Site, plot_id, genus = Genus, species = Species) %>%
    # convert al sp. to species
    mutate(species = (gsub("sp\\.", "species", species))) %>%
    # Remove all uncertain ids
    filter(!grepl("cf", genus) & !grepl("cf", species)) %>%
    mutate(short_name = paste(genus, species)) %>%
    arrange(site_id, plot_id, short_name) %>%
    mutate(site_id = gsub("S", "Site ", site_id)) %>%
    mutate(phylum = "Bryophytes") 

# Combine dataframes
combined_obs <- bind_rows(vascular_obs, bryo_obs) %>%
arrange(site_id, plot_id, short_name, phylum) %>%
relocate(site_id) %>%
group_by(plot_id) %>%
summarise(n_spp = n_distinct(short_name)) %>%
full_join(plot_rgb)

# Correlate the the species richness with the spectral diversity
ggplot() +
geom_point(data = combined_obs, aes(x = n_spp, spectral_cv)) +
theme_classic()

# Point framing 
point_framing <- read_csv("data/botanical_obs/pointframing.csv") %>%
select(-observer, -scribe, -date) %>%
group_by(plot_id, full_name) %>%
summarise(sum_hits = n()) %>%
split(., .$plot_id) %>%
map(function(x){
    x <- mutate(x, p = sum_hits/sum(x$sum_hits),
p_ln_p = p * log(p)) %>%
summarise(shan_div = -1* sum(p_ln_p))
}) %>% bind_rows()

# Combine
combined_obs <- full_join(combined_obs, point_framing)

# Correlate the the cover richness with the spectral diversity
ggplot(data = combined_obs, mapping = aes(x = shan_div, spectral_cv)) +
geom_point() +
geom_smooth(method = "lm") +
theme_classic()
