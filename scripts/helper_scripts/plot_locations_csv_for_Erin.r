# Generate CSV of plot locations for Erin
# Jakob Assmann jakob.assmann@uzh.ch 27 November 2025

# Dependencies
library(tidyverse)
library(sf)

# Read in plot locations 
plots <- read_sf("data/geometries/plot_locations.gpkg")
plots %>% 
select(Name, geom) %>%
mutate(Easting = st_coordinates(geom)[,1],
       Northing = st_coordinates(geom)[,2]) %>%
       st_transform(4326) %>%
mutate(Longitude = st_coordinates(geom)[,1],
       Latitude = st_coordinates(geom)[,2]) %>%
       st_drop_geometry() %>%
       write_csv("data/geometries/plot_locations.csv")


# Load point framing data and summarise
point_frame <- read_csv("data/botanical_obs/pointframing.csv")
cover <- point_frame %>% group_by(plot_id, full_name) %>% 
    summarise(n = n()) %>%
    mutate(perc_plot = n / sum(n)) %>%
    split(.$plot_id) %>%
    map(function(x){
        sum_all <- sum(x$n)
        mutate(x, perc_plot = round(n / sum_all * 100)) %>%
            arrange(desc(n)) %>%
            select(-n)
    }) %>% bind_rows()
write_csv(cover, "data/botanical_obs/point_framing_cover.csv")
