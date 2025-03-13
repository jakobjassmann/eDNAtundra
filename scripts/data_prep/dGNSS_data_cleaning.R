# Data cleaning of the eDNA project dGNSS surveys
# Jakob J. Assmann jakob.assmann@uzh.ch 2023-08-15

# Dependencies
library(tidyverse)
library(sf)

# Load files
ppk_data <- list.files("data/geometries/ppk_results", full.names = T) %>%
    read_csv() %>%
    bind_rows() %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(crs = 32613) %>%
    mutate(site = gsub("(S[1-3]).*", "\\1", Name)) %>%
        mutate(site = case_when(
            site == "S1" ~ "Site 1", 
            site == "S2" ~ "Site 2",
            site == "S3" ~ "Site 3"
        ))
plot(st_geometry(ppk_data))

# Subset plot locations only
plot_locations <- ppk_data %>%
    filter(!(grepl("NW", Name) |
        grepl("NE", Name) |
        grepl("SE", Name) |
        grepl("SW", Name))) 
write_sf(plot_locations, "data/geometries/plot_locations.gpkg")

# Get site boundaries
site1_bounds <- ppk_data %>%
    filter(site == "Site 1") %>%
    summarise() %>%
    st_convex_hull() %>%
    mutate(site = "Site 1")
plot(st_geometry(site1_bounds))
st_area(site1_bounds)
site2_bounds <- ppk_data %>%
    filter(site == "Site 2") %>%
    summarise() %>%
    st_convex_hull() %>%
    mutate(site = "Site 2")
plot(st_geometry(site2_bounds))
st_area(site2_bounds)
site3_bounds <- ppk_data %>%
    filter(site == "Site 3") %>%
    summarise() %>%
    st_convex_hull() %>%
    st_bbox()
site3_bounds <- (site3_bounds - c(1, 0, 0, 0)) %>%
    st_as_sfc() %>%
    st_as_sf() %>%
    mutate(
        site = "Site 3"
    ) %>%
        rename(geometry = x)
plot(st_geometry(site3_bounds))
st_area(site3_bounds)

# Combine and safe
site_boundaries <- bind_rows(site1_bounds, site2_bounds, site3_bounds)
write_sf(site_boundaries, "data/geometries/site_boundaries.gpkg")
