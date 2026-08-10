# Convert GCP coordinates to UTM13N for Pix4D
# Jakob Assmann 10 August 2026

# Dpendencies
library(tidyverse)
library(sf)

# Load coordinates
site_1 <- read_csv(
    "data/geometries/ppk_results/eDNA Site 1 GCPs for Pix4D.csv"
) %>%
    st_as_sf(coords = c(x = "X", y = "Y"), crs = 4326)
site_2 <- read_csv(
    "data/geometries/ppk_results/eDNA Site 2 GCPs for Pix4D.csv"
) %>%
    st_as_sf(coords = c(x = "X", y = "Y"), crs = 4326)
site_3 <- read_csv(
    "data/geometries/ppk_results/eDNA Site 3 GCPs for Pix4D.csv"
) %>%
    st_as_sf(coords = c(x = "X", y = "Y"), crs = 4326)

# reporject to 32613
site_1 <- st_transform(site_1, 32613)
site_2 <- st_transform(site_2, 32613)
site_3 <- st_transform(site_3, 32613)

# Re_export
site_1 %>%
    mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(Name, X, Y, Z, Acc_XY, Acc_Z) %>%
    write_csv(
        "data/geometries/ppk_results/eDNA Site 1 GCPs for Pix4D UTM13N.csv"
    )
site_2 %>%
    mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(Name, X, Y, Z, Acc_XY, Acc_Z) %>%
    write_csv(
        "data/geometries/ppk_results/eDNA Site 2 GCPs for Pix4D UTM13N.csv"
    )
site_3 %>%
    mutate(X = st_coordinates(.)[, 1], Y = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(Name, X, Y, Z, Acc_XY, Acc_Z) %>%
    write_csv(
        "data/geometries/ppk_results/eDNA Site 3 GCPs for Pix4D UTM13N.csv"
    )
