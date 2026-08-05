# Heat maps of summed ASV frequencies for each taxs, repeated for ITS2 and trnl
# Jakob J. Assmann jakob.assmann@gu.se 4 August 2026

# Dependencies
library(tidyverse)
library(cowplot)

## ITS2

# Load asv observations
its2_taxa <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa.csv") %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id))

# Plot site and plot level detection heatmap
its2_taxa <- its2_taxa %>%
    mutate(
        short_name = paste(genus, if_else(is.na(species), "sp.", species))
    ) %>%
    group_by(site_id, plot_id, subsample_id, sample_id, short_name) %>%
    summarise(count = sum(count))

# Rearrange short name for plotting
unique_spp <- sort(unique(its2_taxa$short_name), decreasing = TRUE)
its2_taxa <- mutate(
    its2_taxa,
    short_name = factor(short_name, levels = unique_spp)
)

# Determine global limits for plotting
its2_taxa %>%
    group_by(site_id, sample_id, short_name) %>%
    summarise(count_sum = sum(count)) %>%
    group_by(site_id) %>%
    summarise(count_sum_max = max(count_sum))
#   site_id count_sum_max
#   <chr>           <dbl>
# 1 Site 1         100971
# 2 Site 2          96724
# 3 Site 3         100359
# => 103000 seems a reasonable limit for a global scale when also considering
# trnl below

# Plot each site
map(unique(its2_taxa$site_id), function(site_id_to_plot) {
    # site_id_to_plot <- "Site 1"
    species_count_plot <- ggplot(filter(
        its2_taxa,
        site_id == site_id_to_plot
    )) +
        geom_tile(
            aes(
                x = sample_id,
                y = short_name,
                fill = count
            ),
            width = 0.8,
            height = 0.6,
            color = "white",
            linewidth = 0.1,
            linejoin = "mitre",
            lineend = "square"
        ) +
        geom_vline(xintercept = seq(4.5, 16 * 4.5, 4), colour = "white") +
        scale_x_discrete(position = "top") +
        scale_y_discrete(drop = FALSE) +
        scale_fill_continuous(
            limits = c(1, 103000),
            trans = "log10",
            breaks = c(1, 10, 100, 1000, 10000, 100000),
            labels = c("1", "10", "100", "1000", "10000", "100000")
        ) +
        labs(
            x = "",
            y = "",
            fill = "Total count\nof ASVs",
            title = paste0(
                site_id_to_plot,
                " - ",
                "Sum of ITS2 ASV counts per taxon"
            )
        ) +
        # coord_flip() +
        theme_cowplot() +
        theme(
            panel.background = element_rect(fill = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 0.5)
        )
    save_plot(
        paste0("figures/ASV_counts_ITS2_", site_id_to_plot, ".png"),
        species_count_plot,
        base_height = 12,
        bg = "white"
    )
})

## trnl

# Load asv observations
trnl_taxa <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa.csv") %>%
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(site_id = gsub("S", "Site ", site_id))

# Plot site and plot level detection heatmap
trnl_taxa <- trnl_taxa %>%
    mutate(
        short_name = paste(genus, if_else(is.na(species), "sp.", species))
    ) %>%
    group_by(site_id, plot_id, subsample_id, sample_id, short_name) %>%
    summarise(count = sum(count))

# Rearrange short name for plotting
unique_spp <- sort(unique(trnl_taxa$short_name), decreasing = TRUE)
trnl_taxa <- mutate(
    trnl_taxa,
    short_name = factor(short_name, levels = unique_spp)
)

# Determine global limits for plotting
trnl_taxa %>%
    group_by(site_id, sample_id, short_name) %>%
    summarise(count_sum = sum(count)) %>%
    group_by(site_id) %>%
    summarise(count_sum_max = max(count_sum))
#   site_id count_sum_max
#   <chr>           <dbl>
# 1 Site 1         100971
# 2 Site 2          96724
# 3 Site 3         100359
# => 103000 seems a reasonable limit for a global scale

map(unique(trnl_taxa$site_id), function(site_id_to_plot) {
    # site_id_to_plot <- "Site 1"
    species_count_plot <- ggplot(filter(
        trnl_taxa,
        site_id == site_id_to_plot
    )) +
        geom_tile(
            aes(
                x = sample_id,
                y = short_name,
                fill = count
            ),
            width = 0.8,
            height = 0.6,
            color = "white",
            linewidth = 0.1,
            linejoin = "mitre",
            lineend = "square"
        ) +
        geom_vline(xintercept = seq(4.5, 16 * 4.5, 4), colour = "white") +
        scale_x_discrete(position = "top") +
        scale_y_discrete(drop = FALSE) +
        scale_fill_continuous(
            trans = "log10",
            breaks = c(1, 10, 100, 1000, 10000, 100000),
            limits = c(1, 103000),
            labels = c("1", "10", "100", "1000", "10000", "100000")
        ) +
        labs(
            x = "",
            y = "",
            fill = "Total count\nof ASVs",
            title = paste0(
                site_id_to_plot,
                " - ",
                "Sum of trnL ASV counts per taxon"
            )
        ) +
        # coord_flip() +
        theme_cowplot() +
        theme(
            panel.background = element_rect(fill = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 0.5)
        )
    save_plot(
        paste0("figures/ASV_counts_trnl_", site_id_to_plot, ".png"),
        species_count_plot,
        base_height = 12,
        bg = "white"
    )
})

# End of File
