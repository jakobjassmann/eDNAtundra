# Detection maps for vascular plants and bryophytes
# Site and plot level
# Jakob J. Assmann jakob.assmann@gu.se July 2027

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)

# Load vascular observations
obs_combined <- read_csv("data/analysis/detection_species_genus.csv")

# Site level detection map
(site_level_map <- ggplot(
    obs_combined %>%
        mutate(site_id = substr(plot_id, 1, 2)) %>%
        mutate(
            site_id = case_when(
                site_id == "S1" ~ "Site 1",
                site_id == "S2" ~ "Site 2",
                site_id == "S3" ~ "Site 3"
            )
        ) %>%
        select(-plot_id) %>%
        group_by(site_id, short_name) %>%
        summarise(
            botanical = any(botanical),
            eDNA = any(eDNA),
            eDNA_harm = any(eDNA_harm)
        ) %>%
        ungroup() %>%
        mutate(
            site_id = factor(site_id, levels = sort(unique(site_id))),
            short_name = factor(
                short_name,
                levels = sort(unique(short_name), decreasing = T)
            )
        )
) +
    geom_tile(
        aes(
            x = short_name,
            y = site_id,
            fill = paste(
                " botanical:",
                botanical,
                " eDNA_harm:",
                eDNA_harm,
                "eDNA:",
                eDNA
            ),
            colour = paste(
                " botanical:",
                botanical,
                " eDNA_harm:",
                eDNA_harm,
                "eDNA:",
                eDNA
            )
        ),
        width = 0.6,
        height = 0.3,
        linewidth = 1.25,
        linejoin = "mitre",
        lineend = "square"
    ) +
    scale_y_discrete(position = "right") +
    scale_colour_manual(
        values = c("#F28F16", "#A62121", "#223F59", "#223F59"),
        labels = rev(c(
            "botanical and eDNA",
            "botanical and eDNA harmonised",
            "botanical only",
            "eDNA only"
        ))
    ) +
    scale_fill_manual(
        values = c("#F28F16", "#A62121", "#A62121", "#223F59"),
        labels = rev(c(
            "botanical and eDNA",
            "botanical and eDNA harmonised",
            "botanical only",
            "eDNA only"
        ))
    ) +
    labs(x = "", y = "", colour = "Observed in", fill = "Observed in") +
    guides(
        colour = guide_legend(reverse = TRUE),
        fill = guide_legend(reverse = TRUE)
    ) +
    coord_flip() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90)))
save_plot(
    "figures/detection_map_site.png",
    site_level_map,
    bg = "white",
    base_asp = 0.8,
    base_height = 10
)

# Plot level detection map
(plot_level_map <- ggplot(
    obs_combined %>%
        group_by(plot_id, short_name) %>%
        summarise(
            botanical = any(botanical),
            eDNA = any(eDNA),
            eDNA_harm = any(eDNA_harm)
        ) %>%
        ungroup() %>%
        mutate(
            plot_id = factor(plot_id, levels = sort(unique(plot_id))),
            short_name = factor(
                short_name,
                levels = sort(unique(short_name), decreasing = T)
            )
        )
) +
    geom_tile(
        aes(
            x = short_name,
            y = plot_id,
            fill = paste(
                " botanical:",
                botanical,
                " eDNA_harm:",
                eDNA_harm,
                "eDNA:",
                eDNA
            ),
            colour = paste(
                " botanical:",
                botanical,
                " eDNA_harm:",
                eDNA_harm,
                "eDNA:",
                eDNA
            )
        ),
        width = 0.8,
        height = 0.6,
        linewidth = 1,
        linejoin = "mitre",
        lineend = "square"
    ) +
    scale_y_discrete(position = "right") +
    scale_colour_manual(
        values = c("#F28F16", "#A62121", "#223F59", "#223F59"),
        labels = rev(c(
            "botanical and eDNA",
            "botanical and eDNA harmonised",
            "botanical only",
            "eDNA only"
        ))
    ) +
    scale_fill_manual(
        values = c("#F28F16", "#A62121", "#A62121", "#223F59"),
        labels = rev(c(
            "botanical and eDNA",
            "botanical and eDNA harmonised",
            "botanical only",
            "eDNA only"
        ))
    ) +
    labs(x = "", y = "", colour = "Observed in", fill = "Observed in") +
    guides(
        colour = guide_legend(reverse = TRUE),
        fill = guide_legend(reverse = TRUE)
    ) +
    coord_flip() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90)))
save_plot(
    "figures/detection_map_plot.png",
    plot_level_map,
    bg = "white",
    base_asp = 1.5,
    base_height = 10
)
