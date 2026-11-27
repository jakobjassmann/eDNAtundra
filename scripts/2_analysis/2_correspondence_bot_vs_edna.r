# Correspondence Analysis Vascular Plants
# Jakob J. Assmann jakob.assmann@uzh.ch 19 March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)

# Load and prepare data
obs_combined <- read_csv("data/analysis/detection_species_genus.csv") %>%
    # Add site id
    mutate(site_id = substr(plot_id, 1, 2)) %>%
    mutate(
        site_id = case_when(
            site_id == "S1" ~ "Site 1",
            site_id == "S2" ~ "Site 2",
            site_id == "S3" ~ "Site 3"
        )
    ) %>%
    relocate(site_id)

# Generate contingency tables
bot_table_site <- obs_combined %>%
    select(site_id, short_name, botanical) %>%
    filter(botanical)
bot_table_site <- table(
    bot_table_site$site_id,
    bot_table_site$short_name
)
edna_table_site <- obs_combined %>%
    select(site_id, short_name, eDNA) %>%
    filter(eDNA)
edna_table_site <- table(
    edna_table_site$site_id,
    edna_table_site$short_name
)

# Correspondence analysis site level
(ca_bot_site <- CA(bot_table_site))
# -> significant association between site and species
# chi square =  133.2929 (p-value =  0.02039588)

# Visualise site level correspondence analysis
get_eigenvalue(ca_bot_site)
fviz_screeplot(ca_bot_site)
fviz_ca_biplot(ca_bot_site)

# Simplify objects
ca_site_sp_site_dims <- ca_bot_site$row$coord %>%
    as.data.frame() %>%
    mutate(., site = row.names(.))
ca_site_sp_spp_dims <- ca_bot_site$col$coord %>%
    as.data.frame() %>%
    mutate(., species = row.names(.))

ca_plot_site_bot <- ggplot() +
    geom_point(
        data = ca_site_sp_spp_dims,
        aes(x = `Dim 1`, y = `Dim 2`)
    ) +
    geom_text_repel(
        data = ca_site_sp_spp_dims,
        # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
        aes(
            x = `Dim 1`,
            y = `Dim 2`,
            label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
        ),
        max.overlaps = 20
    ) +
    geom_label(
        data = ca_site_sp_site_dims,
        aes(
            x = `Dim 1`,
            y = `Dim 2`,
            label = site,
            color = site
        ),
        label.size = NA,
        fill = "#FFFFFF88",
        nudge_x = -0.25,
        hjust = 0.5,
        size = 28 / .pt
    ) +
    geom_point(
        data = ca_site_sp_site_dims,
        aes(x = `Dim 1`, y = `Dim 2`, fill = site),
        colour = "white",
        shape = 21,
        size = 6,
    ) +
    scale_x_continuous(
        limits = c(-1.5, 1.5),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    scale_y_continuous(
        limits = c(-1.6, 1.6),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    annotate(
        "text",
        x = Inf,
        y = Inf,
        hjust = 1.2,
        vjust = 1.8,
        label = "Botanical Surveys",
        size = 24 / .pt,
        fontface = "bold",
        colour = "black"
    ) +
    annotate(
        "text",
        x = Inf,
        y = -Inf,
        hjust = 1.2,
        vjust = -1,
        label = "p-value < 0.05",
        size = 16 / .pt
    ) +
    labs(
        x = paste0("Dim 1 (", round(ca_bot_site$eig[1, 2], 2), "%)"),
        y = paste0("Dim 2 (", round(ca_bot_site$eig[2, 2], 2), "%)")
    ) +
    theme_cowplot(24) +
    theme(
        legend.position = "none",
        panel.border = element_rect(color = "black")
    )


# Correspondence analysis site level - eDNA
(ca_edna_site <- CA(edna_table_site))
# -> significant association between site and species
# chi square =  114.497 (p-value =  0.03038606 )

# Visualise site level correspondence analysis
get_eigenvalue(ca_edna_site)
fviz_screeplot(ca_edna_site)
fviz_ca_biplot(ca_edna_site)

# Simplify objects
ca_site_sp_site_dims <- ca_edna_site$row$coord %>%
    as.data.frame() %>%
    mutate(., site = row.names(.)) %>%
    mutate(`Dim 1` = -1 * `Dim 1`, `Dim 2` = -1 * `Dim 2`)
ca_site_sp_spp_dims <- ca_edna_site$col$coord %>%
    as.data.frame() %>%
    mutate(., species = row.names(.)) %>%
    mutate(`Dim 1` = -1 * `Dim 1`, `Dim 2` = -1 * `Dim 2`)

ca_plot_site_edna <- ggplot(
    mapping = aes(x = 1 * `Dim 1`, y = -1 * `Dim 2`)
) +
    geom_point(
        data = ca_site_sp_spp_dims
    ) +
    geom_text_repel(
        data = ca_site_sp_spp_dims,
        # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
        aes(
            label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
        ),
        max.overlaps = 20
    ) +
    geom_label(
        data = ca_site_sp_site_dims,
        aes(
            label = site,
            color = site
        ),
        label.size = NA,
        fill = "#FFFFFF88",
        nudge_x = -0.25,
        hjust = 0.5,
        size = 28 / .pt
    ) +
    geom_point(
        data = ca_site_sp_site_dims,
        aes(fill = site),
        colour = "white",
        shape = 21,
        size = 6,
    ) +
    scale_x_continuous(
        limits = c(-1.5, 1.5),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    scale_y_continuous(
        limits = c(-1.6, 1.6),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    annotate(
        "text",
        x = Inf,
        y = Inf,
        hjust = 1.2,
        vjust = 1.8,
        label = "eDNA",
        size = 24 / .pt,
        fontface = "bold",
        colour = "black"
    ) +
    annotate(
        "text",
        x = Inf,
        y = -Inf,
        hjust = 1.2,
        vjust = -1,
        label = "p-value < 0.05",
        size = 15 / .pt
    ) +
    labs(
        x = paste0("Dim 1 (", round(ca_edna_site$eig[1, 2], 2), "%)"),
        y = paste0("Dim 2 (", round(ca_edna_site$eig[2, 2], 2), "%)")
    ) +
    theme_cowplot(24) +
    theme(
        legend.position = "none",
        panel.border = element_rect(color = "black")
    )
save_plot(
    "figures/figure_S3_ca_site.png",
    plot_grid(
        ca_plot_site_bot,
        ca_plot_site_edna,
        labels = c("a)", "b)"),
        label_size = 24
    ),
    bg = "white",
    base_height = 7,
    base_asp = 1.6 * 2
)

#### End of file
################################################################################
