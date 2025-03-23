# Correspondence Analysis Vascular Plants 
# Jakob J. Assmann jakob.assmann@uzh.ch 19 March 2025

# Dependencies
library(tidyverse)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggrepel)

# Load presence/ absence botanical observations and prepare for comaprison
botanical_obs <- read_csv("data/botanical_obs/presence_absence_cleaned.csv") %>%
    select(-observer, -scribe) %>%
    mutate(
        species = gsub("([A-Za-z]*) .*", "\\1", species),
        short_name = paste(genus, species)
    )

# Load presence / absence from eDNA
trnl_obs <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa_in_sample.csv") %>%
    mutate(short_name = paste(genus, species), sq = "trnl")
its2_obs <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa_in_sample.csv") %>%
    mutate(short_name = paste(genus, species), sq = "its2")
eDNA_obs <- bind_rows(trnl_obs, its2_obs) %>%
    group_by(sample_id, short_name) %>%
    group_split() %>%
    map(function(x) {
        # Arrange by count
        x <- arrange(x, count)
        # Get unique markers
        sqs <- sort(unique(x$sq)) %>% paste(collapse = " ")
        # Slice tor row with most counts
        x <- slice(x, 1)
        # Update markers
        x$sq <- sqs
        return(x)
    }) %>%
    bind_rows() %>%
    filter(!is.na(plot_id))

# Remove IDs that were not to the species level (for eDNA keep Salix sp)
botanical_obs <- filter(botanical_obs, species != "species")
eDNA_obs <- filter(eDNA_obs, species != "species" | genus == "Salix") 

# Add site id
botanical_obs <- mutate(botanical_obs, site_id = substr(plot_id, 1,2)) %>%
mutate(site_id = gsub("S", "Site ", site_id)) %>%
arrange(site_id, plot_id, short_name)
eDNA_obs <- mutate(eDNA_obs, site_id = substr(plot_id, 1,2)) %>%
mutate(site_id = gsub("S", "Site ", site_id)) %>%
arrange(site_id, plot_id, short_name)

# Generate a plot specific contingency table
bot_table_plot <- table(
    botanical_obs$plot_id,
    botanical_obs$short_name
)
bot_table_site <- table(
    botanical_obs$site_id,
    botanical_obs$short_name
)
edna_table_plot <- table(
    eDNA_obs$plot_id,
    eDNA_obs$short_name
)
edna_table_site <- table(
    eDNA_obs$site_id,
    eDNA_obs$short_name
)


# Correspondence analysis site level
(ca_bot_site <- CA(bot_table_site))
# -> significant association between site and species
# chi square =  138.1124 (p-value =  0.01003371 )

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
        data = ca_site_sp_spp_dims
        # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
        ,
        aes(
            x = `Dim 1`, y = `Dim 2`,
            label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
        ),
        max.overlaps = 20
    ) +    
    geom_label(
        data = ca_site_sp_site_dims,
        aes(
            x = `Dim 1`, y = `Dim 2`, label = site,
            color = site
        ),
        label.size = NA,
        fill = "#FFFFFF88",
        nudge_x = -0.25,
        hjust = 0.5,
        size = 28 /.pt
    ) +   
    geom_point(
        data = ca_site_sp_site_dims,
        aes(x = `Dim 1`, y = `Dim 2`,
        fill = site),
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
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.8, 
           label = "Botanical Surveys", size = 24 / .pt,
           fontface = "bold", colour = "black") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -1.8, label = "p-value = 0.01",
           size = 24/.pt) +
    labs(x = paste0("Dim 1 (", round(ca_bot_site$eig[1, 2], 2), "%)"), y = paste0("Dim 2 (", round(ca_bot_site$eig[2, 2], 2), "%)")) +
    theme_cowplot(24) +
        theme(
            legend.position = "none",
            panel.border = element_rect(color = "black")
        )


# Correspondence analysis site level - eDNA
(ca_edna_site <- CA(edna_table_site))
# -> significant association between site and species
# chi square =  138.1124 (p-value =   267.9041 (p-value =  3.330311e-23 ) )

# Visualise site level correspondence analysis
get_eigenvalue(ca_edna_site)
fviz_screeplot(ca_edna_site)
fviz_ca_biplot(ca_edna_site)

# Simplify objects
ca_site_sp_site_dims <- ca_edna_site$row$coord %>%
    as.data.frame() %>%
    mutate(., site = row.names(.)) %>%
    mutate(`Dim 1` = -1 * `Dim 1`,
    `Dim 2` = -1 * `Dim 2`)
ca_site_sp_spp_dims <- ca_edna_site$col$coord %>%
    as.data.frame() %>%
    mutate(., species = row.names(.)) %>%
    mutate(`Dim 1` = -1 * `Dim 1`,
    `Dim 2` = -1 * `Dim 2`)

ca_plot_site_edna <- ggplot() +
    geom_point(
        data = ca_site_sp_spp_dims,
        aes(x = `Dim 1`, y = `Dim 2`)
    ) +
    geom_text_repel(
        data = ca_site_sp_spp_dims
        # filter(abs(`Dim 1`) > 0.5 | abs(`Dim 2`) > 0.5)
        ,
        aes(
            x = `Dim 1`, y = `Dim 2`,
            label = gsub("([A-Z])[a-z]* ([a-z]*).*", "\\1\\. \\2", species)
        ),
        max.overlaps = 20
    ) +   
    geom_label(
        data = ca_site_sp_site_dims,
        aes(
            x = `Dim 1`, y = `Dim 2`, label = site,
            color = site
        ),
        label.size = NA,
        fill = "#FFFFFF88",
        nudge_x = -0.25,
        hjust = 0.5,
        size = 28 /.pt
    ) + 
    geom_point(
        data = ca_site_sp_site_dims,
        aes(x = `Dim 1`, y = `Dim 2`,
        fill = site),
        colour = "white",
        shape = 21,
        size = 6,
    ) + 
    scale_x_continuous(
        limits = c(-1.5, 1.5),
        breaks = seq(-1.5, 1.5, 0.5)
    ) +
    # scale_y_reverse() +
  scale_colour_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
  scale_fill_manual(values = c("#F2D64B", "#5F1773", "#F2507B")) +
    geom_vline(xintercept = 0, alpha = 0.5) +
        geom_hline(yintercept = 0, alpha = 0.5) +
        scale_y_continuous(
            limits = c(-1.6, 1.6),
            breaks = seq(-1.5, 1.5, 0.5)
        ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.8, 
           label = "eDNA", size = 24 / .pt,
           fontface = "bold", colour = "black") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -1.8, label = "p-value < 0.001",
           size = 24/.pt) +
    labs(x = paste0("Dim 1 (", round(ca_edna_site$eig[1, 2], 2), "%)"), y = paste0("Dim 2 (", round(ca_edna_site$eig[2, 2], 2), "%)")) +
    theme_cowplot(24) +
        theme(
            legend.position = "none",
            panel.border = element_rect(color = "black")
        )
set.seed(24)
save_plot("figures/fig_2_ca_site.png", 
    plot_grid(ca_plot_site_bot, ca_plot_site_edna),
    bg = "white",
    base_height = 7,
    base_asp = 1.6*2
)



