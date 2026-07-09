# Detection "heatmap"
# Jakob J. Assmann jakob.assmann@gu.se 8 August 2027

# Library
tidyverse()
library(cowplot)
library(ggplot2)


ggplot(
    obs_combined %>%
        # mutate(plot_id = substr(plot_id, 1, 2)) %>%
        # group_by(plot_id, short_name) %>%
        # summarise(
        #     botanical = any(botanical),
        #     eDNA = any(eDNA),
        #     eDNA_genus = any(eDNA_genus)
        # ) %>%
        # ungroup() %>%
        mutate(
            plot_id = factor(plot_id, levels = sort(unique(plot_id))),
            short_name = factor(
                short_name,
                levels = sort(unique(short_name), decreasing = T)
            )
        )
) +
    geom_point(
        aes(
            x = short_name,
            y = plot_id,
            fill = paste0(botanical, eDNA),
            colour = eDNA_genus & botanical & !eDNA
        ),
        size = 5,
        ,
        stroke = 0,
        shape = 22,
        colour = "white"
    ) +
    geom_point(aes(
        x = short_name,
        y = plot_id,
        colour = eDNA_genus & botanical & !eDNA
    )) +
    scale_y_discrete(position = "right") +
    scale_colour_manual(values = c("#FFFFFF00", "white")) +
    scale_fill_manual(
        values = c("#F28F16", "#A62121", "#223F59"),
        labels = c("eDNA only", "botanical only", "botanical and eDNA")
    ) +
    labs(colour = "Observed in") +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
