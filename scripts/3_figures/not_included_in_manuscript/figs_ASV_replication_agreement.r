# Replicate agreement of ASV-based taxon detection
# Jakob J. Assmann jakob.assmann@gu.se 5 August 2026
# vibe coded with GitHub copilot (manually validated)

# Dependencies
library(tidyverse)
library(cowplot)

compute_agreement <- function(taxa_df, marker_name) {
    taxa_df <- taxa_df %>%
        mutate(site_id = substr(plot_id, 1, 2)) %>%
        mutate(site_id = gsub("S", "Site ", site_id)) %>%
        filter(!is.na(species)) %>%
        mutate(
            short_name = paste(genus, species)
        ) %>%
        group_by(site_id, plot_id, subsample_id, short_name) %>%
        summarise(count = sum(count), .groups = "drop")

    # identify pooled libraries (keep both pooled and non-pooled in the frame)
    taxa_df <- taxa_df %>%
        mutate(
            is_pooled = grepl("P", subsample_id, ignore.case = TRUE)
        )

    # separate subsamples (non-pooled) and pooled libraries
    subsamples_df <- taxa_df %>% filter(!is_pooled)
    pooled_df <- taxa_df %>%
        filter(is_pooled) %>%
        select(plot_id, short_name) %>%
        mutate(pooled_pres = TRUE)

    # for each plot x taxon count how many replicates detected it (0..3)
    agreement_df <- subsamples_df %>%
        group_by(site_id, plot_id, short_name) %>%
        summarise(n_detect = n(), .groups = "drop")

    # join pooled presence to agreement frame
    agreement_df <- agreement_df %>%
        left_join(pooled_df, by = c("plot_id", "short_name")) %>%
        mutate(pooled_pres = if_else(is.na(pooled_pres), 0L, pooled_pres))

    # per-plot summary: fraction of taxa detected in all three replicates
    per_plot <- agreement_df %>%
        group_by(site_id, plot_id) %>%
        summarise(
            n_taxa = n(),
            n_all3 = sum(n_detect == 3),
            prop_all3 = n_all3 / n_taxa,
            n_atleast2 = sum(n_detect >= 2),
            prop_atleast2 = n_atleast2 / n_taxa,
            # pooled performance: fraction of taxa (seen in subsamples) detected in pooled library
            n_pooled = sum(pooled_pres == 1),
            prop_pooled = n_pooled / n_taxa,
            .groups = "drop"
        )

    # overall statistic (mean across plots)
    overall <- per_plot %>%
        summarise(
            mean_prop_all3 = mean(prop_all3, na.rm = TRUE),
            median_prop_all3 = median(prop_all3, na.rm = TRUE),
            mean_prop_atleast2 = mean(prop_atleast2, na.rm = TRUE),
            mean_prop_pooled = mean(prop_pooled, na.rm = TRUE)
        ) %>%
        mutate(marker = marker_name)

    # prepare plotting data: distribution of n_detect per site
    plot_df <- agreement_df %>%
        group_by(site_id, plot_id, n_detect) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(site_id, plot_id) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup() %>%
        mutate(n_detect = factor(n_detect, levels = c("1", "2", "3")))

    # pooled per-plot/frame for plotting pooled sensitivity per site
    pooled_plot_df <- per_plot %>%
        group_by(site_id) %>%
        summarise(
            mean_prop_pooled = mean(prop_pooled, na.rm = TRUE),
            .groups = "drop"
        )

    return(list(
        overall = overall,
        per_plot = per_plot,
        plot_df = plot_df,
        pooled_plot_df = pooled_plot_df
    ))
}

plot_agreement_by_site <- function(plot_df, marker_name) {
    p <- ggplot(
        plot_df,
        aes(x = site_id, y = prop, fill = as.factor(n_detect))
    ) +
        geom_col(position = "fill", width = 0.6) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_fill_manual(
            values = rev(c("#00485C", "#EB7E1C", "#D24141")),
            name = "# replicates\nwith detection"
        ) +
        labs(
            x = "",
            y = "Proportion of taxa",
            title = paste0(marker_name, " — Agreement across subsamples")
        ) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5))
    p
}

# Run for ITS2
its2 <- read_csv("data/sequences/its2/edna/cleaned/its2_taxa.csv")
its2_res <- compute_agreement(its2, "ITS2")
its2_plot <- plot_agreement_by_site(its2_res$plot_df, "ITS2")

save_plot(
    "figures/ASV_replication_agreement_ITS2.png",
    its2_plot,
    base_width = 9,
    base_height = 5,
    bg = "white"
)

# pooled sensitivity plot (ITS2)
its2_pooled_plot <- ggplot(
    its2_res$pooled_plot_df,
    aes(x = site_id, y = mean_prop_pooled)
) +
    geom_col(fill = "#2b8cbe", width = 0.6) +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1)
    ) +
    labs(
        x = "",
        y = "Mean proportion of taxa detected in pooled library",
        title = "ITS2 — Pooled sample sensitivity (per site)"
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))

save_plot(
    "figures/ASV_replication_pooled_ITS2.png",
    its2_pooled_plot,
    base_width = 6,
    base_height = 4,
    bg = "white"
)

# Run for trnL
trnl <- read_csv("data/sequences/trnl/edna/cleaned/trnl_taxa.csv")
trnl_res <- compute_agreement(trnl, "trnL")
trnl_plot <- plot_agreement_by_site(trnl_res$plot_df, "trnL")

save_plot(
    "figures/ASV_replication_agreement_trnl.png",
    trnl_plot,
    base_width = 9,
    base_height = 5,
    bg = "white"
)

# pooled sensitivity plot (trnL)
trnl_pooled_plot <- ggplot(
    trnl_res$pooled_plot_df,
    aes(x = site_id, y = mean_prop_pooled)
) +
    geom_col(fill = "#66c2a5", width = 0.6) +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1)
    ) +
    labs(
        x = "",
        y = "Mean proportion of taxa detected in pooled library",
        title = "trnL — Pooled sample sensitivity (per site)"
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))

save_plot(
    "figures/ASV_replication_pooled_trnl.png",
    trnl_pooled_plot,
    base_width = 6,
    base_height = 4,
    bg = "white"
)

# Combined summary table
combined_overall <- bind_rows(its2_res$overall, trnl_res$overall)
print("Overall replication agreement (per-marker):")
print(combined_overall)
# # A tibble: 2 × 5
#   mean_prop_all3 median_prop_all3 mean_prop_atleast2 mean_prop_pooled marker
#            <dbl>            <dbl>              <dbl>            <dbl> <chr>
# 1          0.439              0.4              0.661            0.712 ITS2
# 2          0.487              0.5              0.727            0.696 trnL

# A single scalar to quote: mean proportion of taxa detected in all three replicates (averaged across markers)
quote_stat <- combined_overall %>%
    summarise(
        mean_prop_all3_across_markers = mean(mean_prop_all3, na.rm = TRUE)
    ) %>%
    pull(mean_prop_all3_across_markers)

cat(sprintf(
    "Mean proportion of taxa detected in all three subsamples (avg across markers): %.1f%%\n",
    quote_stat * 100
))
# Mean proportion of taxa detected in all three subsamples (avg across markers): 46.3%

# Also report pooled sensitivity (fraction of taxa seen in subsamples that pooled captured)
pooled_stat <- combined_overall %>%
    summarise(
        mean_prop_pooled_across_markers = mean(mean_prop_pooled, na.rm = TRUE)
    ) %>%
    pull(mean_prop_pooled_across_markers)

cat(sprintf(
    "Mean proportion of taxa detected in pooled library (avg across markers): %.1f%%\n",
    pooled_stat * 100
))
# Mean proportion of taxa detected in pooled library (avg across markers): 70.4%

#### End of file
################################################################################
