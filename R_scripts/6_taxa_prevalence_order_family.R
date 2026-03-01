###############################
# Taxa prevalence and abundance
# Levels: Order and Family
###############################

library(microeco)
library(file2meco)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)

# Master dataset paths
abund_path <- "data/MASTER_abundance_table_no_unclassified.txt"
sample_path <- "data/sample_info_all.csv"

# Output directory
out_dir <- "microeco_output/taxa_prevalence_order_family"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Match table required by mpa2meco
match_table <- read.table("data/samples_all.txt")
match_table$V2 <- match_table$V1

# Read data as microeco object (relative abundance)
mt <- mpa2meco(
  abund_path,
  sample_table = sample_path,
  match_table = match_table,
  rel = TRUE
)

# Same extraction style as RC_3.R
mt_order <- trans_abund$new(dataset = mt, taxrank = "Order")
mt_family <- trans_abund$new(dataset = mt, taxrank = "Family")

n_samples <- nrow(mt$sample_table)
prevalence_cutoff <- 75

compute_taxa_stats <- function(trans_obj, rank_name) {
  long_df <- trans_obj$data_abund %>%
    dplyr::select(Taxonomy, Sample, Abundance) %>%
    mutate(Taxonomy = as.character(Taxonomy),
           Sample   = as.character(Sample)) %>%
    filter(!grepl("Bacteria_unclassified", Taxonomy, fixed = TRUE))

  stats <- long_df %>%
    group_by(Taxonomy, Sample) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Taxonomy) %>%
    summarise(samples_present = sum(Abundance > 0, na.rm = TRUE),
              prevalence_pct  = 100 * samples_present / n_samples,
              mean_abundance  = mean(Abundance, na.rm = TRUE),
              mean_abundance_present = ifelse(samples_present > 0, mean(Abundance[Abundance > 0], na.rm = TRUE), 0),
              .groups = "drop") %>%
    mutate(rank = rank_name,
           common_75 = prevalence_pct > prevalence_cutoff) %>%
    rename(taxon = Taxonomy) %>%
    arrange(desc(samples_present), desc(mean_abundance))

  write_csv(stats, file.path(out_dir, paste0(rank_name, "_taxa_stats.csv")))
  write_csv(
    stats %>% filter(common_75),
    file.path(out_dir, paste0(rank_name, "_taxa_prevalence_gt75.csv"))
  )

  top_n <- min(30, nrow(stats))
  top_stats <- stats %>% slice_head(n = top_n)

  p_prev <- ggplot(
    top_stats %>% mutate(taxon = fct_reorder(taxon, samples_present)),
    aes(x = taxon, y = samples_present, fill = common_75)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#C7AEFF", "FALSE" = "#5BA3D6"),
      labels = c("TRUE" = ">75% samples", "FALSE" = "<=75% samples"),
      name = "Prevalence class"
    ) +
    labs(title = paste0(rank_name, ": taxa prevalence across samples"),
         subtitle = paste0("Top ", top_n, " taxa by number of samples present"),
         x = rank_name,
         y = "Number of samples present") +
    theme_bw(base_size = 12)

  p_mean <- ggplot(
    top_stats %>% mutate(taxon = fct_reorder(taxon, mean_abundance)),
    aes(x = taxon, y = mean_abundance, fill = common_75)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#C7AEFF", "FALSE" = "#4FCFAF"),
      labels = c("TRUE" = ">75% samples", "FALSE" = "<=75% samples"),
      name = "Prevalence class"
    ) +
    labs(
      title = paste0(rank_name, ": mean abundance across samples"),
      subtitle = paste0("Top ", top_n, " taxa by prevalence"),
      x = rank_name,
      y = "Mean abundance"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    filename = file.path(out_dir, paste0(rank_name, "_prevalence_top", top_n, ".png")),
    plot = p_prev,
    width = 10,
    height = 8,
    dpi = 300
  )

  ggsave(
    filename = file.path(out_dir, paste0(rank_name, "_mean_abundance_top", top_n, ".png")),
    plot = p_mean,
    width = 10,
    height = 8,
    dpi = 300
  )

  return(top_stats)

  message(paste0("Finished: ", rank_name))
}

order_top <- compute_taxa_stats(mt_order, rank_name = "Order")
family_top <- compute_taxa_stats(mt_family, rank_name = "Family")


# Four-panel figure (Order/Family x Prevalence/Mean abundance)
plot_df <- bind_rows(order_top, family_top) %>%
  select(rank, taxon, samples_present, mean_abundance, common_75) %>%
  pivot_longer(
    cols = c(samples_present, mean_abundance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      samples_present = "Prevalence (n samples)",
      mean_abundance = "Mean abundance"
    ),
    panel = paste(rank, metric, sep = " | "),
    color_group = ifelse(
      common_75,
      ">75% samples",
      "<=75% samples"
    ),
    taxon_panel = paste0(taxon, "___", panel)
  )

p_four <- ggplot(
  plot_df %>% mutate(taxon_panel = fct_reorder(taxon_panel, value)),
  aes(x = taxon_panel, y = value, fill = color_group)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_fill_manual(
    values = c(">75% samples" = "#C7AEFF", "<=75% samples" = "#5BA3D6")
  ) +
  labs(
    #title = "Order and Family: prevalence and mean abundance",
    #subtitle = "Taxa present in >75% of samples are highlighted",
    x = NULL,
    y = NULL,
    fill = "Prevalence class"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(out_dir, "Order_Family_prevalence_mean_4panel_top30.png"),
  plot = p_four,
  width = 7.2,
  height = 8.8,
  dpi = 300
)

