###############################################
# Core taxa vs differential taxa mean abundance
# Levels: Order and Family
###############################################

# 1) Get core taxa (>75% prevalence) from script 6 outputs.
# 2) Get significant differential taxa from script 2 LEfSe results.
# 3) Compute mean abundance at Order and Family levels with microeco.
# 4) Label each taxon as Core/Differential/Both/Other and plot.

library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(stringr)
library(microeco)

# Output directory
out_dir <- "microeco_output/core_vs_diff_mean_abundance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source required scripts in an isolated environment
e <- new.env(parent = globalenv())

# 1) Core taxa prevalence script (>75% taxa)
source("R_scripts/6_taxa_prevalence_order_family.R", local = e)

# 2) Build mt object
source("R_scripts/1_microeco_alpha_div.R", local = e)

# 3) Differential abundance script
try(source("R_scripts/2_microeco_composition.R", local = e), silent = TRUE)

if (!exists("t0", envir = e) || !exists("t1", envir = e)) {
  stop("Could not find differential abundance objects t0/t1 from R_scripts/2_microeco_composition.R")
}

if (!exists("mt", envir = e)) {
  stop("Could not find microeco dataset object mt from sourced scripts.")
}

prevalence_cutoff <- 75
alpha_cutoff <- 0.05

# -----------------------------
# A) Core taxa from script 6
# -----------------------------
core_order <- read_csv(
  "microeco_output/taxa_prevalence_order_family/Order_taxa_stats.csv",
  show_col_types = FALSE
) %>%
  filter(prevalence_pct > prevalence_cutoff) %>%
  transmute(rank = "Order", taxon, is_core = TRUE)

core_family <- read_csv(
  "microeco_output/taxa_prevalence_order_family/Family_taxa_stats.csv",
  show_col_types = FALSE
) %>%
  filter(prevalence_pct > prevalence_cutoff) %>%
  transmute(rank = "Family", taxon, is_core = TRUE)

core_taxa <- bind_rows(core_order, core_family)

# -----------------------------------------
# B) Differential taxa from script 2 LEfSe
# -----------------------------------------
# Keep significant taxa and extract terminal taxon names:

ord_diff_raw <- e$t0$res_diff
fam_diff_raw <- e$t1$res_diff

if ("P.adj" %in% colnames(ord_diff_raw)) {
  ord_diff_raw <- ord_diff_raw %>% mutate(p_sig = P.adj)
} else {
  ord_diff_raw <- ord_diff_raw %>% mutate(p_sig = P.unadj)
}

if ("P.adj" %in% colnames(fam_diff_raw)) {
  fam_diff_raw <- fam_diff_raw %>% mutate(p_sig = P.adj)
} else {
  fam_diff_raw <- fam_diff_raw %>% mutate(p_sig = P.unadj)
}

ord_diff <- ord_diff_raw %>%
  mutate(
    taxon = str_split(Taxa, "\\|") %>%
      sapply(function(x) tail(x, 1)) %>%
      str_replace("^[a-z]__", "")
  ) %>%
  filter(p_sig < alpha_cutoff) %>%
  distinct(taxon) %>%
  transmute(rank = "Order", taxon, is_diff = TRUE)

fam_diff <- fam_diff_raw %>%
  mutate(
    taxon = str_split(Taxa, "\\|") %>%
      sapply(function(x) tail(x, 1)) %>%
      str_replace("^[a-z]__", "")
  ) %>%
  filter(p_sig < alpha_cutoff) %>%
  distinct(taxon) %>%
  transmute(rank = "Family", taxon, is_diff = TRUE)

diff_taxa <- bind_rows(ord_diff, fam_diff)

# ------------------------------------------
# C) Mean abundance from microeco trans_abund
# ------------------------------------------
ord_abund <- trans_abund$new(dataset = e$mt, taxrank = "Order")$data_abund %>%
  dplyr::select(Taxonomy, Sample, Abundance) %>%
  group_by(Taxonomy, Sample) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Taxonomy) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  transmute(rank = "Order", taxon = as.character(Taxonomy), mean_abundance)

fam_abund <- trans_abund$new(dataset = e$mt, taxrank = "Family")$data_abund %>%
  dplyr::select(Taxonomy, Sample, Abundance) %>%
  group_by(Taxonomy, Sample) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Taxonomy) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  transmute(rank = "Family", taxon = as.character(Taxonomy), mean_abundance)

mean_abund <- bind_rows(ord_abund, fam_abund) %>%
  filter(!grepl("Bacteria_unclassified", taxon, fixed = TRUE))

# -----------------------------------------
# D) Merge labels and define final category
# -----------------------------------------
status_tbl <- mean_abund %>%
  left_join(core_taxa, by = c("rank", "taxon")) %>%
  left_join(diff_taxa, by = c("rank", "taxon")) %>%
  mutate(
    is_core = ifelse(is.na(is_core), FALSE, is_core),
    is_diff = ifelse(is.na(is_diff), FALSE, is_diff),
    status = case_when(
      is_core & is_diff ~ "Core + Differential",
      is_core & !is_diff ~ "Core only",
      !is_core & is_diff ~ "Differential only",
      TRUE ~ "None"
    )
  )

write_csv(status_tbl, file.path(out_dir, "Order_Family_core_diff_mean_abundance.csv"))
write_csv(status_tbl %>% filter(is_core), file.path(out_dir, "core_taxa_order_family.csv"))
write_csv(status_tbl %>% filter(is_diff), file.path(out_dir, "differential_taxa_order_family.csv"))

# -----------------------------------------
# E) Plot top taxa by mean abundance per rank
# -----------------------------------------
top_n <- 50
plot_df <- status_tbl %>%
  group_by(rank) %>%
  arrange(desc(mean_abundance), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(rank_taxon = paste0(rank, "___", taxon))

p <- ggplot(
  plot_df %>% mutate(rank_taxon = fct_reorder(rank_taxon, mean_abundance)),
  aes(x = rank_taxon, y = mean_abundance, fill = status)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ rank, scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(
    values = c(
      "Core + Differential" = "#440154FF",
      "Core only" = "#21908CFF",
      "Differential only" = "#6DCD59FF",
      "None" = "#BDBDBD"
    )
  ) +
  labs(
    #title = "Mean abundance by taxon: core and differential status",
    #subtitle = "Order and Family levels; core taxa are prevalence >75%; differential taxa are LEfSe significant (P.adj < 0.05)",
    x = NULL,
    y = "Mean abundance across all samples",
    fill = "Taxon status"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

p

ggsave(
  filename = file.path(out_dir, "Order_Family_mean_abundance_core_vs_diff_top30.png"),
  plot = p,
  width = 7.2,
  height = 6.8,
  dpi = 300
)


