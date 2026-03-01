############################################################
# Microbiome beta-diversity + PERMANOVA across taxonomic levels
#
# Workflow:
# 1. Read MetaPhlAn-style abundance table into microeco
# 2. Collapse abundances to Order, Family, and Genus
# 3. Compute Bray–Curtis dissimilarities
# 4. Run PERMANOVA (adonis2) for:
#       - Geographic_region
#       - year
#       - month
# 5. Extract R² and significance
# 6. Plot R² across taxonomic levels
############################################################


############################
# 1. Load required packages
############################

library(vegan)      
library(tidyr)      
library(dplyr)      
library(ggplot2)    


############################
# 2. Define input file paths
############################

# MetaPhlAn abundance table 
abund_path  <- "data/MASTER_abundance_table_no_unclassified.txt"

# Sample metadata table
sample_path <- "data/sample_info_all.csv"

# Matching table required by mpa2meco
match_table <- read.table("data/samples_all.txt")
match_table$V2 <- match_table$V1   

############################
# 3. Read data into microeco
############################

# Convert MetaPhlAn output into microeco object
# rel = TRUE → use relative abundances
mt <- mpa2meco(
    abund_path,
    sample_table = sample_path,
    match_table  = match_table,
    rel = TRUE
)

# Extract metadata
meta <- mt$sample_table

# Ensure temporal variables are treated as categorical factors
meta$year  <- as.factor(meta$year)
meta$month <- as.factor(meta$month)


###########################################################
# 4. Collapse abundances to different taxonomic resolutions
###########################################################

# Aggregate relative abundances to desired taxonomic rank
mt_order  <- trans_abund$new(dataset = mt, taxrank = "Order")
mt_family <- trans_abund$new(dataset = mt, taxrank = "Family")
mt_genus  <- trans_abund$new(dataset = mt, taxrank = "Genus")


############################################################
# 5. Convert long-format microeco output to wide matrices
#    (Samples as rows, taxa as columns)
############################################################

# ----- ORDER LEVEL -----
order_wide <- mt_order$data_abund %>%
    dplyr::select(Taxonomy, Sample, Abundance) %>%
    pivot_wider(
        names_from  = Taxonomy,
        values_from = Abundance,
        values_fill = 0
    )

# ----- FAMILY LEVEL -----
family_wide <- mt_family$data_abund %>%
    dplyr::select(Taxonomy, Sample, Abundance) %>%
    pivot_wider(
        names_from  = Taxonomy,
        values_from = Abundance,
        values_fill = 0
    )

# ----- GENUS LEVEL -----
genus_wide <- mt_genus$data_abund %>%
    dplyr::select(Taxonomy, Sample, Abundance) %>%
    pivot_wider(
        names_from  = Taxonomy,
        values_from = Abundance,
        values_fill = 0
    )


############################################################
# 6. Prepare matrices for distance calculation
#    - Set Sample IDs as rownames
#    - Remove Sample column
############################################################

order_mat <- as.data.frame(order_wide)
rownames(order_mat) <- order_mat$Sample
order_mat$Sample <- NULL

family_mat <- as.data.frame(family_wide)
rownames(family_mat) <- family_mat$Sample
family_mat$Sample <- NULL

genus_mat <- as.data.frame(genus_wide)
rownames(genus_mat) <- genus_mat$Sample
genus_mat$Sample <- NULL


###############################################
# 7. Compute Bray–Curtis dissimilarity matrices
###############################################

bray_order  <- vegdist(order_mat,  method = "bray")
bray_family <- vegdist(family_mat, method = "bray")
bray_genus  <- vegdist(genus_mat,  method = "bray")


############################################################
# 8. Run PERMANOVA (adonis2)
#
# by = "margin" → tests each factor marginally,
#                  controlling for other predictors
############################################################

perm_order <- adonis2(
    bray_order ~ Geographic_region + year + month,
    data = meta,
    permutations = 999,
    by = "margin"
)

perm_family <- adonis2(
    bray_family ~ Geographic_region + year + month,
    data = meta,
    permutations = 999,
    by = "margin"
)

perm_genus <- adonis2(
    bray_genus ~ Geographic_region + year + month,
    data = meta,
    permutations = 999,
    by = "margin"
)


############################################################
# 9. Extract R² and p-values from PERMANOVA output
############################################################

extract_r2 <- function(perm, level) {
    df <- as.data.frame(perm)
    df$Factor <- rownames(df)
    
    # Keep only predictors of interest
    df <- df[df$Factor %in% c("Geographic_region", "year", "month"), ]
    
    df$Level <- level
    
    # Return clean summary table
    df[, c("Level", "Factor", "R2", "Pr(>F)")]
}

r2_order  <- extract_r2(perm_order,  "Order")
r2_family <- extract_r2(perm_family, "Family")
r2_genus  <- extract_r2(perm_genus,  "Genus")

# Combine into single summary table
r2_all <- bind_rows(r2_order, r2_family, r2_genus)

# Export this table
write.csv(r2_all, "microeco_output/permanova_braycurtis_taxonomiclevels.csv")


############################################################
# 10. Add significance codes and set factor ordering
############################################################

r2_all$Signif <- cut(
    r2_all$`Pr(>F)`,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
)

r2_all$Level <- factor(
    r2_all$Level,
    levels = c("Order", "Family", "Genus")
)


############################################################
# 11. Plot PERMANOVA R² across taxonomic levels
############################################################

permanova_plot <- ggplot(
    r2_all,
    aes(x = Level, y = R2, color = Factor, group = Factor)
) +
    geom_point(size = 4) +
    geom_line(size = 1) +
    geom_text(
        aes(label = Signif),
        vjust = -0.7,
        size = 6,
        show.legend = FALSE
    ) +
    ylab("PERMANOVA R²") +
    xlab("") +
    scale_y_continuous(
        limits = c(0, 0.135),
        expand = expansion(mult = c(0, 0.05))
    ) +
    scale_color_manual(
        values = c(
            "Geographic_region" = "#D55E00",
            "year"              = "#440154FF",
            "month"             = "#009E73"
        )
    ) +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.92),
        legend.box.background = element_rect(color = "black")
    )

permanova_plot

