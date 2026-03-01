############################################################
# Assess potential effect of time (month + year)
#
# This script evaluates temporal effects on:
#   1) Alpha diversity (Shannon, Simpson, InvSimpson, Pielou)
#
# Predictors tested:
#   - Geographic_region (primary biological variable)
#   - year
#   - month
#
############################################################


############################
# 1. Load libraries
############################

library(vegan)       # beta diversity + PERMANOVA
library(car)         # Type II ANOVA
library(dplyr)       # data manipulation
library(effectsize)  # partial eta squared


############################################################
# PART I — Alpha diversity ~ region + time
############################################################

############################
# 2. Load data via microeco
############################

# Input files
abund_path  <- "data/MASTER_abundance_table_no_unclassified.txt"
sample_path <- "data/sample_info_all.csv"

# Required mapping for mpa2meco
match_table <- read.table("data/samples_all.txt")
match_table$V2 <- match_table$V1

# Convert MetaPhlAn output to microeco object (relative abundances)
mt <- mpa2meco(
    abund_path,
    sample_table = sample_path,
    match_table  = match_table,
    rel = TRUE
)


############################
# 3. Calculate alpha diversity
############################

# Compute alpha diversity metrics
# (group argument does not restrict model; it structures internal output)
t_time <- trans_alpha$new(dataset = mt, group = "Geographic_region")

# Extract long-format alpha table and enforce factor encoding
alpha_long <- t_time$data_alpha %>%
    mutate(
        Geographic_region = factor(Geographic_region),
        year  = factor(year),
        month = factor(month)
    )


############################################################
# 4. Linear model per diversity metric
#
# Model:
#   Value ~ Geographic_region + year + month
#
# - Type II ANOVA (car::Anova)
# - Partial eta² for effect size
############################################################

run_alpha_long <- function(measure_name) {
    
    # Subset for one diversity measure
    df <- alpha_long %>% filter(Measure == measure_name)
    
    # Linear model
    m <- lm(Value ~ Geographic_region + year + month, data = df)
    
    # Type II ANOVA (appropriate when no interaction terms)
    a <- car::Anova(m, type = 2)
    
    # Partial eta squared (effect size)
    e <- effectsize::eta_squared(m, partial = TRUE)
    
    list(model = m, anova = a, eta2 = e)
}


############################
# 5. Run models for all metrics
############################

res_shannon  <- run_alpha_long("Shannon")
res_simpson  <- run_alpha_long("Simpson")
res_invsimp  <- run_alpha_long("InvSimpson")
res_pielou   <- run_alpha_long("Pielou")


############################################################
# 6. Tidy ANOVA + effect size results
############################################################

tidy_alpha <- function(res, measure_name) {
    
    # Extract ANOVA table
    an <- as.data.frame(res$anova)
    an$Factor <- rownames(an)
    
    # Extract eta²
    e2 <- as.data.frame(res$eta2)
    e2$Factor <- e2$Parameter
    
    # Join and keep predictors of interest
    out <- dplyr::left_join(
        dplyr::select(an, Factor, Df, `F value`, `Pr(>F)`),
        dplyr::select(e2, Factor, Eta2_partial),
        by = "Factor"
    ) %>%
        dplyr::mutate(Measure = measure_name) %>%
        dplyr::filter(Factor %in% c("Geographic_region", "year", "month"))
    
    out
}


############################
# 7. Combine into summary table
############################

alpha_table <- bind_rows(
    tidy_alpha(res_shannon, "Shannon"),
    tidy_alpha(res_simpson, "Simpson"),
    tidy_alpha(res_invsimp, "InvSimpson"),
    tidy_alpha(res_pielou,  "Pielou")
)

alpha_table

# Adjust p-values across all tests (FDR correction)
alpha_table$P_adj_FDR <- p.adjust(
    alpha_table$`Pr(>F)`,
    method = "fdr"
)

alpha_table

# Optional export
# write.csv(alpha_table, "microeco_output/alpha_anova_time_eta2.csv", row.names = FALSE)


