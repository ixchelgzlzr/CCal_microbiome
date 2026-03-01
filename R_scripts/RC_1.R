############################################################
# Reviewer-response analysis: "Unclassified" MetaPhlAn reads
#
# Purpose:
#   Quantify whether the proportion of reads assigned to the
#   "unclassified" category varies systematically with:
#     - Geographic region
#     - Sequencing batch
#     - Latitude / longitude
#     - Plant sex
#     - Collection month / year
#
# Rationale:
#   Reviewers often worry that patterns in community composition
#   could be driven by technical artifacts or uneven taxonomic
#   resolution. If "unclassified" varies strongly with geography,
#   batch, or metadata, it could indicate bias in downstream results.
#
# Inputs:
#   - MetaPhlAn abundance table INCLUDING unclassified
#   - Sample metadata CSV
#
# Outputs:
#   - Diagnostics (boxplots, scatterplots)
#   - ANOVA / Welch tests / correlations
#   - Multi-panel figure summarizing key comparisons
############################################################


############################
# 1. Load required libraries
############################

library(microeco)    # (not directly used in this script, but often used elsewhere in pipeline)
library(file2meco)   # (not directly used here; used for MetaPhlAn→microeco conversions)
library(dplyr)       # data manipulation
library(tidyr)       # separate()
library(car)         # leveneTest()
library(ggplot2)     # plotting (implicitly needed; used via ggplot)
library(ggpubr)      # ggarrange()
library(ggplotify)   # (not directly used in this script)
library(ggalluvial)  # (not directly used in this script)
library(ggradar)     # (not directly used in this script)
# Note: ggpubr is loaded twice in the original; once is sufficient.


############################
# 2. Define input file paths
############################

# MetaPhlAn abundance table INCLUDING unclassified taxa
abund_path  <- "output/metaphlan_out_with_unclassified/profiles/abundance_table_with_unclassified.txt"

# Sample metadata (must include at least Sample, Geographic_region, Batch,
# latitude, longitude, sex, collection_date)
sample_path <- "data/sample_info_all.csv"


############################
# 3. Read abundance + metadata
############################

# MetaPhlAn abundance tables typically have a header line we skip (skip = 1)
data    <- read.table(abund_path, skip = 1)

# Sample metadata
samples <- read.csv(sample_path)


############################################################
# 4. Extract ONLY the "unclassified" signal
#
# Assumption encoded in this section:
#   The first two rows of the abundance table correspond to
#   the lines you need to recover "unclassified" abundance.
#
# What this code does:
#   - Keep first two rows only
#   - Transpose so samples become rows
#   - Rename columns to: Sample + unclassified
#   - Drop the first row after transpose (often contains headers)
#   - Convert unclassified column to numeric
############################################################

# keep only info about unclassified (as represented in this file)
data <- data[1:2, ]

# transpose so each sample is a row
data <- t(data)

# label columns after transpose
colnames(data) <- c("Sample", "unclassified")

# drop the first row (often not a real sample row after transpose)
data <- data.frame(data[2:nrow(data), ])

# ensure numeric for downstream stats
data$unclassified <- as.numeric(data$unclassified)


############################
# 5. Merge with metadata
############################

# Join unclassified values to sample metadata by Sample ID
df <- merge.data.frame(data, samples, by = "Sample")


############################################################
# 6. Hypothesis tests / diagnostics
############################################################

############################
# 6a. Geographic region effect
############################

# Visual check: distribution by region
geog_region <- boxplot(unclassified ~ Geographic_region, data = df)

# Check homogeneity of variance across regions (ANOVA assumption)
leveneTest(unclassified ~ Geographic_region, data = df)

# If variance is not significantly different, proceed with standard ANOVA
anova_region <- aov(unclassified ~ Geographic_region, data = df)
summary(anova_region)


############################
# 6b. Sequencing batch effect
############################

# Visual check by sequencing batch
boxplot(unclassified ~ Batch, data = df)

# ANOVA for batch effect
anova_sequencing <- aov(unclassified ~ Batch, data = df)
summary(anova_sequencing)


############################
# 6c. Spatial gradients: latitude and longitude
############################

# Latitude association (scatterplot + Pearson correlation)
plot(unclassified ~ latitude, data = df)
cor.test(df$unclassified, df$latitude, method = "pearson")

# Longitude association (scatterplot + Pearson correlation)
plot(unclassified ~ longitude, data = df)
cor.test(df$unclassified, df$longitude, method = "pearson")


############################
# 6d. Plant sex effect
############################

# Visual check by sex
boxplot(unclassified ~ sex, data = df)

# ANOVA for sex effect
anova_sex <- aov(unclassified ~ sex, data = df)
summary(anova_sex)


############################
# 6e. Collection date effects: month + year
############################

# Split collection_date into month/day/year (assumes "MM/DD/YYYY" format)
df <- separate(df, col = collection_date, into = c("month", "day", "year"), sep = "/")

# Month effect:
boxplot(unclassified ~ month, data = df)
leveneTest(unclassified ~ month, data = df)

anova_month <- aov(unclassified ~ month, data = df)
summary(anova_month)

# Year effect:
# If Levene’s test indicates unequal variance across years, use Welch’s test
leveneTest(unclassified ~ year, data = df)
boxplot(unclassified ~ year, data = df)

# Welch one-way test (robust to unequal variances)
oneway.test(unclassified ~ year, data = df, var.equal = FALSE)


############################################################
# 7. Build a multi-panel figure for the manuscript/rebuttal
#
# Panels:
#   A: unclassified by Geographic region (boxplot)
#   B: unclassified by Sequencing batch (boxplot)
#   C: unclassified by Plant sex (boxplot)
#   D: unclassified vs Latitude (scatter)
#   E: unclassified vs Longitude (scatter)
############################################################

# Panel A: geography
geog <- ggplot(df, aes(x = Geographic_region, y = unclassified)) +
    geom_boxplot(fill = "#35B779FF") +
    theme_bw() +
    xlab("Geographic region") +
    ylab("Proportion of unclassified reads")

# Panel C: sex
sex <- ggplot(df, aes(x = sex, y = unclassified)) +
    geom_boxplot(fill = "#3E4A89FF") +
    theme_bw() +
    xlab("Plant sex") +
    ylab("Proportion of unclassified reads")

# Panel B: sequencing batch (cast to factor so ggplot treats it categorically)
seq_batch <- ggplot(df, aes(x = as.factor(Batch), y = unclassified)) +
    geom_boxplot(fill = "#3E4A89FF") +
    theme_bw() +
    xlab("Sequencing batch") +
    ylab("Proportion of unclassified reads")

# Panel D: latitude
lat <- ggplot(df, aes(x = latitude, y = unclassified)) +
    geom_point(color = "#31688EFF", size = 2) +
    theme_bw() +
    xlab("Latitude") +
    ylab("Proportion of unclassified reads")

# Panel E: longitude
long <- ggplot(df, aes(x = longitude, y = unclassified)) +
    geom_point(color = "#31688EFF", size = 2) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Proportion of unclassified reads")


############################
# 8. Arrange panels
############################

# Top row: B (batch) and C (sex)
top_row <- ggarrange(
    seq_batch, sex,
    ncol = 2,
    labels = c("B", "C")
)

# Bottom row: D (lat) and E (long)
bottom_row <- ggarrange(
    lat, long,
    ncol = 2,
    labels = c("D", "E")
)

# Full figure: A on top, then (B,C), then (D,E)
ggarrange(
    geog,
    top_row,
    bottom_row,
    nrow = 3,
    labels = c("A", "", "")
)