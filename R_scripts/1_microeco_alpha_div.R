##########################
# Using package microeco #
##########################

#install.packages("microeco")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("Rhdf5lib")
# if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# install.packages("file2meco", repos = BiocManager::repositories())
# install.packages("ggalluvial")
# install.packages("ggradar")
# devtools::install_github("gmteunisse/ggnested")
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)


# libraries
library(microeco)
library(file2meco)
library(dplyr)
library(ggalluvial)
library(ggradar)
library(ggpubr)
library(ggplotify)


# using the package file2meco for converting metaphlan files into the format required for microeco

# define the files path
abund_path  <- "data/MASTER_abundance_table_no_unclassified.txt"
sample_path <- "data/sample_info_all.csv"

# make match table required for mpa2meco
match_table <- read.table("data/samples_all.txt")
match_table$V2 <- match_table$V1

# read data
mt <- mpa2meco(abund_path, sample_table = sample_path, match_table = match_table,  rel = T)

# Some standard cleaning just in case 

# organelles
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))

# To clean taxonomy issues
mt$tidy_dataset()

# calculate abundances, using rel = true to indicate relative abundances are used
mt$cal_abund(rel = TRUE)

#-----------------#
# Alpha diversity #
#-----------------#

# compute alpha diversity 
mt$cal_alphadiv()

# return alpha_diversity in the object
mt$alpha_diversity

# save alpha_diversity to a directory
mt$save_alphadiv(dirpath = "microeco_output/alpha_diversity")

# alpha diversity plot
alpha <-mt$alpha_diversity
alpha$sample <- rownames(alpha)
ggplot(alpha, aes (x = sample , y = Shannon)) + geom_col()

# save the taxonomy table to file
mt$tax_table -> taxonomy_table
#write.csv(taxonomy_table, "microeco_output/taxonomy_table.csv")

#----------------#
# Beta diversity #
#----------------#

# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
mt$cal_betadiv(unifrac = F)

# save beta_diversity to a directory
# this will write down matrices of pairwise distances between samples
mt$save_betadiv(dirpath = "microeco_output/beta_diversity")


########################################
# COMPOSITION OF MICROBIAL COMMUNITIES #
########################################

# create trans_abund object at the taxonomic Order level
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7)

# plot of samples by geographic region
p1 <- t1$plot_bar(others_color = "grey70",
            facet = "Geographic_region",
            xtext_keep = T,
            legend_text_italic = FALSE,
            ytitle_size = 14, xtext_size = 8,
            xtitle_keep = T, xtext_angle = 90)
p1


# --------- save ---------
# jpeg("plots/composition_rel_ab.png", units = "cm", width = 36, height = 22, res = 200)
# p1
# dev.off()
# --------- save ---------

# VISUALIZE SUMMARY PER GROUP
# let's try ans summarize composition by geographic region
# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7, groupmean = "Geographic_region")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 <- g1 + theme_bw() + theme(axis.title.y = element_text(size = 12), axis.text.x = element_text(angle = 90))

g1

p1 <- p1 + theme(legend.position = "none")


#### clustering
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7, groupmean = "Geographic_region")
g1 <- t1$plot_bar(coord_flip = T, others_color = "grey70")
g1 <- t1$plot_bar(clustering_plot = TRUE, use_alluvium = F, coord_flip = F, others_color = "grey70")
g1
# In this case, g1 (aplot object) is the combination of different ggplot objects
# to adjust the main plot, please select g1[[1]]
g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1

# make gg object
g1_gg <-as.ggplot(g1) 

# second version of composition plot
ggarrange(p1, g1_gg,
          nrow = 2, heights = c(1.3, 1),
          labels = "auto")


# --------- save ---------
# jpeg("plots/relative_ab_composite.jpg", res = 300, width = 26, height = 22, units = "cm" )
# ggarrange(p1, g1_gg,
#           nrow = 2, heights = c(1.4, 1),
#           labels = "auto")
# dev.off()

# #save as tiff
# tiff("plots/TIFFs/relative_ab_composite.tiff", res = 300, width = 26, height = 22, units = "cm" )
# ggarrange(p1, g1_gg,
#           nrow = 2, heights = c(1.4, 1),
#           labels = "auto")
# dev.off()
# --------- save ---------


