##########################
# Using package microeco #
##########################

# install.packages("microeco")

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



# using the package file2meco for converting metaphlan files into the format required fo microeco

# define the files path
abund_path  <- "data/all_merged/calasterella_merged_all_2.tabular.txt"
sample_path <- "data/sample_info_2.txt"
match_table <- "data/match_table.csv"

# read data
mt <- mpa2meco(abund_path, sample_table = sample_path, match_table = match_table,  rel = T)


# This will remove the lines containing the taxa word regardless of taxonomic ranks and ignoring word case in the tax_table.
# So if you want to filter some taxa not considerd pollutions, please use subset like the previous operation to filter tax_table.
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))

# To clean taxonomy issues
mt$tidy_dataset()

# calculate abundances, using rel = true to indicate relative abundances are used
mt$cal_abund(rel = TRUE)


# make a rarefaction to correct for sequencing depth
# first clone the data
mt_rarefied <- clone(mt)
# use sample_sums to check the sequence numbers in each sample
mt_rarefied$sample_sums() %>% range

# As an example, use 10000 sequences in each sample
mt_rarefied$rarefy_samples(sample.size = 100)

mt_rarefied$sample_sums() %>% range

#--------------
# Alpha diversity
#-------------

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
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


#--------------
# Beta diversity
#-------------

# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
mt$cal_betadiv(unifrac = F)

# save beta_diversity to a directory
mt$save_betadiv(dirpath = "microeco_output/beta_diversity")



########################################
# COMPOSITION OF MICROBIAL COMMUNITIES #
########################################

# create trans_abund object
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7)

# plot of samples by geographic region
p1 <- t1$plot_bar(others_color = "grey60",
            facet = "Geographic.region",
            xtext_keep = T,
            legend_text_italic = FALSE,
            xtitle_keep = T, xtext_angle = 90)


# VISULIZE SUMMARY PER GROUP

# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7, groupmean = "Geographic.region")
g1 <- t1$plot_bar(others_color = "grey80", legend_text_italic = FALSE)
g1 <- g1 + theme_bw() + theme(axis.title.y = element_text(size = 18), axis.text.x = element_text(angle = 90))

g1

p1 <- p1 + theme(legend.position = "none")

p1





# Radar plots
t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 7, groupmean = "Geographic.region")
radar <- t1$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)
radar <- radar + theme(legend.position = "none",
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12)) +
                guides(color = guide_legend(override.aes = list(size = 5)))

radar


#--------------------
# Potential figure 1
#--------------------

# get legend from p1 for combined plot

leg <- get_legend(p1, legend = "right")


# combined plot
ggarrange( p1, radar,
           labels = "auto",
           common.legend = T,
           legend = "right",
           legend.grob = leg,
           widths = (2:1.2))



# show 15 taxa at family level
t1 <- trans_abund$new(dataset = mt, taxrank = "Family", ntaxa = 10)
t1$plot_box(group = "Geographic.region", xtext_angle = 90, plot_flip = T) + theme_light()


# show 40 taxa at Genus level in a heat map
# maybe a supp mat
t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 50)
g1 <- t1$plot_heatmap(facet = "Geographic.region", xtext_keep = F, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
g1
g1 + theme(axis.text.y = element_text(face = 'italic'))



#### clustering

t1 <- trans_abund$new(dataset = mt, taxrank = "Order", ntaxa = 8, groupmean = "Geographic.region")
g1 <- t1$plot_bar(coord_flip = TRUE)
g1 <- g1 + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
g1 <- t1$plot_bar(clustering_plot = TRUE)
# In this case, g1 (aplot object) is the combination of different ggplot objects
# to adjust the main plot, please select g1[[1]]
g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1



