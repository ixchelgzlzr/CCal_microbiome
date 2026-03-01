#############################
# microeco-diversity based  #
#############################

# calculate alpha diversity 
t1 <- trans_alpha$new(dataset = mt, group = "Geographic_region")

# Anova to test for geographic region effect on diversity
t1$cal_diff(method = "anova", formula = "Geographic_region")


# return t1$res_diff
t1$res_diff

# save anova results
anova_diversity <- t1$res_diff
#write.csv(anova_diversity, "microeco_output/anova_diversity.csv")


#----------------#
# BETA diversity #
#----------------#

# create trans_beta object
# For PCoA and NMDS, measure parameter must be provided.
# measure parameter should be either one of names(mt_rarefied$beta_diversity) or a customized symmetric matrix
t1 <- trans_beta$new(dataset = mt, group = "Geographic_region", measure = "bray")


# Make an ordination and visualize
t1$cal_ordination(method = "NMDS", ncomp = 3)

# set up a custome palette
palette <- c("Northcal"          = "#2DB27DFF",
             "Sierras"           = "#2E6E8EFF",
             "Centralcal"        = "#EB5760FF",
             "Peninsular ranges" = "#982D80FF",
             "SCI"               = "#66A61E",
             "Transverse ranges" = "#410F75FF",
             "Desert"            ="#E6AB02")

# plot the PCoA result with confidence ellipse
# in a biplot
biplot <- t1$plot_ordination(plot_color = "Geographic_region", plot_type = c("point", "ellipse"), color_values = palette, ellipse_chull_alpha = 0.05, NMDS_stress_pos = c(1.3, 1)) +
  theme_bw() +
  theme(legend.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.12, 0.8))
biplot

# --------- save ---------
# Save the NMDS
# jpeg("plots/NMDS.jpg", res = 300, width= 26, height= 18, units = "cm")
# biplot
# dev.off()
# --------- save ---------


#---------tryng a 3D approach
# uncomment lines below for a 3D visualization

# library(plotly)
# 
# df <- t1$res_ordination$scores
# 
# plot_ly(
#     data = df,
#     x = ~MDS1, y = ~MDS2, z = ~MDS3,
#     type = "scatter3d", mode = "markers",
#     color = ~Geographic_region,
#     colors = palette,
#     marker = list(size = 5, opacity = 0.85)
# ) %>%
#     layout(
#         scene = list(
#             xaxis = list(title = "NMDS1"),
#             yaxis = list(title = "NMDS2"),
#             zaxis = list(title = "NMDS3")
#         )
#     )


# Calculate perMANOVA (Permutational Multivariate Analysis of Vari-ance) based on the adonis2 function of vegan packagePerMANOVA(Anderson 2001) can be applied to the differential test of distances among groups
t1$cal_manova(group = "Geographic_region")
t1$res_manova

#write.csv(t1$res_manova, "microeco_output/perMANOVA.csv")



###############
# MODELING PART
###############

# Differential abundance test
# It can find the taxa that vary across groups
# We use linear discriminant analysis effect size (LefSE) to determine the taxa that differ between groups
# We transform the data using transformation = 'AST' represents the arc sine square root transformation, to accountfor using relative abundances

palette <- c("Northcal"          = "#2DB27DFF",
             "Sierras"           = "#2E6E8EFF",
             "Centralcal"        = "#EB5760FF",
             "Peninsular ranges" = "#982D80FF",
             "SCI"               = "#66A61E",
             "Transverse ranges" = "#410F75FF",
             "Desert"            ="#E6AB02")


# For geographic region:

# order level
t0 <- trans_diff$new(dataset = mt,
                     method = "lefse",
                     group = "Geographic_region",
                     alpha = 0.05,
                     lefse_subgroup = NULL,
                     p_adjust_method = "none",
                     taxa_level = "Order", 
                     transformation = 'AST')


# family level
t1 <- trans_diff$new(dataset = mt,
                     method = "lefse",
                     group = "Geographic_region",
                     alpha = 0.05,
                     lefse_subgroup = NULL,
                     p_adjust_method = "none",
                     taxa_level = "Family", 
                     transformation = 'AST')

# genus level
t2 <- trans_diff$new(dataset = mt,
                     method = "lefse",
                     group = "Geographic_region",
                     alpha = 0.05,
                     lefse_subgroup = NULL,
                     p_adjust_method = "none",
                     taxa_level = "Genus",
                     transformation = 'AST')

#------------------------------------------------
# Groups with significant differential abundance
#-------------------------------------------------

# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
lda_ord <- t0$plot_diff_bar(keep_prefix = F, 
                            threshold = 3,
                            add_sig = T,
                            color_values = palette) +
    theme(axis.title.x = element_text(size = 12)) + theme_bw() 
lda_ord

lda_fam <- t1$plot_diff_bar(keep_prefix = F, 
                            threshold = 3,
                            add_sig = T,
                            color_values = palette) +
  theme(axis.title.x = element_text(size = 12)) + theme_bw() 
lda_fam


lda_gen <- t2$plot_diff_bar(keep_prefix = F,
                            add_sig = T,
                            threshold = 3,
                            color_values = palette) +
  theme(axis.title.x = element_text(size = 12)) + theme_bw() +
    labs(fill = "Geographic region",
         color = "Geographic region")


lda_gen

# side by side

###########################################################################################
# Manuscript figure: Significant differential abundance at the family and the genus level #
###########################################################################################

ggarrange(
          ncol = 2, labels = "auto", 
          common.legend = T, legend = "bottom")

# --------- save ---------
# jpeg("plots/sig_differential_abundance.jpg", res = 300, units = "cm", width = 20, height = 21)
# ggarrange(lda_gen, lda_fam,
#           ncol = 2, labels = "auto",
#           common.legend = T, legend = "bottom")
# dev.off()
# --------- save ---------



# --------- save ---------
#save as tiff
# tiff("plots/TIFFs/sig_differential_abundance.tiff", res = 300, units = "cm", width = 20, height = 21)
# ggarrange(lda_gen, lda_fam, 
#           ncol = 2, labels = "auto",
#           common.legend = T, legend = "bottom")
# dev.off()
# --------- save ---------


#------------------------------------------------
# Visualize the comparison with bar plot of all the groups
#-------------------------------------------------

# At the family level
fam <- t1$plot_diff_abund(plot_type = "errorbar",
                          use_number = 1:10,
                          xtext_size = 10,
                          ytext_size = 10, ytitle_size = 12,
                          color_values = palette,
                           coord_flip = T,
                           #errorbar_color_black = TRUE,
                           errorbar_addpoint = T,
                           add_sig = F,
                            keep_prefix = F) + theme_bw()

fam

# little hack to add some bars in the back for easier reading
taxa_lev_fam <- levels(fam@data$Taxa)

bg_df_fam <- data.frame(xmin = seq_along(taxa_lev_fam) - 0.5,
                    xmax = seq_along(taxa_lev_fam) + 0.5,
                    fill = rep(c("white", "grey80"), 
                               length.out = length(taxa_lev_fam)))

band_layer_fam <- geom_rect( data = bg_df_fam,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                         inherit.aes = FALSE,
                         alpha = 0.6)

fam3 <- fam
fam3$layers <- c(list(band_layer_fam), fam$layers)

fam_banded <- fam3 + scale_fill_identity()
fam_banded


# differential abundance at the genus level
gen <- t2$plot_diff_abund(plot_type = "errorbar",
                           coord_flip = T,
                           #errorbar_color_black = TRUE,
                           errorbar_addpoint = T,
                           xtext_size = 10,
                           ytext_size = 10,
                           ytitle_size = 12,
                           color_values = palette,
                           add_sig = F,
                           keep_prefix = F) + theme_bw()

gen

# some manual adhoc code to add some grey bands to improve visibility
# get taxa levels
taxa_levels <- levels(gen@data$Taxa)

bg_df <- data.frame(xmin = seq_along(taxa_levels) - 0.5,
                    xmax = seq_along(taxa_levels) + 0.5,
                    fill = rep(c("white", "grey80"), 
                    length.out = length(taxa_levels)))

band_layer <- geom_rect( data = bg_df,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                         inherit.aes = FALSE,
                         alpha = 0.6)

gen3 <- gen
gen3$layers <- c(list(band_layer), gen$layers)

gen_banded <- gen3 + scale_fill_identity()
gen_banded

legend <- get_legend(fam)




# maybe for sup mat
ggarrange(fam_banded, gen_banded,
          ncol = 2,
          labels = "auto",
          widths = c(1, 1),
          align = "hv",
          common.legend = T,
          legend.grob = legend,
          legend = "right")



# --------- save ---------
# jpeg("plots/SUP_diff_abundance_detailed.jpg", res=300, width = 32, height = 18, units = "cm")
# 
# ggarrange(fam_banded, gen_banded,
#           ncol = 2,
#           labels = "auto",
#           widths = c(1, 1),
#           align = "hv",
#           common.legend = T,
#           legend.grob = legend,
#           legend = "right")
# dev.off()
# --------- save ---------


# --------- save ---------
#save tiff
# tiff("plots/TIFFs/SUP_diff_abundance_detailed.tiff", res=300, width = 32, height = 18, units = "cm")
# 
# ggarrange(fam_banded, gen_banded,
#           ncol = 2,
#           labels = "auto",
#           widths = c(1, 1),
#           align = "hv",
#           common.legend = T,
#           legend.grob = legend,
#           legend = "right")
# dev.off()
# --------- save ---------
