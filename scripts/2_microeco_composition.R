################################
# microeco-diversity based
##############################

# calculate alpha diversity metrics
t1 <- trans_alpha$new(dataset = mt, group = "Geographic.region")
# return t1$data_stat
head(t1$data_stat)


# Krukal wallis test
t1$cal_diff(method = "anova")
# return t1$res_diff
head(t1$res_diff)


# calculate alpha diversity
t1 <- trans_alpha$new(dataset = mt, group = "Geographic.region")


t1$cal_diff(method = "anova", formula = "Geographic.region")

 #anova_diversity <- t1$res_diff
#write.csv(anova_diversity, "microeco_output/anova_diversity.csv")

t1$cal_diff(method = "anova")

# Make a plot
t1$plot_alpha(measure = "Shannon",
              add_sig_text_size = 6,
              plot_type = "ggviolin",
              add = "dotplot",
              xtext_size = 15,
              y_increase = 0.9,
              order_x_mean = TRUE,
              fill = "Geographic.region",
              alpha = 0.5)



################
# BETA diversity

# create trans_beta object
# For PCoA and NMDS, measure parameter must be provided.
# measure parameter should be either one of names(mt_rarefied$beta_diversity) or a customized symmetric matrix
t1 <- trans_beta$new(dataset = mt, group = "Geographic.region", measure = "bray")


# Make an ordination and visualize
t1$cal_ordination(method = "PCoA")

# plot the PCoA result with confidence ellipse
# in a biplot
biplot <- t1$plot_ordination(plot_color = "Geographic.region", plot_type = c("point")) +
  theme_bw() +
  theme(legend.background = element_rect(fill = NA, colour = "black"),
        legend.position = c(0.12, 0.8))
biplot

# in a tree
# use replace_name to set the label name, group parameter used to set the color
cluster <- t1$plot_clustering(group = "Geographic.region", replace_name = c("Geographic.region"))

ggarrange(biplot, cluster,
          widths = c(2,1),
          labels = "auto")

t1 <- trans_beta$new(dataset = mt, group = "Geographic.region", measure = "bray")

# Calculate perMANOVA (Permutational Multivariate Analysis of Vari-ance) based on the adonis2 function of vegan packagePerMANOVA(Anderson 2001) can be applied to the differential test of distances among groups
t1$cal_manova(manova_set = "Geographic.region")
t1$res_manova

write.csv(t1$res_manova, "microeco_output/manova.csv")

t1$cal_anosim()
t1$res_anosim


# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = F)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "KW")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(add = c("mean"),
                       plot_type = "ggviolin",
                       xtext_size = 15,
                       y_increase = 0.05,
                       order_x_mean = TRUE,
                       fill = "Geographic.region",
                       alpha = 0.5,
                       add_sig = T,
                       add_sig_text_size = 3)





###############
# MODELING PART
###############


# Differential abundance test
# It can find the taxa that vary across groups
# We use linear discriminant analysis effect size (LefSE) to determine the taxa that differ between groups
# We transform the data using transformation = 'AST' represents the arc sine square root transformation, to accountfor using relative abundances

palette <- c("northcal"   = "#1B9E77",
  "southcal"   = "#D95F02",
  "SC_island"  = "#7570B3",
  "centralcal" = "#E7298A",
  "sierras"    = "#66A61E",
  "desert"     = "#E6AB02"
  )



# For geographic region:

  #family level
t1 <- trans_diff$new(dataset = mt,
                     method = "lefse",
                     group = "Geographic.region",
                     alpha = 0.05,
                     lefse_subgroup = NULL,
                     p_adjust_method = "none",
                     taxa_level = "Family",
                     transformation = 'AST')

  # genus level
t2 <- trans_diff$new(dataset = mt,
                     method = "lefse",
                     group = "Geographic.region",
                     alpha = 0.05,
                     lefse_subgroup = NULL,
                     p_adjust_method = "none",
                     taxa_level = "Genus",
                     transformation = 'AST')

# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
lda_fam <- t1$plot_diff_bar(keep_prefix = F,
                            add_sig = T,
                            color_values = palette) +
  theme(axis.title.x = element_text(size = 12))
lda_fam


lda_gen <- t2$plot_diff_bar(keep_prefix = F,
                            add_sig = T,
                            color_values = palette) +
  theme(axis.title.x = element_text(size = 12))

lda_gen

ggarrange(lda_gen, lda_fam, common.legend = T, legend = "right")


t1$plot_diff_abund(fill = "Geographic.region", alpha = 0.5, add_sig = T)
t1$plot_diff_abund(coord_flip = FALSE)
t1$plot_diff_abund(plot_type = "errorbar")

fam <- t1$plot_diff_abund(#plot_type = "errorbar",
                          xtext_size = 10,
                          ytext_size = 10, ytitle_size = 12,
                          color_values = palette,
                   coord_flip = T,
                   #errorbar_color_black = TRUE,
                   errorbar_addpoint = T,
                   add_sig = F,
                   keep_prefix = F) + theme_bw()

fam


gen <- t2$plot_diff_abund(#plot_type = "errorbar",
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

ggarrange(fam, gen,
          labels = c("a) Family", "b) Genus"),
          legend = "right",
          common.legend = T)

legend <- get_legend(fam)


# all together plot
ggarrange(lda_fam, fam, lda_gen, gen,
          ncol = 2,
          nrow = 2,
          labels = "auto",
          widths = c(1.2, 3),
          align = "hv",
          common.legend = T,
          legend.grob = legend,
          legend = "right")

