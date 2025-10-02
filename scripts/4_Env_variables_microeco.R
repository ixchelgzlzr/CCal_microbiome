#########################################################
# Test for the effect of macroclimatic variables
# in community composition
#########################################################

# add the environmental data
mt$sample_table <- data.frame(mt$sample_table, df_env_points[rownames(mt$sample_table), ])


#selected variables :
#bio4  = Temperature Seasonality
#alt   = altitude
#bio16 = precipitation of wettest quarter
#bio12 = annual precipitation
t1 <- trans_env$new(dataset = mt, env_cols = c("bio4", "alt",  "bio12"))
colnames(t1$data_env) <- c("T_seasonality", "altitude", "annual_PP")

# use Wilcoxon Rank Sum Test as an example
t1$cal_diff(group = "Geographic.region", method = "wilcox")
head(t1$res_diff)


t1$cal_diff(method = "anova", group = "Geographic.region")
t1$res_diff

# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}


plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 2))


t1$cal_autocor()

#------------------------------
# This should also be a figure:
#------------------------------
t1$cal_autocor(group = "Geographic.region")


# use bray-curtis distance for CCA
t1$cal_ordination(method = "CCA",  taxa_level = "Family")

# show the orginal results
t1$trans_ordination()
t1$plot_ordination(plot_color = "Geographic.region")
# the main results of RDA are related with the projection and angles between arrows
# adjust the length of the arrows to show them better
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "Geographic.region",
                   plot_type = c("point"),
                   #centroid_segment_linetype = 1,
                   ellipse_chull_alpha = 0.05,
                   ellipse_chull_fill = T,
                   env_text_color = "blue",
                   color_values = palette,
                   env_arrow_color = "blue",
                   taxa_arrow_color = "black",
                   taxa_text_color = "black"
                   )

# test significance of ordination
t1$cal_ordination_anova()
t1$cal_ordination_envfit()

# see anova results
t1$res_ordination_axis
t1$res_ordination_terms


# Mantel test can be used to check whether there is significant correlations between environmental variables and distance matrix.
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel




