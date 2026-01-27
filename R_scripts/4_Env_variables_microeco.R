#########################################################
# Test for the effect of macroclimatic variables
# in community composition
#########################################################

# add the environmental data
mt$sample_table <- data.frame(mt$sample_table, df_env_points[rownames(mt$sample_table), ])

#selected variables :
#bio4  = Temperature Seasonality
#alt   = altitude
#bio12 = annual precipitation
t1 <- trans_env$new(dataset = mt, env_cols = c("bio4", "alt",  "bio12"))
colnames(t1$data_env) <- c("T_seasonality", "Altitude", "Annual_PP")

# check for significant differences in the environmental variables across greographic regions 
t1$cal_diff(method = "anova", group = "Geographic_region")
t1$res_diff

# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 11) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm")) +
      theme_bw()
}


# a complete figure with the three variables
ggarrange(tmp[[1]], tmp[[2]], tmp[[3]], nrow= 3, common.legend = T, legend = "right")


# --------save-------------#
# jpeg("plots/SUP_climatic_vars_per_region.jpg", res = 300, units = "cm", width = 20, height = 20)
# ggarrange(tmp[[1]], tmp[[2]], tmp[[3]], nrow= 3, common.legend = T, legend = "right")
# dev.off()
# --------save-------------#

# --------save-------------#
# save as tiff
# tiff("plots/TIFFs/SUP_climatic_vars_per_region.tiff", res = 300, units = "cm", width = 20, height = 20)
# ggarrange(tmp[[1]], tmp[[2]], tmp[[3]], nrow= 3, common.legend = T, legend = "right")
# dev.off()
# --------save-------------#


# #------------------------------
# # This should also be a figure:
# #------------------------------
t1$cal_autocor(group = "Geographic_region")

# use bray-curtis distance for CCA
t1$cal_ordination(method = "CCA",  taxa_level = "Family")

# show the orginal results
t1$trans_ordination(show_taxa = 20)
t1$plot_ordination(plot_color = "Geographic_region")

# adjust the length of the arrows to show them better
t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)

# t1$res_rda_trans is the transformed result for plotting
t1$plot_ordination(plot_color = "Geographic_region",
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



# --------save-------------#
#save ordination plot for CCA
# jpeg("plots/SUP_CCA.jpg", res = 300, units = "cm", width = 24, height = 16)
# 
# t1$plot_ordination(plot_color = "Geographic_region",
#                    plot_type = c("point", "ellipse"),
#                    #centroid_segment_linetype = 1,
#                    ellipse_chull_alpha = 0.05,
#                    ellipse_chull_fill = T,
#                    env_text_color = "blue",
#                    color_values = palette,
#                    env_arrow_color = "blue",
#                    taxa_arrow_color = "black",
#                    taxa_text_color = "black")
# dev.off()
# --------save-------------#


# --------save-------------#
#save ordination plot for CCA as TIFF
# tiff("plots/TIFFs/SUP_CCA.tiff", res = 300, units = "cm", width = 24, height = 16)
# 
# t1$plot_ordination(plot_color = "Geographic_region",
#                    plot_type = c("point", "ellipse"),
#                    #centroid_segment_linetype = 1,
#                    ellipse_chull_alpha = 0.05,
#                    ellipse_chull_fill = T,
#                    env_text_color = "blue",
#                    color_values = palette,
#                    env_arrow_color = "blue",
#                    taxa_arrow_color = "black",
#                    taxa_text_color = "black"
# )
# dev.off()
# --------save-------------#


# test significance of ordination
t1$cal_ordination_anova()
t1$cal_ordination_envfit()

# see anova results
t1$res_ordination_axis
t1$res_ordination_terms


# Save outputs
write.csv(t1$res_ordination_axis, "microeco_output/environmental_ordination_sig.csv")
write.csv(t1$res_ordination_terms, "microeco_output/env_correlation_with_CCA.csv")

# Mantel test can be used to check whether there is significant correlations between environmental variables and distance matrix.
t1$cal_mantel(use_measure = "bray")
# return t1$res_mantel
t1$res_mantel

write.csv(t1$res_mantel, "microeco_output/mantel_test_env_variables_on_microbiome.csv")



