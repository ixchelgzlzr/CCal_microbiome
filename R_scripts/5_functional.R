
#########################################################
# FUNCTIONAL COMPOSITION
#########################################################

library(cowplot)

# create object of trans_func
t2 <- trans_func$new(mt)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database
t2$cal_spe_func(prok_database = "FAPROTAX")

t2$res_spe_func[1:5, 1:2]

# calculate the percentages for communities
# here do not consider the abundance
t2$cal_func_FR(abundance_weighted = F, perc = T)


t2$trans_func_FR()

func <- t2$plot_func_FR() +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 9),
        legend.position = "right")

func

# --------- save ---------
# save
# jpeg("plots/funcional.jpg", res = 300, units = "cm", width = 30, height = 18)
# func 
# dev.off()
# --------- save ---------

# --------- save ---------
# save as tiff
# tiff("plots/TIFFs/funcional.tiff", res = 300, units = "cm", width = 30, height = 18)
# func 
# dev.off()
# --------- save ---------


#differential abundance of functional groups

t4 <- trans_func$new(mt)
t4$cal_spe_func(prok_database = "FAPROTAX")
t4$cal_spe_func_perc(abundance_weighted = T)
# it is better to clone a dataset
tmp_mt <- clone(mt)
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp <- as.data.frame(t(t4$res_spe_func_perc), check.names = FALSE)
# assign the table back to taxa_abund list for further analysis
tmp_mt$taxa_abund$func <- tmp
# select the "func" in taxa_abund list in trans_diff
t5 <- trans_diff$new(dataset = tmp_mt, method = "lefse",
                     group = "Geographic_region",
                     taxa_level = "func",
                     p_adjust_method = "none")


palette <- c("Northcal"          = "#2DB27DFF",
             "Sierras"           = "#2E6E8EFF",
             "Centralcal"        = "#EB5760FF",
             "Peninsular ranges" = "#982D80FF",
             "SCI"               = "#66A61E",
             "Transverse ranges" = "#410F75FF",
             "Desert"            ="#E6AB02")


# plot differential abundance
diff_abund <- t5$plot_diff_abund(plot_type = "errorbar", use_number = 1:20,
                                 add_sig = T, alpha = 0.5, color_values = palette) +
  ggplot2::ylab("Relative abundance (%)") +
  theme(axis.text.y =  element_text(size = 10),
        axis.title.x = element_text(size = 10))

diff_abund


# again little hack to add lines behind
taxa_levels <- levels(diff_abund@data$Taxa)

bg_df <- data.frame(xmin = seq_along(taxa_levels) - 0.5,
                    xmax = seq_along(taxa_levels) + 0.5,
                    fill = rep(c("white", "grey80"), 
                               length.out = length(taxa_levels)))

band_layer <- geom_rect( data = bg_df,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                         inherit.aes = FALSE,
                         alpha = 0.6)

diff_abund_2 <- diff_abund
diff_abund_2$layers <- c(list(band_layer), diff_abund_2$layers)

diff_abund_banded <- diff_abund_2 + scale_fill_identity()
diff_abund_banded


# --------- save ---------
# save this figure 
# jpeg("plots/SUP_funcional_diff_ab.jpg", res = 300, units = "cm", width = 26, height = 20)
# diff_abund_banded 
# dev.off()
# 
# tiff("plots/TIFFs/SUP_funcional_diff_ab.tiff", res = 300, units = "cm", width = 26, height = 20)
# diff_abund_banded 
# dev.off()
# --------- save ---------






