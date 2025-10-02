
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
t2$cal_spe_func_perc(abundance_weighted = F)


t2$trans_spe_func_perc()

func <- t2$plot_spe_func_perc() +
  theme(axis.text.x = element_text(angle = 90, size = 11),
        legend.position = "right")

func


t3 <- trans_env$new(dataset = mt, env_cols = c("bio4", "alt", "bio12"))
colnames(t3$data_env) <- c("T_seasonality", "altitude", "annual_PP")

t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
t3$plot_cor()



# calculate correlations for different groups using parameter by_group
t3$cal_cor(add_abund_table = t2$res_spe_func_perc,
           cor_method = "spearman",
           by_group = "Geographic.region",
           p_adjust_method = "fdr")
# return t1$res_cor
t3$plot_cor()



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
                     group = "Geographic.region",
                     taxa_level = "func",
                     p_adjust_method = "none")


# plot differential abundance
diff_abund <- t5$plot_diff_abund(plot_type = "errorbar",
                                 add_sig = T, alpha = 0.5) +
  ggplot2::ylab("Relative abundance (%)") +
  theme(axis.text.y =  element_text(size = 10),
        axis.title.x = element_text(size = 10))



#combined plots
plot_grid(func, diff_abund,
          labels = "auto",
          ncol = 1,
          align = "v",
          axis = "lr",
          rel_heights = c(2,1.5))


