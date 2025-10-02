# Off-target metagenomics: Leveraging whole genome resequencing data to characterize the bacteriome of *Calasterella californica* across California.
Public repository hosting the data and code for the microbiome analyses of the liverwort *Calasterella californica*.


This repository contains all the code and data necessary to replicate our results. All the data is in the folder `data`, and the code is in the folder `scripts`. The scripts are numbered by sequence of analyses, for example, running the second script requires running the script 1. 

* `scripts/1_microeco_alpha_div.R` - Pre-processing and alpha diversity analyses.
* `scripts/2_microeco_composition.R` - Composition visualization and differential abundance analyses. 
* `scripts/3_processing_env_var_microeco.R` - Obtaining macroclimatic variables from bioclimatic layers.
* `scripts/4_Env_variables_microeco.R` - Testing for the effect of macroclimatic variables on microbial composition.
* `scripts/5_functional.R` - Profiling and visualizing functional composition of the microbiome, including differential abundance analyses. 

All the scripts are commented. This code heavily relies on the R package `microeco`, which is throuroughly documented at [https://chiliubio.github.io/microeco_tutorial/](https://chiliubio.github.io/microeco_tutorial/).


A few other files that might be of interest are:
* `microeco_output/taxonomy_table.csv` - the taxonomic table of all the bacterial taxa recorded in _C. californica_ samples.
* `data/sample_info_2.xlsx` - Information of the _C. californica_ samples used in this study.