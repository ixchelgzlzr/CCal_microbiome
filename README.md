# Off-target metagenomics: Leveraging whole genome resequencing data to characterize the bacteriome of *Calasterella californica* across California.
Public repository hosting the data and code for the microbiome analyses of the liverwort *Calasterella californica*.

This repository all the data necessary to replicate our results. In the folder `data` we provide the MASTER_bacteriome results (`data/MASTER_abundance_table_no_unclassified.txt`), as well as metadata information of the samples (`data/sample_info_all.csv`), and a few helper files useful to run the analyses. 

There are three folders with scripts:

### `scripts/`
Contains bash scripts to perform all the different steps of this pipeline, from cleaning reads to running metaphlan locally. The provided scripts rely on mantaining the folder structure of this repository. It also assumes that the raw reads are stored in `data/raw_reads/`.

* step 1: `scripts/1_clean_reads.sh` Clean reads with trimmomatic

* step 2: `scripts/2_rename_files.sh` Intermediate step to rename files with shorter names.

* Step 3: `scripts/3_filter_out_liverwort_reads.sh` Map reads to reference genome and keep only the ones that did not map out.

* Step 4: `scripts/4_run_metaphlan_locally.sh` Run metaphlan locally. Alternatively, if you can run this step in a computer cluster, check the instructions below on `savio_scripts`.

* Step 5, `scripts/5_combine_results.sh` Combine individual samples results into a single file for analysis. 

Additionally, there are two scripts that can be useful:

* `scripts/install_software.sh` commands to install MetaPhlan locally.

* `scripts/3.1_get_stats.sh` get the number of reads before and after cleaning.


### `savio_scripts/` 
Contains scripts to run metaphlan on a computer cluster. These scripts will be cluster dependent but they can serve as a guide. 
**Notes about Running metaphlan on savio:** This required a lot of troubleshooting. Eventually what worked was to make a conda enviornment, then install metaphlan using pip 
`pip install metaphlan`
Then downloading the database manually using the script savio_scripts/install_metaphlan.sh
The best approach was to run one sample at the time using multiple cores, instead of paralleliza multiple samples, because you quickly run into memory issues. There are a lot of details in the way metaphlan is set up that shouldn't be changed for SAVIO, e.g. --offline prevents metaphlan from trying to download newer version of the database on a htpc node, which would be super slow. Also, i point directly to the folder where my database lives. So much troubleshooting!! But it worled. The script savio_scripts/savio_sequential.sh == scripts/4_run_metaphlan_on_savio_sequential.sh, just trying to keep teh folder scripts self contained and with all the code. 


### `R_scripts/` 
Contains all the code necesarry to analyze the output of metaphlan. All the scripts are commented. This code heavily relies on the R package `microeco`, which is throuroughly documented at [https://chiliubio.github.io/microeco_tutorial/](https://chiliubio.github.io/microeco_tutorial/).

* `R_scripts/1_microeco_alpha_div.R` - Pre-processing and alpha diversity analyses.
* `R_scripts/2_microeco_composition.R` - Composition visualization and differential abundance analyses. 
* `R_scripts/3_processing_env_var_microeco.R` - Obtaining macroclimatic variables from bioclimatic layers.
* `R_scripts/4_Env_variables_microeco.R` - Testing for the effect of macroclimatic variables on microbial composition.
* `R_scripts/5_functional.R` - Profiling and visualizing functional composition of the microbiome, including differential abundance analyses. 






