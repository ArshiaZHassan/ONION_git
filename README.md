# ONION_git
 Code repository for preprint 'Dimensionality reduction methods for extracting functional networks from large-scale CRISPR screens'
 https://doi.org/10.1101/2023.02.22.529573
 Repository for autoencoder-normalization code: https://github.com/csbio/ae-norm 
 Repository for input data and output files : https://zenodo.org/record/7671685#.Y_gi9nbMK5c
 
 ## Directory Set Up
 -Download ONION-git repo.
 
 -All source codes (R scripts) are in the ONION-git/src directory.
 
 -Create 'data' directory in ONION-git.
 
 -Input data required for this are uploaded here - https://zenodo.org/record/7671685#.Y_gi9nbMK5c
 
 -Create 'output' directory in ONION-git. Create 'ae','pca','rpca' and 'onion' directories in the 'output' directory.
 
 -Download FLEX and put in the same directory as the ONION-git.
 
 -Required directory and file organization
 
    .
     ├── ONION-git
     │   ├── src                     
     │   │   ├── ...                                                    # provided R scripts
     │   ├── data
     │   │   ├── 20Q2_GLS_fdr.npy                                       
     │   │   ├── 20Q2_GLS_sign.npy
     │   │   ├── Achilles_gene_effect.csv
     │   │   ├── depmap_q2_2020_nona_mean.tsv
     │   │   ├── depmap_q2_2020_nona_mean_rst_clp_mms.tsv
     │   │   ├── all_genes_20q2.txt
     │   │   ├── Mitochondial_genelist_1_26_2021_genes.tsv
     │   │   ├── olfactory_receptors.csv
     │   │   ├── ae_tanh_e1_depmap_20q2
     │   │   │   ├── 20q2_epochs_1_latent_1_normalized_ae.tsv
     │   │   │   ├── 20q2_epochs_1_latent_2_normalized_ae.tsv
     │   │   │   ├── 20q2_epochs_1_latent_3_normalized_ae.tsv
     │   │   │   ├── 20q2_epochs_1_latent_4_normalized_ae.tsv
     │   │   │   ├── 20q2_epochs_1_latent_5_normalized_ae.tsv
     │   │   │   ├── 20q2_epochs_1_latent_10_normalized_ae.tsv
     │   ├── output
     │   │   ├── ae
     │   │   │   ├── ...                                              # outputs are created after running R scripts
     │   │   ├── pca
     │   │   │   ├── ...                                              # outputs are created after running R scripts
     │   │   ├── rpca
     │   │   │   ├── ...                                              # outputs are created after running R scripts
     │   │   ├── onion
     │   │   │   ├── ...                                              # outputs are created after running R scripts
     ├── FLEX
     │   ├── ...                                                      # Directories from FLEX package
 
## Run instruction
-Required packages and dependencies should be installed prior to running the scripts. Required packages are listed in the 'R packages'.

-Run from inside src directory in the given order. Some scripts are dependent on outputs from other scripts. Please refer to script_dependency_flow_.pdf for script dependency.

-Please refer to data_flow_.pdf and data_flow_2_.pdf for script input/output flow.

     Rscript pre_process.R
     Rscript pca_normalization_pipeline.R
     Rscript rpca_normalization_pipeline.R
     Rscript ae_normalization_pipeline.R
     Rscript network_correlation_pipeline.R
     Rscript pca_onion_pipeline.R
     Rscript rpca_onion_pipeline.R
     Rscript ae_onion_pipeline.R
     Rscript pre_process_gls.R
     Rscript onion_evaluation_pipeline.R
     Rscript AUPRC_barplot_layer_evaluation_pipeline.R
     Rscript AUPRC_barplot_onion_evaluation_pipeline.R
 

## Script list
Pre-process DepMap data:

   pre_process.R

PCA-normalization:

   pca_normalization_pipeline.R

RPCA-normalization:

   rpca_normalization_pipeline.R

AE-normalization:

   ae_normalization_pipeline.R

PC-Onion:

   pca_onion_pipeline.R

RPC-Onion:

   rpca_onion_pipeline.R

AE-Onion:

   ae_onion_pipeline.R

Evaluation:

   network_correlation_pipeline.R

   pre_process_gls.R

   onion_evaluation_pipeline.R

   AUPRC_barplot_layer_evaluation_pipeline.R

   AUPRC_barplot_onion_evaluation_pipeline.R
