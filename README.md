# ONION_git
 onion_manuscript local rep
 
 ## Directory Set Up
 -Download ONION-git repo.
 -All source codes (R scripts) are in the ONION-git/src directory.
 -Create 'input' directory in ONOION-git. The following files should be in the 'input' directory to successfully run all the scripts.
 -Create 'output' directory in ONOION-git. Create 'ae','pca','rpca' and 'onion' directories in the 'output' directory.
 -Download FLEX and put in the same directory as the ONION-git.
 -Directory and file organization
 	-ONION-git
   -src
    -<provided R scripts>
   -input
    -20Q2_GLS_fdr.npy
    -20Q2_GLS_sign.npy
    -Achilles_gene_effect.csv
    -all_genes_20q2.txt
    -depmap_q2_2020_nona_mean.tsv
    -depmap_q2_2020_nona_mean_rst_clp_mms.tsv
    -Mitochondial_genelist_1_26_2021_genes.tsv
    -olfactory_receptors.csv
    -ae_tanh_e1_depmap_20q2
     -20q2_epochs_1_latent_1_normalized_ae.tsv
     -20q2_epochs_1_latent_2_normalized_ae.tsv
     -20q2_epochs_1_latent_3_normalized_ae.tsv
     -20q2_epochs_1_latent_4_normalized_ae.tsv
     -20q2_epochs_1_latent_5_normalized_ae.tsv
     -20q2_epochs_1_latent_10_normalized_ae.tsv
   -output
    -ae
    -pca
    -rpca
    -onion
  -FLEX


## Workflow
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
