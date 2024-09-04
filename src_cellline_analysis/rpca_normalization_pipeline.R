library("rpca", lib.loc="../../../local/R/library/")
library("SNFtool", lib.loc="../../../local/R/library/")

#Load DepMap 20q2 Achilles_gene_effect data (NA values replaced with gene-level mean)
input_file <- "../data/depmap_q2_2020_nona_mean_t.tsv"
data <- read.csv(input_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

#Set output directory
output_folder <- "../output/rpca"

#Weights to regulate RPCA lambda hyperparameter values in creating normalized data sets
pc_list<-c(.7,.8,.9,1,1.1,1.2,1.3)

print('Normalizing ...')
for(pc in pc_list)
{
  print(pc)
  # set RPCA lambda hyperparameter (suqre-root of maximum between row-number(genes) and column-number(cell-lines))
  lam  = pc/(sqrt(max(nrow(data),ncol(data)))) # 1/(sqrt(max(# of genes,# of cell lines)))
  print(lam)
  #scale data prior to applying rpca
  data <- scale(data)
  #apply rpca to extract sparse component from original data
  pca <- rpca(data, lambda = lam)
  
  norm_data <- pca$S #Get RPCA-normalized data (variable S of rpca output contain sparse component or normlaized data)
  projected <- pca$L #Get RPCA-reconstructed data (variable L of rpca output contain low-rank component or reconstructed data)
  dimnames(norm_data) <- dimnames(data)
  dimnames(projected) <- dimnames(data)

#save normalized and reconstructed data
  write.table(norm_data, file.path(output_folder, paste('normalized_pca_', pc,".tsv", sep = "")), sep = "\t",
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(projected, file.path(output_folder, paste('projected_pca_', pc,".tsv", sep = "")), sep = "\t",
              quote = FALSE, row.names = TRUE, col.names = TRUE)
			}
