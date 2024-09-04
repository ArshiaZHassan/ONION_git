library("rpca", lib.loc="../../../local/R/library/")
library("SNFtool", lib.loc="../../../local/R/library/")

#Load DepMap 20q2 Achilles_gene_effect data (NA values replaced with gene-level mean)
input_file <- "../data/depmap_q2_2020_nona_mean_t.tsv"
data <- read.csv(input_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

#Set output directory
output_folder <- "../output/pca"

#apply Principal Component Analysis on data
pca <- prcomp(data, center = TRUE, scale. = TRUE)

pc_list<-c(1,3,5,7,9,11,13,15,17,19)

print('Normalizing ...')
for(pc in pc_list)
{
  print(pc)
  #Get Normalized data removing first 'pc' number of principal components
  v <- pca$rotation[,1:pc]
  projected <- as.matrix(data) %*% v %*% t(v)
  norm_data <- data - projected
  dimnames(norm_data) <- dimnames(data)
  dimnames(projected) <- dimnames(data)

#save normalized and reconstructed data
  write.table(norm_data, file.path(output_folder, paste('normalized_pca_', pc,".tsv", sep = "")), sep = "\t",
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(projected, file.path(output_folder, paste('projected_pca_', pc,".tsv", sep = "")), sep = "\t",
              quote = FALSE, row.names = TRUE, col.names = TRUE)
			}
