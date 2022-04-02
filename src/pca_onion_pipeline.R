
#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)
#ONION normalization with PCA: 
#Apply Principal Component Analysis(PCA) to create normalized data and integrate with SNF(Similarity Network Fusion)
#Input: DepMap 20q2 Achilles_gene_effect data preprocessed to replace NA values with gene-level mean (rows:genes, columns:cell lines)
#Output: .Rdata matrix w (gene-level integrated similarity-network)

#SNFtool - http://compbio.cs.toronto.edu/SNF/SNF/Software.html - SNFtool_v2.1.tar.gz

library(SNFtool)

#Load DepMap gene effect-size data (NA values replaced with gene-level mean)
input_file <- "../data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(input_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

#Set output directory
output_folder <- "../output/onion"

#Set SNF hyper-parameter values 	
K= 5 #number of neighbors
alpha= 3*.1 # hyperparameter
T = 10; 	# Number of Iterations

#apply Principal Component Analysis on data
pca <- prcomp(data, center = TRUE, scale. = TRUE)

#number of top principal components to remove in creating normalized data sets
pc_list<-c(1,3,5,7,9,11,13,15,17,19)

s_list <- list()

print('Normalizing and creating similarity networks ...')
for(pc in pc_list)
{
  print(pc)
  #Get Normalized data removing first 'pc' number of principal components
  v <- pca$rotation[,1:pc]
  projected <- as.matrix(data) %*% v %*% t(v)
  norm_data <- data - projected
  dimnames(norm_data) <- dimnames(data)
  
  #create gene-level similarity-network with PCC
  S_t <- data.frame(t(norm_data))
  sim_net <- cor(S_t, use = "all.obs",method="pearson")
  
  #create affinity-network from similarity-network (required step for integration with SNF)
  s = affinityMatrix(1-sim_net, K, alpha)
  
  #append affinity-network to integration list
  s_list <- append(s_list,list(s))

}

#Integrate affinity-networks with SNF
print('Integrating ...')
w = SNF(s_list, K, T)
print('end run SNF')

#Store integrated network as .Rdata (matrix variable name is w)
#Network is gene-level with genes as rows and columns
w <- as.data.frame(w)
row.names(w)<-row.names(data)
colnames(w)<-row.names(data)
save(w,file=file.path(output_folder, paste("snf_run_pca_10_5_3.Rdata",sep = '')))


