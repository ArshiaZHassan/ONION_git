
#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)
#ONION normalization with AE: 
#Integrate auto-encoder(AE) generated normalized data using SNF(Similarity Network Fusion)
#Input:
	#DepMap 20q2 Achilles_gene_effect data preprocessed to replace NA values with gene-level mean (rows:genes, columns:cell lines) 
	#Auto-encoder generated normalized data from DepMap 20q2 Achilles_gene_effect data(rows:genes, columns:cell lines)
#Output: .Rdata matrix w (gene-level integrated similarity-network)

#SNFtool - http://compbio.cs.toronto.edu/SNF/SNF/Software.html - SNFtool_v2.1.tar.gz

library("SNFtool", lib.loc="../../../local/R/library/")

#Directory path to AE-normalized data generated from DepMap 20q2 Achilles_gene_effect data
input_folder <- "../data/ae_tanh_e1_depmap_20q2/"

#Set output directory
output_folder <- "../output/onion"

#Set SNF hyper-parameter values  	
K= 5 #number of neighbors
alpha= 5*.1 # hyperparameter
T = 10; 	# Number of Iterations

s_list <- list()

#Load DepMap 20q2 Achilles_gene_effect data (NA values replaced with gene-level mean)
original_file <- "../data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(original_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)
#Create gene-level similarity-network with PCC
data_t <- data.frame(t(data))
sim_net <- cor(data_t, use = "all.obs",method="pearson")
#Create affinity-network from similarity-network (required step for integration with SNF)
s = affinityMatrix(1-sim_net, K, alpha)
#Append affinity-network to integration list
s_list <- append(s_list,list(s))

#Latent-space sizes used in AE-normalization
pc_list<-c(1,2,3,4,5,10)

print('Normalizing and creating similarity networks ...')
for(pc in pc_list)
{
  print(pc)
  #Load AE-normalized data for latent-space size 'pc'
  input_file_name <- paste(input_folder,'20q2_epochs_1_latent_',pc,'_normalized_ae.tsv',sep = '')
  norm_data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1)
  
  #create gene-level similarity-network with PCC
  norm_data_t <- data.frame(t(norm_data))
  sim_net <- cor(norm_data_t, use = "all.obs",method="pearson")
  
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
save(w,file=file.path(output_folder, paste("snf_run_ae_7_5_5.Rdata",sep = '')))


