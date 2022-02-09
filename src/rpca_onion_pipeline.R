
#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)
#ONION normalization with rpca: 
#Apply Robust Principal Component Analysis(RPCA) to create normalized data and integrate with SNF(Similarity Network Fusion)
#Input: DepMap 20q2 Achilles_gene_effect data preprocessed to replace NA values with gene-level mean
#Output: .Rdata matrix w (gene-level similarity network)

#rpca tool - https://cran.r-project.org/web/packages/rpca/index.html - rpca_0.2.3.tar.gz
#SNFtool - http://compbio.cs.toronto.edu/SNF/SNF/Software.html - SNFtool_v2.1.tar.gz

library("rpca", lib.loc="../../../local/R/library/")
library("SNFtool", lib.loc="../../../local/R/library/")

#Load DepMap gene effect-size file (NA values replaced with gene-level mean)
input_file <- "../data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(input_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

#Set output directory
output_folder <- "../output/onion"
	
#Set SNF hyper-parameter values
K= 5 #number of neighbors
alpha= 5*.1 # hyperparameter
T = 10; 	# Number of Iterations

#Set weight list to regulate RPCA lambda hyperparameter values in creating normalized data sets
pc_list<-c(.7,.8,.9,1,1.1,1.2,1.3)

s_list <- list()

print('Normalizing and creating similarity networks ...')
for(pc in pc_list)
{
  print(pc)
  # set RPCA lambda hyperparameter (suqre-root of maximum between row-number(genes) and column-number(cell-lines))
  lam  = pc/(sqrt(max(nrow(data),ncol(data))))
  print(lam)
  #scale data prior to applying rpca
  data <- scale(data)
  #apply rpca to extract sparse component from original data
  pca <- rpca(data, lambda = lam)
  
  #Get Normalized data (variable S of rpca output contain sparse component or normlaized data)
  S_ <- pca$S
  dimnames(S_) <- dimnames(data)
  
  #create gene-level similarity-network with PCC
  S_t <- data.frame(t(S_))
  sim_net <- cor(S_t, use = "all.obs",method="pearson")
  
  #create affinity matrix from similarity-network (required step for integration with SNF)
  s = affinityMatrix(1-sim_net, K, alpha)
  
  #append network to integration list
  s_list <- append(s_list,list(s))

}

#Integrate affinity networks with SNF
print('Integrating ...')
w = SNF(s_list, K, T)
print('end run SNF')

#Store integrated network as .Rdata (matrix variable name is w)
#Network is gene-level with genes as rows and columns
w <- as.data.frame(w)
row.names(w)<-row.names(data)
colnames(w)<-row.names(data)
save(w,file=file.path(output_folder, paste("snf_run_rpca_7_5_5.Rdata",sep = '')))


