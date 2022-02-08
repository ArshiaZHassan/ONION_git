
library("SNFtool", lib.loc="../../../local/R/library/")

start_time <- Sys.time()
#args = commandArgs(trailingOnly=TRUE)

input_file <- "../data/depmap_q2_2020_nona_mean.tsv"
output_folder <- "../output/onion"
data <- read.csv(input_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)

## set all the parameters:
#K = as.numeric(args[1]);		
#alpha = as.numeric(args[2])*.1;  	
K= 5 #number of neighbors, usually (10~30)
alpha= 3*.1 # hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

s_list <- list()

pca <- prcomp(data, center = TRUE, scale. = TRUE)

pc_list<-c(1,3,5,7,9,11,13,15,17,19)

print('Normalizing and creating similarity networks ...')
for(pc in pc_list)
{
  print(pc)
  v <- pca$rotation[,1:pc]
  projected <- as.matrix(data) %*% v %*% t(v)
  norm_data <- data - projected
  dimnames(norm_data) <- dimnames(data)
  
  S_t <- data.frame(t(norm_data))
  sim_net <- cor(S_t, use = "all.obs",method="pearson")
  
  s = affinityMatrix(1-sim_net, K, alpha)
  s_list <- append(s_list,list(s))

}

print('Integrating ...')
w = SNF(s_list, K, T)
print('end run SNF')
w <- as.data.frame(w)
row.names(w)<-row.names(data)
colnames(w)<-row.names(data)
save(w,file=file.path(output_folder, paste("snf_run_pca_10_5_3.Rdata",sep = '')))

end_time <- Sys.time()
print("Time: ")
print(end_time - start_time)

print("Memory profile: ")
print(memory.profile())

print("Object-wise Memory profile: ")

sapply(ls(), function(x) {print(object.size(get(x)),standard = "legacy", units="Gb") })

gc()

