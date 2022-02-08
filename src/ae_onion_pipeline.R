
library("SNFtool", lib.loc="../../../local/R/library/")

start_time <- Sys.time()
#args = commandArgs(trailingOnly=TRUE)

start_time <- Sys.time()
input_folder <- "../data/ae_tanh_e1_depmap_20q2/"
output_folder <- "../output/onion"

## set all the parameters:
#K = as.numeric(args[1]);		
#alpha = as.numeric(args[2])*.1;  	
K= 5 #number of neighbors, usually (10~30)
alpha= 5*.1 # hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

s_list <- list()

#include original data in the integration process
original_file <- "../data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(original_file, sep = "\t", header = TRUE, row.names = 1,stringsAsFactors=F)
data_t <- data.frame(t(data))
sim_net <- cor(data_t, use = "all.obs",method="pearson")
s = affinityMatrix(1-sim_net, K, alpha)
s_list <- append(s_list,list(s))

pc_list<-c(1,2,3,4,5,10)

print('Normalizing and creating similarity networks ...')
for(pc in pc_list)
{
  print(pc)
  input_file_name <- paste(input_folder,'20q2_epochs_1_latent_',pc,'_normalized_ae.tsv',sep = '')
  norm_data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1)
  norm_data_t <- data.frame(t(norm_data))
  sim_net <- cor(norm_data_t, use = "all.obs",method="pearson")
  
  s = affinityMatrix(1-sim_net, K, alpha)
  s_list <- append(s_list,list(s))

}

print('Integrating ...')
w = SNF(s_list, K, T)
print('end run SNF')
w <- as.data.frame(w)
row.names(w)<-row.names(data)
colnames(w)<-row.names(data)
save(w,file=file.path(output_folder, paste("snf_run_ae_7_5_5.Rdata",sep = '')))

end_time <- Sys.time()
print("Time: ")
print(end_time - start_time)

print("Memory profile: ")
print(memory.profile())

print("Object-wise Memory profile: ")

sapply(ls(), function(x) {print(object.size(get(x)),standard = "legacy", units="Gb") })

gc()

