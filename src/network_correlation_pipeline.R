#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)
#Create network correlation scatterplots
#Input:
	#DepMap 20q2 Achilles_gene_effect data preprocessed to replace NA values with gene-level mean (rows:genes, columns:cell lines) 
	#Auto-encoder generated normalized data from DepMap 20q2 Achilles_gene_effect data(rows:genes, columns:cell lines)
	#DepMap 20q2 Achilles_gene_effect data preprocessed 
		# replace NA values with gene-level mean (rows:genes, columns:cell lines)
		# row-standardized, clipped, min-max scaled [-1,1]
#Output: pdf file (scatter plot of PCC between un-normalized and processed networks)


#rpca tool - https://cran.r-project.org/web/packages/rpca/index.html - rpca_0.2.3.tar.gz

library("ramify", lib.loc="../../../local/R/library/")
library("crayon", lib.loc="../../../local/R/library/")
library("pillar", lib.loc="../../../local/R/library/")
library("withr", lib.loc="../../../local/R/library/")
library("labeling", lib.loc="../../../local/R/library/")
library("farver", lib.loc="../../../local/R/library/")
library("digest", lib.loc="../../../local/R/library/")
library("ggplot2", lib.loc="../../../local/R/library/")
library("RColorBrewer", lib.loc="../../../local/R/library/")

library("rpca", lib.loc="../../../local/R/library/")

# Function:
    # create_network_correlation_scatter_plot
# Arguments:
    # cor_all : PCC Data; type: data.frame (required); columns: PCC norm. net., PCC recon. net, method parameter value
	# cols1 : color list
	# output_folder : output file path
	# outfile_name : output file name
# Output: 
    # pdf file (scatter plot of PCC between un-normalized and processed networks)
# Description:
    # Produce scatter plot of PCC between un-normalized and processed networks
create_network_correlation_scatter_plot<-function(cor_all, cols1, output_folder, outfile_name){
print("In plot")
	ggplot(data =cor_all, aes(x=as.character(pc))) +
	geom_point(data = cor_all, aes(y=`s`								 
								 ,shape = "Normalized"),fill=cols1[1]
			 , colour=cols1[4],alpha = 1,size=15) + 
	geom_point(data = cor_all, aes(y=`l`								 
								 ,shape = "Reconstructed"),fill=cols1[1]
			 , colour=cols1[4],alpha = 1,size=15) + 
	xlab("Layers") +
	ylab("PCC with un-norm. network") +
	scale_y_continuous(limits= seq(0,1))+
	scale_x_discrete(limits = as.character(cor_all$pc), labels = as.character(cor_all$pc))+
	scale_color_manual(
					 values = c(cols1[3], cols1[3])) +
	scale_shape_manual(limits = c("Normalized","Reconstructed"),
					 labels = c("Normalized","Reconstructed"), 
					 #values=c(16,17)) +
					  values=c(21,24)) +
	theme_bw() +
	labs(
	shape = ''#"Pearson Correlation\n between networks"
	) +
	theme(panel.border = element_blank(), #panel.grid.major = element_blank(), #axis.text.x = element_blank(),
		panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		axis.line = element_line(colour = "black",size = 2),
		legend.title = element_text(family="ArialMT",size = 30,colour = "black"),
		legend.text = element_text(family="ArialMT",size = 30,colour = "black"),
		axis.title = element_text(family="ArialMT",size = 30,colour = "black"),
		axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
		axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
		#axis.text = element_blank(),
		axis.text = element_text(family="ArialMT",size = 30,colour = "black"),
		axis.text.x = element_text(angle = 90)
		#,axis.ticks.length.y = unit(1, "cm") #, legend.position = "none"
	)
	aspect_ratio <- 1.5
	ggsave(file.path(output_folder, outfile_name),height = 7 , width = 7*1.5, useDingbats=FALSE)

}

#########################################################################################################################
#load depmap data
input_file_name <- "../data/depmap_q2_2020_nona_mean.tsv"
data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1)
#create similarity network from normalized data
data_t <- data.frame(t(data)) #transpose data - genes as rows
sim_net <- cor(data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
data_flat <- flatten(as.matrix(sim_net))
#########################################################################################################################
#Create network-correlation scatter plots for AE-normalized and AE-reconstructed data
print("AE")
ae_folder <- "../data/ae_tanh_e1_depmap_20q2/"

#load pre-processed data - row-standardized, clipped, min-max scaled [-1,1]
pre_pro_file_name <- "../data/depmap_q2_2020_nona_mean_rst_clp_mms.tsv"
pre_pro_data <- read.csv(pre_pro_file_name, sep = "\t", header = TRUE, row.names = 1)

cor_list_data_norm <- c()
cor_list_data_proj <- c()

#Set output directory
output_folder <- "../output/ae/"
#Latent-space sizes used in AE-normalization
layer_list<-c(1,2,3,4,5,10)
for(layer in layer_list)
{
	print(layer)
  #generate normalized and reconstructed data
  #Load AE-normalized data for latent-space size 'layer'
  input_file_name <- paste(ae_folder,'20q2_epochs_1_latent_',layer,'_normalized_ae.tsv',sep = '')
  norm_data <- read.csv(input_file_name, sep = "\t", header = TRUE, row.names = 1)
	#use pre-processed data to extract reconstructed data
  projected <- pre_pro_data - norm_data
  dimnames(norm_data) <- dimnames(data)
  dimnames(projected) <- dimnames(data)
  
  print("Normalized data")
  #create similarity network from normalized data
  norm_data_t <- data.frame(t(norm_data)) #transpose data - genes as columns
  sim_net <- cor(norm_data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between normalized-network and un-normalized-network
  norm_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_norm <-c(cor_list_data_norm,cor(data_flat, norm_data_flat, method="pearson"))
  
  print("Projected data")
  #create similarity network from reconstructed data
  projected_t <- data.frame(t(projected)) #transpose data - genes as columns
  sim_net <- cor(projected_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between reconstructed-network and un-normalized-network
  proj_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_proj <-c(cor_list_data_proj,cor(data_flat, proj_data_flat, method="pearson"))
  
}

#plot network PCC values
cor_all <- data.frame(cor_list_data_norm)
names(cor_all)<-c('s')
cor_all['l'] <- cor_list_data_proj
cor_all['pc'] <- layer_list

cols1<- brewer.pal(n = 9, "BuPu")[6:9]

create_network_correlation_scatter_plot(cor_all, cols1, output_folder, "ae_pcc_real_nw_vs_norm_nw_6_scatter.pdf")
#########################################################################################################################
#Create network-correlation scatter plots for PCA-normalized and PCA-reconstructed data
print("PCA")

cor_list_data_norm <- c()
cor_list_data_proj <- c()
#Set output directory
output_folder <- "../output/pca/"

#apply Principal Component Analysis on data
pca <- prcomp(data, center = TRUE, scale. = TRUE)

#number of top principal components to remove in creating normalized and reconstructed data
layer_list<-c(1,3,5,7,9,11,13,15,17,19)
for(layer in layer_list)
{
	print(layer)
  #generate normalized and reconstructed data
  v <- pca$rotation[,1:layer]
  projected <- as.matrix(data) %*% v %*% t(v) #reconstructed data
  norm_data <- data - projected #normalized data
  dimnames(norm_data) <- dimnames(data)
  dimnames(projected) <- dimnames(data)
  
  print("Normalized data")
  #create similarity network from normalized data
  norm_data_t <- data.frame(t(norm_data)) #transpose data - genes as columns
  sim_net <- cor(norm_data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between normalized-network and un-normalized-network
  norm_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_norm <-c(cor_list_data_norm,cor(data_flat, norm_data_flat, method="pearson"))
  
  print("Projected data")
  #create similarity network from reconstructed data
  projected_t <- data.frame(t(projected)) #transpose data - genes as columns
  sim_net <- cor(projected_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between reconstructed-network and un-normalized-network
  proj_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_proj <-c(cor_list_data_proj,cor(data_flat, proj_data_flat, method="pearson"))
  
}

#plot network PCC values
cor_all <- data.frame(cor_list_data_norm)
names(cor_all)<-c('s')
cor_all['l'] <- cor_list_data_proj
cor_all['pc'] <- layer_list

cols1<- brewer.pal(n = 9, "GnBu")[6:9]

create_network_correlation_scatter_plot(cor_all, cols1, output_folder, "pca_pcc_real_nw_vs_norm_nw_10_scatter.pdf")
#########################################################################################################################
#Create network-correlation scatter plots for RPCA-normalized and RPCA-reconstructed data

print("RPCA")
cor_list_data_norm <- c()
cor_list_data_proj <- c()

#Set output directory
output_folder <- "../output/rpca/"

#Weights to regulate RPCA lambda hyperparameter values in creating normalized data sets
layer_list<-c(.7,.8,.9,1,1.1,1.2,1.3)
for(layer in layer_list)
{
	print(layer)
  #generate normalized and reconstructed data
  # set RPCA lambda hyperparameter (suqre-root of maximum between row-number(genes) and column-number(cell-lines))
  lam  = layer/(sqrt(max(nrow(data),ncol(data)))) # 1/(sqrt(max(# of genes,# of cell lines)))
  print(lam)
  #scale data prior to applying rpca
  data <- scale(data)
  #apply rpca to original data
  pca <- rpca(data, lambda = lam)
  
  #Get RPCA-normalized data (variable S of rpca output contain sparse component or normlaized data)
  norm_data <- pca$S
  #Get RPCA-reconstructed data (variable L of rpca output contain low-rank component or reconstructed data)
  projected <- pca$L
  dimnames(norm_data) <- dimnames(data)
  dimnames(projected) <- dimnames(data)
  
  print("Normalized data")
  #create similarity network from normalized data
  norm_data_t <- data.frame(t(norm_data)) #transpose data - genes as columns
  sim_net <- cor(norm_data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between normalized-network and un-normalized-network
  norm_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_norm <-c(cor_list_data_norm,cor(data_flat, norm_data_flat, method="pearson"))
  
  print("Reconstructed data")
  #create similarity network from reconstructed data
  projected_t <- data.frame(t(projected)) #transpose data - genes as columns
  sim_net <- cor(projected_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
  #Calculate network-network PCC between reconstructed-network and un-normalized-network
  proj_data_flat <- flatten(as.matrix(sim_net))
  cor_list_data_proj <-c(cor_list_data_proj,cor(data_flat, proj_data_flat, method="pearson"))
  
}

#plot network PCC values
cor_all <- data.frame(cor_list_data_norm)
names(cor_all)<-c('s')
cor_all['l'] <- cor_list_data_proj
cor_all['pc'] <- layer_list
cor_all['pc'] <- paste(format(as.numeric(cor_all$pc)*.007),sep = '') # .007*layer 0.007 is the approximate value of 1/(sqrt(max(# genes,# cell lines)))

cols1<- brewer.pal(n = 9, "BuGn")[4:7]

create_network_correlation_scatter_plot(cor_all, cols1, output_folder, "rpca_pcc_real_nw_vs_norm_nw_7_scatter.pdf")
#########################################################################################################################