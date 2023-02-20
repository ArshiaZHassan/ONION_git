library("stringi", lib.loc="../../local/R/library/")
library("stringr", lib.loc="../../local/R/library/")
library("rstudioapi", lib.loc="../../local/R/library/")
library("desc", lib.loc="../../local/R/library/")
library("withr", lib.loc="../../local/R/library/")
library("ps", lib.loc="../../local/R/library/")
library("usethis", lib.loc="../../local/R/library/")
library("devtools", lib.loc="../../local/R/library/")
library("RColorBrewer", lib.loc="../../local/R/library/")
library("rpca", lib.loc="../../local/R/library/")
library("ramify", lib.loc="../../local/R/library/")

generate_corum_similarity_plot<-function(corum_sim_list, output_folder, output_file_name_noext, curve_names, cols, title)
{
  
  profile_similarity_plot <- paste(output_file_name_noext, 'CORUM_complex_PR_',  sep = "")
  labs <- c('TP', 'Precision')
  PlotPRSimilarity(corum_sim_list, outfile.name = profile_similarity_plot, 
                   outfile.type = 'pdf', fig.labs = labs, fig.title = title,
                   subsample = T,
                   legend.names = curve_names, 
                   legend.color = cols, save.figure = TRUE)
    PlotPRSimilarity(corum_sim_list, outfile.name = profile_similarity_plot, 
                   outfile.type = 'png', fig.labs = labs, fig.title = title,
                   subsample = T,
                   legend.names = curve_names, 
                   legend.color = cols, save.figure = TRUE)
  files <- c( 
    paste(profile_similarity_plot, '.pdf', sep = ""), paste(profile_similarity_plot, '.png', sep = "")
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    file.remove(f)
  }
}

generate_stepwise_contribution_data<-function(complex_info,data_complex_info, output_folder, output_file_name_noext)
{
  pairs_in <- data.frame(true = complex_info$true, 
                         predicted = complex_info$predicted, 
                         ID =  as.character(complex_info$ID), 
                         stringsAsFactors = FALSE)
  stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = complex_info$predicted, 
                                               value.true = complex_info$true, 
                                               x.axis = 'TP', y.axis = 'precision')
  
  precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.1, max(stepwise_contrib$y), 0.025)) #.1 0.025
  precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
  
  stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex_info)
  stepwise_contrib_file <- paste('Stepwise_cont_',output_file_name_noext,'.tsv', sep = "")
  write.table(stepwise_contrib, stepwise_contrib_file, sep="\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #cols<-brewer.pal(10, "Paired")
  stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
  structure_plot_file <- paste('Cont_struct_', output_file_name_noext, sep = "")
  PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, #ccol =cols,
                            min.pairs = 10, outfile.name = structure_plot_file, outfile.type = 'pdf',
                            save.figure = TRUE)
  files <- c( 
    stepwise_contrib_file, paste(structure_plot_file, '.pdf', sep = "")
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    file.remove(f)
  }
}

start_time <- Sys.time()
print("Loading data...")
#path to flex
flex_folder <- "../FLEX_R-master"
#load flex package
load_all(flex_folder)
#generate CORUM standard 
corum <- MakeCoAnnotationFromGeneSymbols(data_standard = data_complex, overlap_length = 1,file_location = "corum")

input_file_name <- "./data/5k_pbmc_v3_filtered_subset_7_norm.tsv"
output_folder <- "./output/5k_pbmc_v3_filtered_subset_7/"
if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
file_ext <- "5k_pbmc_v3_filtered_subset_7_"

#load ribosomal gene list
ribo_gene_list_file_genes <- "./data/ribosomal_gene_list.txt"
ribo_data_genes <- read.table(ribo_gene_list_file_genes,stringsAsFactors = F)

print("Evaluating original data...")
#load data
data <- read.table(input_file_name,sep='\t',header = T, row.names=1, stringsAsFactors = F )#genes as rows

#create similarity network from normalized data
data_t <- data.frame(t(data)) #transpose data - make genes as columns
sim_net <- cor(data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
sim_net[is.na(sim_net)] <- 0 #NA imputation with 0

print("Generating corum-complex standard evaluation data...")
#generate corum-complex standard evaluation results for original data
complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net)
corum_sim <- list(list(true = complex$true, predicted = complex$predicted))

layer = 0
stepwise_contribution_file_name_noext <- paste(file_ext,layer,sep = '')
generate_stepwise_contribution_data(complex, data_complex, output_folder, stepwise_contribution_file_name_noext)

print("Generating corum-complex standard evaluation data no ribo. ...")
#remove ribosomal gene pairs
complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(ribo_data_genes$V1), replace = F)
corum_sim_no_ribo <- list(list(true = complex$true, predicted = complex$predicted))

#curve name and color setting for original data PR curve
curve_labels <- c('Un-normalized')
col_set <- c()
cols0<-brewer.pal(4, "Dark2")
col_set <- c(cols0[4])

#color setting for PR curves
cols1<- brewer.pal(n = 9, "BuGn")[2:8] #rpca
col_set <- c(col_set, cols1)

layer_list <-c(.7,.8,.9,1,1.1,1.2,1.3)
s_list <- list()
for(layer in layer_list)
{
	print(layer)
	#generate normalized data
	lam  = layer/(sqrt(max(nrow(data),ncol(data))))
	print(lam)
	#data is already scaled 
	data <- scale(data)
	pca <- rpca(as.matrix(data), lambda = lam)
	
	print('RPCA applied.')
	norm_data <- pca$S
	dimnames(norm_data) <- dimnames(data)
	print('Components extracted.')
	#save normalized data
	write.table(norm_data, file.path(output_folder, paste('normalized_rpca_', layer,".tsv", sep = "")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
	
	print("Normalized data")
	print("Generating corum-complex standard evaluation data...")
	#create similarity network from normalized data
	norm_data_t <- data.frame(t(norm_data)) #transpose data - make genes as columns
	sim_net <- cor(norm_data_t, use = "all.obs",method="pearson") #use pearson correlation to create similarity network
	sim_net[is.na(sim_net)] <- 0 #NA imputation with 0
	
	#generate corum-complex standard evaluation results
	complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net)	
	corum_sim <- append(corum_sim, list(list(true = complex$true, predicted = complex$predicted)))
	
	#remove ribosomal gene pairs
	complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(ribo_data_genes$V1), replace = F)
	corum_sim_no_ribo <- append(corum_sim_no_ribo, list(list(true = complex$true, predicted = complex$predicted)))	
		
	val = as.double(layer)*.009 #0.00897664621#0.00546048753
	label = val#bquote(paste(lambda,'=',.(val),sep = ''))
	curve_labels <- c(curve_labels, label)
	
	#add similarity network to a list
	s_list <- append(s_list,list(flatten(as.matrix(sim_net))))
	
}
print("RPCA Normalization complete.")

print("Integrating with maximum.")
#Integrate similarity networks
pred_mat<-matrix(unlist(s_list),c(length(s_list),length(s_list[[1]])),byrow=T)
pred_max<- apply(pred_mat,2,max)
max_sim_net =matrix(pred_max, c(dim(sim_net)[1], dim(sim_net)[2]),byrow=T)
sim_net <- as.data.frame(max_sim_net)
row.names(sim_net)<-row.names(data)
colnames(sim_net)<-row.names(data)

print("Generating corum-complex standard evaluation data...")
#generate corum-complex standard evaluation results from integrated network
complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, sim_net)
corum_sim <- append(corum_sim, list(list(true = complex$true, predicted = complex$predicted)))
#remove ribosomal gene pairs
complex <- getSubsetOfCoAnnRemovePairs(corum, data.frame(complex), list(ribo_data_genes$V1), replace = F)
corum_sim_no_ribo <- append(corum_sim_no_ribo, list(list(true = complex$true, predicted = complex$predicted)))

cols2<- brewer.pal(n = 9, "Set2")[4] #rpca
col_set <- c(col_set, cols2)
curve_labels <- c(curve_labels, 'RPCO')

print("Generating PR curves.")
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "All gene pairs" )
generate_corum_similarity_plot(corum_sim_no_ribo, output_folder, paste(file_ext,'no_ribo_fp_pairs',sep=''), curve_labels, col_set, "No ribo. gene pairs" )

end_time <- Sys.time()
print(end_time - start_time)

print("Memory profile: ")
print(memory.profile())

print("Object-wise Memory profile: ")
sapply(ls(), function(x) {print(object.size(get(x)),standard = "legacy", units="Gb") })