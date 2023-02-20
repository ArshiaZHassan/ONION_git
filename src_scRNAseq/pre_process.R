library(dplyr)
library(Seurat)

pbmc.data <- Read10X(data.dir = "../data/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc5k") #create Seurat object

#extract genes for which the number of cells with non-zero values is smaller than or equal to 50
data_matrix = pbmc[["RNA"]]@data 
data_matrix_df_ = as.data.frame(as.matrix(data_matrix))
data_matrix_df_ = data_matrix_df_[ rowSums(data_matrix_df_ > 0) > 50, ] 
gene_list = row.names(data_matrix_df_)

#Filter data
subset.matrix <- pbmc[["RNA"]]@counts[gene_list, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
pbmc2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
pbmc2[["percent.mt"]] <- PercentageFeatureSet(pbmc2, pattern = "^MT-")
pbmc2 <- subset(pbmc2, subset = nFeature_RNA > 100 & nFeature_RNA < 4500 & percent.mt < 7) #10

#normalize and scale
pbmc2 <- NormalizeData(pbmc2)
data_matrix = pbmc2[["RNA"]]@data
data_matrix_df_ = as.data.frame(as.matrix(data_matrix))
write.table(data_matrix_df_, "../output/5k_pbmc_v3_filtered_subset_7_norm.tsv",quote = F,sep='\t')

all.genes <- rownames(pbmc2)
pbmc2 <- ScaleData(pbmc2, features = all.genes)
data_matrix_2 = pbmc2[["RNA"]]@scale.data
data_matrix_df = as.data.frame(as.matrix(data_matrix_2))

write.table(data_matrix_df, "../output/5k_pbmc_v3_filtered_subset_7_norm_scaled.tsv",quote = F,sep='\t')
