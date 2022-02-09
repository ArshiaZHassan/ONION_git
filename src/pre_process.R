#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)
#Pre-process DepMap 2020 Q2 Achilles_gene_effect data (rows:cell lines, columns:genes) 
#Input: DepMap 2020 Q2 Achilles_gene_effect data https://ndownloader.figshare.com/files/22629068
#Output: Pre-processed DepMap 20q2 data (rows:genes, columns:cell lines)
# - replace NA values with gene-level mean 
# - extract gene names from row labels


# Function:
    # split_gene_name
# Arguments:
    # x : gene label. type: str (required)
# Return: 
    # gene name (str)
# Description:
    # extract gene name from gene label

split_gene_name<-function(x)
{
  x1<-strsplit(x,"\\.")[[1]][1] #split x with delimiter \. and retrieve first sub-string
  return( x1)
}

#Load DepMap 20q2 Achilles_gene_effect data (rows:cell lines, columns:genes)
input_file_name <- "../data/Achilles_gene_effect.csv"
Achilles_gene_effect <- read.csv(file = input_file_name,row.names = 1, header = TRUE)

#Extract gene names from gene-labels(column names) and set column names to gene names
cols<-c(colnames(Achilles_gene_effect))
cols2<-unlist(lapply(cols,split_gene_name))
colnames(Achilles_gene_effect) <- c(cols2)

#transpose matrix (rows:genes, columns:cell lines)
Achilles_gene_effect_t <- data.frame(t(Achilles_gene_effect))

#replace NA values with row-wise(gene-wise) mean values
k <- which(is.na(Achilles_gene_effect_t), arr.ind=TRUE)
Achilles_gene_effect_t[k] <- rowMeans(Achilles_gene_effect_t, na.rm=TRUE)[k[,1]]# replace na with row average

#save pre-processed DepMap 20q2 gene_effect data
write.table(Achilles_gene_effect_t,file="../data/depmap_q2_2020_nona_mean.tsv",sep='\t',row.names = TRUE,col.names = TRUE, quote=FALSE)

