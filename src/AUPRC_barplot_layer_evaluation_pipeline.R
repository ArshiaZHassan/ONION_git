#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)

library(crayon)
library(pillar)
library(withr)
library(labeling)
library(farver)
library(digest)

library(dplyr)
library(tidyr)
library(tools)
library(stringi)
library(stringr)

library(rstudioapi)
library(desc)
library(withr)
library(ps)
library(usethis)
library(devtools)

library(RColorBrewer)
library(ggthemes)
library(ggplot2)

library(gplots)
library(PRROC)
library(ramify)

get_complex_subset_up <- function(com_df, meth, col_n, up_t, low_t)
{
  temp_row <-as.data.frame(matrix(0,0,6))
  colnames(temp_row) <- col_n
  temp_row[1,1] <- meth
  com_df_sub<-com_df[which(com_df$AUPRC<=low_t & com_df[[meth]]>up_t),]
  temp_row[1,2] <- nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>1 & com_df$Length<4) & (com_df$AUPRC<=low_t &
                                                                    ((com_df[[meth]]>up_t )
                                                                    ))
  ),]
  temp_row[1,3] <- nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>3 & com_df$Length<6) & (com_df$AUPRC<=low_t &
                                                                    ((com_df[[meth]]>up_t )
                                                                    ))
  ),]
  temp_row[1,4] <- nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>5 & com_df$Length<10) & (com_df$AUPRC<=low_t &
                                                                     ((com_df[[meth]]>up_t )
                                                                     ))
  ),]
  temp_row[1,5] <- nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>9 & com_df$Length<10000) & (com_df$AUPRC<=low_t &
                                                                        ((com_df[[meth]]>up_t )
                                                                        ))
  ),]
  temp_row[1,6] <- nrow(com_df_sub)
  #temp_row <- data.frame(paste0('pc_',pc), nrow(com_df_sub))
  #colnames(temp_row) <- c('method','count')
  return(temp_row)
}

get_complex_subset_down <- function(com_df, meth, col_n, up_t, low_t)
{
  temp_row <-as.data.frame(matrix(0,0,6))
  colnames(temp_row) <- col_n
  temp_row[1,1] <- meth
  
  com_df_sub<-com_df[which(com_df$AUPRC>up_t & com_df[[meth]]<=low_t),]
  temp_row[1,2] <- -nrow(com_df_sub)
  
  com_df_sub<-com_df[which((com_df$Length>1 & com_df$Length<4) & (com_df$AUPRC>up_t &
                                                                    ((com_df[[meth]]<=low_t )
                                                                    ))
  ),]
  temp_row[1,3] <- -nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>3 & com_df$Length<6) & (com_df$AUPRC>up_t &
                                                                    ((com_df[[meth]]<=low_t )
                                                                    ))
  ),]
  temp_row[1,4] <- -nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>5 & com_df$Length<10) & (com_df$AUPRC>up_t &
                                                                     ((com_df[[meth]]<=low_t )
                                                                     ))
  ),]
  temp_row[1,5] <- -nrow(com_df_sub)
  com_df_sub<-com_df[which((com_df$Length>9 & com_df$Length<10000) & (com_df$AUPRC>up_t &
                                                                        ((com_df[[meth]]<=low_t )
                                                                        ))
  ),]
  temp_row[1,6] <- -nrow(com_df_sub)
  #temp_row <- data.frame(paste0('pc_',pc), nrow(com_df_sub))
  #colnames(temp_row) <- c('method','count')
  return(temp_row)
}

###################################################################################################
generate_auprc_barplot <- function(up,low)
{
	
	pco_flag=1
	aeo_flag=1
	rpco_flag=1

	orig_folder <- "../output/onion/"
	orig_file <- 'Complex_AUPRC_depmap_q2_2020_onion_0.txt'

	pc_folder_ae =  "../output/ae/"
	pc_folder_rpca =  "../output/rpca/"
	pc_folder_pca =  "../output/pca/"

	snf_folder_aeo = "../output/onion/"
	snf_file_aeo = "Complex_AUPRC_aeo.txt"

	snf_folder_pco =  "../output/onion/"
	snf_file_pco = "Complex_AUPRC_pco.txt"

	snf_folder_rpco =  "../output/onion/"
	snf_file_rpco = "Complex_AUPRC_rpco.txt"

	output_folder_2 =  "../output/onion/"
	output_plot_name = paste("auprc_bar_plot_",low,"_",up,"_layers_onion",sep='')

	output_plot_2 = paste(output_plot_name,".pdf",sep='')

	pc_list_rpca<-c('0.7','0.8','0.9','1','1.1','1.2','1.3')
	pc_list_pca<-c(1,3,5,7,9,11,13,15,17,19)
	pc_list_ae<-c(1:5,10)
	#pc_list<-c()
	###################################################################################################

	AUPRC_unnorm <- read.table(file.path(orig_folder,orig_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	AUPRC_unnorm<-AUPRC_unnorm[order(AUPRC_unnorm$ID),]
	method_labels = c()
	method_names = c()

	count = 0

	for(pc in pc_list_rpca)
	{
	  
	  AUPRC_norm <- read.table(file.path(pc_folder_rpca, paste("Complex_AUPRC_normalized_depmap_q2_2020_rpca_",pc,".txt",sep = '')), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  #label = paste(pc,sep = '')
	  val = as.double(pc)*.007
	  label = bquote(paste(lambda,'=',.(val),sep = ''))
	  if(pc>9)
	  {
		name = paste('rpc_d_',pc,sep = '')
	  }else
	  {
		name = paste('rpc_',pc,sep = '')
	  }
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1

	}

	for(pc in pc_list_pca)
	{
	  
	  AUPRC_norm <- read.table(file.path(pc_folder_pca, paste("Complex_AUPRC_normalized_depmap_q2_2020_pca_",pc,".txt",sep = '')), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  #label = paste(pc,sep = '')
	  label = paste('PC=',pc,sep = '')
	  if(pc>9)
	  {
		name = paste('pc_d_',pc,sep = '')
	  }else
	  {
		name = paste('pc_',pc,sep = '')
	  }
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	  

	}

	for(pc in pc_list_ae)
	{
	  
	  AUPRC_norm <- read.table(file.path(pc_folder_ae, paste("Complex_AUPRC_normalized_depmap_q2_2020_ae_",pc,".txt",sep = '')), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  #label = paste(pc,sep = '')
	  label = paste('LS=',pc,sep = '')
	  if(pc>9)
	  {
		name = paste('ae_d_',pc,sep = '')
	  }else
	  {
		name = paste('ae_',pc,sep = '')
	  }
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1

	}


	if(aeo_flag==1){
	  AUPRC_norm <- read.table(file.path(snf_folder_aeo, snf_file_aeo), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  label = 'AEO'
	  name = 'x_AEO'
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	}

	if(pco_flag==1){
	  AUPRC_norm <- read.table(file.path(snf_folder_pco, snf_file_pco), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  label = 'PCO'
	  name = 'y_PCO'
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	}

	if(rpco_flag==1){
	  AUPRC_norm <- read.table(file.path(snf_folder_rpco, snf_file_rpco), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  #adjust legnth of un-normalized list
	  if(count == 0)
	  {
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		#data frame to load all AUPRC data
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  
	  label = 'RPCO'
	  name = 'z_RPCO'
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	}


	###################################################################################################

	col_n <- c('method','count','1_len_2_3','2_len_4_5','3_len_6_9','4_len_10_')
	jump_count <- as.data.frame(matrix(0,0,6))
	colnames(jump_count) <- col_n



	for(pc in pc_list_rpca){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('rpc_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('rpc_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	for(pc in pc_list_pca){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('pc_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('pc_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	for(pc in pc_list_ae){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('ae_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('ae_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	if(aeo_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('x_AEO',sep = ''), col_n,up,low))
	}
	if(pco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('y_PCO',sep = ''), col_n,up,low))
	}
	if(rpco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('z_RPCO',sep = ''), col_n,up,low))
	}



	for(pc in pc_list_rpca){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('rpc_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('rpc_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	for(pc in pc_list_pca){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('pc_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('pc_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	for(pc in pc_list_ae){
	  if(pc>9)
	  {
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('ae_d_',pc,sep = ''), col_n,up,low))
		
	  }else{
		jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('ae_',pc,sep = ''), col_n,up,low))
		
	  }
	}

	if(aeo_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('x_AEO',sep = ''), col_n,up,low))
	}
	if(pco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('y_PCO',sep = ''), col_n,up,low))
	}
	if(rpco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('z_RPCO',sep = ''), col_n,up,low))
	}

	gc()

	len_labels<-col_n[3:6]
	com_df_2 <- gather(jump_count,
					   key = "lens",
					   value = "counts",
					   c(all_of(len_labels))
					   
	)

	###################################################################################################
	font_family = "ArialMT"
	my_blue = brewer.pal(n = 9, "GnBu")[6:9]
	my_purple = brewer.pal(n = 9, "BuPu")[6:9]
	my_green = brewer.pal(n = 9, "BuGn")[4:7]
	#my_orange = brewer.pal(n = 9, "YlOrBr")[3:6]

	color_map = c(rep(my_purple, times = 6), 
				  rep(my_blue, times = 10),
				  rep(my_green, times = 7),
				  my_purple,
				  my_blue, 				  
				  my_green
				  )
	text_y = -20#.7 -15#.1 -145 #.5 -20 #.3 -40
	line_y = -text_y + 1
	ggplot(com_df_2, aes(fill=interaction(lens, method), y=counts, x=method)) + 
	  geom_bar(position="stack", stat="identity") +
	  geom_hline(yintercept = 0, linetype = 1, size = 1) +
	  geom_vline(xintercept = 23.5, linetype = 2, size = 1) +
	  xlab("") +
	  ylab("# of complexes") +
	  #ylab("") +
	  scale_fill_manual(
		values= color_map,
		#breaks = len_labels,
		#labels=c("2-3", "4-5", '6-9','10-'),
		aesthetics = "fill"
	  ) +
	  theme_bw() +
	  labs(fill = "Complex \n Size") +
	  theme(legend.key.size = unit(2, 'cm'), 
			title = element_text(family=font_family,size = 20,colour = "black"),
			axis.title = element_text(family=font_family,size = 40,colour = "black"), 
			axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
			axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
			axis.text = element_text(family=font_family,size = 30,colour = "black"),
			legend.title = element_text(family=font_family,size = 30,colour = "black"),
			legend.text = element_text(family=font_family,size = 30,colour = "black"),
			panel.border = element_blank(), panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.length=unit(.5, "cm"),
			#axis.text.x = element_blank(),
			#axis.text.y = element_blank(),
			#axis.ticks.x = element_blank(),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
			#panel.grid.minor = element_line(colour = "gray",size = 1),
			axis.line.y = element_line(colour = "black",size = 1)
			, axis.line.x = element_line(colour = "black",size = 1)
			, legend.position = "none"
			
	  ) +
	  labs(subtitle  =paste("t=",low))+
	  scale_x_discrete(
		breaks=method_names,
		labels=method_labels 
	  )

	###################################################################################################
	ggsave(file.path(output_folder_2 , output_plot_2),height = 9 , width = 9*1.5 )
}


t_list=c(.1,.3,.5,.7)
for(t in t_list)
{
	high <- t
	low <- t
	generate_auprc_barplot(high,low)
	}
rm(list=ls())
