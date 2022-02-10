#@author:  Arshia Zernab Hassan (hassa418@d.umn.edu)

library("crayon", lib.loc="../../../local/R/library/")

library("pillar", lib.loc="../../../local/R/library/")
library("withr", lib.loc="../../../local/R/library/")
library("labeling", lib.loc="../../../local/R/library/")
library("farver", lib.loc="../../../local/R/library/")
library("digest", lib.loc="../../../local/R/library/")
library("dplyr", lib.loc="../../../local/R/library/")


library("tidyr", lib.loc="../../../local/R/library/")
library("tools", lib.loc="../../../local/R/library/")
library("stringi", lib.loc="../../../local/R/library/")
library("stringr", lib.loc="../../../local/R/library/")

library("rstudioapi", lib.loc="../../../local/R/library/")
library("desc", lib.loc="../../../local/R/library/")
library("withr", lib.loc="../../../local/R/library/")
library("ps", lib.loc="../../../local/R/library/")
library("usethis", lib.loc="../../../local/R/library/")
library("devtools", lib.loc="../../../local/R/library/")

library("RColorBrewer", lib.loc="../../../local/R/library/")
library("ggthemes", lib.loc="../../../local/R/library/")
library("ggplot2", lib.loc="../../../local/R/library/")
library("gplots", lib.loc="../../../local/R/library/")
library("PRROC", lib.loc="../../../local/R/library/")
library("ramify", lib.loc="../../../local/R/library/")


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
	gls_flag=1
	pco_flag=0
	aeo_flag=1
	rpco_flag=1
	olf_flag=1
	flip = 0


	orig_folder <- "../output/onion/"
	orig_file <- 'Complex_AUPRC_depmap_q2_2020_onion_0.txt'

	snf_folder_aeo = "../output/onion/"
	snf_file_aeo = "Complex_AUPRC_aeo.txt"

	snf_folder_pco = "../output/onion/"
	snf_file_pco = "Complex_AUPRC_pco.txt"

	snf_folder_rpco = "../output/onion/"
	snf_file_rpco = "Complex_AUPRC_rpco.txt"

	gls_folder =  "../output/onion/"
	gls_file = paste("Complex_AUPRC_gls.txt",sep = '')

	olf_folder =  "../output/onion/"
	olf_file = paste("Complex_AUPRC_olf.txt",sep = '')

	output_folder = "../output/onion/"
	output_plot = paste("auprc_barplot_onion_",low,"_",up,"_benchmark.pdf",sep='')

	#l <- str_subset(list.dirs(snf_folder, full.names = F, recursive = F), "^snf_k_")


	###################################################################################################

	AUPRC_unnorm <- read.table(file.path(orig_folder,orig_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	AUPRC_unnorm<-AUPRC_unnorm[order(AUPRC_unnorm$ID),]
	method_labels = c()
	method_names = c()

	count = 0
	
	if(olf_flag==1){
	  
	  AUPRC_norm <- read.table(file.path(olf_folder, olf_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  if(count == 0){
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  name = paste('z_olf',sep = '')
	  label = paste('OLF',sep = '')
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	}

	if(gls_flag==1){
	  
	  AUPRC_norm <- read.table(file.path(gls_folder, gls_file), stringsAsFactors=FALSE, sep = "\t", header = T, quote = '')
	  AUPRC_norm<-AUPRC_norm[order(AUPRC_norm$ID),]
	  if(count == 0){
		AUPRC_unnorm <- AUPRC_unnorm[which(AUPRC_unnorm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
		com_df <- AUPRC_unnorm
	  }
	  AUPRC_norm <- AUPRC_norm[which(AUPRC_norm$ID %in% intersect(AUPRC_unnorm$ID, AUPRC_norm$ID)),]
	  name = paste('z_gls',sep = '')
	  label = paste('GLS',sep = '')
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
	  name = 'y_RPCO'
	  
	  com_df[[name]] <- AUPRC_norm$AUPRC[match(AUPRC_norm$ID,com_df$ID)]
	  method_labels <- c(method_labels, label)
	  method_names <- c(method_names, name)
	  count = 1
	}
	
	###################################################################################################

	col_n <- c('method','count','1_len_2_3','2_len_4_5','3_len_6_9','4_len_10_')
	jump_count <- as.data.frame(matrix(0,0,6))
	colnames(jump_count) <- col_n


	if(aeo_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('x_AEO',sep = ''), col_n,up,low))
	}
	if(pco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('y_PCO',sep = ''), col_n,up,low))
	}
	if(rpco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('y_RPCO',sep = ''), col_n,up,low))
	}

	if(gls_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('z_gls',sep = ''), col_n,up,low))
	  
	}
	if(olf_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_up(com_df, paste('z_olf',sep = ''), col_n,up,low))
	  
	}

	if(aeo_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('x_AEO',sep = ''), col_n,up,low))
	}
	if(pco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('y_PCO',sep = ''), col_n,up,low))
	}
	if(rpco_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('y_RPCO',sep = ''), col_n,up,low))
	}
	if(gls_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('z_gls',sep = ''), col_n,up,low))
	}
	if(olf_flag==1){
	  jump_count <- rbind(jump_count, get_complex_subset_down(com_df, paste('z_olf',sep = ''), col_n,up,low))
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
	my_orange = brewer.pal(n = 9, "YlOrBr")[2:5]
	my_red = brewer.pal(n = 9, "Reds")[2:5]
	color_map = c(my_purple,
				  #my_blue, 
				  my_green,
				  my_red,
				  my_orange
				  
	)

	ggplot(com_df_2, aes(fill=interaction(lens, method), y=counts, x=method)) + 
	  geom_bar(position="stack", stat="identity") +
	  geom_hline(yintercept = 0, linetype = 1, size = 1) +
	  xlab("Normalization") +
	  #xlab("") +
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
			axis.title.x = element_text(margin = margin(t = 40, r = 0, b = 0, l = 0)),
			axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
			axis.text = element_text(family=font_family,size = 30,colour = "black"),
			legend.title = element_text(family=font_family,size = 40,colour = "black"),
			legend.text = element_text(family=font_family,size = 40,colour = "black"),
			panel.border = element_blank(), panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.ticks.length=unit(.5, "cm"),
			#axis.text.x = element_blank(),
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

	ggsave(file.path(output_folder , output_plot),height = 9, width = 9)
}

t_list=c(.1,.3,.5,.7)
for(t in t_list)
{
	high <- t
	low <- t
	generate_auprc_barplot(high,low)
	}
rm(list=ls())