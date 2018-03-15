
library(stringr)
library(ggplot2)
library(dplyr)
require(reshape)
library(cowplot)
library(RWebLogo)
library(VennDiagram)

#############################################
#create a matrix with only the significant GoF positions
gof_matrix_melt_only_GOFpos <- gof_matrix_melt %>%
  filter(X2 == 34 |
           X2 == 35 |
           X2 == 64 |
           X2 == 68 |
           X2 == 69 |
           X2 == 103 |
           X2 == 134 |
           X2 == 135)

#reassign their X axis position
gof_matrix_melt_only_GOFpos$mutposnum <- 0
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==34)] <- 1
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==35)] <- 2
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==64)] <- 3
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==68)] <- 4
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==69)] <- 5
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==103)] <- 6
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==134)] <- 7
gof_matrix_melt_only_GOFpos$mutposnum[which(gof_matrix_melt_only_GOFpos$X2==135)] <- 8

gof_matrix_melt_only_GOFpos$aanum <- 0
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="A")] <- 12
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="C")] <- 10
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="D")] <- 5
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="E")] <- 4
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="F")] <- 19
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="G")] <- 11
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="H")] <- 3
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="I")] <- 15
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="K")] <- 1
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="L")] <- 14
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="M")] <- 16
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="N")] <- 6
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="P")] <- 17
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="Q")] <- 7
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="R")] <- 2
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="S")] <- 9
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="T")] <- 8
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="V")] <- 13
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="W")] <- 20
gof_matrix_melt_only_GOFpos$aanum[which(gof_matrix_melt_only_GOFpos$X1=="Y")] <- 18

gof_matrix_melt_only_GOFpos_wnum <- gof_matrix_melt_only_GOFpos %>%
  inner_join(gof_matrix_num_melt,by=c("X1","X2")) %>%
  dplyr::rename(mutnum=value.y,value=value.x)

########################################
#plot the data from all mutants:
ggplot(gof_matrix_melt_only_GOFpos_wnum, 
       aes(x=mutposnum, y=aanum,
           fill=value,
           label=mutnum)) +
  geom_tile() +
  geom_text() +
  labs(x = "Position (aa)",
       y ="Amino acid",color="") +
  scale_fill_gradient(low = "blue", 
                      high = "red",
                      name="Fitness",
                      na.value="grey",
                      limit = c(0,5)) +
  theme_minimal()+
  scale_x_continuous(name="Position (aa)", 
                     breaks=c(1,2,3,4,5,6,7,8),
                     labels=c("34","35","64","68","69","103","134","135"))+
  scale_y_continuous(name="Amino acid", 
                     breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                     labels=c("K","R","H","E","D","N","Q","T","S","C","G","A","V","L","I","M","P","Y","F","W"))
ggsave(file=paste("plots/GOF/GOF_fit_nummut.pdf",sep=""),
       width=4.5, height=4.5, units="in")
ggsave(file=paste("plots/GOF/GOF_fit_nummut.png",sep=""),
       width=4.5, height=4.5, units="in")

############################
# this shows the same GoF data as above plot determined using the BMS data
BMS_matrix_all_melt_GOFonly <- BMS_matrix_all_melt %>%
  filter(X1 != "X") %>%
  filter(X2 == 34 | X2 == 35 | X2 == 64 | X2 == 68 | X2 == 69 | X2 == 103 | X2 == 134 | X2 == 135)

BMS_matrix_all_melt_GOFonly$mutposnum <- 0
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==34)] <- 1
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==35)] <- 2
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==64)] <- 3
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==68)] <- 4
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==69)] <- 5
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==103)] <- 6
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==134)] <- 7
BMS_matrix_all_melt_GOFonly$mutposnum[which(BMS_matrix_all_melt_GOFonly$X2==135)] <- 8

names(BMS_matrix_all_melt_GOFonly)

BMS_matrix_all_melt_GOFonly_wnum <- BMS_matrix_all_melt_GOFonly %>%
  inner_join(BMS_matrix_all_num_melt,by=c("X1","X2")) %>%
  dplyr::rename(mutnum=value.y,value=value.x)

ggplot(BMS_matrix_all_melt_GOFonly_wnum, 
       aes(x=mutposnum, y=aanum,
           fill=value,
           label=mutnum)) +
  geom_tile() +
  geom_text() +
  labs(x = "Position (aa)",
       y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", high = "red", mid="gold",
                       name="Fitness",na.value="grey", 
                       limit = c(-5,1.1*max(BMS_matrix_all_melt_GOFonly$value))) +
  theme_minimal()+
  scale_x_continuous(name="Position (aa)", 
                     breaks=c(1,2,3,4,5,6,7,8),
                     labels=c("34","35","64","68","69","103","134","135"))+
  scale_y_continuous(name="Amino acid", 
                     breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                     labels=c("K","R","H","E","D","N","Q","T","S","C","G","A","V","L","I","M","P","Y","F","W"))
ggsave(file=paste("plots/GOF/BMS_GOF_fit_nummut.pdf",sep=""),
       width=4.5, height=4.5, units="in")
ggsave(file=paste("plots/GOF/BMS_GOF_fit_nummut.png",sep=""),
       width=4.5, height=4.5, units="in")


#######################################################
# web logo plots for GOF positions vs complementing

GOF_dropout_IDs_unique <- unique(dropout_mutants_GOF$ID)
ecoli_map_dropouts_GOF <- read.csv(file=paste("alignments/MSAdropouts/NP_418091.csv",sep=""),head=TRUE,sep=",")
GOF_MSAdropout_positions <- ecoli_map_dropouts_GOF$msa_aanum[which(ecoli_map_dropouts_GOF$orth_aanum %in% GOF_fitness_collapsed_by_pos_2sigma$position)]

mutants_to_plot_index_of_mut_to_drop = numeric()
for (i in 1:length(GOF_dropout_IDs_unique)){
  mutant_current_temp <- GOF_dropout_IDs_unique[i]
  if(!file.exists(paste("alignments/MSAdropouts/",mutant_current_temp,".csv",sep=""))){
    print(mutant_current_temp)
    mutants_to_plot_index_of_mut_to_drop[length(mutants_to_plot_index_of_mut_to_drop)+1] <- which(GOF_dropout_IDs_unique==mutant_current_temp)
  }
}

subchar <- function(string, pos, char) { 
  for(i in pos) { 
    string <- gsub(paste("^(.{", i-1, "}).", sep=""), 
                   paste("\\1",char,sep=""), 
                   string) 
  } 
  string 
} 

GOF_weblogo_dropouts = rep("--------",length(GOF_dropout_IDs_unique))

#loop through all dropouts and find residue at GOF positions
for (i in 1:length(GOF_dropout_IDs_unique)){
  mutant_current <- GOF_dropout_IDs_unique[i]
  
  #get the MSA mapping
  mutant_map <- read.csv(file=paste("alignments/MSAdropouts/",mutant_current,".csv",sep=""),head=TRUE,sep=",")
  
  #load the full seq and add an M at start
  GOF_seq_temp <- paste("M",
                        as.character(orcollapse3$seq[which(orcollapse3$ID==mutant_current)]),
                        sep="")[1]
  
  for (j in 1:length(GOF_MSAdropout_positions)){
    GOF_temp_ortho_res_num <- mutant_map$orth_aanum[which(mutant_map$msa_aanum==GOF_MSAdropout_positions[j])]
    if (length(GOF_temp_ortho_res_num) != 0){
      GOF_weblogo_dropouts[i] <- subchar(GOF_weblogo_dropouts[i],
                                         j,
                                         substr(GOF_seq_temp,
                                                GOF_temp_ortho_res_num,
                                                GOF_temp_ortho_res_num))
    }
  }
}

ecoli_map_complement <- read.csv(file=paste("alignments/MSA_BMS_1/NP_418091.csv",sep=""),head=TRUE,sep=",")
GOF_MSAcomplement_positions <- ecoli_map_complement$msa_aanum[which(ecoli_map_complement$orth_aanum %in% GOF_fitness_collapsed_by_pos_2sigma$position)]
GOF_weblogo_complement = rep("--------",length(BMS_ortho_ID_list$ID))

#loop through all dropouts and find residue at GOF positions
for (i in 1:length(BMS_ortho_ID_list$ID)){
  mutant_current <- BMS_ortho_ID_list$ID[i]
  
  #get the MSA mapping
  mutant_map <- read.csv(file=paste("alignments/MSA_BMS_1/",mutant_current,".csv",sep=""),head=TRUE,sep=",")
  
  #load the full seq and add an M at start
  GOF_seq_temp <- paste("M",
                        as.character(orcollapse3$seq[which(orcollapse3$ID==mutant_current)]),
                        sep="")[1]
  
  for (j in 1:length(GOF_MSAcomplement_positions)){
    GOF_temp_ortho_res_num <- mutant_map$orth_aanum[which(mutant_map$msa_aanum==GOF_MSAcomplement_positions[j])]
    if (length(GOF_temp_ortho_res_num) != 0){
      GOF_weblogo_complement[i] <- subchar(GOF_weblogo_complement[i],
                                           j,
                                           substr(GOF_seq_temp,
                                                  GOF_temp_ortho_res_num,
                                                  GOF_temp_ortho_res_num))
    }
  }
}

# Now for an example using an alignment as an R character vector

# WebLogo GOF
weblogo(seqs=GOF_weblogo_dropouts,
        errorbars=FALSE,
        file.out="plots/GOF/GOF_dropout_weblogo.pdf",
        title='Dropout orthologs',
        fineprint='',
        label='',
        units="probability",
        xlabel='Residue position',
        annotate=GOF_fitness_collapsed_by_pos_2sigma$position)

# WebLogo complement
weblogo(seqs=GOF_weblogo_complement,
        errorbars=FALSE,
        file.out="plots/GOF/GOF_complement_weblogo.pdf",
        title='Complementing orthologs',
        fineprint='',
        label='',
        units="probability",
        xlabel='Residue position',
        annotate=GOF_fitness_collapsed_by_pos_2sigma$position)
