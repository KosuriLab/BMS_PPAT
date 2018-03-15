# Broad Mutational Scanning (BMS) of PPAT
# Calin Plesa
#

#install any missing packages
list.of.packages <- c("stringr","ggplot2", "dplyr","reshape", "tidyr","cowplot", "RWebLogo", "VennDiagram","ggExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")

library(stringr)
library(ggExtra)
library(ggplot2)
library(dplyr)
require(reshape)
library(tidyr)
library(cowplot)

#set current directory to current file directory
setwd("/Users/calin/GoogleDrive/postdoc/PPAT/BMScode")

#load data file with:
#perfects - info on all homologs
#mutants - info on all homolog mutants
load("PPATdata.RData")

#load EVmutations data (generated using E. coli PPAT on evmutation.org)
EVmutPPAT <- read.csv(file="ev_couplings/PPAT_EVmut.csv",head=TRUE,sep=",")
EVmutPPAT$subs <- as.character(EVmutPPAT$subs)
EVmutPPAT$pos <- as.numeric(EVmutPPAT$pos)

#directory containing the parsed MSA mappings for each homolog
BMS_MSA_directory = "alignments/MSA_BMS_1"

#mapping of amino acids + X (deletion) to numerical values
BMS_aa_list<-data.frame(aa=c('G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T','X'),aanum=c(1:21))
BMS_aa_dim <- length(BMS_aa_list$aa)

#amino acids grouped by class
AAclass <- read.csv(file="AAclass.csv",head=TRUE,sep=",")

#length of E. coli PPAT homolog which collapsing data on
BMS_ref_len <- 159

#minimum number of BCs to use a mutant
BMS_min_BCs = 1

#largest number of mutations to use
#currently limited to 5 due to mutation naming scheme
BMS_max_mutations = 5

#minimum fitness to include the mutation
BMS_min_fitness <- -1

#grab the perfects (homologs) with fitness above minimum fitness
BMS_ortho_ID_list <- perfects %>%
  filter(globalfit14 > BMS_min_fitness) %>%
  dplyr::select(ID)

#save csv file with homolog IDs
#write.table(BMS_ortho_ID_list, file = paste("data/BMS_ortho_ID_list_min_fitness_",as.character(BMS_min_fitness),".csv",sep=""), sep = ",", row.names = F,quote=F,col.names = F)

#make sure all of the perfects are in the MSA
BMS_index_of_mut_to_drop = numeric()
for (i in 1:length(BMS_ortho_ID_list$ID)){
  mutant_current_temp <- BMS_ortho_ID_list$ID[i]
  if(!file.exists(paste(BMS_MSA_directory,"/",mutant_current_temp,".csv",sep=""))){
    print(mutant_current_temp)
    BMS_index_of_mut_to_drop[length(BMS_index_of_mut_to_drop)+1] <- which(BMS_ortho_ID_list$ID==mutant_current_temp)
  }
}
#remove IDs not in MSA
if (length(BMS_index_of_mut_to_drop) > 0){
  BMS_ortho_ID_list <- BMS_ortho_ID_list[-BMS_index_of_mut_to_drop,]
}
rm(BMS_index_of_mut_to_drop,
   mutant_current_temp)

#load ecoli MSA mapping
BMS_ecoli_map <- read.csv(file=paste(BMS_MSA_directory,"/NP_418091.csv",sep=""),head=TRUE,sep=",")

#E.coli PPAT seq:
BMS_ecoli_seq <- "MQKRAIYPGTFDPITNGHIDIVTRATQMFDHVILAIAASPSKKPMFTLEERVALAQQATAHLGNVEVVGFSDLMANFARNQHATVLIRGLRAVADFEYEMQLAHMNRHLMPELESVFLMPSKEWSFISSSLVKEVARHQGDVTHFLPENVHQALMAKLA"

#grab the perfect homologs (no mutations)
BMS_mutants_temp <- mutants %>%
  filter(mutations == 0 &
           numprunedBCs >= BMS_min_BCs) %>%
  dplyr::rename(ID=IDalign) %>%
  semi_join(BMS_ortho_ID_list,by="ID") %>%#filtering join
  ungroup() %>%
  dplyr::rename(IDalign=ID)

#make a new data frame which will keep info
BMS_fitness_map <- data.frame(position=numeric(),
                              aa=character(),
                              mutations=numeric(),
                              fitness=numeric(),
                              posortho=numeric(),
                              ingap=character(),
                              mutID=character(),
                              ID=character())

#loop over all perfects
#what this does is create a table of each amino acid at 
#each aligned position and it's fitness
#THIS IS SLOW
for (i in 1:nrow(BMS_mutants_temp)){
  
  #current homolog ID:
  mutant_current_temp <- BMS_mutants_temp$IDalign[i]
  
  #get the MSA mapping:
  mutant_map_temp <- read.csv(file=paste(BMS_MSA_directory,"/",mutant_current_temp,".csv",sep=""),head=TRUE,sep=",")
  
  #load the full homolog seq and add an M at start
  BMS_seq_temp <- paste("M",as.character(BMS_mutants_temp$seq[i]),sep="")
  
  #get fitness for this homolog
  BMS_fit_temp <- BMS_mutants_temp$globalfit14[i]
  
  #loop over all residues in homolog seq
  for (j in 1:nchar(BMS_seq_temp)){
    
    #find the corresponding residue in E.coli using MSA
    BMS_cons_aanum <- mutant_map_temp$msa_aanum[which(mutant_map_temp$orth_aanum == j)]
    
    #does this map to a non-gap residue in E.coli?
    if (BMS_cons_aanum %in% BMS_ecoli_map$msa_aanum){
      
      #get the corresponding E.coli residue
      e_coli_residue <- BMS_ecoli_map$orth_aanum[which(BMS_ecoli_map$msa_aanum == BMS_cons_aanum)]
      
      #aa at this residue
      BMS_aa_temp <- substr(BMS_seq_temp, j, j)
      
      #add info for this residue to dataframe as a new row
      BMS_fitness_map <- rbind(BMS_fitness_map,
                               data.frame(position=e_coli_residue,
                                          aa=BMS_aa_temp,
                                          mutations=0,
                                          fitness=BMS_fit_temp,
                                          posortho=j,
                                          ingap="No",
                                          mutID=mutant_current_temp,
                                          ID=mutant_current_temp))
    } else {
      #if it's here it maps to a gap
      #these are labeled with position -1 and ingap = "Yes"
      BMS_fitness_map <- rbind(BMS_fitness_map,
                               data.frame(position=-1,
                                          aa=BMS_aa_temp,
                                          mutations=0,
                                          fitness=BMS_fit_temp,
                                          posortho=j,
                                          ingap="Yes",
                                          mutID=mutant_current_temp,
                                          ID=mutant_current_temp))
      
    }
  }
  
}

#save all residues from homologs as a row with position and fitness values
write.table(BMS_fitness_map, file = paste("output/BMS_fitness_map_0.csv",sep=""), sep = ",", row.names = F,quote=F,col.names = T)

#this is where the table is collapsed by position and aa
#for each E.coli position and aa, determine the median, 
#number of points, and std
BMS_fitness_collapsed <- BMS_fitness_map %>%
  filter(position > 0) %>%
  group_by(position, aa) %>%
  summarise(fitval=median(fitness),
            numpoints=n(),
            stdfit=sd(fitness))

#these matrices have the fitness for each aa at each position:
BMS_matrix_perfects = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)
BMS_matrix_perfects_num = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)
BMS_matrix_perfects_sd = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)

#populate matrices with the data
for (i in 1:nrow(BMS_fitness_collapsed)){
  BMS_matrix_perfects[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed$aa[i])),BMS_fitness_collapsed$position[i]] <- as.numeric(BMS_fitness_collapsed$fitval[i])
  BMS_matrix_perfects_num[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed$aa[i])),BMS_fitness_collapsed$position[i]] <- as.numeric(BMS_fitness_collapsed$numpoints[i])
  BMS_matrix_perfects_sd[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed$aa[i])),BMS_fitness_collapsed$position[i]] <- as.numeric(BMS_fitness_collapsed$stdfit[i])
}

#assign aa names and position numbers to matrices
rownames(BMS_matrix_perfects)<-BMS_aa_list$aa
colnames(BMS_matrix_perfects)<-c(1:BMS_ref_len)
rownames(BMS_matrix_perfects_num)<-BMS_aa_list$aa
colnames(BMS_matrix_perfects_num)<-c(1:BMS_ref_len)
rownames(BMS_matrix_perfects_sd)<-BMS_aa_list$aa
colnames(BMS_matrix_perfects_sd)<-c(1:BMS_ref_len)

BMS_matrix_perfects_melt <- melt(BMS_matrix_perfects)
BMS_matrix_perfects_num_melt <- melt(BMS_matrix_perfects_num)
BMS_matrix_perfects_sd_melt <- melt(BMS_matrix_perfects_sd)

#plot the median fitness data using only homolog info:
ggplot(data = BMS_matrix_perfects_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="Fitness",na.value="grey", limit = c(-5,5)) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,155,5))
ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_0_fit.pdf",sep=""))

#plot the number of points at each position:
ggplot(data = BMS_matrix_perfects_num_melt, aes(x=X2, y=X1, fill=log10(value))) +
  geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient(low = "blue", high = "red",name="log(#points)",na.value="grey", limit = c(0,max(BMS_matrix_perfects_num_melt$value))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,155,5))
ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_0_num.pdf",sep=""))

#plot the std data:
ggplot(data = BMS_matrix_perfects_sd_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="std(Fitness)",na.value="grey", limit = c(0,1.1*max(BMS_matrix_perfects_sd_melt$value))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,155,5))
ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_0_sd.pdf",sep=""))


##############################
# add in informtion from the mutants

BMS_fitness_map_all <- BMS_fitness_map

#in this loop only the values for the mutates residues are added in
#not every single aa in the sequence, since the mutations is what
#impacts the change in fitness
#loop over mutants at some distance and collect the fitness values
for (BMS_cur_mut_num in 1:BMS_max_mutations){
  
  #grab the mutants
  BMS_mutants_temp <- mutants %>%
    filter(mutations == BMS_cur_mut_num &
             numprunedBCs >= BMS_min_BCs) %>%
    dplyr::rename(ID=IDalign) %>%
    semi_join(BMS_ortho_ID_list,by="ID") %>%#filtering join
    ungroup() %>%
    dplyr::rename(IDalign=ID)
  
  #initialize new df
  BMS_fitness_map1 <- data.frame(position=numeric(),
                                 aa=character(),
                                 mutations=numeric(),
                                 fitness=numeric(),
                                 posortho=numeric(),
                                 ingap=character(),
                                 mutID=character(),
                                 ID=character())
  
  #loop over mutants
  for (i in 1:nrow(BMS_mutants_temp)){
    
    #mutant base name
    mutant_current_temp <- BMS_mutants_temp$IDalign[i]
    
    #length of the name
    name_size = nchar(paste(mutant_current_temp,"_",sep=""))
    
    #residue mappings
    mutant_map_temp <- read.csv(file=paste(BMS_MSA_directory,"/",mutant_current_temp,".csv",sep=""),head=TRUE,sep=",")
    
    #this mutants fitness
    BMS_fit_temp <- BMS_mutants_temp$globalfit14[i]
    
    #grab the mut name
    mutations_names <- as.character(BMS_mutants_temp$mutID[i])
    
    #grab only the relevant portion of the name
    mutations_names <- substr(mutations_names, name_size+1, nchar(mutations_names))
    
    ## split mutation string at non-digits
    s <- strsplit(mutations_names, "_")
    
    for (mutnum in 1:BMS_cur_mut_num){
      
      #grab the corresponding mutation string
      mutcurr<-s[[1]][mutnum]
      
      #get the position
      mutpos <- as.numeric(str_extract(mutcurr, "[0-9]+"))
      
      #get ending aa
      to_aa <- substr(mutcurr, nchar(mutpos)+2, nchar(mutcurr))
      
      #find the number in the consensus seq
      BMS_cons_aanum <- mutant_map_temp$msa_aanum[which(mutant_map_temp$orth_aanum == mutpos)]
      
      #does this map to a non-gap
      if (BMS_cons_aanum %in% BMS_ecoli_map$msa_aanum){
        
        #the corresponding e.coli residue
        e_coli_residue <- BMS_ecoli_map$orth_aanum[which(BMS_ecoli_map$msa_aanum == BMS_cons_aanum)]
        
        #add this point to the data
        BMS_fitness_map1 <- rbind(BMS_fitness_map1,
                                  data.frame(position=e_coli_residue,
                                             aa=to_aa,
                                             mutations=BMS_cur_mut_num,
                                             fitness=BMS_fit_temp,
                                             posortho=mutpos,
                                             ingap="No",
                                             mutID=as.character(BMS_mutants_temp$mutID[i]),
                                             ID=mutant_current_temp))
        
      } else {
        #if it's here it maps to a gap
        
        #add this point to the data
        BMS_fitness_map1 <- rbind(BMS_fitness_map1,
                                  data.frame(position=-1,
                                             aa=to_aa,
                                             mutations=BMS_cur_mut_num,
                                             fitness=BMS_fit_temp,
                                             posortho=mutpos,
                                             ingap="Yes",
                                             mutID=as.character(BMS_mutants_temp$mutID[i]),
                                             ID=mutant_current_temp))
        
      }
      
    }
    
  }
  
  #add these mutants onto the existing data
  BMS_fitness_map_all <- rbind(BMS_fitness_map_all,BMS_fitness_map1)
  
  write.table(BMS_fitness_map1, file = paste("output/BMS_fitness_map_",as.character(BMS_cur_mut_num),".csv",sep=""), sep = ",", row.names = F,quote=F,col.names = T)
  
  #collapse values onto matrix
  #calc all data and stats
  BMS_fitness_collapsed_all <- BMS_fitness_map_all %>%
    filter(position > 0) %>%
    group_by(position, aa) %>%
    summarise(fitval=median(fitness),
              numpoints=n(),
              stdfit=sd(fitness))
  
  #these matrices have the fitness/num/sd for each aa at each position:
  BMS_matrix_all = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)
  BMS_matrix_all_num = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)
  BMS_matrix_all_sd = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)
  
  #populate matrix
  for (i in 1:nrow(BMS_fitness_collapsed_all)){
    
    BMS_matrix_all[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed_all$aa[i])),BMS_fitness_collapsed_all$position[i]] <- as.numeric(BMS_fitness_collapsed_all$fitval[i])
    BMS_matrix_all_num[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed_all$aa[i])),BMS_fitness_collapsed_all$position[i]] <- as.numeric(BMS_fitness_collapsed_all$numpoints[i])
    BMS_matrix_all_sd[which(BMS_aa_list$aa==as.character(BMS_fitness_collapsed_all$aa[i])),BMS_fitness_collapsed_all$position[i]] <- as.numeric(BMS_fitness_collapsed_all$stdfit[i])
  }
  
  #assign row (aa) and col (pos) names
  rownames(BMS_matrix_all)<-BMS_aa_list$aa
  colnames(BMS_matrix_all)<-c(1:BMS_ref_len)
  rownames(BMS_matrix_all_num)<-BMS_aa_list$aa
  colnames(BMS_matrix_all_num)<-c(1:BMS_ref_len)
  rownames(BMS_matrix_all_sd)<-BMS_aa_list$aa
  colnames(BMS_matrix_all_sd)<-c(1:BMS_ref_len)
  
  BMS_matrix_all_melt <- melt(BMS_matrix_all)
  BMS_matrix_all_num_melt <- melt(BMS_matrix_all_num)
  BMS_matrix_all_sd_melt <- melt(BMS_matrix_all_sd)
  
  #plot the data from these mutants:
  #median fitness data
  ggplot(data = BMS_matrix_all_melt, aes(x=X2, y=X1, fill=value)) +
    geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
    scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="Fitness",na.value="grey", limit = c(-5,5)) +
    theme_minimal() + 
    scale_x_continuous(breaks=seq(0,155,5))
  ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_max",toString(BMS_cur_mut_num),"_fit.pdf",sep=""))
  
  #number of points
  ggplot(data = BMS_matrix_all_num_melt, aes(x=X2, y=X1, fill=log10(value))) +
    geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
    scale_fill_gradient(low = "blue", high = "red",name="log(#points)",na.value="grey", limit = c(0,max(BMS_matrix_all_num_melt$value))) +
    theme_minimal() + 
    scale_x_continuous(breaks=seq(0,155,5))
  ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_max",toString(BMS_cur_mut_num),"_num.pdf",sep=""))
  
  #STD plot
  ggplot(data = BMS_matrix_all_sd_melt, aes(x=X2, y=X1, fill=value)) +
    geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
    scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="std(Fitness)",na.value="grey", limit = c(0,1.1*max(BMS_matrix_all_sd_melt$value))) +
    theme_minimal() + 
    scale_x_continuous(breaks=seq(0,155,5))
  ggsave(file=paste("plots/BMS_allpos_min",toString(BMS_min_BCs),"_max",toString(BMS_cur_mut_num),"_sd.pdf",sep=""))
}

#save fitness values
write.table(BMS_matrix_all_melt, 
            file = paste("output/BMS_matrix_all_melt.csv",sep=""),
            sep = ",", row.names = F,quote=F,col.names = T)
#save number of datapoints
write.table(BMS_matrix_all_num_melt, 
            file = paste("output/BMS_matrix_all_num_melt.csv",sep=""),
            sep = ",", row.names = F,quote=F,col.names = T)
#save std
write.table(BMS_matrix_all_sd_melt, 
            file = paste("output/BMS_matrix_all_sd_melt.csv",sep=""),
            sep = ",", row.names = F,quote=F,col.names = T)

######################################
# build custom plot (FIG 4)

####
# add new column for WT residues
BMS_matrix_all_melt$WTcolor <- NA
#X1 is aa, X2 is position

BMS_matrix_all_melt$aanum <- 0

#label the E. coli PPAT WT residues at each position
for (i in 1:nrow(BMS_matrix_all_melt)){
  
  #find the corresponding residue in E.coli using MSA
  BMS_cons_aanum <- BMS_matrix_all_melt$X2[i]
  
  #is this the WT E.coli residue?
  if (BMS_matrix_all_melt$X1[i] == substr(BMS_ecoli_seq, BMS_cons_aanum, BMS_cons_aanum)){
    
    #assign WT color
    BMS_matrix_all_melt$WTcolor[i] <- "red"
    
  }
  
  BMS_matrix_all_melt$aanum[i] <- which(BMS_aa_list$aa == BMS_matrix_all_melt$X1[i])
  
}

#relative solvent accesibility data
RSA_1H1T <- read.table( "structures/1H1T.RSAs.txt" , skip = 2 , header = FALSE )
#secondary structure data
SS_1H1T <- read.table( "structures/1H1T.SSs.txt" , skip = 2 , header = FALSE )
#conservation data
cons_1H1T <- read.table( "structures/PPAT_JSD_cons.txt" , skip = 5 , header = FALSE )

colnames(cons_1H1T) <- c("pos","cons","col")
cons_1H1T$pos <- cons_1H1T$pos+1

protein_info_1H1T <- RSA_1H1T %>%
  right_join(SS_1H1T, by="V1")

colnames(protein_info_1H1T) <- c("pos","RSA","SS")

protein_info_1H1T$cons <- 0

#add conservation value to data
for (i in 1:nrow(cons_1H1T)){
  #get MSA position: cons_1H1T$pos[i]
  #find if this position exists in the e.coli MSA
  if (cons_1H1T$pos[i] %in% BMS_ecoli_map$msa_aanum){
    #get corresponding residue
    e_coli_residue <- BMS_ecoli_map$orth_aanum[which(BMS_ecoli_map$msa_aanum==cons_1H1T$pos[i])]
    protein_info_1H1T$cons[which(protein_info_1H1T$pos==e_coli_residue)] <- cons_1H1T$cons[i]
  }
}

protein_info_1H1T <- BMS_matrix_all_melt %>%
  filter(X1 != "X") %>% #remove deletions
  group_by(X2) %>% #X2 is position
  dplyr::rename(pos=X2) %>%
  summarise(avgfit=mean(value,na.rm=TRUE),
            numcov=20-sum(is.na(value))) %>%
  right_join(protein_info_1H1T, by="pos")

protein_info_1H1T_melt <- melt(as.data.frame(protein_info_1H1T),
                               id=c("pos"))

#yval is the y-axis position to be used by each set of data in the plot
protein_info_1H1T_melt$yval <- 0
protein_info_1H1T_melt$yval[which(protein_info_1H1T_melt$variable=="avgfit")] <- 23
protein_info_1H1T_melt$yval[which(protein_info_1H1T_melt$variable=="RSA")] <- 24
protein_info_1H1T_melt$yval[which(protein_info_1H1T_melt$variable=="SS")] <- 25
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="X")] <- 21.3

#determine the coverage
(nrow(BMS_matrix_all_melt %>%
            filter(X1 != "X") %>%
            dplyr::select(value)) -
  sum(is.na(BMS_matrix_all_melt %>%
              filter(X1 != "X") %>%
              dplyr::select(value))))/nrow(BMS_matrix_all_melt %>%
                                             filter(X1 != "X") %>%
                                             dplyr::select(value))



BMS_Fig4_yminlim <- -2
BMS_Fig4_ymaxlim <- 30
BMS_Fig4_xminlim <- -1
BMS_Fig4_xmaxlim <- 160

#get the data determined from only 1 aa mutants
#this will be used to calculate the impact of single aa deletions
BMS_1mut_del <- read.table( "output/BMS_fitness_map_1.csv" , skip = 0 , sep=",",header = TRUE )
#calc all data and stats
BMS_1mut_del_collapsed_all <- BMS_1mut_del %>%
  filter(position > 0) %>%
  group_by(position, aa) %>%
  summarise(fitval=median(fitness),
            numpoints=n(),
            stdfit=sd(fitness))

#these matrices have the fitness/num/sd for each aa at each position:
BMS_1mut_del_matrix = matrix(rep(NA, BMS_ref_len*BMS_aa_dim),nrow=BMS_aa_dim,ncol=BMS_ref_len)

#populate matrix
for (i in 1:nrow(BMS_1mut_del_collapsed_all)){
  
  BMS_1mut_del_matrix[which(BMS_aa_list$aa==as.character(BMS_1mut_del_collapsed_all$aa[i])),BMS_1mut_del_collapsed_all$position[i]] <- as.numeric(BMS_1mut_del_collapsed_all$fitval[i])
}

rownames(BMS_1mut_del_matrix)<-BMS_aa_list$aa
colnames(BMS_1mut_del_matrix)<-c(1:BMS_ref_len)

BMS_1mut_del_matrix_melt <- melt(BMS_1mut_del_matrix)
#take only the deletions
BMS_1mut_del_matrix_melt <- BMS_1mut_del_matrix_melt %>%
  filter(X1=="X")
BMS_1mut_del_matrix_melt$WTcolor <- NA
BMS_1mut_del_matrix_melt$aanum <- 21.3

#this is the rest of the data
BMS_matrix_all_melt <- BMS_matrix_all_melt %>%
  filter(X1!="X")

#join 1aa deletions to rest of data
BMS_matrix_all_melt <- rbind(BMS_matrix_all_melt,BMS_1mut_del_matrix_melt)

BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="A")] <- 12
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="C")] <- 10
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="D")] <- 5
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="E")] <- 4
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="F")] <- 19
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="G")] <- 11
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="H")] <- 3
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="I")] <- 15
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="K")] <- 1
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="L")] <- 14
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="M")] <- 16
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="N")] <- 6
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="P")] <- 17
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="Q")] <- 7
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="R")] <- 2
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="S")] <- 9
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="T")] <- 8
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="V")] <- 13
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="W")] <- 20
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="Y")] <- 18
BMS_matrix_all_melt$aanum[which(BMS_matrix_all_melt$X1=="X")] <- 22.6

BMS_matrix_all_melt_WT <- BMS_matrix_all_melt %>%
  filter(WTcolor=="red")

###############
##FIG 4
p4 <- ggplot() +
  #BMS matrix
  geom_rect(data=BMS_matrix_all_melt,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1, fill=value))+
  #WT seq
  geom_rect(data=BMS_matrix_all_melt_WT,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1), color="green",alpha=0)+
  #avg fit
  geom_rect(data=protein_info_1H1T,aes(xmin=pos, xmax=pos+1, ymin=21.3, ymax=22.3,fill=avgfit))+
  labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="Fitness",na.value="grey", limit = c(-5,1.1*max(BMS_matrix_all_melt$value))) +
  geom_text(data=BMS_matrix_all_melt[1:21,], aes(x=-1, y=aanum+0.5, label=X1), size=3.5)+
  geom_text(data=data.frame(pos=seq(0,150,30)),aes(x=pos+0.5,y=0,label=pos))+
  #TxGH motif
  geom_segment(aes(x = 15, y = 1, xend = 15, yend = 0), colour = "red")+
  geom_segment(aes(x = 19, y = 1, xend = 19, yend = 0), colour = "red")+
  #ATP interacting residues
  geom_segment(aes(x = 7.5, y = 1, xend = 7.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 18.5, y = 1, xend = 18.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 89.5, y = 1, xend = 89.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 91.5, y = 1, xend = 91.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 124.5, y = 1, xend = 124.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 129.5, y = 1, xend = 129.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 130.5, y = 1, xend = 130.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 17.5, y = 1, xend = 17.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 21.5, y = 1, xend = 21.5, yend = 0), colour = "blue")+
  geom_segment(aes(x = 120.5, y = 1, xend = 120.5, yend = 0), colour = "blue")+
  
  #4'-phosphopantetheine interacting residues
  geom_segment(aes(x = 10.5, y = 1, xend = 10.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 42.5, y = 1, xend = 42.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 74.5, y = 1, xend = 74.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 88.5, y = 1, xend = 88.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 8.5, y = 1, xend = 8.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 9.5, y = 1, xend = 9.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 37.5, y = 1, xend = 37.5, yend = 0), colour = "magenta")+
  geom_segment(aes(x = 102.5, y = 1, xend = 102.5, yend = 0), colour = "magenta")+
  
  geom_text(aes(x=18,y=-1.5,label="TxGH\nmotif"),size=3)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())+
  xlim(BMS_Fig4_xminlim,BMS_Fig4_xmaxlim)+
  ylim(BMS_Fig4_yminlim,BMS_Fig4_ymaxlim)

p4
ggsave(file="plots/Fig4_p4_matrix_legend.pdf",
       width=13.3, height=6.5, units="in")
p4 + theme(legend.position="none")
ggsave(file="plots/Fig4_p4_matrix_nolegend.pdf",
       width=13.3, height=6.5, units="in")
#Saving 13.3 x 6.51 in image

#fitness limits (average for position)
max(protein_info_1H1T$avgfit)
min(protein_info_1H1T$avgfit)

#take only beta-strand and alpha-helix SS info
protein_info_1H1T_no_loop <- protein_info_1H1T %>%
  filter(SS != "loop")
#plot SS
p5 <- ggplot()+
  geom_segment(aes(x = 1, y = 24.5, xend = 160, yend = 24.5), colour = "black")+
  geom_rect(data=protein_info_1H1T_no_loop,aes(xmin=pos, xmax=pos+1, ymin=24, ymax=25, fill=SS))+
  xlim(BMS_Fig4_xminlim,BMS_Fig4_xmaxlim)+
  ylim(BMS_Fig4_yminlim,BMS_Fig4_ymaxlim)+
  labs(x = "", y ="",color="") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())
p5
ggsave(file="plots/Fig4_p5_SS_legend.pdf",
       width=13.3, height=6.5, units="in")
p5 + theme(legend.position="none")
ggsave(file="plots/Fig4_p5_SS_nolegend.pdf",
       width=13.3, height=6.5, units="in")

#plot RSA
p6 <- ggplot(protein_info_1H1T)+
  geom_rect(aes(xmin=pos, xmax=pos+1, ymin=25.3, ymax=26.3, fill=RSA))+
  xlim(BMS_Fig4_xminlim,BMS_Fig4_xmaxlim)+
  ylim(BMS_Fig4_yminlim,BMS_Fig4_ymaxlim)+
  labs(x = "", y ="",color="") +
  scale_fill_gradient(low = "white", high = "red",name="RSA",na.value="grey") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())
p6
ggsave(file="plots/Fig4_p6_RSA_legend.pdf",
       width=13.3, height=6.5, units="in")
p6 + theme(legend.position="none")
ggsave(file="plots/Fig4_p6_RSA_nolegend.pdf",
       width=13.3, height=6.5, units="in")

#conservation
p7 <- ggplot(protein_info_1H1T)+
  geom_rect(aes(xmin=pos, xmax=pos+1, ymin=26.6, ymax=27.6, fill=cons))+
  xlim(BMS_Fig4_xminlim,BMS_Fig4_xmaxlim)+
  ylim(BMS_Fig4_yminlim,BMS_Fig4_ymaxlim)+
  labs(x = "", y ="",color="") +
  scale_fill_gradient(low = "white", high = "red",name="Cons",na.value="grey") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())
p7
#scale_x_continuous(breaks=seq(0,155,5)) +
ggsave(file="plots/Fig4_p7_cons_legend.pdf",
       width=13.3, height=6.5, units="in")
p7 + theme(legend.position="none")
ggsave(file="plots/Fig4_p7_cons_nolegend.pdf",
       width=13.3, height=6.5, units="in")

#calc sitewise EVmutations
EVmutsitewise <- EVmutPPAT %>%
  group_by(pos) %>%
  summarise(siteEV=mean(effect_prediction_epistatic))

protein_info_1H1T <- EVmutsitewise %>%
  inner_join(protein_info_1H1T, by="pos")
  
p8 <- ggplot(EVmutsitewise)+
  geom_rect(aes(xmin=pos, xmax=pos+1, ymin=27.9, ymax=28.9, fill=siteEV))+
  xlim(BMS_Fig4_xminlim,BMS_Fig4_xmaxlim)+
  ylim(BMS_Fig4_yminlim,BMS_Fig4_ymaxlim)+
  labs(x = "", y ="",color="") +
  scale_fill_gradient(low = "red", high = "white",name="Cons",na.value="grey") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())
p8
ggsave(file="plots/Fig4_p8_EVmut_legend.pdf",
       width=13.3, height=6.5, units="in")
p8 + theme(legend.position="none")
ggsave(file="plots/Fig4_p8_EVmut_nolegend.pdf",
       width=13.3, height=6.5, units="in")

#the saved pdfs should be merged manually in Illustrator
#End of FIG4 code

########################
# more plots

# avg fitness vs conservation
p <- ggplot(protein_info_1H1T,aes(x=cons,y=avgfit))+
  geom_smooth()+
  geom_point()+
  labs(x = "Position conservation", y ="Average fitness at position",color="") +
  ylim(-2.7,2.7)
p <- ggExtra::ggMarginal(p,type = "histogram",
                         xparams = list(bins=30),
                         yparams = list(bins=20),
                         col = 'black',
                         fill = 'gray')
p
ggsave(file="plots/avg_fit_vs_cons.pdf")
rm(p)

cor.test(protein_info_1H1T$cons,protein_info_1H1T$avgfit)

####################################
# site coverage vs avg fitness

p <- ggplot(protein_info_1H1T,aes(x=avgfit,y=numcov/20*100))+
  geom_smooth()+
  geom_point()+
  labs(x = "Average fitness at position", y ="Position mutational coverage (%)",color="")
p
ggsave(file="plots/coverage_vs_avg_fit.pdf",
       width=6, height=5, units="in")
rm(p)

cor.test(protein_info_1H1T$avgfit,protein_info_1H1T$numcov/20*100)

####################################
# avg fitness vs secondary structure
p <- ggplot(protein_info_1H1T,aes(x=SS,y=avgfit))+
  geom_boxplot()+
  geom_jitter()+
  labs(x = "Secondary structure", y ="Average fitness at position",color="")
p
ggsave(file="plots/avg_fit_vs_ss.pdf")
rm(p)

####################################
# avg fitness vs RSA
p2 <- ggplot(protein_info_1H1T,aes(x=RSA,y=avgfit))+
  geom_smooth()+
  geom_point(alpha=0.7)+
  labs(x = "Relative solvent accesibility", y ="Average fitness at position",color="")+
  ylim(-2.7,2.7)
p2 <- ggExtra::ggMarginal(p2,type = "histogram",
                         xparams = list(bins=30),
                         yparams = list(bins=20),
                         col = 'black',
                         fill = 'gray')
p2
ggsave(file="plots/avg_fit_vs_rsa.pdf")
rm(p2)
cor.test(protein_info_1H1T$RSA,protein_info_1H1T$avgfit)

####################################
# avg fitness vs EVmutations

p20 <- ggplot(protein_info_1H1T,aes(x=siteEV,y=avgfit))+
  labs(x = "Average site EVMutations statistical energy", 
       y ="Average fitness at position",color="")+
  geom_smooth()+
  geom_point()
p20
ggsave(file="plots/avg_fit_vs_evmut.pdf")
rm(p20)
cor.test(protein_info_1H1T$avgfit,protein_info_1H1T$siteEV)
############################################
#save the avf fit to be used to color structural models

#for Pymol:
#save the avg fit to a file to be used by pymol to color using b-factors:
#it's rescale to be between 0 and 40
write(40/max(protein_info_1H1T$avgfit+-1*min(protein_info_1H1T$avgfit))*(protein_info_1H1T$avgfit+-1*min(protein_info_1H1T$avgfit)), file = "structures/avg_fit_as_new_BFactors.txt",
      ncolumns = 1,append = FALSE, sep = "")

#for Chimera:
#https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/defineattrib.html
sink("structures/avg_fit_1H1T_chimera_attribute.txt")
cat("# BMS data for 1H1T\n")
cat("attribute: percentExposed\n")
cat("match mode: 1-to-1\n")
cat("recipient: residues\n")
for (i in 1:nrow(protein_info_1H1T)){
  cat(paste("\t:",
            as.character(protein_info_1H1T$pos[i]),
            "\t",
            as.character(protein_info_1H1T$avgfit[i]),
            "\n",
            sep=""))
}

sink()

#####################################
#plot of WT seq fitness vs everything else

ggplot(BMS_matrix_all_melt,aes(x=WTcolor,y=value))+
  geom_boxplot()+
  labs(x = "Residue", y ="Fitness",color="")+ 
  scale_x_discrete(name ="Residue type",
                   labels=c("Wildtype E. coli PPAT","Other"))
ggsave(file="plots/WT_vs_other_fitness.pdf",
       width=5.5, height=4.5, units="in")

#values for WT
mean(as.numeric(unlist(BMS_matrix_all_melt %>% filter(WTcolor=="red") %>% dplyr::select(value))))
median(as.numeric(unlist(BMS_matrix_all_melt %>% filter(WTcolor=="red") %>% dplyr::select(value))))
sd(as.numeric(unlist(BMS_matrix_all_melt %>% filter(WTcolor=="red") %>% dplyr::select(value))))

#values for non-WT
mean(as.numeric(unlist(BMS_matrix_all_melt %>% filter(is.na(WTcolor) & !is.na(value)) %>% dplyr::select(value))))
median(as.numeric(unlist(BMS_matrix_all_melt %>% filter(is.na(WTcolor) & !is.na(value)) %>% dplyr::select(value))))
sd(as.numeric(unlist(BMS_matrix_all_melt %>% filter(is.na(WTcolor) & !is.na(value)) %>% dplyr::select(value))))

#################################
#plot STD of BMS

BMS_matrix_all_sd_melt <- BMS_matrix_all_sd_melt %>%
  filter(X1!="X")

BMS_matrix_all_sd_melt$aanum <- 0

BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="A")] <- 12
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="C")] <- 10
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="D")] <- 5
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="E")] <- 4
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="F")] <- 19
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="G")] <- 11
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="H")] <- 3
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="I")] <- 15
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="K")] <- 1
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="L")] <- 14
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="M")] <- 16
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="N")] <- 6
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="P")] <- 17
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="Q")] <- 7
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="R")] <- 2
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="S")] <- 9
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="T")] <- 8
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="V")] <- 13
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="W")] <- 20
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="Y")] <- 18
BMS_matrix_all_sd_melt$aanum[which(BMS_matrix_all_sd_melt$X1=="X")] <- 22.6


p4 <- ggplot() +
  #BMS matrix
  geom_rect(data=BMS_matrix_all_sd_melt,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1, fill=value))+
  #WT seq
  geom_rect(data=BMS_matrix_all_melt_WT,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1), color="green",alpha=0)+
  labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", high = "red", mid="gold",name="Std",na.value="grey", midpoint = 3, limit = c(0,1.1*max(BMS_matrix_all_sd_melt$value))) +#
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())

p4
###############################
#plot STD only for sites with at least 4 points
BMS_epistasis <- BMS_matrix_all_sd_melt %>%
  dplyr::rename(sd=value) %>%
  inner_join(BMS_matrix_all_num_melt, by=c("X1","X2")) %>%
  dplyr::rename(num=value) %>%
  filter(num>4)

p4 <- ggplot() +
  #BMS matrix
  geom_rect(data=BMS_epistasis,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1, fill=sd))+
  #WT seq
  geom_rect(data=BMS_matrix_all_melt_WT,aes(xmin=X2, xmax=X2+1, ymin=aanum, ymax=aanum+1), color="green",alpha=0)+
  geom_text(data=BMS_matrix_all_melt[1:21,], aes(x=-1, y=aanum+0.5, label=X1), size=3.5)+
  geom_text(data=data.frame(pos=seq(0,150,30)),aes(x=pos+0.5,y=0,label=pos))+
  labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient(low = "gold", high = "red", name="Std",na.value="grey", limit = c(3,1.1*max(BMS_epistasis$sd))) +#
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line=element_blank())

p4

BMS_epistasis_pos <- BMS_epistasis %>%
  group_by(X2) %>%
  summarise(numaa=n(),meansd=mean(sd),stdsd=sd(sd))
  
#by position:
#number of points
ggplot(BMS_epistasis_pos,aes(x=X2,y=numaa)) +
  geom_point()
#STD
ggplot(BMS_epistasis_pos,aes(x=X2,y=meansd)) +
  geom_point()
#STD STD
ggplot(BMS_epistasis_pos,aes(x=X2,y=stdsd)) +
  geom_point()

########################
#output the data used to generate FIG 4
#add number of points to fitness values
BMS_info <- right_join(BMS_matrix_all_melt %>%
                         select(X1,X2,value) %>%
                         dplyr::rename(AA=X1,Pos=X2,avgfitness=value),
                       BMS_matrix_all_num_melt %>%
                         select(X1,X2,value) %>%
                         dplyr::rename(AA=X1,Pos=X2,numpoints=value),
                       by=c("AA", "Pos"))
#add STD values
BMS_info <- right_join(BMS_info,
                       BMS_matrix_all_sd_melt%>%
                         select(X1,X2,value) %>%
                         dplyr::rename(AA=X1,Pos=X2,sd=value),
                       by=c("AA", "Pos"))
#save CSV
write.table(BMS_info,
            file = paste("output/BMS_info.csv",sep=""), 
            sep = ",", row.names = F,quote=F,col.names = T)

#################
# output a file with all mutants info

names(mutants)

write.table(mutants %>%
              dplyr::select(mutID,
                            IDalign,
                            numprunedBCs,
                            mutations,
                            globalfit14,
                            seq) %>%
              filter(mutations < 21 & mutations > -1),
            file = paste("output/all_mutants_minimal_db_mut20aa.csv",sep=""), 
            sep = ",", 
            row.names = F,
            quote=F,
            col.names = F)

###########################

source("GOF_paper.R")

save.image(file = "BMSanalysis.RData")
#save(orcollapse3info, orcollapse3, mutants, perfects, perfects_tree, file = "PPATdata.RData")
