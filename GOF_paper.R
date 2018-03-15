# Gain of Function Analysis
# Calin Plesa
#

library(stringr)
library(ggplot2)
library(ggExtra)
library(dplyr)
require(reshape)
library(cowplot)
library(RWebLogo)
library(VennDiagram)

#load data
#load("PPATdata.RData")
load("dropout_mutants_GOF.RData")

###################################
# find GOF for dropouts:

mutants_to_plot <- data.frame(ID=unique(dropout_mutants_GOF$ID))

#make sure all of the chosen dropouts are in the MSA
mutants_to_plot_index_of_mut_to_drop = numeric()
for (i in 1:nrow(mutants_to_plot)){
  mutant_current_temp <- mutants_to_plot$ID[i]
  if(!file.exists(paste("alignments/MSAdropouts/",mutant_current_temp,".csv",sep=""))){
    print(mutant_current_temp)
    mutants_to_plot_index_of_mut_to_drop[length(mutants_to_plot_index_of_mut_to_drop)+1] <- which(mutants_to_plot$ID==mutant_current_temp)
  }
}
#any not in MSA?
if (length(mutants_to_plot_index_of_mut_to_drop) > 0){
  mutants_to_plot <- mutants_to_plot[-mutants_to_plot_index_of_mut_to_drop,]
}

#read in the E.coli map
ecoli_map <- read.csv(file=paste("alignments/MSAdropouts/NP_418091.csv",sep=""),head=TRUE,sep=",")
ecoli_map_all <- read.csv(file=paste("alignments/MSAmap_all/NP_418091.csv",sep=""),head=TRUE,sep=",")

#make a new data frame which will keep all info
GOF_fitness_map <- data.frame(position=numeric(),
                              aa=character(),
                              mutations=numeric(),
                              fitness=numeric(),
                              posortho=numeric(),
                              ingap=character(),
                              mutID=character(),
                              ID=character())

aminoacids<-data.frame(aa=c('G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T'),aanum=c(1:20))

#loop over all low-fitness homologs
for (iii in 1:nrow(mutants_to_plot)){
  
  #current homolog:
  mutant_current <- as.character(mutants_to_plot$ID[iii])
  
  #length of name
  name_size = nchar(paste(mutant_current,"_",sep=""))
  
  #get the MSA mapping
  mutant_map <- read.csv(file=paste("alignments/MSAdropouts/",mutant_current,".csv",sep=""),head=TRUE,sep=",")
  
  #grab the mutants which are within a certain distance
  GOFmutIDinfo_temp <- dropout_mutants_GOF %>%
    filter(ID == mutant_current)
  
  #loop over all mutants for this construct:
  for (mn in 1:nrow(GOFmutIDinfo_temp)) {
    
    #this mutants fitness
    gof_fit_temp <- GOFmutIDinfo_temp$globalfit14[mn]
    
    #grab the mut name
    mutations_names <- as.character(GOFmutIDinfo_temp$mutID[mn])
    
    #grab only the relevant portion of the name
    mutations_names <- substr(mutations_names, name_size+1, nchar(mutations_names))
    
    ## split mutation string at non-digits
    s <- strsplit(mutations_names, "_")
    
    for (mutnum in 1:GOFmutIDinfo_temp$mutations[mn]){
      
      #grab the corresponding mutation string
      mutcurr<-s[[1]][mutnum]
      
      #get the position
      mutpos <- as.numeric(str_extract(mutcurr, "[0-9]+"))
      
      #get ending aa
      to_aa <- substr(mutcurr, nchar(mutpos)+2, nchar(mutcurr))
      
      #find the number in the consensus seq
      gof_cons_aanum <- mutant_map$msa_aanum[which(mutant_map$orth_aanum == mutpos)]
      
      #does this map to a non-gap
      if (gof_cons_aanum %in% ecoli_map$msa_aanum){
        
        #the corresponding e.coli residue
        e_coli_residue <- ecoli_map$orth_aanum[which(ecoli_map$msa_aanum == gof_cons_aanum)]
        
        #add this point to the data
        GOF_fitness_map <- rbind(GOF_fitness_map,
                                 data.frame(position=e_coli_residue,
                                            aa=to_aa,
                                            mutations=GOFmutIDinfo_temp$mutations[mn],
                                            fitness=gof_fit_temp,
                                            posortho=mutpos,
                                            ingap="No",
                                            mutID=GOFmutIDinfo_temp$mutID[mn],
                                            ID=GOFmutIDinfo_temp$ID[mn]))
        
      } else {
        #if it's here it maps to a gap
        
        #add this point to the data
        GOF_fitness_map <- rbind(GOF_fitness_map,
                                 data.frame(position=-1,
                                            aa=to_aa,
                                            mutations=GOFmutIDinfo_temp$mutations[mn],
                                            fitness=gof_fit_temp,
                                            posortho=mutpos,
                                            ingap="Yes",
                                            mutID=GOFmutIDinfo_temp$mutID[mn],
                                            ID=GOFmutIDinfo_temp$ID[mn]))
        
      }
    }
  }
}

rm(iii,
   mn,
   mutnum,
   gof_cons_aanum,
   mutant_current,
   name_size,
   mutant_map,
   GOFmutIDinfo_temp,
   mutants_to_plot,
   gof_fit_temp,
   mutations_names,
   mutpos,
   mutcurr,
   to_aa,
   s)

write.table(GOF_fitness_map, file = "output/GOF_fitness_map.csv", sep = ",", row.names = F,quote=F,col.names = T)

#calc all data and stats
GOF_fitness_collapsed_all <- GOF_fitness_map %>%
  filter(position > 0) %>%
  group_by(position, aa) %>%
  summarise(fitval=median(fitness),
            numpoints=n(),
            stdfit=sd(fitness),
            numortho=length(unique(ID)))

gof_aa_dim <- nrow(aminoacids)
gof_ref_len <- nrow(ecoli_map)

#these matrices have the median fitness, number of mutants, sd, 
#and number of homologs for each aa at each position:
gof_matrix = matrix(rep(NA, gof_ref_len*gof_aa_dim),nrow=gof_aa_dim,ncol=gof_ref_len)
gof_matrix_num = matrix(rep(NA, gof_ref_len*gof_aa_dim),nrow=gof_aa_dim,ncol=gof_ref_len)
gof_matrix_sd = matrix(rep(NA, gof_ref_len*gof_aa_dim),nrow=gof_aa_dim,ncol=gof_ref_len)
gof_matrix_numortho = matrix(rep(NA, gof_ref_len*gof_aa_dim),nrow=gof_aa_dim,ncol=gof_ref_len)

#populate matrix
for (i in 1:nrow(GOF_fitness_collapsed_all)){
  gof_matrix[which(aminoacids$aa==as.character(GOF_fitness_collapsed_all$aa[i])),GOF_fitness_collapsed_all$position[i]] <- as.numeric(GOF_fitness_collapsed_all$fitval[i])
  gof_matrix_num[which(aminoacids$aa==as.character(GOF_fitness_collapsed_all$aa[i])),GOF_fitness_collapsed_all$position[i]] <- as.numeric(GOF_fitness_collapsed_all$numpoints[i])
  gof_matrix_sd[which(aminoacids$aa==as.character(GOF_fitness_collapsed_all$aa[i])),GOF_fitness_collapsed_all$position[i]] <- as.numeric(GOF_fitness_collapsed_all$stdfit[i])
  gof_matrix_numortho[which(aminoacids$aa==as.character(GOF_fitness_collapsed_all$aa[i])),GOF_fitness_collapsed_all$position[i]] <- as.numeric(GOF_fitness_collapsed_all$numortho[i])
}

#assign row and col names
rownames(gof_matrix)<-aminoacids$aa
colnames(gof_matrix)<-c(1:gof_ref_len)
rownames(gof_matrix_num)<-aminoacids$aa
colnames(gof_matrix_num)<-c(1:gof_ref_len)
rownames(gof_matrix_sd)<-aminoacids$aa
colnames(gof_matrix_sd)<-c(1:gof_ref_len)
rownames(gof_matrix_numortho)<-aminoacids$aa
colnames(gof_matrix_numortho)<-c(1:gof_ref_len)

gof_matrix_melt <- melt(gof_matrix)
gof_matrix_num_melt <- melt(gof_matrix_num)
gof_matrix_sd_melt <- melt(gof_matrix_sd)
gof_matrix_numortho_melt <- melt(gof_matrix_numortho)

############################
#plot the data from all mutants:
ggplot(data = gof_matrix_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile() +
  labs(x = "Position (aa)",
       y ="Amino acid",color="") +
  scale_fill_gradient(low = "blue", 
                       high = "red",
                       name="Fitness",
                       na.value="grey",
                       limit = c(0,5)) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,150,10))
ggsave(file=paste("plots/GOF/GOF_fit.pdf",sep=""))
ggsave(file=paste("plots/GOF/GOF_fit.png",sep=""))

ggplot(data = gof_matrix_num_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile()+ 
  labs(x = "Position (aa)",
       y ="Amino acid",color="") +
  scale_fill_gradient(low = "blue", 
                      high = "red",
                      name="#points",
                      na.value="grey", 
                      limit = c(0,max(gof_matrix_num_melt$value))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,150,10))
ggsave(file=paste("plots/GOF/GOF_num.pdf",sep=""))
ggsave(file=paste("plots/GOF/GOF_num.png",sep=""))

ggplot(data = gof_matrix_numortho_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile()+ 
  labs(x = "Position (aa)",
       y ="Amino acid",color="") +
  scale_fill_gradient(low = "blue", 
                      high = "red",
                      name="#homologs",
                      na.value="grey", 
                      limit = c(0,max(gof_matrix_numortho_melt$value))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,150,10))
ggsave(file=paste("plots/GOF/GOF_numhomologs.pdf",sep=""))
ggsave(file=paste("plots/GOF/GOF_numhomologs.png",sep=""))

ggplot(data = gof_matrix_sd_melt, aes(x=X2, y=X1, fill=value)) +
  geom_tile()+ labs(x = "Position (aa)", y ="Amino acid",color="") +
  scale_fill_gradient2(low = "blue", 
                       high = "red", 
                       mid="gold",
                       name="std(Fitness)",
                       na.value="grey", 
                       limit = c(0,1.1*max(gof_matrix_sd_melt$value))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,150,10))
ggsave(file=paste("plots/GOF/GOF_sd.pdf",sep=""))
ggsave(file=paste("plots/GOF/GOF_sd.pdf",sep=""))

GOF_fitness_collapsed_by_pos <- GOF_fitness_map %>%
  filter(position > 0) %>%
  group_by(position) %>%
  summarise(fitval=median(fitness),
            numpoints=n(),
            stdfit=sd(fitness),
            numortho=length(unique(ID)))

##############################################
#FIG 4C GOF
p <- ggplot(GOF_fitness_collapsed_by_pos, aes(x=position, y=numpoints, color=numortho)) +
  geom_segment(aes(x = 0, y = mean(numpoints)+2*sd(numpoints), xend = 160, yend = mean(numpoints)+2*sd(numpoints)),linetype=2,colour = "blue")+
  geom_segment(aes(x = 0, y = mean(numpoints), xend = 160, yend = mean(numpoints)),linetype=2,colour = "red")+
  geom_point(size=1.8)+
  labs(x = "Position (aa)", y ="Number of gain-of-function mutants",color="") +
  scale_color_gradient(low = "blue", 
                       high = "red",
                       name="Num.\nUniq.\nOrtho.",
                       na.value="grey", 
                       limit = c(0,1.1*max(GOF_fitness_collapsed_by_pos$numortho))) +
  scale_x_continuous(breaks=seq(0,160,20))+
  theme(legend.position="left")
p <- ggExtra::ggMarginal(p,type = "histogram",
                    margins = "y",
                    bins=21,
                    col = 'black',
                    fill = 'red')
p
ggsave(p, 
       file=paste("plots/GOF/Fig_5_GOF_hist.pdf",
                  sep=""),
       width = 6.7, height = 5, units = "in")
ggsave(p, 
       file=paste("plots/GOF/Fig_5_GOF_hist.png",
                  sep=""),
       width = 6.7, height = 5, units = "in")

GOF_fitness_collapsed_by_pos_2sigma <- GOF_fitness_collapsed_by_pos %>%
  filter(numpoints >= (mean(GOF_fitness_collapsed_by_pos$numpoints)+2*sd(GOF_fitness_collapsed_by_pos$numpoints)))

GOF_fitness_collapsed_by_pos_2sigma <- protein_info_1H1T %>%
  dplyr::rename(position=pos) %>%
  filter(position %in% GOF_fitness_collapsed_by_pos_2sigma$position) %>%
  right_join(GOF_fitness_collapsed_by_pos_2sigma,by="position") %>%
  arrange(position)

p2 <- ggplot(GOF_fitness_collapsed_by_pos_2sigma,aes(x=cons,y=RSA,color=as.factor(position),shape=as.factor(position)))+
  geom_point(alpha=0.9,size=3)+
  labs(x = "Site Conservation", y ="Relative Solvent Accessibility",color="Residue") +
  scale_color_manual(name = "Residue",
                     values = c("red", 
                                "red", 
                                "blue", 
                                "red",
                                "red",
                                "green",
                                "magenta",
                                "magenta",
                                "magenta"))+
  scale_shape_manual(name = "Residue",
                     values = c(0,1,2,3,4,5,6,7,8,9))
p2
save_plot("plots/GOF/GOF_hit_stats.pdf", p2,
          base_aspect_ratio = 1.3)
save_plot("plots/GOF/GOF_hit_stats.png", p2,
          base_aspect_ratio = 1.3)


ggplot(GOF_fitness_collapsed_by_pos, aes(numpoints,fill="red")) +
  geom_histogram(bins=max(GOF_fitness_collapsed_by_pos$numpoints))+
  labs(x = "GOF hits per residue", y ="Counts",color="") +
  theme_minimal() + 
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(0,160,10))
ggsave(file=paste("plots/GOF/GOF_bypos_hist.pdf",sep=""))

ggplot(GOF_fitness_collapsed_by_pos, aes(x=position, y=numortho, color=fitval)) +
  geom_point()+
  labs(x = "Position (aa)", y ="Number of homologs with GOF mutants",color="") +
  scale_color_gradient(low = "blue", 
                       high = "red",
                       name="Fitness",
                       na.value="grey", 
                       limit = c(0,1.1*max(GOF_fitness_collapsed_by_pos$fitval))) +
  theme_minimal() + 
  scale_x_continuous(breaks=seq(0,160,10))
ggsave(file=paste("plots/GOF/GOF_bypos_numhomologs.pdf",sep=""))

sd(GOF_fitness_collapsed_by_pos$numpoints)
mean(GOF_fitness_collapsed_by_pos$numpoints)

floor(mean(GOF_fitness_collapsed_by_pos$numpoints)+2*sd(GOF_fitness_collapsed_by_pos$numpoints))
mean(GOF_fitness_collapsed_by_pos$numpoints)+2*sd(GOF_fitness_collapsed_by_pos$numpoints)

GOF_cutoff <- mean(GOF_fitness_collapsed_by_pos$numpoints)+2*sd(GOF_fitness_collapsed_by_pos$numpoints)

GOF_top_hits <- GOF_fitness_collapsed_by_pos %>%
  filter(numpoints >= GOF_cutoff) %>%
  arrange(desc(numpoints))

EVCinfo = read.csv("ev_couplings/PPAT_EC_residues_summary_table.csv", skip=1, head=TRUE)  # read csv file
SCAinfo = read.csv("SCA/SCAsectors.csv", head=TRUE)  # read csv file

EVCinfo$EVCrank <- numeric(nrow(EVCinfo))
EVCinfo$EVCrank <- c(1:nrow(EVCinfo))
GOF_top_hits$SCAsector <- NA
GOF_top_hits$EVCrank <- NA
GOF_top_hits$cumulative.strength <- NA
GOF_top_hits$EC.strength <- NA
GOF_top_hits$ECconservation <- NA
GOF_top_hits$number.of.ECs <- NA
GOF_top_hits$MSAall_res_num <- NA
GOF_top_hits$MSAdrop_res_num <- NA

for (i in 1:nrow(GOF_top_hits)){
  
  EVC_rownum <- which(EVCinfo$residue.index==GOF_top_hits$position[i])
  if (length(EVC_rownum)>0){
    GOF_top_hits$cumulative.strength[i] <- EVCinfo$cumulative.strength[EVC_rownum]
    GOF_top_hits$EC.strength[i] <- EVCinfo$EC.strength[EVC_rownum]
    GOF_top_hits$ECconservation[i] <- EVCinfo$conservation[EVC_rownum]
    GOF_top_hits$number.of.ECs[i] <- EVCinfo$number.of.ECs[EVC_rownum]
    GOF_top_hits$EVCrank[i] <- EVCinfo$EVCrank[EVC_rownum]
    
  }
  SCA_rownum <- which(SCAinfo$Residue==GOF_top_hits$position[i])
  if (length(SCA_rownum)>0){
    GOF_top_hits$SCAsector[i] <- SCAinfo$Sector[SCA_rownum]
  }
  GOF_top_hits$MSAdrop_res_num[i] <- ecoli_map$msa_aanum[which(ecoli_map$orth_aanum == GOF_top_hits$position[i])]
  GOF_top_hits$MSAall_res_num[i] <- ecoli_map_all$msa_aanum[which(ecoli_map_all$orth_aanum == GOF_top_hits$position[i])]
  
}
rm(EVC_rownum,
   SCA_rownum)

GOF_top_hits$color <- NA
GOF_top_hits$color[which(GOF_top_hits$position==68)] <- "red"
GOF_top_hits$color[which(GOF_top_hits$position==69)] <- "red"
GOF_top_hits$color[which(GOF_top_hits$position==35)] <- "red"
GOF_top_hits$color[which(GOF_top_hits$position==34)] <- "red"
GOF_top_hits$color[which(GOF_top_hits$position==134)] <- "magenta"
GOF_top_hits$color[which(GOF_top_hits$position==135)] <- "magenta"
GOF_top_hits$color[which(GOF_top_hits$position==137)] <- "magenta"
GOF_top_hits$color[which(GOF_top_hits$position==103)] <- "green"
GOF_top_hits$color[which(GOF_top_hits$position==64)] <- "aquamarine"

#make pml file for top GOF hits
sink("output/dropout_GOF.pml",append = FALSE)
cat("fetch 1h1t, async=0\n")
cat("hide all\n")
cat("color grey80, chain A\n")
cat("color grey80, chain B\n")
cat("set cartoon_flat_sheets, 0\n")
cat("set cartoon_smooth_loops, 0\n")
cat("set cartoon_fancy_helices, 1\n")
cat("show cartoon, chain A\n")
cat("show cartoon, chain B\n")
for (i in 1:nrow(GOF_top_hits)){
  cat(paste("select resi ",as.character(GOF_top_hits$position[i])," and chain A\n"),sep="")
  cat(paste("color ",as.character(GOF_top_hits$color[i]),", sele\n"))
  cat(paste("select resi ",as.character(GOF_top_hits$position[i])," and chain B\n"),sep="")
  cat(paste("color ",as.character(GOF_top_hits$color[i]),", sele\n"))
}
cat("select none\n")
sink()

write.table(GOF_top_hits, file = "output/GOF_top_hits.csv", sep = ",", row.names = F,quote=F,col.names = T)

########################################
# Do GOF mutations cluster?

#calc all IDs which drop out:
#dropout_mutants_GOF$ID

#IDs with GOF in each patch:
#red pos: 34, 35, 68, 69
GOF_red_IDs <- GOF_fitness_map %>%
  filter(position == 34 |
           position == 35 |
           position == 68 |
           position == 69) %>%
  select(ID) %>%
  unique(.)
GOF_red_IDs$set <- "red"


#magenta: 134, 135, 137
GOF_magenta_IDs <- GOF_fitness_map %>%
  filter(position == 134 |
           position == 135) %>%
  select(ID) %>%
  unique(.)
GOF_magenta_IDs$set <- "magenta"

#green: 103
GOF_green_IDs <- GOF_fitness_map %>%
  filter(position == 103) %>%
  select(ID) %>%
  unique(.)
GOF_green_IDs$set <- "green"

#blue: 64
GOF_blue_IDs <- GOF_fitness_map %>%
  filter(position == 64) %>%
  select(ID) %>%
  unique(.)
GOF_blue_IDs$set <- "blue"

GOF_patch_data_for_Venn <- rbind(GOF_red_IDs, GOF_magenta_IDs)
GOF_patch_data_for_Venn <- rbind(GOF_patch_data_for_Venn, GOF_green_IDs)
GOF_patch_data_for_Venn <- rbind(GOF_patch_data_for_Venn, GOF_blue_IDs)


dev.off()
#plot(venneuler(GOF_patch_data_for_Venn))
#dev.new()
#quartz()
# Reference four-set diagram
venn.plot <- draw.quad.venn(
  area1 = length(GOF_red_IDs$ID),
  area2 = length(GOF_magenta_IDs$ID),
  area3 = length(GOF_green_IDs$ID),
  area4 = length(GOF_blue_IDs$ID),
  n12 = length(intersect(GOF_red_IDs$ID,GOF_magenta_IDs$ID)),
  n13 = length(intersect(GOF_red_IDs$ID,GOF_green_IDs$ID)),
  n14 = length(intersect(GOF_red_IDs$ID,GOF_blue_IDs$ID)),
  n23 = length(intersect(GOF_magenta_IDs$ID,GOF_green_IDs$ID)),
  n24 = length(intersect(GOF_magenta_IDs$ID,GOF_blue_IDs$ID)),
  n34 = length(intersect(GOF_green_IDs$ID,GOF_blue_IDs$ID)),
  n123 = length(intersect(intersect(GOF_red_IDs$ID,GOF_magenta_IDs$ID),GOF_green_IDs$ID)),
  n124 = length(intersect(intersect(GOF_red_IDs$ID,GOF_magenta_IDs$ID),GOF_blue_IDs$ID)),
  n134 = length(intersect(intersect(GOF_red_IDs$ID,GOF_green_IDs$ID),GOF_blue_IDs$ID)),
  n234 = length(intersect(intersect(GOF_magenta_IDs$ID,GOF_green_IDs$ID),GOF_blue_IDs$ID)),
  n1234 = length(intersect(intersect(intersect(GOF_red_IDs$ID,GOF_magenta_IDs$ID),GOF_green_IDs$ID),GOF_blue_IDs$ID)),
  category = c("34,35,68,69", "134,135", "103", "64"),
  fill = c("red", "magenta", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("red", "magenta", "green", "blue"))
# Writing to file
#tiff(filename = "Quad_Venn_diagram.tiff", compression = "lzw");
#grid.draw(venn.plot);
ggsave(file=paste("plots/GOF/GOF_Venn.pdf",sep=""))
ggsave(file=paste("plots/GOF/GOF_Venn.png",sep=""))

###############################
# dropouts with no GoFs

dropout_noGOF <- orcollapse3info %>%
  filter(globalfit14 < -2.5 & !is.na(globalfit14)) %>%
  dplyr::select(ID) %>%
  filter(!(ID %in% dropout_mutants_GOF$ID))

dropout_noGOF <- perfects_tree %>%
  inner_join(dropout_noGOF,by="ID")
  
mean(dropout_noGOF$numprunedBCs)

dropout_GOF <- orcollapse3info %>%
  filter(globalfit14 < -2.5 & !is.na(globalfit14)) %>%
  dplyr::select(ID) %>%
  filter(ID %in% dropout_mutants_GOF$ID)

dropout_GOF <- perfects_tree %>%
  inner_join(dropout_GOF,by="ID")

mean(dropout_GOF$numprunedBCs)

write.table(dropout_GOF %>%
              dplyr::select(ID,PctIdentEcoli,TaxID,Source,Taxa1,Taxa2,Taxa3,Taxa4,globalfit14),
            file = "output/dropout_orgs_with_GOF.csv", sep = ",", row.names = F,quote=F,col.names = T)

#################

source("GOFsignificant.R")