#Prepare FASTA files for alignments

library(dplyr)

#set current directory to current file directory
#setwd("/Users/calin/GoogleDrive/postdoc/PPAT/BMScode")

#load data file with:
#load data file with:
#perfects - info on all homologs
#mutants - info on all homolog mutants

#orcollapse3info
#orcollapse3
#variables from BMS analysis

load("PPATdata.RData")

#minimum fitness to include the mutation
BMS_min_fitness <- -1

#grab the perfects (homologs) with fitness above minimum fitness
BMS_ortho_ID_list <- perfects %>%
  filter(globalfit14 > BMS_min_fitness) %>%
  dplyr::select(ID)

#save csv file with homolog IDs
write.table(BMS_ortho_ID_list, file = paste("data/BMS_ortho_ID_list_min_fitness_",as.character(BMS_min_fitness),".csv",sep=""), sep = ",", row.names = F,quote=F,col.names = F)

#get low-fitness homologs (dropouts)
dropout_mutants_GOF <- orcollapse3info %>%
  filter(globalfit14 < -2.5) %>%
  dplyr::select(ID, numMutants)

#make sure none of the collapsed info has corresponding perfects with high fitness
dropout_fake_index = numeric()
for (i in 1:nrow(dropout_mutants_GOF)){
  if (dropout_mutants_GOF$ID[i] %in% perfects$ID){
    if (perfects$globalfit14[which(perfects$ID==dropout_mutants_GOF$ID[i])] > -2.5){#was -1
      dropout_fake_index[length(dropout_fake_index)+1] <- i
    }
  }
}
if (length(dropout_fake_index) > 0){
  dropout_mutants_GOF <- dropout_mutants_GOF[-dropout_fake_index,]
}
rm(dropout_fake_index)

#remove negative controls and seq which don't align well
#bad MSA or missing TxGH motif:
bad_MSA_IDs = c("WP_013999525",
                "KEH27566",
                "KHN36098",
                "XP_003381324",
                "EMT24624",
                "KMZ74788",
                "XP_008352456",
                "KHG26015",
                "WP_051012154",
                "AGN24346",
                "WP_026791263",
                "WP_042225014",
                "EQA34904",
                "WP_052606875",
                "WP_016361007",
                "WP_051506186",
                "WP_051282410",
                "WP_050330521",
                "EAU53429",
                "WP_050806102",
                "WP_051264928",
                "WP_051412112",
                "AHB12988")

#total number of mutants for dropouts
deltemp <- mutants %>%
  dplyr::rename(ID=IDalign) %>%
  filter(mutations < 6 &
           mutations > 0 &
           numprunedBCs >0 &
           !(ID %in% bad_MSA_IDs) &
           (ID %in% dropout_mutants_GOF$ID))

rm(deltemp)

#get dropout mutants and select those with positive fitness
#mutants dataf is generated in R_plot_all_mutants.R
dropout_mutants_GOF <- mutants %>%
  dplyr::rename(ID=IDalign) %>%
  filter(globalfit14 > 0 &
           mutations < 6 &
           mutations > 0 &
           numprunedBCs >0 &
           !(ID %in% bad_MSA_IDs)) %>%
  inner_join(dropout_mutants_GOF, by="ID") %>%
  ungroup()

#########
#notes
#using perfects low-fitness homologs with >=5 BarCodes:
#268 mutants, with BCs >0, across 26 homolog IDs
#25 mutants, with BCs >1, across 13 homolog IDs
#########
#using perfects low-fitness homologs with  >1 BarCodes:
#552 mutants, with BCs >0, across 58 homolog IDs
#########
#using collapse3:
#385 mutants, with BCs >0, across 55 IDs
#########

#number of dropouts with GOF mutants
length(unique(dropout_mutants_GOF$ID))

#mutants with no perfect info
#dropout_mutants_GOF_temp <- dropout_mutants_GOF[which(dropout_mutants_GOF$ID %in% unique(dropout_mutants_GOF$ID[which(!(dropout_mutants_GOF$ID %in% perfects$ID))])),]
#BCs_temp <- BCs %>% filter(IDalign %in% dropout_mutants_GOF_temp$ID)

write.table(dropout_mutants_GOF,
            file = "output/dropout_mutants_GOF.csv",
            sep = ",",
            row.names = F,
            quote=F,
            col.names = T)

write.table(unique(dropout_mutants_GOF$ID),
            file = "data/dropout_mutants_GOF_ID_for_MSA.csv",
            sep = ",",
            row.names = F,
            quote=F,
            col.names = F)

save(dropout_mutants_GOF, file = "dropout_mutants_GOF.RData")
