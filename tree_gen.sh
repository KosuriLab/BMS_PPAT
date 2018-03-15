#################################
# Make MSAs:

cd ../alignments/

t_coffee PPAT_aa_all_cleanMSA_wEcoli.fasta -n_core=30 -run_name=PPAT_aa_all_cleanMSA_wEcoli_tcoffee

#low-fitness dropouts
t_coffee PPAT_aa_dropouts_wEcoli.fasta -n_core=30 -run_name=PPAT_aa_dropouts_wEcoli_tcoffee

#BMS ortho min fitness -1
t_coffee BMS_ortho_min_fitness_-1.fasta -n_core=30 -run_name=PPAT_ortho_min_fitness_-1_wEcoli_tcoffee
