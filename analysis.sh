
echo "BMS analysis scripts, Calin Plesa.\n"

#conda create --name py3k --file py3k_env.txt

echo "Setting python 3 environment in anaconda.\n"
#source activate py3k

echo "Generate RSA and SS files from PPAT dssp file.\n"
#this uses the Jesse Bloom's dssp module from mapmuts
#https://github.com/jbloom/mapmuts/blob/master/src/dssp.py
python parse_dssp.py

echo "Make csv files needed for alignments.\n"
rscript prepFASTAforALN.R

echo "Make FASTA files using the csv files.\n"
python gen_tree_files.py

echo "Generate alignments using t_coffee.\n"
#./tree_gen.sh

echo "Parsing alignments.\n"
python map_aligned_residues.py

echo "BMS analysis in R.\n"
#carry out the BMS analysis (complementing homologs and their mutants)
rscript --verbose BMS_paper.R
#this will also call GOF_paper.R to carry out the GoF analysis (low-fitness homologs and their GoF mutants)
#this will also call GOFsignificant.R to compare the significant residues found in the GoF analysis to those same positions in the BMS data
