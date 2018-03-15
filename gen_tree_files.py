import time
import Bio
import csv
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#generate Fasta files

##################################
#FUNCTIONS:

def getFastaSeqs(filename):
    fastaseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        #add start codon
        tempseq = str("M"+seqrec.seq.translate())
        stop_index = tempseq.find('GT*')
        newseq = Seq(str(tempseq[0:stop_index]),generic_protein)
        record = SeqRecord(newseq,id=seqrec.id,description="")
        fastaseqs.append(record)
    handle.close()
    return fastaseqs

##################################
#INPUTS:

base_path = ""
analysis_dir = base_path+'data/'
trees_dir = base_path+"alignments/"

protein_seq_file = 'lib-vars.fasta'

# csv file from with bad MSA alignments:
MSA_info = analysis_dir+'MSA_bad_ID_list.csv'

#csv files with IDs
input_files_with_IDs = [analysis_dir+'dropout_mutants_GOF_ID_for_MSA.csv',
                        analysis_dir+'BMS_ortho_ID_list_min_fitness_-1.csv']

##################################
#OUTPUTS:

#these are fasta files:
output_files = [trees_dir+'PPAT_aa_dropouts_wEcoli.fasta',
                trees_dir+'BMS_ortho_min_fitness_-1.fasta']

FASTA_cleanMSA_wEcoli = trees_dir+'PPAT_aa_all_cleanMSA_wEcoli.fasta'

##################################
#VARIABLES:

num_files = len(input_files_with_IDs)
#list of lists with IDs
synthlist = [[] for i in range(num_files)]

#list of IDs with bad MSA
MSAlist = []

##################################
#CODE:

print("load reference protein seqs and translate")
#this returns dict with ID as key and protein seq as value
buildseqs = getFastaSeqs(protein_seq_file)

for i in range(num_files):
    #use the csv files from R:
    for line in open(input_files_with_IDs[i]):
        synthlist[i].append(line.strip())

for line in open(MSA_info):
    MSAlist.append(line.strip())

#######################################################
#make new FATA files:

FASTAaa_AllSynth = [[] for i in range(num_files)]
FASTAaa_AllGoodMSA = []
synth_counter = [0 for i in range(num_files)]
MSA_counter = 0
#loop over all orthologs:
for const in buildseqs:
    for i in range(num_files):
        if const.id in synthlist[i]:
            FASTAaa_AllSynth[i].append(const)
            synth_counter[i] += 1
    if not (const.id in MSAlist):
        FASTAaa_AllGoodMSA.append(const)
        MSA_counter += 1

# #write FASTA files without E.coli added:
# print("Wrote "+str(synth_counter)+" records to "+FASTA_synth)
# output_handle_synth = open(FASTA_synth, "w")
# SeqIO.write(FASTAaa_AllSynth, output_handle_synth, "fasta")
# output_handle_synth.close()
# 
# print("Wrote "+str(collapse_counter)+" records to "+FASTA_mutant_collapse)
# output_handle_collapse = open(FASTA_mutant_collapse, "w")
# SeqIO.write(FASTAaa_Collapse, output_handle_collapse, "fasta")
# output_handle_collapse.close()

#add in E.coli to buildseqs
newseq = Seq("MQKRAIYPGTFDPITNGHIDIVTRATQMFDHVILAIAASPSKKPMFTLEERVALAQQATAHLGNVEVVGFSDLMANFARNQHATVLIRGLRAVADFEYEMQLAHMNRHLMPELESVFLMPSKEWSFISSSLVKEVARHQGDVTHFLPENVHQALMAKLA",generic_protein)
record = SeqRecord(newseq,id="NP_418091",description="")

#add ecoli to records
for i in range(num_files):
    FASTAaa_AllSynth[i].append(record)
FASTAaa_AllGoodMSA.append(record)

#write FASTA files with E.coli added:
for i in range(num_files):
    print("Wrote "+str(synth_counter[i]+1)+" records to "+output_files[i])
    output_handle_synth = open(output_files[i], "w")
    SeqIO.write(FASTAaa_AllSynth[i], output_handle_synth, "fasta")
    output_handle_synth.close()

print("Wrote "+str(MSA_counter+1)+" records to "+FASTA_cleanMSA_wEcoli)
output_handle_MSA = open(FASTA_cleanMSA_wEcoli, "w")
SeqIO.write(FASTAaa_AllGoodMSA, output_handle_MSA, "fasta")
output_handle_MSA.close()