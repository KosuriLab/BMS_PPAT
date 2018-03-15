import time
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

##################################
#FUNCTIONS:
def getFastaSeqs(filename):
    fastaseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastaseqs.append(seqrec)
    handle.close()
    return fastaseqs

##################################
#INPUTS:

base_path = ""
trees_path_prefix = base_path+"alignments/"

#clustal format alignment file
align_file_in = [trees_path_prefix+"PPAT_aa_all_cleanMSA_wEcoli_tcoffee.aln",
                trees_path_prefix+"PPAT_aa_dropouts_wEcoli_tcoffee.aln",
                trees_path_prefix+"PPAT_ortho_min_fitness_-1_wEcoli_tcoffee.aln"]

#number of seqs in each alignment file
num_samples_in_file = [1142,56,498]

##################################
#OUTPUTS:

msa_map_out_path = [trees_path_prefix+"MSAmap_all/",
                    trees_path_prefix+"MSAdropouts/",
                    trees_path_prefix+"MSA_BMS_1/"]


for alni in range(len(align_file_in)):
    
    ##################################
    #VARIABLES:
    
    #ID as key, align as value
    align_dict = dict()
    
    #num_samples = 454
    num_samples = num_samples_in_file[alni]
    
    #pos key, consensus pos val
    IDaadictlist = [dict() for x in range(num_samples)]
    
    IDtoindexdict = dict()
    indexdtoIDict = dict()
    
    ##################################
    #CODE:
    
    line_count = 0
    #loop over all alignments:
    print(align_file_in[alni])
    for line in open(align_file_in[alni]):
        #skip header
        if line_count > 1:
            listWords = line.split('    ')
            ID = listWords[0]
            align = line[16:].rstrip()
            if ID.strip() != "":
                align_dict[ID] = align_dict.get(ID, "") + align
        line_count += 1
    
    counter = 0
    for ID in align_dict:
        #print(ID)
        #print(align_dict[ID])
        IDtoindexdict[ID] = counter
        indexdtoIDict[counter]=ID
        align = align_dict[ID]
        
        aacounter = 1
        
        
        for i in range(len(align)):
            if align[i] != "-":
                
                #print(str(counter)+" "+str(aacounter))
                IDaadictlist[counter][aacounter]=i+1
                aacounter += 1
        counter += 1
        
    for i in range(len(IDaadictlist)):
        #print(indexdtoIDict[i])
        csvfile = open(str(msa_map_out_path[alni]+indexdtoIDict[i]+".csv"), 'w')
        fieldnames = ['orth_aanum','msa_aanum']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for j in IDaadictlist[i]:
            #print(str(j)+" "+str(IDaadictlist[i][j]))
            #save all data:
            writer.writerow({'orth_aanum':str(j),'msa_aanum':str(IDaadictlist[i][j])})
        csvfile.close()
    
    
    #######################################################
