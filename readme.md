## Broad Mutational Scanning (BMS) analysis script

This set of scripts provides an example broad mutational scanning analysis where data from many related protein homologs and their mutants is combined and collapsed for further analysis and visualization. This script demonstrates the approach used in: _Plesa C, Sidore AM, Lubock N, Zhang D, Kosuri S. Multiplexed Gene Synthesis in Emulsions for Exploring Protein Functional Landscapes. Science 359, 343–347 (2018), DOI: 10.1126/science.aao5167_. 

For updated versions of this code please check [https://github.com/KosuriLab](https://github.com/KosuriLab)

#### Requirements:
* anaconda
 * python 3.5 environment
  * biopython
 * R 3.3.2
* t_coffee (for alignments)

#### Data input

Download the PPAT dataset from: 
[https://doi.org/10.6084/m9.figshare.5990896.v1](https://doi.org/10.6084/m9.figshare.5990896.v1)

PPATdata.RData contains all data necessary for the analysis. Here is a brief descirption of each dataframe:

| Variable  |  Description |
|---|---|
| *mutants*  |  Each row here in this dataframe is a unique mutant with a corresponding fitness. |
| orcollapse3  | Similar to *mutants* but each row is a mutant within a distance of 3 a.a., without the requirement that the mutant appears in both replicates.  |
| orcollapse3info  | For each homolog the median fitness of all mutants within distance of 3 a.a.  |
| perfects  | The fitness of homologs using only the perfect a.a. sequence (synonymous mutations allowed).  |
| perfects_tree  | The information for all designed homologs.  |

Here is a brief descirption of the most important variables in each dataframe:

***_mutants and orcollapse3_***

| Variable  |  Description |
|---|---|
| mutID  |  The NCBI Accession ID of the closest homolog plus the mutations annotated as underscore + old amino acid + residue position + new amino acid. In case of more than 5 mutations the SHA256 hash of the mutant's sequence is used. |
| IDalign  | The NCBI Accession ID of the closest aligned homolog.  |
| mutations  | The number of a.a. mutations from the closest homolog.  |
| seq  | The a.a. sequence without the starting Met.  |
| pct_ident  | The percentage a.a. identity relative to the closest homolog.  |
| globalfit14  | The fitness of this mutant determined using both replicates.   |
| fitSA14  | The fitness of this mutant in replicate A.  |
| fitSB14  | The fitness of this mutant in replicate B.  |
| numprunedBCs  | The number of barcodes for this mutant (after low-read BCs are pruned).  |


***_perfects\_tree_***

| Variable  |  Description |
|---|---|
| ID  | The NCBI Accession ID of the closest aligned homolog.  |
| PctIdentEcoli  | The percentage a.a. identity relative to E. coli PPAT.  |
| TaxID  | NCBI Taxonomy ID for source organism.  |
| Source  | The name of the source organism.  |
| Taxa1, Taxa2,...  | Taxonomy levels with 1 as top.   |
| numBCs_all  | The number of barcodes for this homolog before pruning.  |

#### Procedure

The analysis can be carried out by running the shell script _analysis.sh_. This file calls a number of scripts:

1. **_parse\_dssp.py_**  
This is used to generate Relative-Solvent-Accessibility and Secondary-Structure information files from the PPAT dssp file. This uses the Jesse Bloom's dssp module from mapmuts.

2. **_prepFASTAforALN.R_**  
This is used to generate csv files with homolog sequences which will be used in the alignments.

3. **_gen\_tree\_files.py_**  
This will generate FASTA files from the csv files (for the alignments) and add in the E. coli PPAT sequence.

4. **_map\_aligned\_residues.py_**  
This will parse the alignemtns and generate tables of which homolog's residue corresponds to which position in the alignment table so that co-aligned residues can be determined.

5. **_BMS\_paper.R_**  
This will create a table of all residues and their fitness for all complementing homologs and their mutants. This data is then collapsed onto a reference sequence, in this case E. coli. Once complete it will call **_GOF\_paper.R_** to carry out the analysis of gain-of-funtion (GoF) mutants for low-fitness homologs. Finally **_GOFsignificant.R_** will compare the significant residues found in the GoF analysis to those same positions in the BMS data.

