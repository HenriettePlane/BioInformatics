# main file for our amino acid analyses
import re

# FINAL PROJECT
standard_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
    "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
    "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
    "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
    "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
    "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

hscale={"R":-4.5,"K":-3.9,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"H":-3.2,
    "P":-1.6,"Y":-1.3,"W":-0.9,"S":-0.8,"T":-0.7,"G":-0.4,"A":1.8,"M":1.9,
    "C":2.5,"F":2.8,"L":3.8,"V":4.2,"I":4.5}

#### FUNCTION 1: ##############################################################
# TRANSLATE Reads in a DNA sequence fasta file, 
# and outputs a Fasta file of protein translations and it
# should be able to do standard (eukaryotic) translation only.
# The first function reads in a fasta file of DNA sequences and:
# (1) Writes a file with the protein sequences (default name ->
# translated_fasta.txt)
# (2) Returns a dictionary of the sequences with key = title, sequence = value
# If a fasta like called "dna.txt" looks like this:
# >seq1
# ATGGAG
# >seq2
# ATGAGA
#
# Then the dictionary will be: {'seq1':'ME','seq2':'MR'}
#
## TEST THIS FUNCTION WITH: fasta_aa_version.txt

### Solution for function 1 ###
# Step 1: Write a simple helper function so that the individual functionality is easier to test
# findprots is a simple function that takes in a DNA anti sense string, cleans it, breaks it into codons,
# translates it into RNA by changing T's into U's and then uses the standard_code dictionairy to find
# the amino acid that matches the codon. It returns a string of amino acids.
def find_aa(dna,aa_dict=standard_code):
    dna = dna.upper().strip() # clean the DNA
    dna = re.sub('T','U',dna) # change T's into U's
    # create a string of amino acids
    aa_string = ''
    for i in range(0,len(dna),3):
        if i+3 <= len(dna):
            if aa_dict[dna[i:i+3]] == '*':
                break
            else:
                aa_string += aa_dict[dna[i:i+3]]
        else: pass
    return aa_string

# testing the function
print('-------- FUNCTION 1 TESTS ------------------')
print(find_aa('ATGTCAAAGT'))

# Step 2: write the main function that calls the helper function.
# function dna2prot has two arguments: f1 is the name of the input fasta text file; f2 is the default name
# of the translated protein fasta file
def dna2prot(f1, f2="translated_fasta.txt"):
    # read lines of dna fasta file into a list
    fasta = open(f1,'r')
    dna_list = []
    for line in fasta:
        dna_list.append(line.strip())
    # create a dictionary of proteins and write the titles and protein strings to a fasta file
    prot_dict = {}
    fout = open(f2,'w')
    for i in range(0,len(dna_list),2):
        key = dna_list[i]
        value = find_aa(dna_list[i+1])
        prot_dict[key] = value
        fout.write(key + '\n' + value + '\n')
    return(prot_dict)

# testing the complete function
test = dna2prot('fasta_aa_version.txt')
print(test.items())

######## FUNCTION 2:#######################################################################
# Reads in a dictionary of protein sequences (see dna2prot) and
# (1) Writes a table of amino acids counts for each protein
# to a tab-separated text file.
# (2) Returns a dictionary of the data that looks like this:
# {'seq1': [(A,7),(C,3), ... (Y,1)], 'seq2':[(A,2),(C,0), ... (Y,2)]}
#
# All the amino acids should be in the dictionary and if there are no examples
# of a particular amino acid, the value would be zero
#The output file should look like this (with different values):
#(Put the AAs alphabetical too.)
# """
# AminoAcid A C D ... Total
# Seq1 7 3 9 100
# Seq2 ...."""
# def aa_counts(prot_dict,f2="aatable.txt"):
#     return

#Function 2 Solution#

def aa_counts(prot_dict, f2="aatable.txt"):
    
    f2=open("aatable.txt", "w")
    
    results={}
    amino_acids=hscale.keys()
    amino_acids=sorted(amino_acids)
    header=f"AminoAcid {' '.join(amino_acids)} TOTAL \n"
    f2.write(header)
    
    j=0
    

    keys_list=list(prot_dict.keys())
    values_list=list(prot_dict.values())
        
    for prot_seq in values_list:
        counts_dict = {key: 0 for key in amino_acids}
        for aa in prot_seq:
            if aa in counts_dict:
                counts_dict[aa] +=1
                 
        #print(counts_dict)
        counts_tuple=list(counts_dict.items())
        #print(counts_tuple)
        results[keys_list[j]] = counts_tuple
        
        j=j+1
            
    for key, value in results.items():    
        TOTAL=0
        second_values=[]
        results[key]= counts_tuple
        for items in value:
            TOTAL=TOTAL+int(items[1])
            second_values.append(str(items[1]))
          
        line = f"{key}      {' '.join(second_values)}  {TOTAL}\n"
        f2.write(line)
        
            
    f2.close()
    return results
 
   
#result=aa_counts(protein_dict)
#print(result)

######### FUNCTION 3: ####################################################################
# Reads in a dictionary protein sequences, and finds
# all instances of the motif.
# (1) Writes a table of motifs and hits for each sequence (see below)
# to a tab-separated text file.
# (2) Returns a dictionary with a list of all the motifs. For example,
# the motif "V[A-Z]R[ML]" might return
# {'seq1': ['VVRL','VWRM'] 'seq2':[], 'seq3':['VFRL','VQRM']}
# The output file should look like this:
# """
# SeqName Motif Hits
# Seq1 MN[A-Z] 2
# Seq2 MN[A-Z] 0
# Seq3 MN[A-Z] 1
# """


### Solution for Function 3 ####
def motif_finder(prot_dict, motif, f2="motifs.txt"):
    motif_dict = {}
    file = open(f2,'w')
    file.write('SeqName' + '\t' + 'Motif' + '\t' +  'Hits\n')
    for key,value in prot_dict.items():
       m = re.findall(motif,value)
       hits = len(m)
       line = (key + '\t' + motif + '\t' + str(hits) + '\n')
       print(line)
       file.write(line)
       motif_dict[key] = m
    return motif_dict

test3 = motif_finder(test,'A.')
print('-------- FUNCTION 3 TESTS ----------------------')
print(test3.items())


############ FUNCTION 4: ##########################################################################
# Reads in a dictionary protein sequences, scans all overlapping
# windows and calculates averge hydrophobicity
# (1) Writes a table of the avg scores (see below)
# (2) Returns a dictionary with a list hydrophobicity scores
# {'seq1': [1.2, 0.5, -2.3] 'seq2':[-0.3, -0.1, 0.5]}
#The output file should look like this:
# """
# Window
# SeqName 1 2 3
# Seq1 1.2 0.5 -2.3
# Seq2 -0.3 -0.1 0.5
# """
# def hydrophobicity_analysis(prot_dict, window=5, f2="hydro.txt"):
#     return

#Function 4 Solution#

def hydrophobicity_analysis(prot_dict, window_size, f2="hydro.txt"):        
    f2=open("hydro.txt", "w")
    f2.write(f'\t Window \nSeqName 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t \n')
    results=prot_dict
    j=0
    
    keys_list=list(prot_dict.keys())
    values_list = list(prot_dict.values())
    #print(values_list)
    for prot_seq in values_list:
        #print(prot_seq)
        windows=[]

        for i in range(0,len(prot_seq),1):
            window=prot_seq[i:i+window_size]
            if len(window)==window_size:
                windows.append(window)
            
        avg_values=[]
        for win in windows:
            hydro_counter=0
            for aa in win: 
                hydro_counter += hscale.get(aa)
            average=hydro_counter/window_size
            avg_values.append(f'{average:.2f}')
            results[keys_list[j]] = avg_values
            
        j=j+1

    for key, value in results.items():    
        line= key + "\t"+"\t".join(value) + "\n"
        f2.write(line)
        
    f2.close()
     
    return results
