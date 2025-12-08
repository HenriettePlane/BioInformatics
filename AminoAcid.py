# main file for our amino acid analyses
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

# Function 1: TRANSLATE Reads in a DNA sequence fasta file,
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
# def dna2prot(f1, f2="translated_fasta.txt"):
# """f1 is the name of the input fasta text file; f2 is the default name
# of the translated protein fasta file"""
#     return


#Function 2: Reads in a dictionary of protein sequences (see dna2prot) and
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

#Function 3: Reads in a dictionary protein sequences, and finds
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

# def motif_finder(prot_dict, motif, f2="motifs.txt"):
#     return

#Function 4: Reads in a dictionary protein sequences, scans all overlapping
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
