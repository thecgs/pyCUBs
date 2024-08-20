# pyCUBs: Codon Usage Bias，CUB base on python

## List of features that have been developed and are planned to be developed

[Done] Observed number of occurrences of codon (Obs) analysis


[Done] Stat frequency of cusp software. cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp


[Done] Relative synonymous codon usage (RSCU) analysis


[Done] Draw codon barplot


[Done] Parity rule 2 (PR2) analysis


[Done] Draw Neutrality curve


[Done] Draw ENC plot


[Done] Cups sortware anslysis


[Done part] codonW sortware anslysis


[Plan] CAI analysis


[Plan] Codon Bias Index (CBI)


[Plan] Frequency of OPtimal codons


[Plan] Corresponding analysis (COA)


## Dependencies and test environment

scipy >= v1.11.4


numpy >= v1.26.3


seaborn >= v0.13.1


matplotlib >= v3.8.2


python >= v3.11.5

## test data

Min.mt.fasta data from https://plantscience.cn/cn/article/doi/10.11913/PSJ.2095-0837.2022.20229


```python
import sys
sys.path.append('/mnt/nfs1/jupyter/pyCUBs/pyCUBs/') # import library path
import pycubcore

# Switch the test work path and import test data
import os
os.chdir('/mnt/nfs1/jupyter/pyCUBs/')
inputfile = "./test_data/Min.mt.fasta"

# Module Introduction
"""
├── codontables.py  # Codon table module
├── fastaio.py      # Fasta IO module
└── pycubcore.py    # CUB Core computing module
"""

import codontables
print("\nTable of available genetic codons:\n",codontables.CodonTables())
import fastaio
help(fastaio.FastaIO)

# Some information about the library
print(pycubcore.__author__)
print("pycubcore function: ", pycubcore.__all__)
print("version: ",pycubcore.__version__)
```

    
    Table of available genetic codons:
     Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
    
    Translate Tables/Genetic Codes:
     1: Standard
     2: Vertebrate Mitochondrial
     3: YeastMitochondrial
     4: Mold Mitochondrial, Protozoan Mitochondrial, Coelenterate Mitochondrial, Mycoplasma, Spiroplasma
     5: Invertebrate Mitochondrial
     6: Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear
     9: Echinoderm Mitochondrial, Flatworm Mitochondrial
    10: Euplotid Nuclear
    11: Bacterial, Archaeal, Plant Plastid
    12: Alternative Yeast Nuclear
    13: Ascidian Mitochondrial
    14: Alternative Flatworm Mitochondrial
    16: Chlorophycean Mitochondrial
    21: Trematode Mitochondrial
    22: Scenedesmus obliquus Mitochondrial
    23: Thraustochytrium Mitochondrial
    24: Rhabdopleuridae Mitochondrial
    25: Candidate Division SR1, Gracilibacteria
    26: Pachysolen tannophilus Nuclear
    27: Karyorelict Nuclear
    28: Condylostoma Nuclear
    29: Mesodinium Nuclear
    30: Peritrich Nuclear
    31: Blastocrithidia Nuclear
    33: Cephalodiscidae Mitochondrial UAA-Tyr
    
    Help on function FastaIO in module fastaio:
    
    FastaIO(file: str)
        Read a fasta file path and support compressed files ending in ".gz", 
        or accept a handle of "_io.TextIOWrapper" class.
    
    Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2024/05/17
    pycubcore function:  ['GetObs', 'GetFranction', 'GetFrequency', 'GetRSCU', 'DrawCodonBarplot', 'GetCusp', 'GetcodonW', 'NPA', 'DrawNPA', 'GetNC', 'GetGC3s', 'ENC', 'DrawENC', 'Find4Dtv', 'GetPR2', 'PR2', 'DrawPR2']
    version:  v0.01



```python
## example1: Effective number of codons (ENC) analysis
#help(pycubcore.ENC)
ENCResult = pycubcore.ENC(inputfile, genetic_codes=16) #Genetic_Codes=16, Selected genetic codon table 16 Plant chloroplasts
#help(pycubcore.DrawENC)
pycubcore.DrawENC(ENCResult, show_gene_name=False)
```


    
![png](output_2_0.png)
    



```python
# example2: Parity rule 2 (PR2) analysis.
# help(pycubcore.PR2)
PR2Result = pycubcore.PR2(inputfile, genetic_codes=16)
# help(pycubcore.DrawPR2)
pycubcore.DrawPR2(PR2Result, show_gene_name=True) #Show gene name
```


    
![png](output_3_0.png)
    



```python
# example3: Neutral plot analysis.
NPAResult = pycubcore.NPA(inputfile, genetic_codes=16)
pycubcore.DrawNPA(NPAResult, show_gene_name=True)
```


    
![png](output_4_0.png)
    



```python
# example4: compute Obs, Franction, Frequency, RSCU
GeneName, Seqence = next(fastaio.FastaIO(inputfile))
print("input seqence:\n")
print(">"+GeneName+"\n", Seqence, sep="")

# Calculate
Obs = pycubcore.GetObs(seqences=Seqence, genetic_codes=1) #The seqences parameter can be a string of gene sequence or a list. For example, a list of all CDS sequences of a chloroplast, so as to calculate the Obs in the whole organism
Franction=pycubcore.GetFranction(Obs)
Frequency=pycubcore.GetFrequency(Obs)
RSCU=pycubcore.GetRSCU(Obs)

# Visualization
print("\nVisualization result：\n")
pycubcore.DrawCodonBarplot(Obs, ylabel="Number")
pycubcore.DrawCodonBarplot(Franction, ylabel="Franction")
pycubcore.DrawCodonBarplot(Frequency, ylabel="Frequency")
pycubcore.DrawCodonBarplot(RSCU, ylabel="RSCU")
```

    input seqence:
    
    >psbA
    ATGACTGCAATTTTAGAGAGACGCGAAAGCGAAAGCCTATGGGGTCGCTTCTGTAACTGGATAACCAGCACTGAGAACCGTCTTTACATTGGATGGTTTGGTGTTTTGATGATCCCTACCTTATTGACCGCAACTTCTGTATTTATTATCGCCTTCATTGCTGCTCCTCCAGTAGATATTGATGGTATTCGTGAACCTGTTTCTGGGTCTCTACTTTACGGAAACAATATTATCTCTGGTGCCATTATTCCTACTTCTGCAGCTATAGGATTGCACTTTTACCCGATATGGGAAGCGGCATCCGTTGATGAATGGTTATACAATGGTGGTCCTTATGAATTGATTGTTCTACACTTCTTACTTGGTGTAGCTTCTTACATGGGTCGTGAGTGGGAACTAAGTTTCCGTCTGGGTATGCGCCCTTGGATTGCTGTTGCATATTCAGCTCCTGTTGCAGCTGCAACTGCTGTTTTCTTGATCTACCCAATCGGTCAAGGAAGCTTCTCTGATGGTATGCCCCTAGGAATCTCTGGTACTTTCAACTTCATGATTGTATTCCAGGCTGAGCACAACATTCTTATGCACCCATTTCACATGTTAGGTGTGGCTGGTGTATTCGGCGGCTCCCTATTCAGTGCTATGCATGGTTCCTTGGTAACTTCAAGTTTGATCAGGGAAACCACTGAAAATGAATCTGCTAATGAAGGTTACAGATTCGGTCAAGAGGAAGAAACTTATAATATCGTAGCTGCTCATGGTTATTTTGGCCGATTGATCTTCCAATATGCTAGTTTCAACAATTCTCGTTCTTTACATTTCTTCCTAGCTGCTTGGCCTGTAGTAGGTATCTGGTTCACTGCTTTAGGTATTAGTACCATGGCTTTCAACCTAAATGGTTTCAATTTCAACCAATCCGTAGTTGACAGTCAAGGTCGTGTAATTAACACTTGGGCTGATATCATCAACCGTGCTAACCTTGGTATGGAAGTTATGCATGAACGTAATGCTCACAACTTCCCTCTAGACCTAGCTGCTGTTGAAGCTCCATCCATAAATGGATAA
    
    Visualization result：
    



    
![png](output_5_1.png)
    



    
![png](output_5_2.png)
    



    
![png](output_5_3.png)
    



    
![png](output_5_4.png)
    



```python
# example5: Get similar results to Cusp software
# Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp

GeneName, Seqence = next(fastaio.FastaIO(inputfile))
print("Input seqence: \n")
print(">"+GeneName+"\n", Seqence, sep="")
print("\n")

Obs = pycubcore.GetObs(seqences=Seqence, genetic_codes=1) 
print("Cusp result: \n")
print(pycubcore.GetCusp(Obs, human_format=True)) #human_format=True, Human-readable output, otherwise machine-readable.
```

    Input seqence: 
    
    >psbA
    ATGACTGCAATTTTAGAGAGACGCGAAAGCGAAAGCCTATGGGGTCGCTTCTGTAACTGGATAACCAGCACTGAGAACCGTCTTTACATTGGATGGTTTGGTGTTTTGATGATCCCTACCTTATTGACCGCAACTTCTGTATTTATTATCGCCTTCATTGCTGCTCCTCCAGTAGATATTGATGGTATTCGTGAACCTGTTTCTGGGTCTCTACTTTACGGAAACAATATTATCTCTGGTGCCATTATTCCTACTTCTGCAGCTATAGGATTGCACTTTTACCCGATATGGGAAGCGGCATCCGTTGATGAATGGTTATACAATGGTGGTCCTTATGAATTGATTGTTCTACACTTCTTACTTGGTGTAGCTTCTTACATGGGTCGTGAGTGGGAACTAAGTTTCCGTCTGGGTATGCGCCCTTGGATTGCTGTTGCATATTCAGCTCCTGTTGCAGCTGCAACTGCTGTTTTCTTGATCTACCCAATCGGTCAAGGAAGCTTCTCTGATGGTATGCCCCTAGGAATCTCTGGTACTTTCAACTTCATGATTGTATTCCAGGCTGAGCACAACATTCTTATGCACCCATTTCACATGTTAGGTGTGGCTGGTGTATTCGGCGGCTCCCTATTCAGTGCTATGCATGGTTCCTTGGTAACTTCAAGTTTGATCAGGGAAACCACTGAAAATGAATCTGCTAATGAAGGTTACAGATTCGGTCAAGAGGAAGAAACTTATAATATCGTAGCTGCTCATGGTTATTTTGGCCGATTGATCTTCCAATATGCTAGTTTCAACAATTCTCGTTCTTTACATTTCTTCCTAGCTGCTTGGCCTGTAGTAGGTATCTGGTTCACTGCTTTAGGTATTAGTACCATGGCTTTCAACCTAAATGGTTTCAATTTCAACCAATCCGTAGTTGACAGTCAAGGTCGTGTAATTAACACTTGGGCTGATATCATCAACCGTGCTAACCTTGGTATGGAAGTTATGCATGAACGTAATGCTCACAACTTCCCTCTAGACCTAGCTGCTGTTGAAGCTCCATCCATAAATGGATAA
    
    
    Cusp result: 
    
    #Coding GC 42.84%
    #1st letter GC 50.0%
    #2nd letter GC 43.22%
    #3rd letter GC 35.31%
    
    #Codon AA Fraction Frequency Number
    TTT	Phe	0.192	14.124	5
    TTC	Phe	0.808	59.322	21
    TCT	Ser	0.393	31.073	11
    TCC	Ser	0.179	14.124	5
    TCA	Ser	0.071	5.65	2
    TCG	Ser	0.0	0.0	0
    AGT	Ser	0.214	16.949	6
    AGC	Ser	0.143	11.299	4
    TAT	Tyr	0.417	14.124	5
    TAC	Tyr	0.583	19.774	7
    TGT	Cys	1.0	2.825	1
    TGC	Cys	0.0	0.0	0
    TTA	Leu	0.226	19.774	7
    TTG	Leu	0.258	22.599	8
    CTT	Leu	0.161	14.124	5
    CTC	Leu	0.0	0.0	0
    CTA	Leu	0.323	28.249	10
    CTG	Leu	0.032	2.825	1
    TAA	*	1.0	2.825	1
    TGA	*	0.0	0.0	0
    TAG	*	0.0	0.0	0
    TGG	Trp	1.0	28.249	10
    CCT	Pro	0.6	25.424	9
    CCC	Pro	0.067	2.825	1
    CCA	Pro	0.267	11.299	4
    CCG	Pro	0.067	2.825	1
    CAT	His	0.4	11.299	4
    CAC	His	0.6	16.949	6
    CGT	Arg	0.533	22.599	8
    CGC	Arg	0.2	8.475	3
    CGA	Arg	0.067	2.825	1
    CGG	Arg	0.0	0.0	0
    AGA	Arg	0.133	5.65	2
    AGG	Arg	0.067	2.825	1
    CAA	Gln	0.833	14.124	5
    CAG	Gln	0.167	2.825	1
    ATT	Ile	0.484	42.373	15
    ATC	Ile	0.387	33.898	12
    ATA	Ile	0.129	11.299	4
    ACT	Thr	0.688	31.073	11
    ACC	Thr	0.312	14.124	5
    ACA	Thr	0.0	0.0	0
    ACG	Thr	0.0	0.0	0
    AAT	Asn	0.455	28.249	10
    AAC	Asn	0.545	33.898	12
    ATG	Met	1.0	33.898	12
    GTT	Val	0.455	28.249	10
    GTC	Val	0.0	0.0	0
    GTA	Val	0.5	31.073	11
    GTG	Val	0.045	2.825	1
    GCT	Ala	0.714	70.621	25
    GCC	Ala	0.057	5.65	2
    GCA	Ala	0.2	19.774	7
    GCG	Ala	0.029	2.825	1
    GAT	Asp	0.714	14.124	5
    GAC	Asp	0.286	5.65	2
    GGT	Gly	0.697	64.972	23
    GGC	Gly	0.091	8.475	3
    GGA	Gly	0.182	16.949	6
    GGG	Gly	0.03	2.825	1
    GAA	Glu	0.762	45.198	16
    GAG	Glu	0.238	14.124	5
    



```python
# example6: Get similar results to codonW software.
GeneName, Seqence = next(fastaio.FastaIO(inputfile))
print("Input seqence: \n")
print(">"+GeneName+"\n", Seqence, sep="")
print("\n")

Obs = pycubcore.GetObs(seqences=Seqence, genetic_codes=16)
print("codonW result: \n")
print(pycubcore.GetcodonW(Obs, human_format=True)) #human_format=True, Human-readable output, otherwise machine-readable.
```

    Input seqence: 
    
    >psbA
    ATGACTGCAATTTTAGAGAGACGCGAAAGCGAAAGCCTATGGGGTCGCTTCTGTAACTGGATAACCAGCACTGAGAACCGTCTTTACATTGGATGGTTTGGTGTTTTGATGATCCCTACCTTATTGACCGCAACTTCTGTATTTATTATCGCCTTCATTGCTGCTCCTCCAGTAGATATTGATGGTATTCGTGAACCTGTTTCTGGGTCTCTACTTTACGGAAACAATATTATCTCTGGTGCCATTATTCCTACTTCTGCAGCTATAGGATTGCACTTTTACCCGATATGGGAAGCGGCATCCGTTGATGAATGGTTATACAATGGTGGTCCTTATGAATTGATTGTTCTACACTTCTTACTTGGTGTAGCTTCTTACATGGGTCGTGAGTGGGAACTAAGTTTCCGTCTGGGTATGCGCCCTTGGATTGCTGTTGCATATTCAGCTCCTGTTGCAGCTGCAACTGCTGTTTTCTTGATCTACCCAATCGGTCAAGGAAGCTTCTCTGATGGTATGCCCCTAGGAATCTCTGGTACTTTCAACTTCATGATTGTATTCCAGGCTGAGCACAACATTCTTATGCACCCATTTCACATGTTAGGTGTGGCTGGTGTATTCGGCGGCTCCCTATTCAGTGCTATGCATGGTTCCTTGGTAACTTCAAGTTTGATCAGGGAAACCACTGAAAATGAATCTGCTAATGAAGGTTACAGATTCGGTCAAGAGGAAGAAACTTATAATATCGTAGCTGCTCATGGTTATTTTGGCCGATTGATCTTCCAATATGCTAGTTTCAACAATTCTCGTTCTTTACATTTCTTCCTAGCTGCTTGGCCTGTAGTAGGTATCTGGTTCACTGCTTTAGGTATTAGTACCATGGCTTTCAACCTAAATGGTTTCAATTTCAACCAATCCGTAGTTGACAGTCAAGGTCGTGTAATTAACACTTGGGCTGATATCATCAACCGTGCTAACCTTGGTATGGAAGTTATGCATGAACGTAATGCTCACAACTTCCCTCTAGACCTAGCTGCTGTTGAAGCTCCATCCATAAATGGATAA
    
    
    codonW result: 
    
    T3s    C3s    A3s    G3s    Nc     GC3s   GC     L_sym  L_aa  
    0.5033 0.2730 0.2964 0.0901 43.01  0.311  0.430  331    353   
    

