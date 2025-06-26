#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2025/06/26"
__all__ = ["translate",
           "get_Obs",
           "get_Fraction",
           "get_Frequency", 
           "get_Fop",
           "get_CBI",
           "get_CAI",
           "get_Gravy",
           "get_Aromo",
           "get_RSCU",
           "get_NC", 
           "get_PR2",
           "get_GC123", 
           "get_Relative_Adaptiveness", 
           "get_emboss_cutfile_from_Obs", 
           "get_Obs_from_emboss_cutfile",
           "get_codonw_cai_file_from_Obs",
           "get_optimal_codons_from_codonw_coafile",
           "draw_codon_barplot", 
           "get_cusp_like",
           "get_codonW_like",
           "find_four_codon_AA",
           "TreePlotter",
           "NPA_analysis",
           "ENC_analysis",
           "PR2_analysis", 
           "RSCU_analysis",
           "AAComposition_analysis",
           "StopCodon_analysis"]

__version__ = "v2.00"

import os
import io
import sys
import math
import prince
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as ss
from fastaio import FastaIO
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from collections import defaultdict
from skbio import DistanceMatrix
from skbio.tree import TreeNode, nj, upgma, gme, bme
from scipy.spatial.distance import pdist, squareform
from codontables import CodonTables, Seq3toSeq1, CBI_and_Fop_preset

def translate(cds, genetic_code):
    """
    Description 
    ----------
    CDS translate protein.
    
    Parameters
    ----------
    cds: CDS seqence string.
    genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
    """
    
    protein = ""
    codotable = CodonTables().get(genetic_code, aaseq3=False)
    for i in range(0, len(cds), 3):
        protein += codotable.get(cds[i:i+3], "X")
    return protein
    
def get_Obs(seqences, genetic_code, aaseq3=True):
    """
    Description 
    ----------
    Calculate Observed number of occurrences of codon.
    
    Parameters
    ----------
    seqences: {str, list, file path} a seqence string, or a seqences list, or a fasta or fasta.gz format file path.
    genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
    aaseq3: if value is True, amino acid three letter code.
    """
    
    codontable = CodonTables().get(genetic_code, aaseq3)
    
    def _FastaIO(inputfile):
        for ID, Seq in FastaIO(inputfile):
            yield Seq
    
    if isinstance(seqences, str):
        if os.path.exists(seqences):
            seqences = _FastaIO(seqences)
        else:
            seqences = [seqences.upper()]
            
    elif isinstance(seqences, list) or isinstance(seqences, tuple):
        pass
    else:
        seqences = [str(seqences).upper()]
        
    Codons = {}
    for Seq in seqences:
        Seq = Seq.upper()
        
        # check seqence length
        if len(Seq)%3 == 1:
            Seq = Seq[:-1]
        elif len(Seq)%3 == 2:
            Seq = Seq[:-2]
        else:
            Seq = Seq
            
        for i, j in zip(range(0,len(Seq), 3), range(3, len(Seq)+3, 3)):
            if Seq[i:j] in Codons:
                Codons[Seq[i:j]] += 1
            else:
                Codons.setdefault(Seq[i:j], 1)
    Obs = {}       
    for Codon in codontable:
        AA = codontable[Codon]
        if AA in Obs:
            Obs[AA][Codon] = Codons.get(Codon, 0)
        else:
            Obs.setdefault(AA, {Codon: Codons.get(Codon, 0)})
    return Obs

def get_codonw_caifile_from_Obs(Obs, outfile=None):
    """
    Description
    ----------
    To create a condow cai file from the get_Obs return value of the constructed high expression gene.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    """
    
    RA = get_Relative_Adaptiveness(Obs)
    a1 = {}
    for aa in RA.Relative_Adaptiveness_dict:
        for c in RA.Relative_Adaptiveness_dict[aa]:
            a1.setdefault(c,RA.Relative_Adaptiveness_dict[aa][c])
    a2 =  ["TTT","TCT","TAT","TGT","TTC","TCC","TAC","TGC","TTA","TCA","TAA","TGA","TTG","TCG","TAG","TGG",
           "CTT","CCT","CAT","CGT","CTC","CCC","CAC","CGC","CTA","CCA","CAA","CGA","CTG","CCG","CAG","CGG",
           "ATT","ACT","AAT","AGT","ATC","ACC","AAC","AGC","ATA","ACA","AAA","AGA","ATG","ACG","AAG","AGG",
           "GTT","GCT","GAT","GGT", "GTC","GCC","GAC","GGC", "GTA","GCA","GAA","GGA","GTG","GCG","GAG","GGG"]
    
    if outfile==None:
        out = sys.stdout
    else:
        out = open(outfile, 'w')
    
    for c in a2:
        print(a1[c], file=out)
    if outfile!=None:
        out.close()
    return None

def get_emboss_cutfile_from_Obs(Obs, outfile=None):
    """
    Description
    ----------
    To create a emboss cutfile from the get_Obs return value of the constructed high expression gene.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    """
    
    Frequency = get_Frequency(Obs)
    Fraction = get_Fraction(Obs)
    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, 'w')
    
    print("# The file building with pycubs (https://github.com/thecgs/pycubs) \n# Codon AA Fraction Frequency Number", file=out)
    for aa in Frequency.Frequency_dict:
        for c in Frequency.Frequency_dict[aa]:
            if len(aa) !=1:
                aa1 = Seq3toSeq1[aa]
            else:
                aa1=aa
            print(c, aa1, 
                  round(Fraction.Fraction_dict[aa][c],3), 
                  round(Frequency.Frequency_dict[aa][c],3), 
                  Obs[aa][c], file=out, sep = "    ")
    if outfile !=None:
        out.close()
    return None

def get_Obs_from_emboss_cutfile(file, aaseq3=True):
    """
    Description
    ----------
    Gets the Obs object from the emboss cut file with the same value as the get_Obs function returns.
    
    Parameters
    ----------
    file : a emboss cut file.
    aaseq3: if value is True, amino acid three letter code.
    """
    Seq1toSeq3 = {Seq3toSeq1[k]:k for k in Seq3toSeq1}
    Obs = {}
    with open(file, 'r') as f:
        for l in f:
            if l.strip() != "" and not l.startswith('#'):
                l = l.strip('\n').split()
                if aaseq3==True:
                    aa = Seq1toSeq3[l[1]]
                else:
                    aa = l[1]
                    
                if aa in Obs:
                    Obs[aa][l[0]] = int(l[4])
                else:
                    Obs.setdefault(aa, {l[0]: int(l[4])})
    return Obs

def get_optimal_codons_from_codonw_coafile(file):
    """
    Description 
    ----------
    Obtain a list of optimal codons from the cbi.coa or fop.coa file generated by the codonw software for downstream get_CBI and get_Fop analysis.

    Parameters
    ----------
    file: {cbi.coa, fop.coa} The file  include of optimal codons generated by codonw software.
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    a1 = []
    with open(file, 'r') as f:
        for l in f:
            a1.extend(l.strip('\n').split(','))
    a2 =  ["TTT","TCT","TAT","TGT","TTC","TCC","TAC","TGC","TTA","TCA","TAA","TGA","TTG","TCG","TAG","TGG",
           "CTT","CCT","CAT","CGT","CTC","CCC","CAC","CGC","CTA","CCA","CAA","CGA","CTG","CCG","CAG","CGG",
           "ATT","ACT","AAT","AGT","ATC","ACC","AAC","AGC","ATA","ACA","AAA","AGA","ATG","ACG","AAG","AGG",
           "GTT","GCT","GAT","GGT", "GTC","GCC","GAC","GGC", "GTA","GCA","GAA","GGA","GTG","GCG","GAG","GGG"
          ]
    cs = []
    for c, v in zip(a2, a1):
        if v=='3':
            cs.append(c)
    return sorted(cs)

def get_Fop(Obs, optimal_codons="Escherichia coli"):
    """
    Description 
    ----------
    Frequency of Optimal codons (Fop) (Ikemura 1981). 
    This index, is the ratio of optimal codons to synonymous codons (genetic 
    code dependent). Optimal codons for several species are in-built.
    By default, the optimal codons of E. coli are assumed.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    optimal_codons: It can be the return value from the get_optimal_codons_from_codonw_coafile function,
                    or it can be a preset value from the codonw software, such as {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}
                    or it can be a custom list containing the best codons.
    Reference
    ----------
    [1] Ikemura, Toshimichi. "Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes." Journal of molecular biology 146.1 (1981): 1-21.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """

    if isinstance(optimal_codons, str):
        if optimal_codons in CBI_and_Fop_preset:
            optimal_codons = CBI_and_Fop_preset[optimal_codons][1]
        else:
            raise TypeError('Preset in {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}')
    elif isinstance(optimal_codons, list):
        pass
    
    aa_opt = set()
    for aa in Obs:
        for c in Obs[aa]:
            if c in optimal_codons:
                aa_opt.add(aa)
    Nopt = 0
    Ntot = 0
    for aa in Obs:
        if aa !="*" and len(Obs[aa]) !=1:
            if aa in aa_opt:
                Ntot += sum(Obs[aa].values())
            for c in Obs[aa]:
                if c in optimal_codons:
                    Nopt += Obs[aa][c]
    Fop = Nopt / Ntot
    return Fop

def get_CBI(Obs, optimal_codons="Escherichia coli"):
    """
    Description 
    ----------
    Codon Bias Index (CBI). 
    Codon bias index is another measure of directional codon bias, it measures 
    the extent to which a gene uses a subset of optimal codons. CBI is similar 
    to Fop as used by Ikemura, with expected usage used as a scaling factor. In a 
    gene with extreme codon bias, CBI will equal 1.0, in a gene with random 
    codon usage CBI will equal 0.0. Note that it is possible for the number of 
    optimal codons to be less than expected by random change. This results in a 
    negative value for CBI.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    optimal_codons: It can be the return value from the get_optimal_codons_from_codonw_coafile function,
                    or it can be a preset value from the codonw software, such as {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}
                    or it can be a custom list containing the best codons.
    Reference
    ----------
    [1] Bennetzen, Jeffrey L., and Benjamin D. Hall. "Codon selection in yeast." Journal of Biological Chemistry 257.6 (1982): 3026-3031.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    if isinstance(optimal_codons, str):
        if optimal_codons in CBI_and_Fop_preset:
            optimal_codons = CBI_and_Fop_preset[optimal_codons][1]
        else:
            raise TypeError('Preset in {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}')
    elif isinstance(optimal_codons, list):
        pass
    
    aa_opt = set()
    for aa in Obs:
        for c in Obs[aa]:
            if c in optimal_codons:
                aa_opt.add(aa)
    Nopt = 0
    Nran = 0
    Ntot = 0
    for aa in Obs:
        if aa !="*" and len(Obs[aa]) !=1:
            if aa in aa_opt:
                optimal_codons_number_aa = 0
                for c in Obs[aa]:
                    Ntot += Obs[aa][c]
                    if c in optimal_codons:
                        Nopt += Obs[aa][c]
                        optimal_codons_number_aa += 1
                Nran += (sum(Obs[aa].values())* optimal_codons_number_aa) /  len(Obs[aa].values())

    CBI = (Nopt - Nran) / (Ntot - Nran)
    return CBI

def get_Relative_Adaptiveness(Obs):
    """
    
    Description
    ----------
    Relative adaptiveness (w) from codonW (by John Peden).

    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    class Relative_Adaptiveness():
        def __init__(self, obj):
            self.Relative_Adaptiveness_dict = obj
            
        def draw_barplot(self,
                         ylabel='Relative adaptiveness',
                         title=None,
                         color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space = 0.16,
                         ax=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            obj: return value of get_Obs, get_Fraction, get_Frequency, get_RSCU or get_Relative_Adaptiveness function.
            figsize: (8,4)
            ylabel: ylabel of plot.
            title: title of plot.
            width: bar spacing width. default=0.9
            color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
            remove_stop_codon: {bool} remove stop codon.
            ax: {None, Aexs}
            codon_space: {0-1} codon spacing.
            """
            draw_codon_barplot(self.Relative_Adaptiveness_dict, ylabel=ylabel, title=title, 
                               color_preset=color_preset, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space)
            return None
    
        def __str__(self):
            STR = ""
            for AA in self.Relative_Adaptiveness_dict:
                STR += AA + ":\n"
                for codon in self.Relative_Adaptiveness_dict[AA]:
                    STR += f"    {codon}: {self.Relative_Adaptiveness_dict[AA][codon]}\n"
            return STR
        
        def __repr__(self):
            return self.__str__()
    
    aa_order = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 
                'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

    Relative_Adaptiveness_dict = {}
    Fraction = get_Fraction(Obs)
    for aa in Fraction.Fraction_dict:
        Relative_Adaptiveness_dict.setdefault(aa, {})
        max_fraction_per_aa = max(Fraction.Fraction_dict[aa].values())
        if max_fraction_per_aa!=0:
            for c in Fraction.Fraction_dict[aa]:
                Relative_Adaptiveness_dict[aa].setdefault(c, Fraction.Fraction_dict[aa][c] / max_fraction_per_aa)
        else:
            for c in Fraction.Fraction_dict[aa]:
                Relative_Adaptiveness_dict[aa].setdefault(c, 0)
            
    Relative_Adaptiveness_dict_order = {}
    for aa in aa_order:
        if aa in Relative_Adaptiveness_dict:
            Relative_Adaptiveness_dict_order.setdefault(aa, {codon: Relative_Adaptiveness_dict[aa][codon] for codon in sorted(Relative_Adaptiveness_dict[aa].keys())})
    return Relative_Adaptiveness(Relative_Adaptiveness_dict_order)

def get_CAI(Obs, ref_Obs, model="codonw", genetic_code=None):
    """
    Description 
    ----------
    Codon Adaptation Index (CAI). 
    CAI is a measurement of the relative adaptiveness of the codon usage of a gene towards the codon usage of highly expressed genes. 
    
    Parameters
    ----------
    Obs : get_Obs function return value. Observed number of occurrences of codon in a query gene.
    ref_Obs : get_Obs function return value. Observed number of occurrences of codon in a reference set of genes.
    model : {emboss, codonw} if model==emboss, Non-synonymous codons and termination codons (dependent on genetic code) are unexcluded. 
                             if model==codonw, Non-synonymous codons and termination codons (dependent on genetic code) are excluded. 
    References
    ----------
    [1] Sharp, Paul M., and Wen-Hsiung Li. "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." Nucleic acids research 15.3 (1987): 1281-1295.
    """
    
    Relative_Adaptiveness = get_Relative_Adaptiveness(ref_Obs)
    w_list = []
    if model == "emboss":
        for aa in Obs:
            for c in Obs[aa]:
                w = Relative_Adaptiveness.Relative_Adaptiveness_dict[aa][c]
                if w != 0:
                    w_list.extend([w]*Obs[aa][c])
                else:
                    w_list.extend([0.01]*Obs[aa][c])
                    
    elif model == "codonw":
        for aa in Obs:
            if aa !="*" and len(Obs[aa].keys()) !=1:
                for c in Obs[aa]:
                    w = Relative_Adaptiveness.Relative_Adaptiveness_dict[aa][c]
                    if w != 0:
                        w_list.extend([w]*Obs[aa][c])
                    else:
                        w_list.extend([0.01]*Obs[aa][c])
                        
    cai = math.exp(sum([math.log(w, math.e) for w in w_list]) /len(w_list))
    return cai

def get_Gravy(Obs):
    """
    Description
    ----------
    The general average hydropathicity or (GRAVY) score, for the hypothetical translated gene product.
    It is calculated as the arithmetic mean of the sum of the hydropathic indices of each amino acid.
    The more the Gravy value is biased towards a negative value, the stronger the hydrophilicity of the protein;
    the more the Gravy value is biased towards a positive value, the stronger the hydrophobicity of the protein
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Reference
    ----------
    [1] Kyte, Jack, and Russell F. Doolittle. "A simple method for displaying the hydropathic character of a protein." Journal of molecular biology 157.1 (1982): 105-132.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    hydropathy_index = {'Ala': 1.8, 'Arg': -4.5, 'Asn': -3.5, 'Asp': -3.5, 'Cys': 2.5, 'Gln': -3.5, 'Glu': -3.5, 'Gly': -0.4, 'His': -3.2, 'Ile': 4.5, 
                        'Leu': 3.8, 'Lys': -3.9, 'Met': 1.9, 'Phe': 2.8, 'Pro': -1.6, 'Ser': -0.8, 'Thr': -0.7, 'Trp': -0.9, 'Tyr': -1.3, 'Val': 4.2, 
                        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    total = 0
    aa_num = 0
    for aa in hydropathy_index:
        if aa in Obs:
            aa_num += sum(Obs[aa].values())
            total += sum(Obs[aa].values()) * hydropathy_index[aa]
    return total/ aa_num

def get_Aromo(Obs):
    """
    Description
    ----------
    Aromaticity score of protein. This is the frequency of aromatic amino  acids (Phe, Tyr, Trp)
    in the hypothetical translated gene product.

    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    aromatic_aa_num = 0
    aromatic_aas= ["Phe", "Tyr", "Trp", "F", "Y", "W"]
    for aromatic_aa in aromatic_aas:
        if aromatic_aa in Obs:
            aromatic_aa_num += sum(Obs[aromatic_aa].values())
    aa_num = sum([sum(Obs[aa].values()) for aa in Obs if aa != "*"])
    return aromatic_aa_num / aa_num

def get_RSCU(Obs):
    """
    Description
    ----------
    Calculate relative synonymous codon usage (RSCU).
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Reference
    ----------
    [1] Sharp, Paul M., Therese MF Tuohy, and Krzysztof R. Mosurski. "Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes." Nucleic acids research 14.13 (1986): 5125-5143.
    """
    
    class RSCU():
        def __init__(self, obj):
            self.RSCU_dict = obj
            
        def draw_barplot(self,
                         ylabel='RSCU',
                         title=None,
                         color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space = 0.16,
                         ax=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            obj: return value of get_Obs, get_Fraction, get_Frequency, get_RSCU or get_Relative_Adaptiveness function.
            figsize: (8,4)
            ylabel: ylabel of plot.
            title: title of plot.
            width: bar spacing width. default=0.9
            color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
            remove_stop_codon: {bool} remove stop codon.
            ax: {None, Aexs}
            codon_space: {0-1} codon spacing.
            """
            draw_codon_barplot(self.RSCU_dict, ylabel=ylabel, title=title, 
                               color_preset=color_preset, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space)
            return None
    
        def __str__(self):
            STR = ""
            for AA in self.RSCU_dict:
                STR += AA + ":\n"
                for codon in self.RSCU_dict[AA]:
                    STR += f"    {codon}: {self.RSCU_dict[AA][codon]}\n"
            return STR
        
        def __repr__(self):
            return self.__str__()
    
    aa_order = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 
                'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
    
    RSCU_dict = {}
    for AA in Obs:
        RSCU_dict.setdefault(AA, {})
        Total = sum(Obs[AA].values())
        n = len(Obs[AA])
        for Codon in Obs[AA]:
            if Total == 0:
                RSCU_dict[AA][Codon] = 0
            else:
                RSCU_dict[AA][Codon] = Obs[AA][Codon]*n/Total
    RSCU_dict_order = {}
    for aa in aa_order:
        if aa in RSCU_dict:
            RSCU_dict_order.setdefault(aa, {codon: RSCU_dict[aa][codon] for codon in sorted(RSCU_dict[aa].keys())})
    return RSCU(RSCU_dict_order)

def get_Fraction(Obs):
    """
    Description 
    ----------
    Calculate Fraction of codon.
    Fraction represents the proportion of each codon in the codon encoding the amino acid, i.e.
    Fraction = the number of occurrences of a codon/the number of occurrences of all codons of the amino acid encoded by the codon.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Note
    ----------
    Cusp software is consistent with the calculated results
    Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    class Fraction():
        def __init__(self, obj):
            self.Fraction_dict = obj
            
        def draw_barplot(self,
                         ylabel='Fraction',
                         title=None,
                         color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space = 0.16,
                         ax=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            obj: return value of get_Obs, get_Fraction, get_Frequency, get_RSCU or get_Relative_Adaptiveness function.
            figsize: (8,4)
            ylabel: ylabel of plot.
            title: title of plot.
            width: bar spacing width. default=0.9
            color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
            remove_stop_codon: {bool} remove stop codon.
            ax: {None, Aexs}
            codon_space: {0-1} codon spacing.
            """
            draw_codon_barplot(self.Fraction_dict, ylabel=ylabel, title=title, 
                               color_preset=color_preset, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space)
            return None
    
        def __str__(self):
            STR = ""
            for AA in self.Fraction_dict:
                STR += AA + ":\n"
                for codon in self.Fraction_dict[AA]:
                    STR += f"    {codon}: {self.Fraction_dict[AA][codon]}\n"
            return STR
        
        def __repr__(self):
            return self.__str__()
        
    Fraction_dict = {}
    for AA in Obs:
        Fraction_dict.setdefault(AA, {})
        Total = sum(Obs[AA].values())
        for Codon in Obs[AA]:
            if Total == 0:
                Fraction_dict[AA][Codon] = 0
            else:
                Fraction_dict[AA][Codon] = Obs[AA][Codon]/Total
            
    aa_order = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 
                'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
    Fraction_dict_order = {}
    for aa in aa_order:
        if aa in Fraction_dict:
            Fraction_dict_order.setdefault(aa, {codon: Fraction_dict[aa][codon] for codon in sorted(Fraction_dict[aa].keys())})
    return Fraction(Fraction_dict_order)

def get_Frequency(Obs):
    """
    Description
    ----------
    Calculate frequency of codon.
    Frequency indicates the Frequency of the codon occurrence in the total gene codon encoding, 
    generally expressed as the number of the codon occurrence in 1000 codons.
    Frequency = the number of the codon occurrence *1000/ the total number of all codons of the gene
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Note
    ----------
    Cusp software is consistent with the calculated results
    Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    class Frequency():
        def __init__(self, obj):
            self.Frequency_dict = obj
            
        def draw_barplot(self,
                         ylabel='Frequency',
                         title=None,
                         color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space=0.16,
                         ax=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            obj: return value of get_Obs, get_Fraction, get_Frequency, get_RSCU or get_Relative_Adaptiveness function.
            figsize: (8,4)
            ylabel: ylabel of plot.
            title: title of plot.
            width: bar spacing width. default=0.9
            color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
            remove_stop_codon: {bool} remove stop codon.
            ax: {None, Aexs}
            codon_space: {0-1} codon spacing.
            """
            draw_codon_barplot(self.Frequency_dict, ylabel=ylabel, title=title, 
                               color_preset=color_preset, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space)
            return None
    
        def __str__(self):
            STR = ""
            for AA in self.Frequency_dict:
                STR += AA + ":\n"
                for codon in self.Frequency_dict[AA]:
                    STR += f"    {codon}: {self.Frequency_dict[AA][codon]}\n"
            return STR
        
        def __repr__(self):
            return self.__str__()
    
    Total = sum([sum(Obs[AA].values()) for AA in Obs])
    Frequency_dict = {}
    for AA in Obs:
        Frequency_dict.setdefault(AA, {})
        for Codon in Obs[AA]:
            Frequency_dict[AA][Codon] = Obs[AA][Codon]/Total*1000
            
    aa_order = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 
                'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
    Frequency_dict_order = {}
    for aa in aa_order:
        if aa in Frequency_dict:
            Frequency_dict_order.setdefault(aa, {codon: Frequency_dict[aa][codon] for codon in sorted(Frequency_dict[aa].keys())})
    return Frequency(Frequency_dict_order)

def draw_codon_barplot(obj, 
                       ylabel=None,
                       title=None,
                       color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                       width=0.9,
                       remove_stop_codon=True,
                       figsize=(8,4),
                       codon_space=0.16,
                       ax=None):
    """
    Description
    ----------
    Draw a codons barplot.
    
    Parameters
    ----------
    obj: return value of get_Obs, get_Fraction, get_Frequency, get_RSCU or get_Relative_Adaptiveness function.
    figsize: (8,4)
    ylabel: ylabel of plot.
    title: title of plot.
    width: bar spacing width. default=0.9
    color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
    remove_stop_codon: {bool} remove stop codon.
    ax: {None, Aexs}
    """
    
    if isinstance(obj, dict):
        obj = obj.copy()
    else:
        if "RSCU_dict" in list(dir(obj)):
            obj = obj.RSCU_dict.copy()
            if ylabel == None:
                ylabel = "RSCU"
        if "Fraction_dict" in list(dir(obj)):
            obj = obj.Fraction_dict.copy()
            if ylabel == None:
                ylabel = "Fraction"
        if "Frequency_dict" in list(dir(obj)):
            obj = obj.Frequency_dict.copy()
            if ylabel == None:
                ylabel = "Frequency"
        if "Relative_Adaptiveness_dict" in list(dir(obj)):
            obj = obj.Relative_Adaptiveness_dict.copy()
            if ylabel == None:
                ylabel = "Relative adaptiveness"
    
    if remove_stop_codon:
        if "*" in obj:
            del obj['*']
    
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    
    if isinstance(color_preset, str):
        cols = plt.colormaps.get_cmap(color_preset).colors
    else:
        cols = color_preset
    
    cex = max([sum(obj[AA].values()) for AA in obj])*codon_space
    for AA in obj:
        value = 0
        values = []
        colors = []
        codons = []
        for codon, color in zip(obj[AA], cols):
            value += obj[AA][codon]
            values.append(value)
            colors.append(color)
            codons.append(codon)
            
        for y, codon, value, color in zip(reversed(range(len(codons))), reversed(codons), reversed(values), reversed(colors)):
            ax.bar(AA, value, width, label=AA, fc=color)
            ax.text(x=AA, y=(-y*cex-2.2*cex)/3, ha="center", va="center", s=codon, fontdict=dict(fontsize=8, color='black', family='monospace'),
                     bbox={'facecolor': color, 'edgecolor':color, 'pad':1})        
    ax.margins(x=0.01)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return None

def get_cusp_like(Obs, human_format=False):
    """
    Description
    ----------
    The calculated result is consistent with result of cusp software.
    
    Parameters
    ----------  
    Obs: get_Obs function return value.
    human_format: {bool} if value is True, return human readable format.
    
    Reference
    ----------
    Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    ATGC1 = {}
    ATGC2 = {}
    ATGC3 = {}
    ATGC3_2d = {}
    ATGC3_4d = {}
    ATGC3_2d_4d = {}
    for AA in Obs:
        for Codon in Obs[AA]:
            if Codon[0] not in ATGC1:
                ATGC1.setdefault(Codon[0], Obs[AA][Codon])
            else:
                ATGC1[Codon[0]] += Obs[AA][Codon]
            if Codon[1] not in ATGC2:
                ATGC2.setdefault(Codon[1], Obs[AA][Codon])
            else:
                ATGC2[Codon[1]] += Obs[AA][Codon]
            if Codon[2] not in ATGC3:
                ATGC3.setdefault(Codon[2], Obs[AA][Codon])
            else:
                ATGC3[Codon[2]] += Obs[AA][Codon]
    Total_base_num = sum(ATGC1.values()) + sum(ATGC2.values())+sum(ATGC3.values())
    GC = (ATGC1.get('C', 0) + ATGC1.get('G', 0) + ATGC2.get('C', 0) + ATGC2.get('G', 0) + ATGC3.get('C', 0) + ATGC3.get('G', 0))/Total_base_num
    GC1 = (ATGC1.get('C', 0) + ATGC1.get('G', 0))/sum(ATGC1.values())
    GC2 = (ATGC2.get('C', 0) + ATGC2.get('G', 0))/sum(ATGC2.values())
    GC3 = (ATGC3.get('C', 0) + ATGC3.get('G', 0))/sum(ATGC3.values())
    Fraction = get_Fraction(Obs)
    Frequency = get_Frequency(Obs)
    CupsResult = [{"Coding GC": GC,
                   "1st letter GC": GC1,
                   "2nd letter GC": GC2,
                   "3rd letter GC": GC3}, 
                  {"Fraction": Fraction,
                   "Frequency":Frequency,
                   "Number": Obs}
                 ]
    
    if human_format:
        out = io.StringIO()
        for k in CupsResult[0]:
            print('#'+k, str(round(CupsResult[0][k]*100, 2))+'%', file=out)
        print('\n#Codon AA Fraction Frequency Number', file=out)
        for AA in CupsResult[1]['Number']:
            for Codon in CupsResult[1]['Fraction'][AA]:
                print(Codon, AA, 
                      round(CupsResult[1]['Fraction'][AA][Codon], 3),
                      round(CupsResult[1]['Frequency'][AA][Codon], 3),
                      CupsResult[1]['Number'][AA][Codon], sep='\t', file=out)
        out.seek(0)
        CupsResult = out.read()
        out.close()
        return CupsResult
    else:
        return CupsResult

def get_codonW_like(Obs, human_format=False):
    """
    Description
    ----------
    Return codonW software calculate result.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    human_format: {bool} if value is True, return human readable format.
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    
    def X3s(Obs, base=["A", "T", "G", "C"]):
        X3s_Obs = {}
        for AA in Obs:
            for Codon in Obs[AA]:
                if Codon[2] == base:
                    X3s_Obs.setdefault(AA, Obs[AA])
                    continue
        X3s_codons = {}
        for AA in X3s_Obs:
                if AA in ['*']:
                    continue
                if len(X3s_Obs[AA]) < 2:
                    continue
                for Codon in Obs[AA]:
                    if Codon[2] not in X3s_codons:
                        X3s_codons.setdefault(Codon[2], Obs[AA][Codon])
                    else:
                        X3s_codons[Codon[2]] += Obs[AA][Codon]
        return X3s_codons.get(base, 0)/sum(X3s_codons.values())
    
    ATGC = {}
    for AA in Obs:
        if AA in ['*']:
            continue
        for Codon in Obs[AA]:
            if Codon[0] not in ATGC:
                ATGC.setdefault(Codon[0], Obs[AA][Codon])
            else:
                ATGC[Codon[0]] += Obs[AA][Codon]
            if Codon[1] not in ATGC:
                ATGC.setdefault(Codon[1], Obs[AA][Codon])
            else:
                ATGC[Codon[1]] += Obs[AA][Codon]
            if Codon[2] not in ATGC:
                ATGC.setdefault(Codon[2], Obs[AA][Codon])
            else:
                ATGC[Codon[2]] += Obs[AA][Codon]
    
    GC = (ATGC.get('C', 0) + ATGC.get('G', 0))/sum(ATGC.values())
    L_aa = int(sum(ATGC.values())/3)
    
    ATGC3s = {}
    for AA in Obs:
        if AA in ['*']:
            continue
        if len(Obs[AA]) < 2:
            continue
        for Codon in Obs[AA]:
            if Codon[2] not in ATGC3s:
                ATGC3s.setdefault(Codon[2], Obs[AA][Codon])
            else:
                ATGC3s[Codon[2]] += Obs[AA][Codon]
    GC3s= ATGC3s.get('G', 0)/sum(ATGC3s.values()) + ATGC3s.get('C', 0)/sum(ATGC3s.values())
    L_sym = sum(ATGC3s.values())
    
    #GCn3 = (GC - G3s -C3s )/(L_aa *3 - L_sym)
    
    A3s = X3s(Obs, base="A")
    T3s = X3s(Obs, base="T")
    G3s = X3s(Obs, base="G")
    C3s = X3s(Obs, base="C")
    
    Nc = get_NC(Obs)
    
    Gravy = get_Gravy(Obs)
    Aromo = get_Aromo(Obs)
    
    codonWResult = {"A3s":A3s, "T3s":T3s, "G3s":G3s, "C3s":C3s, "GC3s": GC3s, "GC": GC,
                    "Nc": Nc, 'L_sym': L_sym, 'L_aa':L_aa, 'Gravy': Gravy, 'Aromo': Aromo}
    if human_format:
        out = io.StringIO()
        print("{:<6s}".format("T3s"),
              "{:<6s}".format("C3s"),
              "{:<6s}".format("A3s"),
              "{:<6s}".format("G3s"),
              "{:<6s}".format("Nc"),
              "{:<6s}".format("GC3s"),
              "{:<6s}".format("GC"),
              "{:<6s}".format("L_sym"), 
              "{:<6s}".format("L_aa"), 
              "{:<6s}".format("Gravy"),
              "{:<6s}".format("Aromo"),
              file=out,
             )
        
        print("{:<6.4f}".format(round(codonWResult["T3s"], 4)), 
              "{:<6.4f}".format(round(codonWResult["C3s"], 4)), 
              "{:<6.4f}".format(round(codonWResult["A3s"], 4)), 
              "{:<6.4f}".format(round(codonWResult["G3s"], 4)),
              "{:<6.2f}".format(round(codonWResult["Nc"], 2)),
              "{:<6.3f}".format(round(codonWResult["GC3s"], 3)), 
              "{:<6.3f}".format(round(codonWResult["GC"], 3)),
              "{:<6}".format(codonWResult["L_sym"]), 
              "{:<6}".format(codonWResult["L_aa"]),
              "{:<6.6f}".format(codonWResult["Gravy"]),
              "{:<6.6f}".format(codonWResult["Aromo"]),
              file=out,
             )
        out.seek(0)
        codonWResult = out.read()
        out.close()
    return codonWResult

def get_GC123(Obs, sym=True):
    """
    Description
    ----------
    Return GC1, GC2, GC3, GC12 value.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    sym: {bool} only synonymous codons model,
         if the value is True, amino acids without synonymous codons and stop codons are deleted.
    """
    
    ATGC1 = {}
    ATGC2 = {}
    ATGC3 = {}
    for AA in Obs:
        if sym:
            if AA in ['*']:
                continue
            if len(Obs[AA]) < 2:
                continue
        for Codon in Obs[AA]:
            if Codon[0] not in ATGC1:
                ATGC1.setdefault(Codon[0], Obs[AA][Codon])
            else:
                ATGC1[Codon[0]] += Obs[AA][Codon]
            if Codon[1] not in ATGC2:
                ATGC2.setdefault(Codon[1], Obs[AA][Codon])
            else:
                ATGC2[Codon[1]] += Obs[AA][Codon]
            if Codon[2] not in ATGC3:
                ATGC3.setdefault(Codon[2], Obs[AA][Codon])
            else:
                ATGC3[Codon[2]] += Obs[AA][Codon]
    if sum(ATGC1.values()) !=0:
        GC1= ATGC1.get('G', 0)/sum(ATGC1.values()) + ATGC1.get('C', 0)/sum(ATGC1.values())
    else:
        GC1 = None
    if sum(ATGC2.values()):
        GC2= ATGC2.get('G', 0)/sum(ATGC2.values()) + ATGC2.get('C', 0)/sum(ATGC2.values())     
    else:
        GC2 = None
    if sum(ATGC3.values()):
        GC3= ATGC3.get('G', 0)/sum(ATGC3.values()) + ATGC3.get('C', 0)/sum(ATGC3.values())
    else:
        GC3 = None
    if GC1 !=None and GC2 !=None:
        GC12 = (GC1 + GC2)/2
    else:
        GC12 = None
    return {"GC1":GC1, "GC2":GC2, "GC12":GC12, "GC3":GC3}

def get_NC(Obs):
    """
    Description
    ----------
    The effective number of codons (NC).
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    
    Reference
    ----------
    [1] Wright, Frank. "The ‘effective number of codons’ used in a gene." Gene 87.1 (1990): 23-29.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    
    def GetF(Obs):
        if len(Obs) == 0:
            return None
        else:
            Fs = []
            for Acid in Obs:
                S = 0
                n = sum(Obs[Acid].values())
                if n==0:
                    continue
                for Codon in Obs[Acid]:
                    S += pow(Obs[Acid][Codon]/n, 2)
                if n-1 != 0:
                    F = (n*S-1)/(n-1)
                    Fs.append(F)
            #Fs = list(filter(lambda x: x!=0, Fs))
            if Fs == []:
                Fs = None
            return Fs
    
    Obs1 = Obs.copy()
    del Obs1["*"]
    
    F1_Obs = {}; F2_Obs = {}; F3_Obs = {}; F4_Obs = {}; F5_Obs = {}; F6_Obs = {}; F7_Obs = {}; F8_Obs = {}
    for Acid in Obs1:
        if len(Obs1[Acid]) == 1:
            F1_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 2:
            F2_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 3:
            F3_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 4:
            F4_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 5:
            F5_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 6:
            F6_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 7:
            F7_Obs.setdefault(Acid, Obs1[Acid])
        if len(Obs1[Acid]) == 8:
            F8_Obs.setdefault(Acid, Obs1[Acid])
    
    Fx_Obs = [F2_Obs, F3_Obs, F4_Obs, F5_Obs, F6_Obs, F7_Obs, F8_Obs]
    
    Fxs = []
    for FxObs in Fx_Obs:
        Fxs.append(GetF(FxObs))
    
    Fxmeans = []
    for Fx in Fxs:
        Fxmeans.append(sum(Fx)/len(list(filter(lambda x: x!=0, Fx))) if (Fx != None and len(list(filter(lambda x: x!=0, Fx)))!=0) else None)
    if (Fxmeans[1] == None) or (Fxmeans[1] == 0):
        try:
            Fxmeans[1] = (Fxmeans[0] + Fxmeans[2])/2
        except:
            return None
    Total = [len(F1_Obs)]
    for FxObs, Fxmean in zip(Fx_Obs, Fxmeans):
        Total.append(len(FxObs)/Fxmean if (Fxmean !=None) else 0)
        
    ENC = sum(Total)
    if ENC > 61 or ENC < 20:
        ENC = None
    return ENC

def find_four_codon_AA(Obs):
    """
    Description:
    ----------
    Find the four-codon amino acids
    
    Parameters
    ----------
    Obs: get_Obs function return value.
    """
    
    four_codon_AA = {}
    for AA in Obs:
        if len(Obs[AA].keys()) > 3:
            m = {}
            for Codon in Obs[AA].keys():
                if Codon[:2] not in m:
                    m.setdefault(Codon[:2], [Codon])
                else:
                    m[Codon[:2]].append(Codon)
            for k in m:
                if len(m[k]) == 4:
                    four_codon_AA.setdefault(AA, m[k])
    return four_codon_AA

def get_PR2(Obs):
    """
    Description:
    ----------
    Parity rule 2 (PR2) analysis.
    
    Parameters
    ----------
    Obs: get_Obs function return value.
        
    Reference
    ----------
    [1] Sueoka, Noboru. "Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of the DNA G+ C content of third codon position." Gene 238.1 (1999): 53-58.
    [2] Nasrullah, Izza, et al. "Genomic analysis of codon usage shows influence of mutation pressure, natural selection, and host features on Marburg virus evolution." BMC evolutionary biology 15 (2015): 1-15.
    """
    
    four_codon_AA = find_four_codon_AA(Obs)
    four_codon = [i for k in four_codon_AA for i in four_codon_AA[k]]
    ATGC3_four_codon = {}
    for AA in Obs:
        for Codon in Obs[AA]:
            if Codon not in four_codon:
                continue
            if Codon[2] not in ATGC3_four_codon:
                ATGC3_four_codon.setdefault(Codon[2], Obs[AA][Codon])
            else:
                ATGC3_four_codon[Codon[2]] += Obs[AA][Codon]
    if sum(ATGC3_four_codon.values()) != 0:                
        A3_four_codon = ATGC3_four_codon.get('A', 0)/sum(ATGC3_four_codon.values())
        T3_four_codon = ATGC3_four_codon.get('T', 0)/sum(ATGC3_four_codon.values())
        G3_four_codon = ATGC3_four_codon.get('G', 0)/sum(ATGC3_four_codon.values())
        C3_four_codon = ATGC3_four_codon.get('C', 0)/sum(ATGC3_four_codon.values())
    else:
        A3_four_codon = 0
        T3_four_codon = 0
        G3_four_codon = 0
        C3_four_codon = 0
    
    if A3_four_codon + T3_four_codon != 0:
        AT3_bias_four_codon = A3_four_codon / (A3_four_codon + T3_four_codon)
    else:
        AT3_bias_four_codon = None
    if G3_four_codon + C3_four_codon != 0:
        GC3_bias_four_codon = G3_four_codon / (G3_four_codon + C3_four_codon)
    else:
        GC3_bias_four_codon = None
    return {"A3/(A3+T3)|4": AT3_bias_four_codon, "G3/(G3+C3)|4": GC3_bias_four_codon}

class NPA_analysis():
    def __init__(self, inputfile, genetic_code, sym=True):
        """
        Description
        ----------
        Neutral plot analysis.
        
        Parameters
        ----------
        inputfile: a fasta or fasta.gz format file include of CDS seqence.
        genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        sym: {bool} only synonymous codons model,
             if the value is True, amino acids without synonymous codons and stop codons are deleted.
        """
        
        GC1 = []
        GC2 = []
        GC3 = []
        GC12 = []
        GeneName = []
        for ID, Seq in FastaIO(inputfile):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            res = get_GC123(Obs, sym)
            if res["GC1"] !=None and res["GC2"] !=None and res["GC12"] !=None and res["GC3"] !=None:
                GC1.append(res["GC1"])
                GC2.append(res["GC2"])
                GC3.append(res["GC3"])
                GC12.append(res["GC12"])
                GeneName.append(ID)
        
        PearsonRResult = ss.pearsonr(GC3, GC12)
        R = PearsonRResult[0]
        P = PearsonRResult[1]
        slope, intercept = np.polyfit(GC3, GC12, 1)
        
        self.sym = sym
        self.R = R
        self.P = P
        self.slope = slope
        self.intercept = intercept
        self.GC1 = GC1
        self.GC2 = GC2
        self.GC3 = GC3
        self.GC12 = GC12
        self.GeneName = GeneName
        
    def __str__(self):
        return str(self.get_df())
    
    def __repr__(self):
        return self.__str__()
    
    def get_df(self):
        return pd.DataFrame({"Gene Name":self.GeneName,
                             "GC3": self.GC3,
                             "GC12": self.GC12, 
                             "GC1": self.GC1,
                             "GC2": self.GC2})
    
    def draw_NPA_plot(self, 
                      show_gene_names=False,
                      figsize = (6,4),
                      gene_names_size=10,
                      gene_names_color="#0A0A0A",
                      point_color="#4F845C",
                      line_color="#C25759", 
                      point_size = 20,
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax = None):
        """
        Description
        ----------
        Draw NPA plot.

        Parameters
        ----------
        NPAResult: NPA function return value.
        show_gene_names: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_names_size: font size of gene name. 
        gene_names_color: font color of gene name. 
        point_color: point color.
        line_color: strand line color.
        point_size: point size.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
        ax: {None, Aexs}
        """
        xs = self.GC3
        ys = self.GC12
        labels = self.GeneName
        
        if xlabel == None:
            if self.sym:
                xlabel = "P$_3$/GC$_3$$_s$"
            else:
                xlabel = "GC$_3$"

        if ylabel == None:
            if self.sym:
                ylabel = "P$_1$$_2$/GC$_1$$_,$$_2$$_s$"
            else:
                ylabel =  "GC$_1$$_,$$_2$"
        
        if title == None:
            title = "Neutral plot analysis"
            formula = '$y$ = {:.4f}$x$ + {:.4f}; $R$$^2$ = {:.4f}; $P$ = {:.4f}'.format(self.slope, self.intercept, pow(self.R, 2), self.P)
            title = title+"\n"+formula
            
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        sns.regplot(x=xs, y=ys, 
                    fit_reg=True,
                    scatter_kws={"color": point_color, 's': point_size, 'alpha':1}, 
                    line_kws={"color":line_color, "linewidth": 2}, 
                    label="Regression Line",
                    ax=ax)
        
        if show_gene_names:
            for x,y,l in zip(xs, ys, labels):
                if isinstance(show_gene_names, bool):
                    ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
                else:
                    if l in show_gene_names:
                        ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.plot(xs, xs, c=line_color)
        return None

class ENC_analysis():
    def __init__(self, inputfile, genetic_code):
        """
        Description
        ----------
        Effective number of codons (ENC) analysis

        Parameters
        ----------
        inputfile: a fasta or fasta.gz format file include of CDS seqence.
        genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        """
        
        ys = []
        xs = []
        GeneName = []
        for ID, Seq in FastaIO(inputfile):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            xs.append(get_GC123(Obs, sym=True)["GC3"])
            #xs.append(get_GC3s(Obs))
            try:
                ys.append(get_NC(Obs))
            except:
                print(ID)
            GeneName.append(ID)
            
        self.GC3s = xs
        self.ENC = ys
        self.GeneName = GeneName
        
    def __str__(self):
        return str(self.get_df())
    
    def __repr__(self):
        return self.__str__()
    
    def get_df(self):
        return pd.DataFrame({"Gene Name":self.GeneName, "GC3s": self.GC3s, "ENC": self.ENC})
        
    def draw_ENC_plot(self, 
                      figsize=(6,4),
                      show_gene_names=False,
                      gene_names_color="#0A0A0A",
                      gene_names_size=10,
                      point_color="#4F845C",
                      point_size = 20,
                      line_color="#C25759", 
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax=None):
        """
        Description
        ----------
        Draw ENC plot.

        Parameters
        ----------
        ENCResult: ENC function return value.
        figsize: {(6,4)}
        show_gene_names: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_names_size: font size of gene name. 
        gene_names_color: font color of gene name. 
        point_size: point size
        point_color: point color.
        line_color: strand line color.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
        ax: {None, Aexs}
        """
        def gc2enc(value):
            return 2 + value + 29/(pow(value, 2) + pow(1-value, 2))

        xs = self.GC3s
        ys = self.ENC
        labels = self.GeneName
        
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        ax.scatter(x=xs, y=ys, c=point_color,s=point_size)
        ax.set_xlim(0,1)
        ax.set_ylim(0,70)
        
        if xlabel == None:
            xlabel="GC$_3$$_s$"
        if ylabel == None:
            ylabel="ENC"
        if title == None:
            title="ENC plot analysis"
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        if show_gene_names:
            for x,y,l in zip(xs, ys, labels):
                if isinstance(show_gene_names, bool):
                    ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
                else:
                    if l in show_gene_names:
                        ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
        ax.plot(np.arange(0, 1, 0.005), gc2enc(np.arange(0, 1, 0.005)), c=line_color)
        return None

class PR2_analysis():
    def __init__(self, inputfile, genetic_code):
        """
        Description
        ----------
        Parity rule 2 (PR2) analysis.
        
        Parameters
        ----------
        inputfile: a fasta or fasta.gz format file.
        genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.

        Reference
        ----------
        [1] Sueoka N. Intrastrand parity rules of DNA base composition and usage biases of synonymous codons.
            J Mol Evol. 1995 Mar;40(3):318-25. doi: 10.1007/BF00163236. PMID: 7723058.
        [2] Sueoka N. Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of
            the DNA G+C content of third codon position. Gene. 1999 Sep 30;238(1):53-8. doi: 10.1016/s0378-1119(99)00320-0. PMID: 10570983.
        """

        ys = []
        xs = []
        GeneName = []
        for ID, Seq in FastaIO(inputfile):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            res = get_PR2(Obs)
            ys.append(res["A3/(A3+T3)|4"])
            xs.append(res["G3/(G3+C3)|4"])
            GeneName.append(ID)
            
        self.GeneName = GeneName
        self.X_axis = xs
        self.Y_axis = ys
    
    def __str__(self):
        return str(self.get_df())
    
    def __repr__(self):
        return self.__str__()
    
    def get_df(self):
        return pd.DataFrame({"Gene Name":self.GeneName, "A3/(A3+T3)|4": self.X_axis, "G3/(G3+C3)|4": self.Y_axis})
    
    def draw_PR2_plot(self, 
                      show_gene_names=False,
                      figsize=(6,4),
                      gene_names_size=10,
                      gene_names_color="#0A0A0A",
                      point_color="#4F845C",
                      point_size = 20,
                      line_color="#C25759", 
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax = None):
        """
        Description
        ----------
        Draw parity rule 2 (PR2) plot.

        Parameters
        ----------
        PR2Result: PR2 function return value.
        figsize: (6,4)
        show_gene_names: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_names_size: font size of gene name. 
        gene_names_color: font color of gene name. 
        point_color: point color.
        point_size: point size.
        line_color: strand line color.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
        ax: {None, Aexs}
        """
        ys = self.X_axis
        xs = self.Y_axis
        labels = self.GeneName
        
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        ax.scatter(xs, ys, c=point_color,s=point_size)
        ax.axhline(0.5, c=line_color)
        ax.axvline(0.5, c=line_color)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        if title == None:
            title="PR2 plot analysis"
        if xlabel == None:
            xlabel="G$_3$/(G$_3$+C$_3$)|4"
        if ylabel == None:
            ylabel="A$_3$/(A$_3$+T$_3$)|4"
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        if show_gene_names:
            for x,y,l in zip(xs, ys, labels):
                if isinstance(show_gene_names, bool):
                    ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
                else:
                    if l in show_gene_names:
                        ax.text(x, y, l, c=gene_names_color, size=gene_names_size)
        return None

class RSCU_analysis():
    def __init__(self, data, genetic_code):
        """
        Description
        ----------
        RSCU analysis.
        
        Parameters
        ----------
        data: [("species name", "cds file path"), ...]
        genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        """
        
        single_codon_amino_acids = set()
        series_list = []
        RSCUs_dict = {}
        
        for prefix, file in data:
            RSCU_Class = get_RSCU(get_Obs(file, genetic_code=genetic_code))
            RSCUs_dict.setdefault(prefix, RSCU_Class)
            RSCU = RSCU_Class.RSCU_dict
            del RSCU["*"]
            aas = []
            codons = []
            values = []
            for aa in RSCU:
                if len(RSCU[aa])==1:
                    single_codon_amino_acids.add(aa)
                for codon in RSCU[aa]:
                    aas.append(aa)
                    codons.append(codon)
                    values.append(RSCU[aa][codon])
            index = [aas, codons]
            series = pd.Series(values, index=index)
            series.name = prefix
            series.index.names = ("Amino Acid", "Codon")
            series_list.append(series)
        RSCU_df = pd.DataFrame(series_list).T
        RSCU_df = RSCU_df.sort_index()
        self.RSCU_raw_df = RSCU_df
        RSCU_df = RSCU_df.drop(list(single_codon_amino_acids), axis=0)
        
        self.RSCU_df = RSCU_df
        self.single_codon_amino_acids = tuple(single_codon_amino_acids)
        
        df = self.RSCU_df.replace(0, pd.NA)
        df = df.dropna()
        df = df.astype("float64")
        ca = prince.CA(n_components=2)
        self.ca = ca.fit(df)
        pca = prince.PCA(n_components=2)
        self.pca = pca.fit(df)
        self.RSCUs_dict = RSCUs_dict
        self.AAColorShapeType={"Gly": ("#4DBBD5FF", "o", "Nonpolar"),
                               "Ala": ("#4DBBD5FF", "v", "Nonpolar"),
                               "Val": ("#4DBBD5FF", "^", "Nonpolar"),
                               "Leu": ("#4DBBD5FF", "<", "Nonpolar"),
                               "Ile": ("#4DBBD5FF", ">", "Nonpolar"),
                               "Phe": ("#4DBBD5FF", "s", "Nonpolar"),
                               "Pro": ("#4DBBD5FF", "p", "Nonpolar"),
                               "Met": ("#4DBBD5FF", "d", "Nonpolar"),
                               "Trp": ("#4DBBD5FF", "+", "Nonpolar"),
                               "Ser": ("#F39B7FFF", "o", "Polar"), 
                               "Thr": ("#F39B7FFF", "s", "Polar"), 
                               "Cys": ("#F39B7FFF", "p", "Polar"), 
                               "Asn": ("#F39B7FFF", "+", "Polar"), 
                               "Gln": ("#F39B7FFF", "d", "Polar"), 
                               "Tyr": ("#F39B7FFF", "^", "Polar"), 
                               "Asp": ("#8491B4FF", "o", "Acidic"),
                               "Glu": ("#8491B4FF", "s", "Acidic"),
                               "Lys": ("#91D1C2FF", "o", "Basic"), 
                               "Arg": ("#91D1C2FF", "s", "Basic"), 
                               "His": ("#91D1C2FF", "p", "Basic")}
        
    def get_PCA_df(self):
        PCA_df = self.pca.column_correlations
        PCA_df.columns = [f"PCA1 ({self.pca.eigenvalues_summary.iloc[0, 1]})",
                          f"PCA2 ({self.pca.eigenvalues_summary.iloc[1, 1]})"]
        PCA_df.index.name = "Species"
        return PCA_df
        
    def get_COA_df(self):
        ca_column = self.ca.column_coordinates(self.RSCU_df)
        ca_column["Type"] = "Species"
        ca_row = self.ca.row_coordinates(self.RSCU_df)
        ca_row["Type"] = "Codon"
        COA_df = pd.concat([ca_column, ca_row])
        COA_df.columns = [f"COA1 ({self.ca.eigenvalues_summary.iloc[0, 1]})",
                          f"COA2 ({self.ca.eigenvalues_summary.iloc[1, 1]})",
                          "Type"]
        return COA_df
    
    def draw_COA_plot(self, 
                      figsize=(8,8), 
                      species_labels_color="black",
                      species_labels_style="italic",
                      species_labels_size = 8,
                      species_shapes_color = None,
                      species_shapes_size = 200,
                      species_labels_ha = "left",
                      species_labels_va = "bottom",
                      codon_labels_color = "black",
                      codon_labels_size = 8,
                      codon_shapes_size = 100,
                      codon_labels_ha = "left",
                      codon_labels_va = "bottom",
                      show_species_labels = True,
                      show_codon_labels = True,
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      ax=None):
        """  
        Description 
        -----------
        
        Draw correspondence analysis of RSCU plot.

        Parameters
        ----------
        figsize: (8,8)
        species_labels_color: "black"
        species_labels_style: "italic"
        species_labels_size: 8,
        species_shapes_color: {None, {"species1": "red", "species2":"blue", ... }}
        species_shapes_size: 200,
        species_labels_ha: {"left", "right", "top", "bottom", "center"}
        species_labels_va: {"left", "right", "top", "bottom", "center"}
        codon_labels_color: "black",
        codon_labels_size: 8,
        codon_shapes_size: 100,
        codon_labels_ha: {"left", "right", "top", "bottom", "center"}
        codon_labels_va: {"left", "right", "top", "bottom", "center"}
        show_species_labels: {bool, ["species1", "species2", ...]} show species name in plot.
        show_codon_labels: {bool, ["codon1", "codon2", ...]} show codon in plot.
        title {None, str}
        xlabel {None, str}
        ylabel {None, str}
        title_size: 12
        xlabel_size: 12
        ylabel_size: 12
        ax: {None, Aexs}
        """
        
        df = self.RSCU_df.replace(0, pd.NA)
        df = df.dropna()
        df = df.astype("float64")
        row_coords = self.ca.row_coordinates(df)
        col_coords = self.ca.column_coordinates(df)
        
        # Plotting the results
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        # Adding labels
        if show_codon_labels !=None and show_codon_labels !=False:
            if isinstance(show_codon_labels, list):
                for i, index in enumerate(df.index):
                    aa = index[0]
                    codon = index[1]
                    if codon in show_codon_labels:
                        ax.annotate(codon,
                                    (row_coords[0][i], row_coords[1][i]), 
                                    color=codon_labels_color,
                                    fontsize=codon_labels_size,
                                    ha=codon_labels_ha,
                                    va=codon_labels_va)
            else:
                for i, index in enumerate(df.index):
                    aa = index[0]
                    codon = index[1]
                    ax.annotate(codon,
                                (row_coords[0][i], row_coords[1][i]),
                                color=codon_labels_color,
                                fontsize=codon_labels_size,
                                ha=codon_labels_ha,
                                va=codon_labels_va)
        
        if show_species_labels != None and show_species_labels != False:
            if isinstance(show_species_labels, list):
                for i, spceies in enumerate(df.columns):
                    if spceies in show_species_labels:
                        ax.annotate(spceies, 
                                    (col_coords[0][i], col_coords[1][i]),
                                    color=species_labels_color, 
                                    fontsize=species_labels_size,
                                    style = species_labels_style,
                                    ha = species_labels_ha,
                                    va = species_labels_va)
            else:
                for i, spceies in enumerate(df.columns):
                    ax.annotate(spceies,
                                (col_coords[0][i], col_coords[1][i]),
                                color=species_labels_color,
                                fontsize=species_labels_size,
                                style = species_labels_style,
                                ha = species_labels_ha,
                                va = species_labels_va)
                    
        # Adding points
        Pa_dict = {}
        for i in range(0,len(row_coords)):
            aa = row_coords.index[i][0]
            codon = row_coords.index[i][1]
            Pa = ax.scatter(row_coords[0][i], row_coords[1][i], 
                            c=self.AAColorShapeType[aa][0],
                            label=aa, 
                            marker=self.AAColorShapeType[aa][1],
                            s=codon_shapes_size)
            Pa_dict[aa] = Pa
            
        Pa_list = [Pa_dict[k] for k in list(self.AAColorShapeType.keys()) if k in Pa_dict]
        
        Ps_list = []
        for i in range(0, len(col_coords)):
            if isinstance(species_shapes_color, dict):
                Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                c=species_shapes_color.get(col_coords.index[i],"#E64B35FF"), 
                                label=col_coords.index[i],
                                marker="*", 
                                s=species_shapes_size)
                Ps_list.append(Ps)
            else:
                Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                c='#E64B35FF', 
                                label='Species', 
                                marker="*",
                                s=species_shapes_size)
                if i == 0:
                    Ps_list.append(Ps)
                    
        # Adding legend
        ax.legend(handles=[*tuple(Pa_list), *tuple(Ps_list)],
                  loc='upper left', 
                  bbox_to_anchor=(1, 1), 
                  ncol=1,
                  frameon=False,
                  shadow=False,
                  title='',
        )
        
        if title == None:
            title = 'Correspondence analysis of RSCU'
        if xlabel == None:
            xlabel = f'Dim 1 ({self.ca.eigenvalues_summary.iloc[0,1]})'
        if ylabel == None:
            ylabel = f'Dim 2 ({self.ca.eigenvalues_summary.iloc[1,1]})'
        ax.set_title(title, size=title_size)
        ax.set_xlabel(xlabel, size=xlabel_size)
        ax.set_ylabel(ylabel, size=ylabel_size)
        return None
    
    def draw_PCA_plot(self,
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="italic", 
                      labels_size = 8, 
                      shapes_color = '#E64B35FF', 
                      labels_ha = "left", 
                      labels_va = "bottom",
                      shapes_type = '*',
                      shapes_size = 100, 
                      show_labels = True,
                      show_legend = False,
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      ax=None):
        """
        Description 
        -----------
        
        Draw Principal component analysis plot.

        Parameters
        ----------
        figsize: (6,6)
        labels_color: black
        labels_style: {'normal', 'italic', 'oblique'}
        labels_size: 8,
        labels_ha: {"left", "right", "top", "bottom", "center"}
        labels_va: {"left", "right", "top", "bottom", "center"}
        show_labels: {bool, ["species1", "species2", ...]} show column names in plot.
        shapes_color: {None, {"species1": "red", "species2":"blue", ... }}
        shapes_type: {None, {"species1": "*", "species2":"<", ... }}
        shapes_size: 100,
        show_legend: {bool}
        title {None, str}
        xlabel {None, str}
        ylabel {None, str}
        title_size: 12
        xlabel_size: 12
        ylabel_size: 12
        ax: {None, Aexs}
         """
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        for z, i in zip(self.get_shape_and_color(self.pca.column_correlations.shape[0]), range(self.pca.column_correlations.shape[0])):
            x = self.pca.column_correlations[0][i]
            y = self.pca.column_correlations[1][i]
            label = self.pca.column_correlations.index[i]
            shape_type=z[0]
            shape_color=z[1]
            
            if shapes_type !=None:
                if isinstance(shapes_type, dict):
                    shape_type = shapes_type.get(label, z[0])
                else:
                    shape_type = shapes_type
                    
            if shapes_color !=None:
                if isinstance(shapes_color, dict):
                    shape_color = shapes_color.get(label, z[1])
                else:
                    shape_color = shapes_color

            ax.scatter(x,y, label=label, 
                       marker=shape_type, 
                       color=shape_color,
                       s=shapes_size)
            
            if isinstance(show_labels, list):
                if label in show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)
            else:
                if show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)

        if title == None:
            title = f'Principal component analysis of RSCU'
        if xlabel == None:
            xlabel = f"Dim 1 ({self.pca.eigenvalues_summary.iloc[0,1]})"
        if ylabel == None:
            ylabel = f"Dim 2 ({self.pca.eigenvalues_summary.iloc[1,1]})"

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if show_legend == True:
            ax.legend(loc='center left',  
                      bbox_to_anchor=(1, 0.5),
                      ncol=1,
                      frameon=False,
                      shadow=False,
                      title='')
        return None
    
    def get_shape_and_color(self, length):
        shapes = ["o", "s", "^", "<", ">", "v", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "8"]
        colors = ["#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF"]
        n = len(shapes) * len(colors)
        l = []
        for x in shapes:
            for y in colors:
                l.append((x, y))
        res = []
        for i in range(0, length):
            res.append(l[i%n])
        return res
    
    def draw_heatmap(self, 
                     figsize=None, 
                     cmap="Blues",
                     ax=None,
                     ylabels_fontstyle='italic',
                     xlabels_fontstyle='normal',
                     ylabels_fontsize=10,
                     xlabels_fontsize=10,
                     show_ylabels=True,
                     show_xlabels=True
                    ):
        """
        Description 
        -----------
        
        Draw heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}
        ylabels_fontsize: 10
        xlabels_fontsize: 10
        show_ylabels: {bool}
        show_xlabels: {bool}
        """
        if figsize == None:
            figsize = (14, self.RSCU_df.shape[1]/2)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(self.RSCU_df.T, cmap=cmap, ax=ax, 
                    yticklabels=show_ylabels,
                    xticklabels=show_xlabels)
        ax.set_yticklabels(ax.get_yticklabels(),fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle)
        ax.set_xticklabels(ax.get_xticklabels(),fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle)
        ax.set_xlabel("")
        return None
    
    def draw_clustermap(self, 
                        figsize=None, 
                        cmap="Blues",
                        ylabels_fontstyle='italic',
                        xlabels_fontstyle='normal',
                        ylabels_fontsize=10,
                        xlabels_fontsize=10,
                        row_cluster=True,
                        col_cluster=True,
                        show_ylabels=True,
                        show_xlabels=True):
        """
        Description 
        -----------
        
        Draw clustermap of RSCU.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}
        ylabels_fontsize: 10
        xlabels_fontsize: 10
        show_ylabels: {bool}
        show_xlabels: {bool}
        """
        if figsize == None:
            figsize = (14, self.RSCU_df.shape[1]/2)
        
        cmp = sns.clustermap(self.RSCU_df.T,
                             figsize=figsize, 
                             cmap=cmap,
                             row_cluster = row_cluster,
                             col_cluster = col_cluster,
                             yticklabels = show_ylabels,
                             xticklabels = show_xlabels
                            )
        cmp.ax_heatmap.set_yticklabels(cmp.ax_heatmap.get_yticklabels(), fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle)
        cmp.ax_heatmap.set_xticklabels(cmp.ax_heatmap.get_xticklabels(), fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle)
        cmp.ax_heatmap.set_xlabel("")
        return None
    
    def draw_boxplot(self,
                     figsize = None,
                     ax = None,
                     fontstyle = 'italic',
                     fontsize = 10):
        """
        Description 
        -----------
        
        Draw boxplot of RSCU.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        fontstyle: {'normal', 'italic', 'oblique'}
        fontsize: 10
        """
        if figsize == None:
            figsize = (4, self.RSCU_df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(self.RSCU_df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        ax.set_xlabel("RSCU")
        return None
    
    def draw_pearson_heatmap(self,
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='italic',
                             labels_fontsize=10,
                             ax = None
                            ):
        """
        Description 
        -----------
        Draw pearson heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        labels_fontstyle: {'normal', 'italic', 'oblique'}
        labels_fontsize: 10
        """
        if figsize == None:
            figsize = (self.RSCU_df.shape[1]/4+1, self.RSCU_df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        corr = self.RSCU_df.corr(method='pearson')
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, fmt=".2f", cmap="Blues", ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        return None
    
    def draw_RSCU_barplot(self):
        fig, axs = plt.subplots(len(self.RSCUs_dict),
                                figsize=(8, len(self.RSCUs_dict)*4))
        fig.subplots_adjust(hspace=0.8)
        for prefix, ax in zip(self.RSCUs_dict, axs):
            self.RSCUs_dict[prefix].draw_barplot(title=prefix, codon_space=0.25, ax=ax)
        return None
    
    def get_tree(self, metric='euclidean', outgroup="midpoint", tree_method="upgma"):
        """
        Description 
        -----------
        Return a tree object from TreeNode class of scikit-bio. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme}  Method of phylogenetic reconstruction. See also skibio.tree.
        metric : {euclidean, cityblock, braycurtis, canberra, chebyshev,
        correlation, cosine, dice, hamming, jaccard, jensenshannon,
        kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
        russellrao, seuclidean, sokalmichener, sokalsneath,
        sqeuclidean, yule} The distance metric to use. 
        Euclidean distance, high sensitivity, suitable for related species;
        braycurtis distance reduces the impact of extreme values, suitable for distant species or highly variable genes, emphasizing compositional differences, such as ecological data
        outgroup : {None, midpoint, list[Species1, Species2, ...]}  if outgroup == None, return unroot tree.
        """
        dist = pdist(np.matrix(self.RSCU_df.T), metric=metric)
        dist_matrix = DistanceMatrix(squareform(dist), ids=list(self.RSCU_df.T.index))
        if tree_method == 'nj':
            tree = nj(dist_matrix)
        elif tree_method == 'upgma':
            tree = upgma(dist_matrix)
        elif tree_method == "gme":
            tree = gme(dist_matrix)
        elif tree_method == 'bme':
            tree = bme(dist_matrix)
            
        if outgroup == None:
            pass
        elif outgroup == "midpoint":
            tree = tree.root_at_midpoint(inplace=False, reset=True)
        else:
            tree = tree.root_by_outgroup(outgroup)
        return tree
    
    def get_tree_newick_string(self, metric='euclidean', outgroup="midpoint", tree_method="upgma"):
        """
        Description 
        -----------
        Return a tree string of newick format. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme}  Method of phylogenetic reconstruction. See also skibio.tree.
        metric : {euclidean, cityblock, braycurtis, canberra, chebyshev,
        correlation, cosine, dice, hamming, jaccard, jensenshannon,
        kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
        russellrao, seuclidean, sokalmichener, sokalsneath,
        sqeuclidean, yule} The distance metric to use. 
        Euclidean distance, high sensitivity, suitable for related species;
        braycurtis distance reduces the impact of extreme values, suitable for distant species or highly variable genes, emphasizing compositional differences, such as ecological data
        outgroup : {None, midpoint, list[Species1, Species2, ...]}  if outgroup == None, return unroot tree.
        """
        tree = self.get_tree(metric, outgroup, tree_method)
        f = io.StringIO()
        tree.write(f)
        f.seek(0)
        tree_str = f.read()
        f.close()
        return tree_str[:-1]
    
    def draw_tree_plot(self,
                       figsize=(8,8), 
                       ax=None,
                       tree_method="upgma",
                       metric='euclidean',
                       outgroup="midpoint",
                       ignore_branch_length=True,
                       ladderize=True,
                       ladderize_by="size",
                       ladderize_direction="right",
                       leaf_label_size=10,
                       linewidth=1.5,
                       width=10,
                       height=0.7):
        """
        Description 
        -----------
        Draw tree plot of RSCU by NJ method.
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme}  Method of phylogenetic reconstruction. See also skibio.tree.
        metric : {euclidean, cityblock, braycurtis, canberra, chebyshev,
        correlation, cosine, dice, hamming, jaccard, jensenshannon,
        kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
        russellrao, seuclidean, sokalmichener, sokalsneath,
        sqeuclidean, yule} The distance metric to use. 
        Euclidean distance, high sensitivity, suitable for related species;
        braycurtis distance reduces the impact of extreme values, suitable for distant species or highly variable genes, emphasizing compositional differences, such as ecological data
        outgroup : {None, midpoint, list[Species1, Species2, ...]}  if outgroup == None, return unroot tree.
        figsize: tuple of (width, height), optional
        ax : matplotlib Axes, Axes in which to draw the plot, otherwise use the currently-active Axes.        
        ignore_branch_length : {bool} Ignore branch lengths for cladogram
        innode_label_size : {float} Font size for internal node labels.
        ladderize : {bool} Enable ladderize tree sorting.
        ladderize_by : {size, branch_length} Sort criterion.
        ladderize_direction : {left, right} Direction for larger subtrees.
        leaf_label_size : {float} Font size for leaf labels.
        linewidth : {float} Branch line width.
        height : Figure height per leaf node.
        width : Figure width.
        """
        tree = self.get_tree(metric=metric, outgroup=outgroup, tree_method=tree_method)
        plotter = TreePlotter(tree, 
                              ignore_branch_length=ignore_branch_length,
                              ladderize=ladderize,
                              ladderize_by=ladderize_by,
                              ladderize_direction=ladderize_direction,
                              leaf_label_size=leaf_label_size,
                              linewidth=linewidth,
                              width=width,
                              height=height)
        plotter.plot(figsize=figsize, ax=ax)
        return None


class AAComposition_analysis():
    def __init__(self, data, genetic_code):
        """
        Description
        ----------
        Amino acid composition analysis.
        
        Parameters
        ----------
        data: [("species name", "cds file path"), ...]
        genetic_code: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        """
        Seq1toSeq3 = {v:k for k,v in Seq3toSeq1.items() if k!='*'}
        
        AA_count_dict = {}
        for prefix, file in data:
            AA_count = defaultdict(int)
            for ID, Seq in FastaIO(file):
                for i in translate(Seq, genetic_code=genetic_code):
                    AA_count[i] = AA_count[i] + 1
                    
            AA_count_dict.setdefault(prefix, {Seq1toSeq3[AA]:AA_count[AA] for AA in AA_count if AA in Seq1toSeq3})
        
        self.AA_count_df = pd.DataFrame(AA_count_dict)
        self.AA_Fraction_df = pd.DataFrame(data=None, dtype="float64", columns=self.AA_count_df.columns)
        for i in range(0, self.AA_count_df.shape[1]):
            self.AA_Fraction_df.iloc[:,i] = self.AA_count_df.iloc[:,i] / self.AA_count_df.iloc[:,i].sum()
        AA_count_ca = prince.CA(n_components=2)
        self.AA_count_ca = AA_count_ca.fit(self.AA_count_df)
        AA_Fraction_ca = prince.CA(n_components=2)
        self.AA_Fraction_ca = AA_Fraction_ca.fit(self.AA_Fraction_df)
        AA_count_pca = prince.PCA(n_components=2)
        self.AA_count_pca = AA_count_pca.fit(self.AA_count_df)
        AA_Fraction_pca = prince.PCA(n_components=2)
        self.AA_Fraction_pca = AA_Fraction_pca.fit(self.AA_Fraction_df)
        self.AAColorShapeType={"Gly": ("#4DBBD5FF", "o", "Nonpolar"),
                               "Ala": ("#4DBBD5FF", "v", "Nonpolar"),
                               "Val": ("#4DBBD5FF", "^", "Nonpolar"),
                               "Leu": ("#4DBBD5FF", "<", "Nonpolar"),
                               "Ile": ("#4DBBD5FF", ">", "Nonpolar"),
                               "Phe": ("#4DBBD5FF", "s", "Nonpolar"),
                               "Pro": ("#4DBBD5FF", "p", "Nonpolar"),
                               "Met": ("#4DBBD5FF", "d", "Nonpolar"),
                               "Trp": ("#4DBBD5FF", "+", "Nonpolar"),
                               "Ser": ("#F39B7FFF", "o", "Polar"), 
                               "Thr": ("#F39B7FFF", "s", "Polar"), 
                               "Cys": ("#F39B7FFF", "p", "Polar"), 
                               "Asn": ("#F39B7FFF", "+", "Polar"), 
                               "Gln": ("#F39B7FFF", "d", "Polar"), 
                               "Tyr": ("#F39B7FFF", "^", "Polar"), 
                               "Asp": ("#8491B4FF", "o", "Acidic"),
                               "Glu": ("#8491B4FF", "s", "Acidic"),
                               "Lys": ("#91D1C2FF", "o", "Basic"), 
                               "Arg": ("#91D1C2FF", "s", "Basic"), 
                               "His": ("#91D1C2FF", "p", "Basic")}
        
    def get_PCA_df(self, dtype="Fraction"):
        """
        Description
        ----------
        get PCA data.
        
        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        """
        if dtype == "Fraction":
            pca = self.AA_Fraction_pca
        else:
            pca = self.AA_count_pca
        PCA_df = pca.column_correlations
        PCA_df.columns = [f"PCA1 ({pca.eigenvalues_summary.iloc[0, 1]})",
                          f"PCA2 ({pca.eigenvalues_summary.iloc[1, 1]})"]
        PCA_df.index.name = "Species"
        return PCA_df
        
    def get_COA_df(self, dtype="Fraction"):
        """
        Description
        ----------
        get COA data.
        
        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        """
        if dtype == "Fraction":
            ca = self.AA_Fraction_ca
            ca_column = ca.column_coordinates(self.AA_Fraction_df)
            ca_row = ca.row_coordinates(self.AA_Fraction_df)
        else:
            ca = self.AA_count_ca
            ca_column = ca.column_coordinates(self.AA_count_df)
            ca_row = ca.row_coordinates(self.AA_count_df)
        
        ca_column["Type"] = "Species"
        ca_row["Type"] = "Codon"
        COA_df = pd.concat([ca_column, ca_row])
        COA_df.columns = [f"COA1 ({ca.eigenvalues_summary.iloc[0, 1]})",
                          f"COA2 ({ca.eigenvalues_summary.iloc[1, 1]})","Type"]
        return COA_df
    
    def draw_heatmap(self, 
                     dtype="Fraction",
                     figsize=None, 
                     cmap="Blues",
                     ax=None,
                     ylabels_fontstyle='italic',
                     xlabels_fontstyle='normal',
                     ylabels_fontsize=10,
                     xlabels_fontsize=10,
                     show_ylabels=True,
                     show_xlabels=True
                    ):
        """
        Description 
        -----------
        
        Draw heatmap of count or Fraction.
        
        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}
        ylabels_fontsize: 10
        xlabels_fontsize: 10
        show_ylabels: {bool}
        show_xlabels: {bool}
        """
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (14, df.shape[1]/2)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(df.T, cmap=cmap, ax=ax, 
                    yticklabels=show_ylabels,
                    xticklabels=show_xlabels)
        ax.set_yticklabels(ax.get_yticklabels(),fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle)
        ax.set_xticklabels(ax.get_xticklabels(),fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle)
        ax.set_xlabel("")
        return None
    
    def draw_clustermap(self, 
                        dtype="Fraction",
                        figsize=None, 
                        cmap="Blues",
                        ylabels_fontstyle='italic',
                        xlabels_fontstyle='normal',
                        ylabels_fontsize=10,
                        xlabels_fontsize=10,
                        row_cluster=True,
                        col_cluster=True,
                        show_ylabels=True,
                        show_xlabels=True):
        """
        Description 
        -----------
        
        Draw clustermap of count or Fraction.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}
        ylabels_fontsize: 10
        xlabels_fontsize: 10
        show_ylabels: {bool}
        show_xlabels: {bool}
        """
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (14, df.shape[1]/2)
        
        cmp = sns.clustermap(df.T,
                             figsize=figsize, 
                             cmap=cmap,
                             row_cluster = row_cluster,
                             col_cluster = col_cluster,
                             yticklabels = show_ylabels,
                             xticklabels = show_xlabels
                            )
        cmp.ax_heatmap.set_yticklabels(cmp.ax_heatmap.get_yticklabels(), fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle)
        cmp.ax_heatmap.set_xticklabels(cmp.ax_heatmap.get_xticklabels(), fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle)
        cmp.ax_heatmap.set_xlabel("")
        return None
    
    def draw_boxplot(self,
                     dtype="Fraction",
                     figsize = None,
                     ax = None,
                     fontstyle = 'italic',
                     xlabel = None,
                     fontsize = 10):
        """
        Description 
        -----------
        
        Draw boxplot of count or Fraction.
        
        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        figsize: tuple of (width, height), optional
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        fontstyle: {'normal', 'italic', 'oblique'}
        fontsize: 10
        """
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (4, df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        if xlabel == None and dtype == "Fraction":
            xlabel = "Fraction"
        if xlabel == None and dtype == "count":
            xlabel = "Count"
        ax.set_xlabel(xlabel)
        return None
    
    def draw_pearson_heatmap(self,
                             dtype="Fraction",
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='italic',
                             labels_fontsize=10,
                             ax = None
                            ):
        """
        Description 
        -----------
        
        Draw pearson heatmap of count or Fraction.
        
        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        labels_fontstyle: {'normal', 'italic', 'oblique'}
        labels_fontsize: 10
        """
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (df.shape[1]/4+1, df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        corr = df.corr(method='pearson')
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, fmt=".2f", cmap="Blues", ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        return None
    
    def draw_COA_plot(self,
                      dtype="Fraction",
                      figsize=(8,8),
                      species_labels_color="black", 
                      species_labels_style="italic", 
                      species_labels_size = 8, 
                      species_shapes_color = None, 
                      species_shapes_size = 200, 
                      species_labels_ha = "left", 
                      species_labels_va = "bottom",
                      amino_acid_labels_color = "black", 
                      amino_acid_labels_size = 8, 
                      amino_acid_shapes_size = 100,
                      amino_acid_labels_ha = "left",
                      amino_acid_labels_va = "bottom",
                      show_species_labels = True,
                      show_amino_acid_labels = True,
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      ax=None):
        """
        Description
        ----------
        Draw correspondence analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        show_species_labels: {bool, ["species1", "species2", ...]} show species name in plot.
        show_amino_acid_labels: {bool, ["codon1", "codon2", ...]} show amino acid in plot.
        figsize: (8,8)
        species_labels_color: "black"
        species_labels_style: {'normal', 'italic', 'oblique'}
        species_labels_size: 8,
        species_shapes_color: {None, {"species1": "red", "species2":"blue", ... }}
        species_shapes_size: 200,
        species_labels_ha: {"left", "right", "top", "bottom", "center"}
        species_labels_va: {"left", "right", "top", "bottom", "center"}
        amino_acid_labels_color: "black",
        amino_acid_labels_size: 8,
        amino_acid_shapes_size: 100,
        amino_acid_labels_ha: {"left", "right", "top", "bottom", "center"}
        amino_acid_labels_va: {"left", "right", "top", "bottom", "center"}
        title {None, str}
        xlabel {None, str}
        ylabel {None, str}
        title_size: 12
        xlabel_size: 12
        ylabel_size: 12
        ax: {None, Aexs}
        """
        if dtype == "Fraction":
            df = self.AA_Fraction_df
            ca = self.AA_Fraction_ca
            row_coords = ca.row_coordinates(df)
            col_coords = ca.column_coordinates(df)
        else:
            df = self.AA_count_df
            ca = self.AA_count_ca
            row_coords = ca.row_coordinates(df)
            col_coords = ca.column_coordinates(df)
            
        row_coords = row_coords.loc[self.AAColorShapeType.keys(), :]
        
        # Plotting the results
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)

        # Adding labels
        if show_amino_acid_labels !=None and show_amino_acid_labels !=False:
            if isinstance(show_amino_acid_labels, list):
                for i, amino_acid in enumerate(df.index):
                    if amino_acid in show_amino_acid_labels:
                        ax.annotate(amino_acid, 
                                    (row_coords[0][i], row_coords[1][i]), 
                                    color=amino_acid_labels_color, 
                                    fontsize=amino_acid_labels_size, 
                                    ha=amino_acid_labels_ha, 
                                    va=amino_acid_labels_va) 
            else:
                for i, amino_acid in enumerate(df.index):
                    ax.annotate(amino_acid, 
                                (row_coords[0][i], row_coords[1][i]), 
                                color=amino_acid_labels_color,
                                fontsize=amino_acid_labels_size, 
                                ha=amino_acid_labels_ha, 
                                va=amino_acid_labels_va)

        if show_species_labels != None and show_species_labels !=False:
            if isinstance(show_species_labels, list):
                for i, spceies in enumerate(df.columns):
                    if spceies in show_species_labels:
                        ax.annotate(spceies, 
                                    (col_coords[0][i], col_coords[1][i]), 
                                    color=species_labels_color, 
                                    fontsize=species_labels_size, 
                                    style = species_labels_style,
                                    ha = species_labels_ha,
                                    va = species_labels_va)
            else:
                for i, spceies in enumerate(df.columns):
                    ax.annotate(spceies, 
                                (col_coords[0][i], col_coords[1][i]), 
                                color=species_labels_color,
                                fontsize=species_labels_size,
                                style = species_labels_style,
                                ha = species_labels_ha,
                                va = species_labels_va)
        # Adding points            
        Pa_list = []
        for i in range(0,len(row_coords)):
            Pa = ax.scatter(row_coords[0][i], row_coords[1][i], 
                            c=self.AAColorShapeType[row_coords.index[i]][0], 
                            label=row_coords.index[i],
                            marker=self.AAColorShapeType[row_coords.index[i]][1],
                            s=amino_acid_shapes_size)
            Pa_list.append(Pa)

        Ps_list = []
        for i in range(0, len(col_coords)):
            if isinstance(species_shapes_color, dict):
                Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                c=species_shapes_color.get(col_coords.index[i],"#E64B35FF"), 
                                label=col_coords.index[i],
                                marker="*", 
                                s=species_shapes_size)
                Ps_list.append(Ps)
            else:
                Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                 c='#E64B35FF', 
                                 label='Species', 
                                 marker="*", 
                                 s=species_shapes_size)
                if i == 0:
                    Ps_list.append(Ps)

        # Adding legend
        ax.legend(handles=[*tuple(Pa_list), *tuple(Ps_list)],
                  loc='upper left', 
                  bbox_to_anchor=(1, 1),
                  ncol=1,
                  frameon=False,
                  shadow=False,
                  title='')

        # Adding title, xlabel, and ylabel
        if title == None:
            title = 'Correspondence analysis of amino acid composition'
        if xlabel == None:
            xlabel = f'Dim 1 ({ca.eigenvalues_summary.iloc[0,1]})'
        if ylabel == None:
            ylabel = f'Dim 2 ({ca.eigenvalues_summary.iloc[1,1]})'
        ax.set_title(title, size=title_size)
        ax.set_xlabel(xlabel, size=xlabel_size)
        ax.set_ylabel(ylabel, size=ylabel_size)
        return None
    
    
    def draw_PCA_plot(self,
                      dtype="Fraction",
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="italic", 
                      labels_size = 8, 
                      shapes_color = '#E64B35FF', 
                      shapes_type = '*',
                      shapes_size = 100, 
                      labels_ha = "left", 
                      labels_va = "bottom",
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      show_labels = True,
                      show_legend = False,
                      ax=None):
        """
        Description 
        ----------
        Draw principal component analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {Fraction, count} Choose to use amino acid composition count or Fraction.
        show_labels: {bool, ["species1", "species2", ...]} show column names in plot.
        show_legend: {bool}
        figsize: (6,6)
        labels_color: black
        labels_style: {'normal', 'italic', 'oblique'}
        labels_size: 8,
        shapes_color: {None, {"species1": "red", "species2":"blue", ... }}
        shapes_type: {None, {"species1": "*", "species2":"<", ... }}
        shapes_size: 100,
        labels_ha: {"left", "right", "top", "bottom", "center"}
        labels_va: {"left", "right", "top", "bottom", "center"}
        title {None, str}
        xlabel {None, str}
        ylabel {None, str}
        ax: {None, Aexs}
        """
        if dtype == "Fraction":
            pca = self.AA_Fraction_pca
        else:
            pca = self.AA_count_pca
            
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        for z, i in zip(self.get_shape_and_color(pca.column_correlations.shape[0]),
                        range(pca.column_correlations.shape[0])):
            x = pca.column_correlations[0][i]
            y = pca.column_correlations[1][i]
            label = pca.column_correlations.index[i]
            shape_type=z[0]
            shape_color=z[1]
            
            if shapes_type !=None:
                if isinstance(shapes_type, dict):
                    shape_type = shapes_type.get(label, z[0])
                else:
                    shape_type = shapes_type
            
            if shapes_color !=None:
                if isinstance(shapes_color, dict):
                    shape_color = shapes_color.get(label, z[1])
                else:
                    shape_color = shapes_color
            
            ax.scatter(x,y, label=label, 
                       marker=shape_type, 
                       color=shape_color,
                       s=shapes_size)
            
            if isinstance(show_labels, list):
                if label in show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)
            else:
                if show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)
        
        if title == None:
            title = f'Principal component analysis'
        if xlabel == None:
            xlabel = f"Dim1 ({pca.eigenvalues_summary.iloc[0,1]})"
        if ylabel == None:
            ylabel = f"Dim2 ({pca.eigenvalues_summary.iloc[1,1]})"
        
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        if show_legend == True:
            ax.legend(loc='center left',  
                      bbox_to_anchor=(1, 0.5),
                      ncol=1,
                      frameon=False,
                      shadow=False,
                      title='')
        return None
    
    def get_shape_and_color(self, length):
        shapes = ["o", "s", "^", "<", ">", "v", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "8"]
        colors = ["#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF"]
        n = len(shapes) * len(colors)
        l = []
        for x in shapes:
            for y in colors:
                l.append((x, y))
        res = []
        for i in range(0, length):
            res.append(l[i%n])
        return res
    
class StopCodon_analysis():
    def __init__(self, data, genetic_code, incomplete_codon=True):
        StopCodon_Count_dict = {}
        GeneName2StopCodon_dict = {}
        for prefix, file in data:
            StopCodon_Count = defaultdict(int)
            GeneName2StopCodon = {}
            for ID, Seq in FastaIO(file):
                if len(Seq) % 3 == 0:
                    codon = Seq[-3:]
                    if translate(codon, genetic_code=genetic_code) == "*":
                        stop_codon = codon
                        StopCodon_Count[stop_codon] += 1
                        GeneName2StopCodon.setdefault(ID, stop_codon)
                        
                elif len(Seq) % 3 == 1:
                    codon = Seq[-1:]
                    if incomplete_codon==True:
                        stop_codon = codon
                        StopCodon_Count[stop_codon] += 1
                        GeneName2StopCodon.setdefault(ID, stop_codon)
                        
                elif len(Seq) % 3 == 2:
                    codon = Seq[-2:]
                    if incomplete_codon==True:
                        stop_codon = codon
                        StopCodon_Count[stop_codon] += 1
                        GeneName2StopCodon.setdefault(ID, stop_codon)

            StopCodon_Count_dict.setdefault(prefix, StopCodon_Count)
            GeneName2StopCodon_dict.setdefault(prefix, GeneName2StopCodon)
        self.GeneName2StopCodon_dict = GeneName2StopCodon_dict
        StopCodon_Count_df = pd.DataFrame(StopCodon_Count_dict)
        StopCodon_Count_df = StopCodon_Count_df.replace(np.NAN, 0)
        self.StopCodon_Count_df = StopCodon_Count_df.astype('int64')
        self.StopCodon_Count_df.index.name = "Stop Codon"
        
        pca = prince.PCA(n_components=2)
        self.pca = pca.fit(self.StopCodon_Count_df)
        
    def get_df(self):
        return pd.melt(self.StopCodon_Count_df.reset_index(), id_vars="Stop Codon", value_name="Count", var_name="Group")
        
    def draw_barplot(self, figsize=None, ax=None, palette="Blues", legend_font_style="italic", legend_font_size=10):
        """
        Description 
        ----------
        Draw barplot of Stop codon count.
        
        Parameters
        ----------
        figsize: {None, (float, float)}
        palette: {"Blues", "Set1", "Set2" ...}
        legend_font_style: {'normal', 'italic', 'oblique'}
        legend_font_size: 10
        ax: {None, Aexs}
        """
        df = self.get_df()
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(x=df["Stop Codon"], y=df["Count"], hue=df["Group"], hue_order=sorted(set(df["Group"])), palette=palette, ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, shadow=False, title='', prop={'style':legend_font_style, 'size':legend_font_size})
        return None
    
    def draw_clustermap(self, 
                        figsize=None, 
                        cmap="Blues",
                        ylabels_fontstyle='italic',
                        xlabels_fontstyle='normal',
                        ylabels_fontsize=10,
                        xlabels_fontsize=10,
                        row_cluster=True,
                        col_cluster=True,
                        show_ylabels=True,
                        show_xlabels=True):
        """
        Description 
        -----------
        
        Draw clustermap of count or Fraction.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        cmap : matplotlib colormap name or object, or list of colors, optional
               The mapping from data values to color space. If not provided, the
               default will depend on whether ``center`` is set.
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}
        ylabels_fontsize: 10
        xlabels_fontsize: 10
        show_ylabels: {bool}
        show_xlabels: {bool}
        """
        df = self.StopCodon_Count_df
        cmp = sns.clustermap(df.T,
                             figsize=figsize, 
                             cmap=cmap,
                             row_cluster = row_cluster,
                             col_cluster = col_cluster,
                             yticklabels = show_ylabels,
                             xticklabels = show_xlabels
                            )
        cmp.ax_heatmap.set_yticklabels(cmp.ax_heatmap.get_yticklabels(), fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle)
        cmp.ax_heatmap.set_xticklabels(cmp.ax_heatmap.get_xticklabels(), fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle)
        cmp.ax_heatmap.set_xlabel("")
        return None
    
    def draw_boxplot(self,
                     figsize = None,
                     ax = None,
                     fontstyle = 'italic',
                     fontsize = 10):
        """
        Description 
        -----------
        
        Draw boxplot of stop codon count.
        
        Parameters
        ----------
        figsize: tuple of (width, height), optional
        ax : matplotlib Axes, optional 
             Axes in which to draw the plot, otherwise use the currently-active Axes.
        fontstyle: {'normal', 'italic', 'oblique'}
        fontsize: 10
        """
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(self.StopCodon_Count_df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        ax.set_xlabel("Count")
        return None
    
    
    def draw_PCA_plot(self,
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="italic", 
                      labels_size = 8, 
                      shapes_color = '#E64B35FF', 
                      shapes_type = '*',
                      shapes_size = 100, 
                      labels_ha = "left", 
                      labels_va = "bottom",
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      show_labels = True,
                      show_legend = False,
                      ax=None):
        """
        Description 
        ----------
        Draw principal component analysis of stop codon count plot.

        Parameters
        ----------
        show_labels: {bool, ["species1", "species2", ...]} show column names in plot.
        show_legend: {bool}
        figsize: (6,6)
        labels_color: black
        labels_style: {'normal', 'italic', 'oblique'}
        labels_size: 8,
        shapes_color: {None, {"species1": "red", "species2":"blue", ... }}
        shapes_type: {None, {"species1": "*", "species2":"<", ... }}
        shapes_size: 100,
        labels_ha: {"left", "right", "top", "bottom", "center"}
        labels_va: {"left", "right", "top", "bottom", "center"}
        title {None, str}
        xlabel {None, str}
        ylabel {None, str}
        ax: {None, Aexs}
        """
        pca = self.pca
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        for z, i in zip(self.get_shape_and_color(pca.column_correlations.shape[0]),
                        range(pca.column_correlations.shape[0])):
            x = pca.column_correlations[0][i]
            y = pca.column_correlations[1][i]
            label = pca.column_correlations.index[i]
            shape_type=z[0]
            shape_color=z[1]
            
            if shapes_type !=None:
                if isinstance(shapes_type, dict):
                    shape_type = shapes_type.get(label, z[0])
                else:
                    shape_type = shapes_type
            
            if shapes_color !=None:
                if isinstance(shapes_color, dict):
                    shape_color = shapes_color.get(label, z[1])
                else:
                    shape_color = shapes_color
            
            ax.scatter(x,y, label=label, 
                       marker=shape_type, 
                       color=shape_color,
                       s=shapes_size)
            
            if isinstance(show_labels, list):
                if label in show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)
            else:
                if show_labels:
                    ax.annotate(label, (x,y), 
                                color=labels_color,
                                fontsize=labels_size,
                                style = labels_style,
                                ha = labels_ha,
                                va = labels_va)
        if title == None:
            title = f'Principal component analysis of stop codon count'
        if xlabel == None:
            xlabel = f"Dim1 ({pca.eigenvalues_summary.iloc[0,1]})"
        if ylabel == None:
            ylabel = f"Dim2 ({pca.eigenvalues_summary.iloc[1,1]})"
        
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        if show_legend == True:
            ax.legend(loc='center left',  
                      bbox_to_anchor=(1, 0.5),
                      ncol=1,
                      frameon=False,
                      shadow=False,
                      title='')
        return None
    
    def get_shape_and_color(self, length):
        shapes = ["o", "s", "^", "<", ">", "v", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d", "8"]
        colors = ["#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF"]
        n = len(shapes) * len(colors)
        l = []
        for x in shapes:
            for y in colors:
                l.append((x, y))
        res = []
        for i in range(0, length):
            res.append(l[i%n])
        return res

class TreePlotter:
    """Phylogenetic Tree Plotter using scikit-bio's TreeNode"""
    
    def __init__(
        self,
        tree,
        height=0.5,
        width=8,
        ignore_branch_length=False,
        leaf_label_size= 12,
        innode_label_size= 0,
        ladderize= True,
        ladderize_by= "size",  # "size" or "branch_length"
        ladderize_direction= "right",  # "left" or "right"
        linewidth= 1
    ):
        """
        Parameters
        ----------
        tree : TreeNode
            Input tree (must be rooted)
        height : float
            Figure height per leaf node
        width : float
            Figure width
        ignore_branch_length : bool
            Ignore branch lengths for cladogram
        leaf_label_size : float
            Font size for leaf labels
        innode_label_size : float
            Font size for internal node labels
        ladderize : bool
            Enable ladderize tree sorting
        ladderize_by : str
            Sort criterion ("size" or "branch_length")
        ladderize_direction : str
            Direction for larger subtrees ("left" or "right")
        linewidth : float
            Branch line width
        """
        self.tree = tree.copy()
        self.ignore_branch_length = ignore_branch_length
        self.linewidth = linewidth
        self.leaf_label_size = leaf_label_size
        self.innode_label_size = innode_label_size
        
        # Ladderize tree if requested
        if ladderize:
            self._ladderize(
                self.tree, 
                by=ladderize_by, 
                direction=ladderize_direction
            )
        
        # Set unique names for internal nodes
        self._set_unique_node_names()
        
        # Calculate node positions
        self.node_positions = self._calc_node_positions()
        
        # Initialize plot containers
        self._plot_patches: List[Patch] = []
        self._plot_funcs: List[Callable[[Axes], None]] = []
        
        # Set figure size based on leaf count
        num_tips = len(list(self.tree.tips()))
        self.figsize = (width, height * num_tips)

    # =====================
    # Tree Processing Methods
    # =====================
    
    def _set_unique_node_names(self):
        """Set unique names for unnamed internal nodes"""
        for idx, node in enumerate(self.tree.non_tips(include_self=True)):
            if not node.is_tip() and not node.name:
                node.name = f"N_{idx+1}"
    
    def _ladderize(
        self, 
        tree, 
        by="size", 
        direction="right"
    ):
        """Recursive ladderization of tree structure
        
        Parameters
        ----------
        by : str
            Sorting criterion: "size" (subtree tip count) or 
            "branch_length" (subtree total length)
        direction : str
            Direction for larger subtrees: "left" or "right"
        """
        if not tree.children:
            return
        
        # Define sort key function
        def sort_key(child):
            if by == "size":
                return child.count(tips=True)  # Tip count
            elif by == "branch_length":
                # Calculate subtree total length
                return sum(n.length for n in child.traverse() if n.length)
            else:
                raise ValueError(f"Invalid 'by' value: {by}")
        
        # Sort children (larger subtrees on left/right)
        reverse = (direction == "left")
        tree.children.sort(key=sort_key, reverse=reverse)
        
        # Recurse into children
        for child in tree.children:
            self._ladderize(child, by, direction)
    
    def _calc_node_positions(self):
        """Calculate node coordinates (x, y) via depth-first traversal
        
        Returns
        -------
        dict
            Mapping of node names to (x, y) coordinates
        """
        positions = {}
        tips = list(self.tree.tips())
        num_tips = len(tips)
        y_positions = np.linspace(0, num_tips, num_tips, endpoint=False)
        y_index = 0
        
        def _layout(node, depth):
            """Recursive layout function"""
            nonlocal y_index
            
            # Calculate x position: depth or cumulative branch length
            if self.ignore_branch_length:
                x_pos = depth
            else:
                # Calculate cumulative distance from root
                if node.is_root():
                    x_pos = 0
                else:
                    # Traverse from node to root to sum branch lengths
                    current = node
                    x_pos = 0
                    while not current.is_root():
                        x_pos += current.length if current.length else 0
                        current = current.parent
                    
            # Handle tip nodes
            if node.is_tip():
                y_pos = y_positions[y_index]
                positions[node.name] = (x_pos, y_pos)
                y_index += 1
                return y_pos
            
            # Internal node: process children
            child_y_positions = []
            for child in node.children:
                child_y = _layout(child, depth + 1)
                child_y_positions.append(child_y)
            
            # Current node position is average of children's y positions
            node_y = np.mean(child_y_positions)
            positions[node.name] = (x_pos, node_y)
            return node_y
        
        # Start layout from root
        _layout(self.tree, depth=0)
        return positions

    # =====================
    # Plotting Methods
    # =====================
    
    def plot(
        self, 
        figsize=(8,8),
        ax= None):
        """Plot rectangular layout phylogenetic tree
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to plot on
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure
        """
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        
        ax.set_axis_off()
        # Draw branches
        for node in self.tree.traverse():
            if node.is_root():
                continue
                
            parent = node.parent
            x_node, y_node = self.node_positions[node.name]
            x_parent, y_parent = self.node_positions[parent.name]
            
            # Vertical branch (from parent level to child y-position)
            ax.plot([x_parent, x_parent], [y_parent, y_node], 
                    color='black', lw=self.linewidth)
            
            # Horizontal branch (to child position)
            ax.plot([x_parent, x_node], [y_node, y_node], 
                    color='black', lw=self.linewidth)
        
        # Add labels
        for node in self.tree.traverse():
            x, y = self.node_positions[node.name]
            if node.is_tip() and self.leaf_label_size > 0:
                ax.text(x + 0.02, y, node.name, 
                        fontsize=self.leaf_label_size, 
                        va='center', ha='left', fontstyle="italic")
            elif not node.is_tip() and self.innode_label_size > 0:
                ax.text(x, y, node.name, 
                        fontsize=self.innode_label_size, 
                        va='center', ha='center')
        
        # Add custom elements
        for patch in self._plot_patches:
            ax.add_patch(patch)
        for func in self._plot_funcs:
            func(ax)
            
        # Set dynamic axis limits
        self._set_axis_limits(ax)
        return None
    
    def _set_axis_limits(self, ax):
        """Set dynamic axis limits with padding"""
        all_x = [pos[0] for pos in self.node_positions.values()]
        all_y = [pos[1] for pos in self.node_positions.values()]
        
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        
        x_padding = (x_max - x_min) * 0.1
        y_padding = (y_max - y_min) * 0.1
        
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
    
    def highlight_clade(
        self, 
        node_names, 
        color, 
        alpha= 0.3,
        **kwargs
    ):
        """Highlight a clade with background rectangle
        
        Parameters
        ----------
        node_names : list of str
            Node names defining the clade
        color : str
            Fill color for the highlight
        alpha : float
            Transparency of the highlight (0-1)
        **kwargs
            Additional rectangle properties
        """
        # Find common ancestor of specified nodes
        mrca = self.tree.lowest_common_ancestor(node_names)
        
        # Get all tips in the clade
        tips = list(mrca.tips())
        x_vals = [self.node_positions[tip.name][0] for tip in tips]
        y_vals = [self.node_positions[tip.name][1] for tip in tips]
        
        # Create background rectangle
        rect = Rectangle(
            xy=(min(x_vals), min(y_vals)),
            width=max(x_vals) - min(x_vals),
            height=max(y_vals) - min(y_vals),
            color=color,
            alpha=alpha,
            zorder=-10,
            **kwargs
        )
        self._plot_patches.append(rect)
    
    def add_node_label(
        self, 
        node_name, 
        label, 
        size= 8,
        **kwargs
    ):
        """Add custom text label to a specific node
        
        Parameters
        ----------
        node_name : str
            Name of the node to label
        label : str
            Text to display
        size : int
            Font size
        **kwargs
            Additional text properties
        """
        def _add_label(ax: plt.Axes):
            if node_name not in self.node_positions:
                raise ValueError(f"Node '{node_name}' not found in tree")
                
            x, y = self.node_positions[node_name]
            ax.text(x, y, label, size=size, **kwargs)
            
        self._plot_funcs.append(_add_label)
    
    def add_scale_bar(
        self, 
        length, 
        label, 
        position= (0.05, 0.05),
        **kwargs):
        """Add evolutionary distance scale bar
        
        Parameters
        ----------
        length : float
            Length to represent (in branch length units)
        label : str
            Text label for the scale
        position : tuple (x, y)
            Position in axes coordinates (0-1)
        **kwargs
            Additional line/text properties
        """
        def _add_scale(ax: plt.Axes):
            # Convert from axes to data coordinates
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()
            
            x_data = position[0] * (x_max - x_min) + x_min
            y_data = position[1] * (y_max - y_min) + y_min
            
            # Draw scale bar
            ax.plot([x_data, x_data + length], [y_data, y_data], 
                    color='black', lw=2, **kwargs)
            
            # Add label
            ax.text(x_data + length/2, y_data - (y_max-y_min)*0.02, 
                    label, ha='center', va='top', **kwargs)
                    
        self._plot_funcs.append(_add_scale)
