#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2025/07/05"
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
           "get_ENC",
           "get_PR2",
           "get_base_phase_synonymous",
           "get_ATGC_Indices",
           "get_Relative_Adaptiveness", 
           "get_emboss_cutfile_from_Obs", 
           "get_Obs_from_emboss_cutfile",
           "get_Obs_from_CUBE_file",
           "get_codonw_caifile_from_Obs",
           "get_optimal_codons_from_codonw_coafile",
           "get_optimal_codons_from_ENC",
           "draw_codon_barplot", 
           "draw_codon_optimization_plot",
           "get_cusp_like",
           "get_codonW_like",
           "find_four_codon_AA",
           "infer_genetic_code_from_Obs",
           "TreePlotter",
           "NPA_Analysis",
           "ENC_Analysis",
           "PR2_Analysis",
           "Sequence_Indices_Analysis",
           "RSCU_Single_Species_Analysis",
           "RSCU_Multiple_Species_Analysis",
           "AA_Composition_Single_Species_Analysis",
           "AA_Composition_Multiple_Species_Analysis",
           "Stop_Codon_Analysis"]

__version__ = "v2.00"

import os
import io
import re
import sys
import math
import prince
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as ss
import urllib.request
from .fastaio import fastaIO
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from collections import defaultdict
from skbio import DistanceMatrix
from skbio.tree import TreeNode, nj, upgma, gme, bme
from scipy.spatial.distance import pdist, squareform
from .codontables import CodonTables, Seq3toSeq1, CBI_and_Fop_preset, CAI_preset
#mpl.rcParams['svg.fonttype'] = 'none'
#mpl.rcParams['svg.hashsalt'] = 'hello'

def translate(cds, genetic_code):
    """
    Description
    ----------
    Translate the CDS sequence into a protein sequence.
    
    Parameters
    ----------
    cds: str
        A CDS sequence.
      
    genetic_code: int
        A genetic code id, use `pycubs.CodonTables()` for more details.
    """
    
    protein = ""
    codotable = CodonTables().get(genetic_code, aaseq3=False)
    for i in range(0, len(cds), 3):
        protein += codotable.get(cds.upper()[i:i+3], "X")
    return protein
    
def get_Obs(seqences, genetic_code, aaseq3=True):
    """
    Description 
    ----------
    Calculate Observed number of occurrences of codon, and return a dict object.
    
    Parameters
    ----------
    seqences: str, list
        A sequence string, or a sequences list, or a fasta or fasta.gz format file path.
    
    genetic_code: int
        A genetic code id, use `pycubs.CodonTables()` for more details.
    
    aaseq3: bool, default=True
        If the value is True, the amino acid uses a three-letter code.
    """
    
    codontable = CodonTables().get(genetic_code, aaseq3)
    
    def _fastaIO(file):
        for ID, Seq in fastaIO(file):
            yield Seq
    
    if isinstance(seqences, str):
        if os.path.exists(seqences):
            seqences = _fastaIO(seqences)
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

def get_Fraction(Obs):
    """
    Description 
    ----------
    Calculate Fraction of codon.
    Fraction represents the proportion of each codon in the codon encoding the amino acid, i.e.
    Fraction = the number of occurrences of a codon/the number of occurrences of all codons of the amino acid encoded by the codon.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
    Reference
    ----------
    [1] Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    [2] Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
    """
    
    class Fraction():
        def __init__(self, obj):
            self.Fraction_dict = obj
            
        def draw_barplot(self,
                         figsize=(8,4),
                         ylabel='Fraction',
                         title=None,
                         width=0.9,
                         palette=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         remove_stop_codon=True,
                         codon_space = 0.16,
                         ax=None,
                         outfile=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            figsize: tuple, default=(8,4)
                Figure size.
            
            ylabel: str, default=None
                Y-axis label of figure.
            
            title: str, default=None
                Title of figure.
            
            width: float, default=0.9
                The bar spacing width.
            
            palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
                Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
                Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
            
            remove_stop_codon: bool, default=None
                If remove stop codon.
            
            ax: matplotlib Axes, default=None
                Axes object to draw the plot onto, otherwise uses the current Axes.   
            
            codon_space: float, default=0.16
                codon spacing.
            
            outfile: str, default=None
                A path of outfile.
            """
            draw_codon_barplot(self.Fraction_dict, ylabel=ylabel, title=title, 
                               palette=palette, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space, outfile=outfile)
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
    Obs: dict
        get_Obs() function return value.
    
    Reference
    ----------
    [1] Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    [2] Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
    """
    
    class Frequency():
        def __init__(self, obj):
            self.Frequency_dict = obj
            
        def draw_barplot(self,
                         ylabel='Frequency',
                         title=None,
                         palette=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space=0.16,
                         ax=None,
                         outfile=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            figsize: tuple, default=(8,4)
                Figure size.
            
            ylabel: str, default=None
                Y-axis label of figure.
            
            title: str, default=None
                Title of figure.
            
            width: float, default=0.9
                The bar spacing width.
            
            palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
                Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
                Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
            
            remove_stop_codon: bool, default=None
                If remove stop codon.
            
            ax: matplotlib Axes, default=None
                Axes object to draw the plot onto, otherwise uses the current Axes.   
            
            codon_space: float, default=0.16
                codon spacing.
            
            outfile: str, default=None
                A path of outfile.
            """
            draw_codon_barplot(self.Frequency_dict, ylabel=ylabel, title=title, 
                               palette=palette, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space, outfile=outfile)
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

def get_Relative_Adaptiveness(Obs):
    """
    Description
    ----------
    The relative adaptiveness (w) of each codon is the ratio of the usage of each codon,
    to that of the most abundant codon for the same amino acid. 
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
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
                         palette=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space = 0.16,
                         ax=None, 
                         outfile=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            figsize: tuple, default=(8,4)
                Figure size.
            
            ylabel: str, default=None
                Y-axis label of figure.
            
            title: str, default=None
                Title of figure.
            
            width: float, default=0.9
                The bar spacing width.
            
            palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
                Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
                Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
            
            remove_stop_codon: bool, default=None
                If remove stop codon.
            
            ax: matplotlib Axes, default=None
                Axes object to draw the plot onto, otherwise uses the current Axes.   
            
            codon_space: float, default=0.16
                Codon spacing.
            
            outfile: str, default=None
                A path of outfile.
            """
            draw_codon_barplot(self.Relative_Adaptiveness_dict, ylabel=ylabel, title=title, 
                               palette=palette, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space, outfile=outfile)
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

def get_RSCU(Obs):
    """
    Description
    ----------
    Calculate relative synonymous codon usage (RSCU).
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
    Reference
    ----------
    [1] Sharp, Paul M., Therese MF Tuohy, and Krzysztof R. Mosurski. "Codon usage in yeast: cluster analysis clearly differentiates highly 
        and lowly expressed genes." Nucleic acids research 14.13 (1986): 5125-5143.
    """
    
    class RSCU():
        def __init__(self, obj):
            self.RSCU_dict = obj
            
        def draw_barplot(self,
                         ylabel='RSCU',
                         title=None,
                         palette=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                         width=0.9,
                         remove_stop_codon=True,
                         figsize=(8,4),
                         codon_space = 0.16,
                         ax=None,
                         outfile=None):
            """
            Description
            ----------
            Draw a codons barplot.

            Parameters
            ----------
            figsize: tuple, default=(8,4)
                Figure size.
            
            ylabel: str, default=None
                Y-axis label of figure.
            
            title: str, default=None
                Title of figure.
            
            width: float, default=0.9
                The bar spacing width.
            
            palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
                Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
                Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
            
            remove_stop_codon: bool, default=None
                If remove stop codon.
            
            ax: matplotlib Axes, default=None
                Axes object to draw the plot onto, otherwise uses the current Axes.   
            
            codon_space: float, default=0.16
                Codon spacing.
            
            outfile: str, default=None
            """
            draw_codon_barplot(self.RSCU_dict, ylabel=ylabel, title=title, 
                               palette=palette, width=width,
                               remove_stop_codon=remove_stop_codon,
                               figsize=figsize, ax=ax, codon_space=codon_space, outfile=outfile)
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

def draw_codon_barplot(obj, 
                       ylabel=None,
                       title=None,
                       palette=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"],
                       width=0.9,
                       remove_stop_codon=True,
                       figsize=(8,4),
                       codon_space=0.16,
                       ax=None,
                       outfile=None):
    """
    Description
    ----------
    Draw a codons barplot.

    Parameters
    ----------
    obj: dict, RSCU, Fraction, Frequency, Relative_Adaptiveness
        Return value of get_Obs(), get_Fraction(), get_Frequency(), get_RSCU() or get_Relative_Adaptiveness() function.
        
    figsize: tuple, default=(8,4)
        Figure size.
    
    ylabel: str, default=None
        Y-axis label of figure.
    
    title: str, default=None
        Title of figure.
        
    width: float, default=0.9
        The bar spacing width.
        
    palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
        Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
        Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
    
    remove_stop_codon: bool, default=None
        If remove stop codon.
    
    ax: matplotlib Axes, default=None
        Axes object to draw the plot onto, otherwise uses the current Axes.   
    
    codon_space: float, default=0.16
        codon spacing.
    
    outfile: str, default=None
        A path of outfile.
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
    
    if isinstance(palette, str):
        cols = plt.colormaps.get_cmap(palette).colors
    else:
        cols = palette
    
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
    if outfile != None:
        fig.savefig(outfile, dpi=600)
    return None

def get_ATGC_Indices(Obs):
    """
    Description
    ----------
    Calculate 58 indices, including the important GC3s, GC12, GC3, etc.
    
    A: A nucleotides content
    T: T nucleotides content
    G: G nucleotides content
    C: C nucleotides content
    GC: G and C nucleotides content
    AT: A and T nucleotides content
    GC-skew: GC skew
    AT-skew: AT skew
    A1: A nucleotides content at first codon postion
    T1: T nucleotides content at first codon postion
    G1: G nucleotides content at first codon postion
    C1: C nucleotides content at first codon postion
    GC1: G and C nucleotides content at first codon postion
    AT1: A and T nucleotides content at first codon postion
    GC1-skew: GC skew at first codon postion
    AT1-skew: AT skew at first codon postion
    A2: A nucleotides content at second codon postion
    T2: T nucleotides content at second codon postion
    G2: G nucleotides content at second codon postion
    C2: C nucleotides content at second codon postion
    GC2: G and C nucleotides content at second codon postion
    AT2: A and T nucleotides content at second codon postion
    GC2-skew: GC skew at second codon postion
    AT2-skew: AT skew at second codon postion
    A3: A nucleotides content at third codon postion
    T3: T nucleotides content at third codon postion
    G3: G nucleotides content at third codon postion
    C3: C nucleotides content at third codon postion
    GC3: G and C nucleotides content at third codon postion
    AT3: A and T nucleotides content at third codon postion
    GC3-skew: GC skew at third codon postion
    AT3-skew: AT skew at third codon postion
    GC12: G and C nucleotides content at first + second codon postions
    A1s: A content at the first synonymous codon position
    T1s: T content at the first synonymous codon position
    G1s: G content at the first synonymous codon position
    C1s: C content at the first synonymous codon position
    GC1s: G and C content at the first synonymous codon position
    AT1s: A and T content at the first synonymous codon position
    A2s: A content at the second synonymous codon position
    T2s: T content at the second synonymous codon position
    G2s: G content at the second synonymous codon position
    C2s: C content at the second synonymous codon position
    GC2s: G and C content at the second synonymous codon position
    AT2s: A and T content at the second synonymous codon position
    A3s: A content at the third synonymous codon position
    T3s: T content at the third synonymous codon position
    G3s: G content at the third synonymous codon position
    C3s: C content at the third synonymous codon position
    GC3s: G and C content at the third synonymous codon position
    AT3s: G and C content at the third synonymous codon position
    GC12s: G and C content at the first + second synonymous codon position
    L_sym: Number of silent sites
    L_aa: Number of amino acids
    
    The following four indices have codow definitions, It quantifies the usage 
    of each nucleotide at synonymous third codon positions as a proportion 
    of the maximum usage of that nucleotide could have without altering the amino acid composition.
    
    A3s codonW: A3s with codonW definition
    T3s codonW: T3s with codonW definition
    G3s codonW: G3s with codonW definition
    C3s codonW: C3s with codonW definition
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    """
    ATGC1 = defaultdict(int)
    ATGC2 = defaultdict(int)
    ATGC3 = defaultdict(int)    
    ATGC1s = defaultdict(int)
    ATGC2s = defaultdict(int)
    ATGC3s = defaultdict(int)
    for aa in Obs:
        for c in Obs[aa]:
            ATGC1[c[0]] += Obs[aa][c]
            ATGC2[c[1]] += Obs[aa][c]
            ATGC3[c[2]] += Obs[aa][c]
        if aa == "*":
            continue
        elif len(Obs[aa]) < 2:
            continue
        else:
            for c in Obs[aa]:
                ATGC1s[c[0]] += Obs[aa][c]
                ATGC2s[c[1]] += Obs[aa][c]
                ATGC3s[c[2]] += Obs[aa][c]
                
    A1, T1, G1, C1 = ATGC1.get("A", 0), ATGC1.get("T", 0), ATGC1.get("G", 0), ATGC1.get("C", 0)
    A2, T2, G2, C2 = ATGC2.get("A", 0), ATGC2.get("T", 0), ATGC2.get("G", 0), ATGC2.get("C", 0)
    A3, T3, G3, C3 = ATGC3.get("A", 0), ATGC3.get("T", 0), ATGC3.get("G", 0), ATGC3.get("C", 0)
    ATGC1_Total,  ATGC2_Total, ATGC3_Total = sum(ATGC1.values()), sum(ATGC2.values()), sum(ATGC3.values())
    GC1_skew, AT1_skew = (G1 - C1) / (G1 + C1), (A1 - T1) / (A1 + T1)
    GC2_skew, AT2_skew = (G2 - C2) / (G2 + C2), (A2 - T2) / (A2 + T2)
    GC3_skew, AT3_skew = (G3 - C3) / (G3 + C3), (A3 - T3) / (A3 + T3)
    
    if ATGC1_Total!=0:
        GC1, AT1 = (G1 + C1)/ATGC1_Total, (A1 + T1)/ATGC1_Total
    else:
        GC1, AT1 = None, None
    if ATGC2_Total!=0:
        GC2, AT2 = (G2 + C2)/ATGC2_Total, (A2 + T2)/ATGC2_Total
    else:
        GC2, AT2 = None, None
    if ATGC3_Total!=0:
        GC3, AT3 = (G3 + C3)/ATGC3_Total, (A3 + T3)/ATGC3_Total
    else:
        GC3, AT3 = None, None
    if GC1 != None and GC2 !=None:
        GC12 = (GC1 + GC2)/2
    else:
        GC12 = None
        
    A1s, T1s, G1s, C1s = ATGC1s.get("A", 0), ATGC1s.get("T", 0), ATGC1s.get("G", 0), ATGC1s.get("C", 0)
    A2s, T2s, G2s, C2s = ATGC2s.get("A", 0), ATGC2s.get("T", 0), ATGC2s.get("G", 0), ATGC2s.get("C", 0)
    A3s, T3s, G3s, C3s = ATGC3s.get("A", 0), ATGC3s.get("T", 0), ATGC3s.get("G", 0), ATGC3s.get("C", 0)
    ATGC1s_Total, ATGC2s_Total, ATGC3s_Total = sum(ATGC1s.values()), sum(ATGC2s.values()), sum(ATGC3s.values())
    
    if ATGC1s_Total!=0:
        GC1s, AT1s = (G1s + C1s)/ATGC1s_Total, (A1s + T1s)/ATGC1s_Total
    else:
        GC1s, AT1s = None, None
    if ATGC2s_Total!=0:
        GC2s, AT2s = (G2s + C2s)/ATGC2s_Total, (A2s + T2s)/ATGC2s_Total
    else:
        GC2s, AT2s = None, None
    if ATGC3s_Total!=0:
        GC3s, AT3s = (G3s + C3s)/ATGC3s_Total, (A3s + T3s)/ATGC3s_Total
    else:
        GC3s, AT3s = None, None
    if GC1s != None and GC2s !=None:
        GC12s = (GC1s + GC2s)/2
    else:
        GC12s = None
    
    L_sym = ATGC3s_Total
    L_aa = sum([Obs[aa][c] for aa in Obs if aa != "*" for c in Obs[aa]])

    A = ATGC1.get('A', 0) + ATGC2.get("A", 0) + ATGC3.get("A", 0)
    T = ATGC1.get('T', 0) + ATGC2.get("T", 0) + ATGC3.get("T", 0)
    G = ATGC1.get('G', 0) + ATGC2.get("G", 0) + ATGC3.get("G", 0)
    C = ATGC1.get('C', 0) + ATGC2.get("C", 0) + ATGC3.get("C", 0)
    Total_Base = A + T + G + C
    
    A3s_codonW = get_base_phase_synonymous(Obs, base="A", phase=2)
    T3s_codonW = get_base_phase_synonymous(Obs, base="T", phase=2)
    G3s_codonW = get_base_phase_synonymous(Obs, base="G", phase=2)
    C3s_codonW = get_base_phase_synonymous(Obs, base="C", phase=2)
    
    return {"A":A / Total_Base,
            "T":T / Total_Base,
            "G":G / Total_Base,
            "C":C / Total_Base,
            "GC":(G + C) / Total_Base,
            "AT":(A + T) / Total_Base,
            "GC-skew":(G - C) / (G + C),
            "AT-skew":(A - T) / (A + T),
            "A1":A1 / ATGC1_Total,
            "T1":T1 / ATGC1_Total,
            "G1":G1 / ATGC1_Total,
            "C1":C1 / ATGC1_Total,
            "GC1":GC1,
            "AT1":AT1,
            "GC1-skew":GC1_skew,
            "AT1-skew":AT1_skew,
            "A2":A2 / ATGC2_Total,
            "T2":T2 / ATGC2_Total,
            "G2":G2 / ATGC2_Total,
            "C2":C2 / ATGC2_Total,
            "GC2":GC2,
            "AT2":AT2,
            "GC2-skew":GC2_skew,
            "AT2-skew":AT2_skew,
            "A3":A3 / ATGC3_Total,
            "T3":T3 / ATGC3_Total,
            "G3":G3 / ATGC3_Total,
            "C3":C3 / ATGC3_Total,
            "GC3":GC3,
            "AT3":AT3,
            "GC3-skew":GC3_skew,
            "AT3-skew":AT3_skew,
            "GC12":GC12,
            "A1s":A1s / ATGC1s_Total,
            "T1s":T1s / ATGC1s_Total,
            "G1s":G1s / ATGC1s_Total,
            "C1s":C1s / ATGC1s_Total,
            "GC1s":GC1s,
            "AT1s":AT1s,
            "A2s":A2s / ATGC2s_Total,
            "T2s":T2s / ATGC2s_Total,
            "G2s":G2s / ATGC2s_Total,
            "C2s":C2s / ATGC2s_Total,
            "GC2s":GC2s,
            "AT2s":AT2s,
            "A3s":A3s / ATGC3s_Total,
            "T3s":T3s / ATGC3s_Total,
            "G3s":G3s / ATGC3s_Total,
            "C3s":C3s / ATGC3s_Total,
            "GC3s":GC3s,
            "AT3s":AT3s,
            "GC12s":GC12s,
            "A3s codonW": A3s_codonW,
            "T3s codonW": T3s_codonW,
            "G3s codonW": G3s_codonW,
            "C3s codonW": C3s_codonW,
            "L_sym": L_sym,
            "L_aa": L_aa
           }

def get_codonw_caifile_from_Obs(Obs, outfile=None):
    """
    Description
    ----------
    A set of high expression genes set calculation of Obs, converted to condow software recognition of the cai file.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
    outfile: str, default=None
        A path of outfile.
        
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
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
    To create a emboss cutfile from the get_Obs() function return value.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
    outfile: str, default=None
        A path of outfile.
    
    Reference
    ----------
    [1] Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    [2] Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
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
    Gets the Obs object from the emboss cut file with the same value as the get_Obs() function returns.
    
    Parameters
    ----------
    file: str
        A emboss cut file.
    
    aaseq3: bool, default=True
        If the value is True, the amino acid uses a three-letter code.
    
    Reference
    ----------
    [1] Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    [2] Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
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

def get_Obs_from_CUBE_file(file, aaseq3=True):
    """
    Description
    ----------
    Gets the Obs object from the CUBE (https://www.codonbias.cn/) file with the same value as the get_Obs() function returns.
    
    Parameters
    ----------
    file: str
        A CUBE file from https://www.codonbias.cn/download or https://github.com/dyyvgug/CUBE/tree/master/resource/RSCU
        
    aaseq3: bool, default=True
        If the value is True, the amino acid uses a three-letter code.
    
    Example
    ----------
    ref_Obs = get_Obs_from_CUBE_file("https://github.com/dyyvgug/CUBE/blob/master/resource/RSCU/Abditibacterium_utsteinense")
    """
    
    Seq1toSeq3 = {Seq3toSeq1[k]:k for k in Seq3toSeq1}
    Obs = {}
    
    if re.search("https://github.com/dyyvgug/CUBE/blob/master/resource/RSCU/", file):
        url = file.replace("https://github.com/dyyvgug/CUBE/blob/master/resource/RSCU/",
                          "https://raw.githubusercontent.com/dyyvgug/CUBE/refs/heads/master/resource/RSCU/")
        response = urllib.request.urlopen(url)
        f = io.StringIO()
        print(response.read().decode('utf'), file=f)
        f.seek(0)
    else:
        f = open(file, 'r')
        
    for l in f:
        if l.strip() != '' and not l.startswith('#'):
            l = l.strip('\n').split()
            if aaseq3==True:
                if l[0] in Seq1toSeq3:
                    aa = Seq1toSeq3[l[0]]
                elif l[0] == 'STOP':
                    aa = '*'
                else:
                    continue
            else:
                if l[0] in Seq1toSeq3:
                    aa = l[0]
                elif l[0] =="STOP":
                    aa = "*"
                else:
                    continue
            if aa not in Obs:
                Obs.setdefault(aa, {l[1]:int(l[2])})
            else:
                Obs[aa][l[1]] = int(l[2])
    f.close()
    return Obs

def get_optimal_codons_from_codonw_coafile(file):
    """
    Description 
    ----------
    Obtain a list of optimal codons from the cbi.coa or fop.coa file generated by the codonw software as downstream get_CBI() and get_Fop() function input.
    
    Parameters
    ----------
    file: str
        The file include of optimal codons generated by codonw software. such as cbi.coa, and fop.coa.
    
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

def get_optimal_codons_from_ENC(file, genetic_code=11, ratio=0.1, rscu=1, delta_rscu=0.08):
    """
    Description 
    ----------
    According to the ENC value, the optimal codon is selected.
    Usually, the codons with RSCU > 1 and delta RSCU (ΔRSCU) > 0.08 were defined as the optimal codons of the gene.
    
    Parameters
    ----------
    file: str
        A fasta or fasta.gz format file path.
    
    genetic_code: int
        A genetic code id, use `pycubs.CodonTables()` for more details.
    
    ratio: float, default=0.1
        How many sequences were selected from the genes with high and low ENC values to establish low and high bias gene groups. 
        The range of values is 0-1.
    
    rscu: float, default=1
        Lower limit of RSCU in low ENC group. The range of values is 0-6.
    
    delta_rscu: float, default=0.08
        Delta RSCU = RSCU in low ENC group - RSCU in high ENC group
    
    Reference
    ----------
    [1] Zhang, Kun, et al. "Codon usage characterization and phylogenetic analysis of the mitochondrial genome in Hemerocallis citrina." BMC Genomic Data 25.1 (2024): 6.
    [2] Wu, Peng, et al. "Comprehensive analysis of codon bias in 13 Ganoderma mitochondrial genomes." Frontiers in Microbiology 14 (2023): 1170790.
    [3] Gao, Wei, et al. "Intraspecific and interspecific variations in the synonymous codon usage in mitochondrial genomes of 8 pleurotus strains." BMC genomics 25.1 (2024): 456.
    """
    
    ENC_dict = {}
    for ID, Seq in fastaIO(file):
        Obs = get_Obs(Seq, genetic_code=genetic_code)
        enc = get_ENC(Obs)
        if enc != None:
            ENC_dict.setdefault(ID, enc)

    low_enc_critical_value =  sorted(list(ENC_dict.values()))[math.floor(len(ENC_dict) * ratio)]
    high_enc_critical_value = sorted(list(ENC_dict.values()), reverse=True)[math.floor(len(ENC_dict) * ratio)]
    
    low_enc_seq_list = []
    high_enc_seq_list = []
    
    for ID, Seq in fastaIO(file):
        if ID in ENC_dict:
            if ENC_dict[ID] >= high_enc_critical_value:
                high_enc_seq_list.append(Seq)
            #if ENC_dict[ID] < low_enc_critical_value:
            if ENC_dict[ID] <= low_enc_critical_value:
                low_enc_seq_list.append(Seq)
    
    low_RSCU = get_RSCU(get_Obs(low_enc_seq_list, genetic_code=genetic_code))
    high_RSCU = get_RSCU(get_Obs(high_enc_seq_list, genetic_code=genetic_code))
    
    aa_list = []
    codon_list = []
    low_enc_RSCU_list = []
    high_enc_RSCU_list = []
    for aa in low_RSCU.RSCU_dict:
        if aa !="*" and len(low_RSCU.RSCU_dict[aa]) >1:
            for c in low_RSCU.RSCU_dict[aa]:
                aa_list.append(aa)
                codon_list.append(c)
                low_enc_RSCU_list.append(low_RSCU.RSCU_dict[aa][c])
                high_enc_RSCU_list.append(high_RSCU.RSCU_dict[aa][c])
                
    df = pd.DataFrame({
        'Amino acid': aa_list,
        'Codon': codon_list,
        'RSCU of low ENC group': low_enc_RSCU_list,
        'RSCU of high ENC group': high_enc_RSCU_list
    })
    df['Delta RSCU'] = df['RSCU of low ENC group'] - df['RSCU of high ENC group']
    df["Optimal codon"] = (df["RSCU of low ENC group"] > rscu) & (df["Delta RSCU"] > delta_rscu)
    return df

def get_Fop(Obs, optimal_codons="Escherichia coli"):
    """
    Description 
    ----------
    Frequency of Optimal codons (Fop).
    This index, is the ratio of optimal codons to synonymous codons (genetic code dependent).
    Optimal codons for several species are in-built. 
    By default, the optimal codons of E. coli are assumed.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
    optimal_codons: list, pd.DateFrame, str, default="Escherichia coli"
        It can be the return value from the get_optimal_codons_from_codonw_coafile() function, 
        or return value from get_optimal_codons_from_ENC() function,
        or it can be a custom list containing the best codons,
        or it can be a preset value from the codonw software, such as "Escherichia coli",
        "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", 
        "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"
        
    Reference
    ----------
    [1] Ikemura, Toshimichi. "Correlation between the abundance of Escherichia coli transfer RNAs 
        and the occurrence of the respective codons in its protein genes."                      
        Journal of molecular biology 146.1 (1981): 1-21.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    
    if isinstance(optimal_codons, str):
        if optimal_codons in CBI_and_Fop_preset:
            optimal_codons = CBI_and_Fop_preset[optimal_codons][1]
        else:
            raise TypeError('Preset in {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}')
    elif isinstance(optimal_codons, pd.core.frame.DataFrame):
        optimal_codons = optimal_codons[optimal_codons["Optimal codon"]]["Codon"].tolist()
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
    Obs: dict
        get_Obs() function return value.
    
    optimal_codons: list, pd.DateFrame, str, default="Escherichia coli"
        It can be the return value from the get_optimal_codons_from_codonw_coafile() function, 
        or return value from get_optimal_codons_from_ENC() function,
        or it can be a custom list containing the best codons,
        or it can be a preset value from the codonw software, such as "Escherichia coli",
        "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", 
        "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"
                     
    Reference
    ----------
    [1] Bennetzen, Jeffrey L., and Benjamin D. Hall. "Codon selection in yeast."
        Journal of Biological Chemistry 257.6 (1982): 3026-3031.
    [2] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    
    if isinstance(optimal_codons, str):
        if optimal_codons in CBI_and_Fop_preset:
            optimal_codons = CBI_and_Fop_preset[optimal_codons][1]
        else:
            raise TypeError('Preset in {"Escherichia coli", "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"}')
    elif isinstance(optimal_codons, pd.core.frame.DataFrame):
        optimal_codons = optimal_codons[optimal_codons["Optimal codon"]]["Codon"].tolist()
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


def get_CAI(Obs, ref_Obs="Escherichia coli", model="codonw"):
    """
    Description
    ----------
    Codon Adaptation Index (CAI). 
    CAI is a measurement of the relative adaptiveness of the codon usage of a gene towards the codon usage of highly expressed genes. 
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value. 
        Observed number of occurrences of codon in a query gene.
    
    ref_Obs: str, dict, default="Escherichia coli"
        get_Obs() function return value.
        Observed number of occurrences of codon in a reference set of genes.
        Or is preset, such as "Escherichia coli", "Bacillus subtilis", "Saccharomyces cerevisiae"
    
    model: str, default="codonw",
        If model==emboss, Non-synonymous codons and termination codons (dependent on genetic code) are unexcluded. 
        If model==codonw, Non-synonymous codons and termination codons (dependent on genetic code) are excluded. 
    
    References
    ----------
    [1] Sharp, Paul M., and Wen-Hsiung Li. "The codon adaptation index-a measure of directional synonymous codon usage bias,
        and its potential applications." Nucleic acids research 15.3 (1987): 1281-1295.
    """
    
    if isinstance(ref_Obs, str):
        Relative_Adaptiveness = CAI_preset[ref_Obs][1]
    elif isinstance(ref_Obs, dict):
        Relative_Adaptiveness = get_Relative_Adaptiveness(ref_Obs).Relative_Adaptiveness_dict
    Relative_Adaptiveness = {c:Relative_Adaptiveness[aa][c] for aa in Relative_Adaptiveness for c in Relative_Adaptiveness[aa]}
    w_list = []
    if model == "emboss":
        for aa in Obs:
            for c in Obs[aa]:
                w = Relative_Adaptiveness[c]
                if w != 0:
                    w_list.extend([w]*Obs[aa][c])
                else:
                    w_list.extend([0.01]*Obs[aa][c])
                    
    elif model == "codonw":
        for aa in Obs:
            if aa !="*" and len(Obs[aa].keys()) !=1:
                for c in Obs[aa]:
                    w = Relative_Adaptiveness[c]
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
    The more the Gravy value is biased towards a negative value, the stronger the hydrophilicity of the protein.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value. 
    
    Reference
    ----------
    [1] Kyte, Jack, and Russell F. Doolittle. "A simple method for displaying the hydropathic character of a protein."
        Journal of molecular biology 157.1 (1982): 105-132.
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
    Aromaticity score of protein. This is the frequency of aromatic amino acids (Phe, Tyr, Trp)
    in the hypothetical translated gene product.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value. 
    
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


def get_cusp_like(Obs, human_format=False, outfile=None):
    """
    Description
    ----------
    The calculated result is consistent with result of cusp software.
    
    Parameters
    ----------  
    Obs: dict
        get_Obs() function return value.
    
    human_format: bool, default=False
        If value is True, return human readable format.
    
    outfile: str, default=None
        A path of outfile.
    
    Reference
    ----------
    [1] Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    [2] Rice, Peter, Ian Longden, and Alan Bleasby. "EMBOSS: the European molecular biology open software suite." Trends in genetics 16.6 (2000): 276-277.
    """
    
    res = get_ATGC_Indices(Obs)
    GC, GC1, GC2, GC3 = res["GC"], res["GC1"], res["GC2"], res["GC3"]
    Fraction = get_Fraction(Obs)
    Frequency = get_Frequency(Obs)
    CupsResult = [{"Coding GC": GC,
                   "1st letter GC": GC1,
                   "2nd letter GC": GC2,
                   "3rd letter GC": GC3}, 
                  {"Fraction": Fraction.Fraction_dict,
                   "Frequency":Frequency.Frequency_dict,
                   "Number": Obs}
                 ]
    
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
    CupsResult_human_format = out.read()
    out.close()
    
    if outfile!=None:
        out = open(outfile, 'w')
        print(CupsResult_human_format, file=out)
        out.close()
    else:   
        if human_format:
            return CupsResult_human_format
        else:
            return CupsResult

def get_codonW_like(file, genetic_code=1, cai_ref_Obs="Escherichia coli", optimal_codons="Escherichia coli", outfile=None):
    """
    Description
    ----------
    Return codonW software calculate result.
    
    Parameters
    ----------
    file: str
        A fasta or fasta.gz format file path.
    
    genetic_code: int
        A genetic code id, use `pycubs.CodonTables()` for more details.
        
    cai_ref_Obs: str, dict, default="Escherichia coli"
        get_Obs() function return value.
        Observed number of occurrences of codon in a reference set of genes.
        Or is preset, such as "Escherichia coli", "Bacillus subtilis", "Saccharomyces cerevisiae"
        
    optimal_codons: list, pd.DateFrame, str, default="Escherichia coli"
        It can be the return value from the get_optimal_codons_from_codonw_coafile() function, 
        or return value from get_optimal_codons_from_ENC() function,
        or it can be a custom list containing the best codons,
        or it can be a preset value from the codonw software, such as "Escherichia coli",
        "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", 
        "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"
    
    outfile: str, default=None
        A path of outfile.
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    
    #GCn3 = (GC - G3s -C3s )/(L_aa *3 - L_sym)
    tmp = []
    for record in fastaIO(file):
        Obs = get_Obs(record[1], genetic_code=genetic_code)
        CAI = get_CAI(Obs, ref_Obs=cai_ref_Obs, model="codonw")
        Fop = get_Fop(Obs, optimal_codons)
        CBI = get_CBI(Obs, optimal_codons)
        res = get_ATGC_Indices(Obs)
        Nc = get_ENC(Obs)
        Gravy = get_Gravy(Obs)
        Aromo = get_Aromo(Obs)
        tmp.append(pd.DataFrame({"T3s":[round(res["T3s codonW"], 4)],
                                 "C3s":[round(res["C3s codonW"], 4)],
                                 "A3s":[round(res["A3s codonW"], 4)],
                                 "G3s":[round(res["G3s codonW"], 4)],
                                 "CAI":[round(CAI, 3)],
                                 "CBI":[round(CBI, 3)],
                                 "Fop":[round(Fop, 3)],
                                 "Nc":[None if Nc == None else round(Nc, 2)],
                                 "GC3s":[round(res["GC3s"], 3)],
                                 "GC":[round(res["GC"], 3)],
                                 "L_sym":[res["L_sym"]],
                                 "L_aa":[res["L_aa"]],
                                 "Gravy":[round(Gravy, 6)],
                                 "Aromo":[round(Aromo, 6)]}, index=[record[0]]))
    df = pd.concat(tmp)
    df.index.name = 'title'
    
    del tmp
    if outfile!=None:
        df.to_csv(outfile, sep='\t', index=True)
    else:
        return df

def get_base_phase_synonymous(Obs, base=["A", "T", "G", "C"][2], phase=[0,1,2][2]):
    """
    Description
    ----------
    Calculate A1s, A2s, A3s, T1s, T2s, T3s, G1s, G2s, G3s, C1s, C2s, C3s.   
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
        
    base: {A, T, G, C}, default="G"
        Four base types.
        
    phase: {0, 1, 2}, default=2
        Codon phase, starting at 0.
    
    Example
    ----------
    A3s = get_base_phase_synonymous(Obs, base=A, phase=2)
    T3s = get_base_phase_synonymous(Obs, base=T, phase=2)
    G3s = get_base_phase_synonymous(Obs, base=G, phase=2)
    C3s = get_base_phase_synonymous(Obs, base=C, phase=2)
    A1s = get_base_phase_synonymous(Obs, base=A, phase=0)
    
    Reference
    ----------
    [1] Peden, John F. Analysis of codon usage. Diss. University of Nottingham, 2000.
    """
    Xphases_Obs = {}
    for aa in Obs:
        for c in Obs[aa]:
            if c[phase] == base:
                Xphases_Obs.setdefault(aa, Obs[aa])
    Xphases_codons = {}
    for aa in Xphases_Obs:
            if aa in ['*']:
                continue
            if len(Xphases_Obs[aa]) < 2:
                continue
            for c in Obs[aa]:
                if c[phase] not in Xphases_codons:
                    Xphases_codons.setdefault(c[phase], Obs[aa][c])
                else:
                    Xphases_codons[c[phase]] += Obs[aa][c]
    return Xphases_codons.get(base, 0)/sum(Xphases_codons.values())    

def get_ENC(Obs):
    """
    Description
    ----------
    The effective number of codons (ENC).
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
    
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

def find_four_codon_AA(Obs, add_two_codon=False):
    """
    Description:
    ----------
    Find the four-codon amino acids or two- and four-codon amino acids.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
        
    add_two_codon: bool
        If add_two_codon=True, add two-codon amino acids to the result.
    
    Reference
    ----------
    [1] Sueoka, Noboru. "Translation-coupled violation of Parity Rule 2 in human genes is not the cause
        of heterogeneity of the DNA G+ C content of third codon position." Gene 238.1 (1999): 53-58.
    """
    
    four_codon_AA = {}
    for AA in Obs:
        if AA == "*":
            continue
        if len(Obs[AA].keys()) > 1:
            m = {}
            for Codon in Obs[AA].keys():
                if Codon[:2] not in m:
                    m.setdefault(Codon[:2], [Codon])
                else:
                    m[Codon[:2]].append(Codon)
            #print(m)
            for k in m:
                if len(m[k]) == 4:
                    if AA not in four_codon_AA:
                        four_codon_AA.setdefault(AA, m[k])
                    else:
                        four_codon_AA[AA].extend(m[k])
                if add_two_codon:
                    if len(m[k]) ==2:
                        if AA not in four_codon_AA:
                            four_codon_AA.setdefault(AA, m[k])
                        else:
                            four_codon_AA[AA].extend(m[k])
    return four_codon_AA

def get_PR2(Obs, add_two_codon=False):
    """
    Description:
    ----------
    get Parity rule 2 (PR2) analysis X-axis and Y-axis.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
        
    add_two_codon: bool
        If add_two_codon=True, add two-codon amino acids to the result.
    
    Reference
    ----------
    [1] Sueoka, Noboru. "Translation-coupled violation of Parity Rule 2 in human genes is not the cause
        of heterogeneity of the DNA G+ C content of third codon position." Gene 238.1 (1999): 53-58.
    [2] Nasrullah, Izza, et al. "Genomic analysis of codon usage shows influence of mutation pressure, 
        natural selection, and host features on Marburg virus evolution." BMC evolutionary biology 15 (2015): 1-15.
    """
    
    four_codon_AA = find_four_codon_AA(Obs, add_two_codon=add_two_codon)
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
        
    if add_two_codon:
        return {"A3/(A3+T3)|2 & 4": AT3_bias_four_codon, "G3/(G3+C3)|2 & 4": GC3_bias_four_codon}
    else:
        return {"A3/(A3+T3)|4": AT3_bias_four_codon, "G3/(G3+C3)|4": GC3_bias_four_codon}


class NPA_Analysis():
    def __init__(self, file, genetic_code, sym=True):
        """
        Description
        ----------
        Neutral plot analysis.
        
        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file include of CDS seqence.
            
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        
        sym: bool, default=True
            Only synonymous codons model,
            If the value is True, amino acids without synonymous codons and stop codons are deleted.
        """
        
        GC1 = []
        GC2 = []
        GC3 = []
        GC12 = []
        GeneName = []
        for ID, Seq in fastaIO(file):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            res = get_ATGC_Indices(Obs)
            if sym == True:
                if res["GC1s"] !=None and res["GC2s"] !=None and res["GC12s"] !=None and res["GC3s"] !=None:
                    GC1.append(res["GC1s"])
                    GC2.append(res["GC2s"])
                    GC3.append(res["GC3s"])
                    GC12.append(res["GC12s"])
            else:
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
    
    def get_df(self, outfile=None):
        if self.sym:
            df = pd.DataFrame({"Gene Name":self.GeneName, "GC3s": self.GC3, "GC12s": self.GC12, "GC1s": self.GC1, "GC2s": self.GC2})
        else:
            df = pd.DataFrame({"Gene Name":self.GeneName, "GC3": self.GC3, "GC12": self.GC12, "GC1": self.GC1, "GC2": self.GC2})
        
        if outfile != None:
            df.to_csv(outfile, index=None, sep='\t')
        else:
            return df
        
    def draw_NPA_plot(self, 
                      figsize=(6,4),
                      show_gene_names=False,
                      gene_names_size=10,
                      gene_names_color="#0A0A0A",
                      point_color="#4F845C",
                      point_size=20,
                      line_color="#C25759", 
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax=None,
                      outfile=None,
                     ):
        """
        Description
        ----------
        Draw NPA plot.

        Parameters
        ----------
        figsize: tuple, default=(6,4)
            Figure size.
        
        show_gene_names: bool, or list, default=False
            Show gene name in plot. A list include of gene name, or bool value.
       
        gene_names_size, float, default=10
            Font size of gene name. 
        
        gene_names_color: str, default="#0A0A0A"
            Font color of gene name. 
            
        point_color: str, default="#4F845C"
            Point color.
            
        point_size: float, default=20
            Point size.
        
        line_color: str, default="#C25759"
            strand line color.
               
        title: str, default=None
          Title of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        #ax.set_xlim((0, max(*ax.get_xlim(), *ax.get_ylim())))
        #ax.set_ylim((0, max(*ax.get_xlim(), *ax.get_ylim())))
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None

class ENC_Analysis():
    def __init__(self, file, genetic_code):
        """
        Description
        ----------
        Effective number of codons (ENC) analysis

        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file include of CDS seqence.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        
        ys = []
        xs = []
        GeneName = []
        for ID, Seq in fastaIO(file):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            xs.append(get_ATGC_Indices(Obs)["GC3s"])
            try:
                ys.append(get_ENC(Obs))
            except:
                pass
                #print(ID)
            GeneName.append(ID)
            
        self.GC3s = xs
        self.ENC = ys
        self.GeneName = GeneName
        
    def __str__(self):
        return str(self.get_df())
    
    def __repr__(self):
        return self.__str__()
    
    def get_df(self, outfile=None):
        df = pd.DataFrame({"Gene Name":self.GeneName, "GC3s": self.GC3s, "ENC": self.ENC})
        if outfile != None:
            df.to_csv(outfile, sep='\t', index=None)
        else:
            return df 
        
    def draw_ENC_plot(self, 
                      figsize=(6,4),
                      show_gene_names=False,
                      gene_names_size=10,
                      gene_names_color="#0A0A0A",
                      point_size = 20,
                      point_color="#4F845C",
                      line_color="#C25759", 
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax=None, 
                      outfile=None):
        """
        Description
        ----------
        Draw ENC plot.

        Parameters
        ----------
        figsize: tuple, default=(6,4)
            Figure size.
        
        show_gene_names: bool, or list, default=False
            Show gene name in plot. A list include of gene name, or bool value.
       
        gene_names_size, float, default=10
            Font size of gene name. 
        
        gene_names_color: str, default="#0A0A0A"
            Font color of gene name. 
            
        point_color: str, default="#4F845C"
            Point color.
            
        point_size: float, default=20
            Point size.
        
        line_color: str, default="#C25759"
            strand line color.
               
        title: str, default=None
          Title of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None

class PR2_Analysis():
    def __init__(self, file, genetic_code):
        """
        Description
        ----------
        Parity rule 2 (PR2) analysis.
        
        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file include of CDS seqence.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.

        Reference
        ----------
        [1] Sueoka N. Intrastrand parity rules of DNA base composition and usage biases of synonymous codons.
            J Mol Evol. 1995 Mar;40(3):318-25. doi: 10.1007/BF00163236. PMID: 7723058.
        [2] Sueoka N. Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of
            the DNA G+C content of third codon position. Gene. 1999 Sep 30;238(1):53-8. doi: 10.1016/s0378-1119(99)00320-0. PMID: 10570983.
        """
        
        ys_4 = []
        xs_4 = []
        ys_2_and_4 = []
        xs_2_and_4 = []
        GeneName = []
        for ID, Seq in fastaIO(file):
            Obs = get_Obs(seqences=Seq, genetic_code=genetic_code)
            res = get_PR2(Obs, add_two_codon=False)
            ys_4.append(res["A3/(A3+T3)|4"])
            xs_4.append(res["G3/(G3+C3)|4"])
            res1 = get_PR2(Obs, add_two_codon=True)
            ys_2_and_4.append(res1["A3/(A3+T3)|2 & 4"])
            xs_2_and_4.append(res1["G3/(G3+C3)|2 & 4"])
            GeneName.append(ID)
            
        self.GeneName = GeneName
        self.X_axis_4 = xs_4
        self.Y_axis_4 = ys_4
        self.X_axis_2_and_4 = xs_2_and_4
        self.Y_axis_2_and_4 = ys_2_and_4
        
    def __str__(self):
        return str(self.get_df())
    
    def __repr__(self):
        return self.__str__()
    
    def get_df(self, outfile=None):
        df = pd.DataFrame({"Gene Name":self.GeneName, "A3/(A3+T3)| 4": self.X_axis_4, "G3/(G3+C3)| 4": self.Y_axis_4, "A3/(A3+T3)| 2 & 4": self.X_axis_2_and_4, "G3/(G3+C3)| 2 & 4": self.Y_axis_2_and_4})
        if outfile != None:
            df.to_csv(outfile, sep='\t', index=None)
        else:
            return df
    
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
                      ax = None,
                      outfile=None,
                      add_two_codon=False
                     ):
        """
        Description
        ----------
        Draw parity rule 2 (PR2) plot.

        Parameters
        ----------
        figsize: tuple, default=(6,4)
            Figure size.
        
        show_gene_names: bool, or list, default=False
            Show gene name in plot. A list include of gene name, or bool value.
       
        gene_names_size, float, default=10
            Font size of gene name. 
        
        gene_names_color: str, default="#0A0A0A"
            Font color of gene name. 
            
        point_color: str, default="#4F845C"
            Point color.
            
        point_size: float, default=20
            Point size.
        
        line_color: str, default="#C25759"
            strand line color.
               
        title: str, default=None
          Title of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        
        add_two_codon: bool
            If add_two_codon=True, add two-codon amino acids to the result.
        """
        
        if add_two_codon:
            ys = self.X_axis_2_and_4
            xs = self.Y_axis_2_and_4
        else:
            ys = self.X_axis_4
            xs = self.Y_axis_4
        
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
        if xlabel == None and add_two_codon==False:
            xlabel="G$_3$/(G$_3$+C$_3$)|4"
        if xlabel == None and add_two_codon==True:
            xlabel="G$_3$/(G$_3$+C$_3$)|2 & 4"
        if ylabel == None and add_two_codon==False:
            ylabel="A$_3$/(A$_3$+T$_3$)|4"
        if ylabel == None and add_two_codon==True:
            ylabel="A$_3$/(A$_3$+T$_3$)|2 & 4"
            
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
                        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None

class RSCU_Single_Species_Analysis():
    def __init__(self, file, genetic_code):
        """
        Description
        ----------
        Analysis of the RSCU of each gene in a single species.
        
        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file include of CDS seqence.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        
        df_list = []
        for name, seq in fastaIO(file):
            RSCU = get_RSCU(get_Obs(seq, genetic_code=genetic_code))
            df_list.append(pd.DataFrame({name:{(aa,c):RSCU.RSCU_dict[aa][c] for aa in RSCU.RSCU_dict for c in RSCU.RSCU_dict[aa]}}))
        df = pd.concat(df_list, axis=1)
        
        tmp = {}
        for aa, c in df.index:
            if aa in tmp:
                tmp[aa].append(c)
            else:
                tmp.setdefault(aa, [c])
        single_codon_amino_acids = []
        for i in tmp:
            if len(tmp[i]) == 1:
                single_codon_amino_acids.append(i)
        self.single_codon_amino_acids = single_codon_amino_acids
        df = df.astype("float64")
        df.index.names = ("Amino Acid", "Codon")
        self.RSCU_raw_df = df
        df = df.drop(single_codon_amino_acids + ['*'], axis=0)
        df = df.replace(0, 0.0000001)
        self.RSCU_df = df
        ca = prince.CA(n_components=4)
        self.ca = ca.fit(self.RSCU_df)
        pca = prince.PCA(n_components=4)
        self.pca = pca.fit(self.RSCU_df)
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
        
    def get_PCA_df(self, outfile=None):
        """
        Description
        ----------
        Return a dataframe of principal component analysis.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        PCA_df = self.pca.column_correlations
        PCA_df.columns = ["Dim 1 ("+self.pca.eigenvalues_summary.iloc[0, 1]+")", "Dim 2 ("+self.pca.eigenvalues_summary.iloc[1, 1]+")", "Dim 3 ("+self.pca.eigenvalues_summary.iloc[2, 1]+")","Dim 4 ("+self.pca.eigenvalues_summary.iloc[3, 1]+")"]
        PCA_df.index.name = "Genes"
        
        if outfile != None:
            PCA_df.to_csv(outfile, sep='\t')
        else:
            return PCA_df
        
    def get_COA_df(self, outfile=None):
        """            
        Description
        ----------
        Return a dataframe of correspondence analysis.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        
        ca_column = self.ca.column_coordinates(self.RSCU_df)
        ca_column["Type"] = "Genes"
        ca_row = self.ca.row_coordinates(self.RSCU_df)
        ca_row["Type"] = "Codons"
        COA_df = pd.concat([ca_column, ca_row])
        COA_df.columns = ["Dim 1 ("+self.ca.eigenvalues_summary.iloc[0, 1]+")",
                          "Dim 2 ("+self.ca.eigenvalues_summary.iloc[1, 1]+")",
                          "Dim 3 ("+self.ca.eigenvalues_summary.iloc[2, 1]+")",
                          "Dim 4 ("+self.ca.eigenvalues_summary.iloc[3, 1]+")", "Type"]
        if outfile != None:
            COA_df.to_csv(outfile, sep='\t', index=None)
        else:
            return COA_df
        
    def draw_COA_plot(self, 
                      figsize=(8,8), 
                      gene_labels_color="black",
                      gene_labels_style="normal",
                      gene_labels_size = 8,
                      gene_shapes_color = None,
                      gene_shapes_size = 200,
                      gene_labels_ha = "left",
                      gene_labels_va = "bottom",
                      codon_labels_color = "black",
                      codon_labels_size = 8,
                      codon_shapes_size = 100,
                      codon_labels_ha = "left",
                      codon_labels_va = "bottom",
                      show_gene_labels = True,
                      show_codon_labels = True,
                      show_codon = True,
                      show_gene = True,
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      ax=None, 
                      outfile=None):
        """  
        Description 
        -----------
        Draw correspondence analysis of RSCU plot.
        
        Parameters
        ----------
        figsize: tuple, default=(8,8)
            Figure size.
            
        gene_labels_color: str, default="black"
            Gene labels color.
        
        gene_labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Gene labels style.
            
        gene_labels_size: float, default=8
            Gene lables size.
        
        gene_shapes_color: dict, default=None,
            Gene shaple color, can be highly customized, such as {"gene1": "red", "gene2":"blue", ... }.
        
        gene_shapes_size: float, default=200
            Gene shapes size.
        
        gene_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
        
        gene_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
        
        codon_labels_color: str, default="black"
            Codon labels color.
            
        codon_labels_size: float, default=8
            Codon labels size.
        
        codon_shapes_size: float, default=100
            Codon shapes size.
            
        codon_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
            
        codon_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
            
        show_gene_labels: bool, or list, default=False
            Show gene labels in figture. A list include of gene names, or bool value.
        
        show_gene: bool, default=True
            Show genes in figture.
        
        show_codon_labels: bool, or list, default=True
            Show codon labels in figture. A list include of codons, or bool value.

        show_codon: bool, default=True
             show codon in figture.
                 
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        row_coords = self.ca.row_coordinates(self.RSCU_df)
        col_coords = self.ca.column_coordinates(self.RSCU_df)
        
        # Plotting the results
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        # Adding labels
        if show_codon:
            if show_codon_labels !=None and show_codon_labels !=False:
                if isinstance(show_codon_labels, list):
                    for i, index in enumerate(self.RSCU_df.index):
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
                    for i, index in enumerate(self.RSCU_df.index):
                        aa = index[0]
                        codon = index[1]
                        ax.annotate(codon,
                                    (row_coords[0][i], row_coords[1][i]),
                                    color=codon_labels_color,
                                    fontsize=codon_labels_size,
                                    ha=codon_labels_ha,
                                    va=codon_labels_va)
        if show_gene:
            if show_gene_labels != None and show_gene_labels != False:
                if isinstance(show_gene_labels, list):
                    for i, gene in enumerate(self.RSCU_df.columns):
                        if gene in show_gene_labels:
                            ax.annotate(gene, 
                                        (col_coords[0][i], col_coords[1][i]),
                                        color=gene_labels_color, 
                                        fontsize=gene_labels_size,
                                        style = gene_labels_style,
                                        ha = gene_labels_ha,
                                        va = gene_labels_va)
                else:
                    for i, gene in enumerate(self.RSCU_df.columns):
                        ax.annotate(gene,
                                    (col_coords[0][i], col_coords[1][i]),
                                    color=gene_labels_color,
                                    fontsize=gene_labels_size,
                                    style = gene_labels_style,
                                    ha = gene_labels_ha,
                                    va = gene_labels_va)
                    
        # Adding points
        Pa_list = []
        if show_codon:
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
        if show_gene:
            for i in range(0, len(col_coords)):
                if isinstance(gene_shapes_color, dict):
                    Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                    c=gene_shapes_color.get(col_coords.index[i],"#E64B35FF"), 
                                    label=col_coords.index[i],
                                    marker="*", 
                                    s=gene_shapes_size)
                    Ps_list.append(Ps)
                else:
                    Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                    c='#E64B35FF', 
                                    label='Genes', 
                                    marker="*",
                                    s=gene_shapes_size)
                    if i == 0:
                        Ps_list.append(Ps)
                    
        # Adding legend
        ax.legend(handles=[*tuple(Pa_list), *tuple(Ps_list)],
                  loc='upper left', 
                  bbox_to_anchor=(1, 1), 
                  ncol=1,
                  frameon=False,
                  shadow=False,
                  title=''
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_PCA_plot(self,
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="normal", 
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
                      ax=None, 
                      outfile=None):
        """
        Description 
        -----------
        Draw Principal component analysis plot.

        Parameters
        ----------
        figsize: tuple, default=(6,6)
            Figure size.
        
        labels_color: str, default="black"
            Gene labels color.
        
        labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Gene labels style.
            
        labels_size: float, default=8
            Gene lables size.
            
        labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
        
        labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
            
        show_labels: bool, or list, default=False
            Show gene labels in figture. A list include of gene names, or bool value.
            
        shapes_color: dict, default=None
            Gene shaple color, can be highly customized, such as {"gene1": "red", "gene2":"blue", ... }.
        
        shapes_size: float, default=200
            Gene shapes size.
        
        shapes_type: dict, default=None
            Gene shaple type, can be highly customized, such as {"gene1": "*", "gene2":"<", ... }.
        
        show_legend: bool, default=True
            Show legend.
        
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
            
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                     ylabels_fontstyle='normal',
                     xlabels_fontstyle='normal',
                     ylabels_fontsize=10,
                     xlabels_fontsize=10,
                     show_ylabels=True,
                     show_xlabels=True,
                     outfile=None):
        """
        Description 
        -----------
        Draw heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of X-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
             
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_clustermap(self, 
                        figsize=None, 
                        cmap="Blues",
                        ylabels_fontstyle='normal',
                        xlabels_fontstyle='normal',
                        ylabels_fontsize=10,
                        xlabels_fontsize=10,
                        row_cluster=True,
                        col_cluster=True,
                        show_ylabels=True,
                        show_xlabels=True,
                        outfile=None):
        """
        Description 
        -----------
        Draw clustermap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
        
        row_cluster: bool, default=True
            Clustering rows.
            
        col_cluster: bool, default=True
            Clustering cols.
        
        outfile: str, default=None
            A path of outfile.
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
        
        if outfile != None:
            cmp.savefig(outfile, dpi=600)
        return None
    
    def draw_boxplot(self,
                     figsize=None,
                     fontstyle='normal',
                     fontsize=10,
                     ax=None, 
                     outfile=None):
        """
        Description 
        -----------
        Draw boxplot of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of labels.
        
        fontsize: float, default=10
            Font size of labels.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        if figsize == None:
            figsize = (4, self.RSCU_df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(self.RSCU_df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        ax.set_xlabel("RSCU")
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_pearson_heatmap(self,
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='normal',
                             labels_fontsize=10,
                             ax=None, 
                             outfile=None):
        """
        Description 
        -----------
        Draw pearson heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        labels_fontstyle: {'normal', 'italic', 'oblique'} default='normal'
            Font style of labels.
            
        labels_fontsize: float, default=10
            Font size of labels.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_RSCU_barplot(self, outfile=None):
        """
        Description 
        -----------
        Draw RSCU barplot.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        
        fig, axs = plt.subplots(len(self.RSCUs_dict),
                                figsize=(8, len(self.RSCUs_dict)*4))
        fig.subplots_adjust(hspace=0.8)
        for prefix, ax in zip(self.RSCUs_dict, axs):
            self.RSCUs_dict[prefix].draw_barplot(title=prefix, codon_space=0.25, ax=ax)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def get_tree(self, metric='euclidean', outgroup="midpoint", tree_method="nj"):
        """
        Description 
        -----------
        Return a tree object from TreeNode class of scikit-bio. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
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
    
    def get_tree_newick_string(self, metric='euclidean', outgroup="midpoint", tree_method="nj"):
        """
        Description 
        -----------
        Return a tree string of newick format. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
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
                       tree_method="nj",
                       metric='euclidean',
                       outgroup="midpoint",
                       ignore_branch_length=True,
                       innode_label_size=0,
                       ladderize=True,
                       ladderize_by="size",
                       ladderize_direction="right",
                       leaf_label_size=10,
                       linewidth=1.5,
                       width=10,
                       height=0.7, 
                       outfile=None):
        """
        Description 
        -----------
        Draw tree plot of RSCU.
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
                    
        figsize: tuple, default=None
            Figure size.
                
        ignore_branch_length: bool, default=True
            Ignore branch lengths for cladogram.
        
        innode_label_size: float, default=0
            Font size for internal node labels.
        
        ladderize: bool, default=True
            Enable ladderize tree sorting.
        
        ladderize_by: {"size", "branch_length"}, default="size"
            Sort criterion.
        
        ladderize_direction: {"left", "right"}, default="right"
            Direction for larger subtrees.
        
        leaf_label_size: float, default=10
            Font size for leaf labels.
        
        linewidth: float, default=1.5
            Branch line width.
        
        height: float, default=0.7
            Figure height per leaf node.
        
        width: float, default=10
            Figure width.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        tree = self.get_tree(metric=metric, outgroup=outgroup, tree_method=tree_method)
        plotter = TreePlotter(tree, 
                              ignore_branch_length=ignore_branch_length,
                              innode_label_size=innode_label_size,
                              ladderize=ladderize,
                              ladderize_by=ladderize_by,
                              ladderize_direction=ladderize_direction,
                              leaf_label_size=leaf_label_size,
                              linewidth=linewidth,
                              width=width,
                              height=height)
        
        #plotter.add_scale_bar(length=0.5, label="Scale", position=(0, 0))
        plotter.plot(figsize=figsize, ax=ax, outfile=outfile)
        return None
    
class RSCU_Multiple_Species_Analysis():
    def __init__(self, data, genetic_code):
        """
        Description
        ----------
        RSCU analysis of multiple species.
        
        Parameters
        ----------
        data: tuple, such as: [("species name1", "cds file1 path"), ("species name2", "cds file2 path"), ...]
            CDS files for a set of species.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        
        single_codon_amino_acids = set()
        series_list = []
        RSCUs_dict = {}
        
        for prefix, file in data:
            RSCU = get_RSCU(get_Obs(file, genetic_code=genetic_code))
            RSCUs_dict.setdefault(prefix, RSCU)
            aas = []
            codons = []
            values = []
            for aa in RSCU.RSCU_dict:
                if len(RSCU.RSCU_dict[aa])==1:
                    single_codon_amino_acids.add(aa)
                for codon in RSCU.RSCU_dict[aa]:
                    aas.append(aa)
                    codons.append(codon)
                    values.append(RSCU.RSCU_dict[aa][codon])
            index = [aas, codons]
            series = pd.Series(values, index=index)
            series.name = prefix
            series.index.names = ("Amino Acid", "Codon")
            series_list.append(series)
        self.single_codon_amino_acids = tuple(single_codon_amino_acids)
        df = pd.DataFrame(series_list).T
        df = df.sort_index()
        self.RSCU_raw_df = df
        df = df.drop(list(single_codon_amino_acids) + ["*"], axis=0)
        df = df.astype("float64")
        df = df.replace(0, 0.0000001)
        self.RSCU_df = df
        ca = prince.CA(n_components=4)
        self.ca = ca.fit(df)
        pca = prince.PCA(n_components=4)
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
        
    def get_PCA_df(self, outfile=None):
        """
        Description
        ----------
        Return a dataframe of principal component analysis.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        PCA_df = self.pca.column_correlations
        PCA_df.columns = ["Dim 1 ("+self.pca.eigenvalues_summary.iloc[0, 1]+")", "Dim 2 ("+self.pca.eigenvalues_summary.iloc[1, 1]+")", "Dim 3 ("+self.pca.eigenvalues_summary.iloc[2, 1]+")","Dim 4 ("+self.pca.eigenvalues_summary.iloc[3, 1]+")"]
        PCA_df.index.name = "Species"
        
        if outfile != None:
            PCA_df.to_csv(outfile, sep='\t')
        else:
            return PCA_df
        
    def get_COA_df(self, outfile=None):
        """            
        Description
        ----------
        Return a dataframe of correspondence analysis.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        
        ca_column = self.ca.column_coordinates(self.RSCU_df)
        ca_column["Type"] = "Species"
        ca_row = self.ca.row_coordinates(self.RSCU_df)
        ca_row["Type"] = "Codons"
        COA_df = pd.concat([ca_column, ca_row])
        COA_df.columns = ["Dim 1 ("+self.ca.eigenvalues_summary.iloc[0, 1]+")",
                          "Dim 2 ("+self.ca.eigenvalues_summary.iloc[1, 1]+")",
                          "Dim 3 ("+self.ca.eigenvalues_summary.iloc[2, 1]+")",
                          "Dim 4 ("+self.ca.eigenvalues_summary.iloc[3, 1]+")", "Type"]
        if outfile != None:
            COA_df.to_csv(outfile, sep='\t', index=None)
        else:
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
                      show_species=True,
                      show_codon=True,
                      title = None,
                      xlabel = None,
                      ylabel = None,
                      title_size = 12,
                      xlabel_size = 12,
                      ylabel_size = 12,
                      ax=None, 
                      outfile=None):
        """  
        Description 
        -----------
        
        Draw correspondence analysis of RSCU plot.
        
        Parameters
        ----------
        figsize: tuple, default=(8,8)
            Figure size.
            
        species_labels_color: str, default="black"
            Species labels color.
        
        species_labels_style: {'normal', 'italic', 'oblique'} default='italic'
            Species labels style.
            
        species_labels_size: float, default=8
            Species lables size.
        
        species_shapes_color: dict, default=None,
            Species shaple color, can be highly customized, such as {"species1": "red", "species2":"blue", ... }.
        
        species_shapes_size: float, default=200
            Species shapes size.
        
        species_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
        
        species_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
        
        codon_labels_color: str, default="black"
            Codon labels color.
            
        codon_labels_size: float, default=8
            Codon labels size.
        
        codon_shapes_size: float, default=100
            Codon shapes size.
            
        codon_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
            
        codon_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
            
        show_species_labels: bool, or list, default=False
            Show species labels in figture. A list include of species names, or bool value.
        
        show_species: bool, default=True
            Show species in figture.
        
        show_codon_labels: bool, or list, default=True
            Show codon labels in figture. A list include of codons, or bool value.

        show_codon: bool, default=True
             show codon in figture.
                 
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        row_coords = self.ca.row_coordinates(self.RSCU_df)
        col_coords = self.ca.column_coordinates(self.RSCU_df)
        
        # Plotting the results
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        # Adding labels
        if show_codon:
            if show_codon_labels !=None and show_codon_labels !=False:
                if isinstance(show_codon_labels, list):
                    for i, index in enumerate(self.RSCU_df.index):
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
                    for i, index in enumerate(self.RSCU_df.index):
                        aa = index[0]
                        codon = index[1]
                        ax.annotate(codon,
                                    (row_coords[0][i], row_coords[1][i]),
                                    color=codon_labels_color,
                                    fontsize=codon_labels_size,
                                    ha=codon_labels_ha,
                                    va=codon_labels_va)
        if show_species:
            if show_species_labels != None and show_species_labels != False:
                if isinstance(show_species_labels, list):
                    for i, spceies in enumerate(self.RSCU_df.columns):
                        if spceies in show_species_labels:
                            ax.annotate(spceies, 
                                        (col_coords[0][i], col_coords[1][i]),
                                        color=species_labels_color, 
                                        fontsize=species_labels_size,
                                        style = species_labels_style,
                                        ha = species_labels_ha,
                                        va = species_labels_va)
                else:
                    for i, spceies in enumerate(self.RSCU_df.columns):
                        ax.annotate(spceies,
                                    (col_coords[0][i], col_coords[1][i]),
                                    color=species_labels_color,
                                    fontsize=species_labels_size,
                                    style = species_labels_style,
                                    ha = species_labels_ha,
                                    va = species_labels_va)
                    
        # Adding points
        Pa_list = []
        if show_codon:
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
        if show_species:
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                      ax=None, 
                      outfile=None):
        """
        Description 
        -----------
        
        Draw Principal component analysis plot.

        Parameters
        ----------
        figsize: tuple, default=(6,6)
            Figure size.
        
        labels_color: str, default="black"
            Species labels color.
        
        labels_style: {'normal', 'italic', 'oblique'} default='italic'
            Species labels style.
            
        labels_size: float, default=8
            Species lables size.
            
        labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
        
        labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
            
        show_labels: bool, or list, default=False
            Show species labels in figture. A list include of species names, or bool value.
            
        shapes_color: dict, default=None
            Species shaple color, can be highly customized, such as {"species1": "red", "species2":"blue", ... }.
        
        shapes_size: float, default=200
            Species shapes size.
        
        shapes_type: dict, default=None
            Species shaple type, can be highly customized, such as {"species1": "*", "species2":"<", ... }.
        
        show_legend: bool, default=True
            Show legend.
        
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
            
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                     show_xlabels=True,
                     outfile=None):
        """
        Description 
        -----------
        
        Draw heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of X-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
             
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                        show_xlabels=True,
                        outfile=None):
        """
        Description 
        -----------
        
        Draw clustermap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
        
        row_cluster: bool, default=True
            Clustering rows.
            
        col_cluster: bool, default=True
            Clustering cols.
        
        outfile: str, default=None
            A path of outfile.
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
        
        if outfile != None:
            cmp.savefig(outfile, dpi=600)
        return None
    
    def draw_boxplot(self,
                     figsize=None,
                     ax=None,
                     fontstyle='italic',
                     fontsize=10, 
                     outfile=None):
        """
        Description 
        -----------
        
        Draw boxplot of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of labels.
        
        fontsize: float, default=10
            Font size of labels.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        if figsize == None:
            figsize = (4, self.RSCU_df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(self.RSCU_df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        ax.set_xlabel("RSCU")
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_pearson_heatmap(self,
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='italic',
                             labels_fontsize=10,
                             ax=None, 
                             outfile=None):
        """
        Description 
        -----------
        Draw pearson heatmap of RSCU.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        labels_fontstyle: {'normal', 'italic', 'oblique'} default='italic'
            Font style of labels.
            
        labels_fontsize: float, default=10
            Font size of labels.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_RSCU_barplot(self, outfile=None):
        """
        Description 
        -----------
        Draw RSCU barplot.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        
        fig, axs = plt.subplots(len(self.RSCUs_dict),
                                figsize=(8, len(self.RSCUs_dict)*4))
        fig.subplots_adjust(hspace=0.8)
        for prefix, ax in zip(self.RSCUs_dict, axs):
            self.RSCUs_dict[prefix].draw_barplot(title=prefix, codon_space=0.25, ax=ax)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def get_tree(self, metric='euclidean', outgroup="midpoint", tree_method="nj"):
        """
        Description 
        -----------
        Return a tree object from TreeNode class of scikit-bio. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
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
    
    def get_tree_newick_string(self, metric='euclidean', outgroup="midpoint", tree_method="nj"):
        """
        Description 
        -----------
        Return a tree string of newick format. 
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
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
                       tree_method="nj",
                       metric='euclidean',
                       outgroup="midpoint",
                       ignore_branch_length=True,
                       innode_label_size=0,
                       ladderize=True,
                       ladderize_by="size",
                       ladderize_direction="right",
                       leaf_label_size=10,
                       linewidth=1.5,
                       width=10,
                       height=0.7,
                       outfile=None):
        """
        Description 
        -----------
        Draw tree plot of RSCU by NJ method.
        
        Parameters
        ----------
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
                    
        figsize: tuple, default=None
            Figure size.
                
        ignore_branch_length: bool, default=True
            Ignore branch lengths for cladogram.
        
        innode_label_size: float, default=0
            Font size for internal node labels.
        
        ladderize: bool, default=True
            Enable ladderize tree sorting.
        
        ladderize_by: {"size", "branch_length"}, default="size"
            Sort criterion.
        
        ladderize_direction: {"left", "right"}, default="right"
            Direction for larger subtrees.
        
        leaf_label_size: float, default=10
            Font size for leaf labels.
        
        linewidth: float, default=1.5
            Branch line width.
        
        height: float, default=0.7
            Figure height per leaf node.
        
        width: float, default=10
            Figure width.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        tree = self.get_tree(metric=metric, outgroup=outgroup, tree_method=tree_method)
        plotter = TreePlotter(tree, 
                              ignore_branch_length=ignore_branch_length,
                              innode_label_size=0,
                              ladderize=ladderize,
                              ladderize_by=ladderize_by,
                              ladderize_direction=ladderize_direction,
                              leaf_label_size=leaf_label_size,
                              linewidth=linewidth,
                              width=width,
                              height=height)
        #plotter.add_scale_bar(length=0.5, label="Scale", position=(0, 0))
        plotter.plot(figsize=figsize, ax=ax, outfile=outfile)
        return None

class AA_Composition_Single_Species_Analysis():
    def __init__(self, file, genetic_code):
        """
        Description
        ----------
        Analysis of the amino acid composition of each gene in a single species.
        
        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file include of CDS seqence.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        
        Seq1toSeq3 = {v:k for k,v in Seq3toSeq1.items() if k!='*'}
        
        AA_count_list = []
        for ID, Seq in fastaIO(file):
            AA_count = {'Gly':0, 'Ala':0, 'Val':0, 'Leu':0, 'Ile':0,
                        'Glu':0, 'Gln':0, 'Asp':0, 'Asn':0, 'Met':0,
                        'Ser':0, 'Thr':0, 'Phe':0, 'Trp':0, 'Tyr':0,
                        'Arg':0, 'His':0, 'Cys':0, 'Pro':0, 'Lys':0}
            for aa in translate(Seq, genetic_code=genetic_code):
                if (aa != "*") and (aa in Seq1toSeq3):
                    AA_count[Seq1toSeq3[aa]] += 1        
            AA_count_list.append(pd.DataFrame({ID:AA_count}))
        self.AA_count_df = pd.concat(AA_count_list, axis=1)
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
        
    def get_PCA_df(self, dtype="Fraction", outfile=None):
        """
        Description
        ----------
        Return a dataframe of principal component analysis.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
        
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            pca = self.AA_Fraction_pca
        else:
            pca = self.AA_count_pca
        PCA_df = pca.column_correlations
        PCA_df.columns = [f"Dim1 ({pca.eigenvalues_summary.iloc[0, 1]})",
                          f"Dim2 ({pca.eigenvalues_summary.iloc[1, 1]})"]
        PCA_df.index.name = "Gene"
        
        if outfile != None:
            PCA_df.to_csv(outfile, sep='\t')
        else:
            return PCA_df
        
    def get_COA_df(self, dtype="Fraction", outfile=None):
        """
        Description
        ----------
        Return a dataframe of correspondence analysis.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
        
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            ca = self.AA_Fraction_ca
            ca_column = ca.column_coordinates(self.AA_Fraction_df)
            ca_row = ca.row_coordinates(self.AA_Fraction_df)
        else:
            ca = self.AA_count_ca
            ca_column = ca.column_coordinates(self.AA_count_df)
            ca_row = ca.row_coordinates(self.AA_count_df)
        
        ca_column["Type"] = "Gene"
        ca_row["Type"] = "Codon"
        COA_df = pd.concat([ca_column, ca_row])
        COA_df.columns = [f"Dim1 ({ca.eigenvalues_summary.iloc[0, 1]})",
                          f"Dim2 ({ca.eigenvalues_summary.iloc[1, 1]})","Type"]
        
        if outfile != None:
            COA_df.to_csv(outfile, sep='\t', index=None)
        else:
            return COA_df
    
    def draw_heatmap(self, 
                     dtype="Fraction",
                     figsize=None, 
                     cmap="Blues",
                     ax=None,
                     ylabels_fontstyle='normal',
                     xlabels_fontstyle='normal',
                     ylabels_fontsize=10,
                     xlabels_fontsize=10,
                     show_ylabels=True,
                     show_xlabels=True,
                     outfile=None):
        """
        Description 
        -----------
        Draw heatmap of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of X-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
             
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (9, df.shape[1]/3)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(df.T, cmap=cmap, ax=ax, 
                    yticklabels=show_ylabels,
                    xticklabels=show_xlabels)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle, rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle, rotation=0)
        ax.set_xlabel("")
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_clustermap(self, 
                        dtype="Fraction",
                        figsize=None, 
                        cmap="Blues",
                        ylabels_fontstyle='normal',
                        xlabels_fontstyle='normal',
                        ylabels_fontsize=10,
                        xlabels_fontsize=10,
                        row_cluster=True,
                        col_cluster=True,
                        show_ylabels=True,
                        show_xlabels=True, 
                        outfile=None):
        """
        Description 
        -----------
        Draw clustermap of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
        
        row_cluster: bool, default=True
            Clustering rows.
            
        col_cluster: bool, default=True
            Clustering cols.
        
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (9, df.shape[1]/3)
        
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
        if outfile != None:
            cmp.savefig(outfile, dpi=600)
        return None
    
    def draw_boxplot(self,
                     dtype="Fraction",
                     figsize=None,
                     ax=None,
                     fontstyle='normal',
                     xlabel=None,
                     fontsize=10,
                     outfile=None):
        """
        Description 
        -----------
        Draw boxplot of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of labels.
        
        fontsize: float, default=10
            Font size of labels.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if xlabel == None and dtype == "Count":
            xlabel = "Count"
        ax.set_xlabel(xlabel)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_pearson_heatmap(self,
                             dtype="Fraction",
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='normal',
                             labels_fontsize=10,
                             ax=None, 
                             outfile=None):
        """
        Description 
        -----------
        Draw pearson heatmap of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        labels_fontstyle: {'normal', 'italic', 'oblique'} default='italic'
            Font style of labels.
            
        labels_fontsize: float, default=10
            Font size of labels.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_COA_plot(self,
                      dtype="Fraction",
                      figsize=(8,8),
                      gene_labels_color="black", 
                      gene_labels_style="italic", 
                      gene_labels_size=8, 
                      gene_shapes_color=None, 
                      gene_shapes_size=200, 
                      gene_labels_ha="left", 
                      gene_labels_va="bottom",
                      amino_acid_labels_color="black", 
                      amino_acid_labels_size=8, 
                      amino_acid_shapes_size=100,
                      amino_acid_labels_ha="left",
                      amino_acid_labels_va="bottom",
                      show_gene_labels=True,
                      show_gene = True,
                      show_amino_acid_labels=True,
                      show_amino_acid = True,
                      title=None,
                      xlabel=None,
                      ylabel=None,
                      title_size=12,
                      xlabel_size=12,
                      ylabel_size=12,
                      ax=None, 
                      outfile=None):
        """
        Description
        ----------
        Draw correspondence analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=(8,8)
            Figure size.
            
        gene_labels_color: str, default="black"
            Gene labels color.
        
        gene_labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Gene labels style.
            
        gene_labels_size: float, default=8
            Gene lables size.
        
        gene_shapes_color: dict, default=None,
            Gene shaple color, can be highly customized, such as {"gene1": "red", "gene2":"blue", ... }.
        
        gene_shapes_size: float, default=200
            Gene shapes size.
        
        gene_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
        
        gene_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
        
        amino_acid_labels_color: str, default="black"
            Amino acid labels color.
            
        amino_acid_labels_size: float, default=8
            Amino acid labels size.
        
        amino_acid_shapes_size: float, default=100
            Amino acid shapes size.
            
        amino_acid_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
            
        amino_acid_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
            
        show_gene_labels: bool, or list, default=False
            Show gene labels in figture. A list include of gene names, or bool value.
        
        show_gene: bool, default=True
            Show gene in figture.
        
        show_amino_acid_labels: bool, or list, default=True
            Show amino acid labels in figture. A list include of amino acids, or bool value.

        show_amino_acid: bool, default=True
             show amino acids in figture.
                 
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if show_amino_acid:
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
        if show_gene:
            if show_gene_labels != None and show_gene_labels !=False:
                if isinstance(show_gene_labels, list):
                    for i, spceies in enumerate(df.columns):
                        if spceies in show_gene_labels:
                            ax.annotate(spceies, 
                                        (col_coords[0][i], col_coords[1][i]), 
                                        color=gene_labels_color, 
                                        fontsize=gene_labels_size, 
                                        style = gene_labels_style,
                                        ha = gene_labels_ha,
                                        va = gene_labels_va)
                else:
                    for i, spceies in enumerate(df.columns):
                        ax.annotate(spceies, 
                                    (col_coords[0][i], col_coords[1][i]), 
                                    color=gene_labels_color,
                                    fontsize=gene_labels_size,
                                    style = gene_labels_style,
                                    ha = gene_labels_ha,
                                    va = gene_labels_va)
        # Adding points
        Pa_list = []
        if show_amino_acid:
            for i in range(0,len(row_coords)):
                Pa = ax.scatter(row_coords[0][i], row_coords[1][i], 
                                c=self.AAColorShapeType[row_coords.index[i]][0], 
                                label=row_coords.index[i],
                                marker=self.AAColorShapeType[row_coords.index[i]][1],
                                s=amino_acid_shapes_size)
                Pa_list.append(Pa)
        
        
        Ps_list = []
        if show_gene:
            for i in range(0, len(col_coords)):
                if isinstance(gene_shapes_color, dict):
                    Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                    c=gene_shapes_color.get(col_coords.index[i],"#E64B35FF"), 
                                    label=col_coords.index[i],
                                    marker="*", 
                                    s=gene_shapes_size)
                    Ps_list.append(Ps)
                else:
                    Ps = ax.scatter(col_coords[0][i], col_coords[1][i], 
                                     c='#E64B35FF', 
                                     label='Genes', 
                                     marker="*", 
                                     s=gene_shapes_size)
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_PCA_plot(self,
                      dtype="Fraction",
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="italic", 
                      labels_size=8, 
                      shapes_color='#E64B35FF', 
                      shapes_type='*',
                      shapes_size=100, 
                      labels_ha="left", 
                      labels_va="bottom",
                      title=None,
                      xlabel=None,
                      ylabel=None,
                      title_size=12,
                      xlabel_size=12,
                      ylabel_size=12,
                      show_labels=True,
                      show_legend=False,
                      ax=None, 
                      outfile=None):
        """
        Description 
        ----------
        Draw principal component analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=(6,6)
            Figure size.
        
        labels_color: str, default="black"
            Gene labels color.
        
        labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Gene labels style.
            
        labels_size: float, default=8
            Gene lables size.
            
        labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of gene labels.
        
        labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of gene labels.
            
        show_labels: bool, or list, default=False
            Show gene labels in figture. A list include of gene names, or bool value.
            
        shapes_color: dict, default=None
            Gene shaple color, can be highly customized, such as {"gene1": "red", "gene2":"blue", ... }.
        
        shapes_size: float, default=200
            Gene shapes size.
        
        shapes_type: dict, default=None
            Gene shaple type, can be highly customized, such as {"gene1": "*", "gene2":"<", ... }.
        
        show_legend: bool, default=True
            Show legend.
        
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
    
    def get_tree(self, metric='euclidean', outgroup="midpoint", tree_method="nj", dtype="Fraction"):
        """
        Description 
        -----------
        Return a tree object from TreeNode class of scikit-bio. 
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
        
        dist = pdist(np.matrix(df.T), metric=metric)
        dist_matrix = DistanceMatrix(squareform(dist), ids=list(df.T.index))
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
    
    def get_tree_newick_string(self, metric='euclidean', outgroup="midpoint", tree_method="nj", dtype="Fraction"):
        """
        Description 
        -----------
        Return a tree string of newick format. 
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
        """
        tree = self.get_tree(metric, outgroup, tree_method, dtype)
        f = io.StringIO()
        tree.write(f)
        f.seek(0)
        tree_str = f.read()
        f.close()
        return tree_str[:-1]
    
    def draw_tree_plot(self,
                       figsize=(8,8), 
                       ax=None,
                       tree_method="nj",
                       metric='euclidean',
                       outgroup="midpoint",
                       ignore_branch_length=True,
                       innode_label_size=0,
                       ladderize=True,
                       ladderize_by="size",
                       ladderize_direction="right",
                       leaf_label_size=10,
                       linewidth=1.5,
                       width=10,
                       height=0.7, 
                       outfile=None, 
                       dtype="Fraction"):
        """
        Description 
        -----------
        Draw tree plot of amino acid composition.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
                    
        figsize: tuple, default=None
            Figure size.
                
        ignore_branch_length: bool, default=True
            Ignore branch lengths for cladogram.
        
        innode_label_size: float, default=0
            Font size for internal node labels.
        
        ladderize: bool, default=True
            Enable ladderize tree sorting.
        
        ladderize_by: {"size", "branch_length"}, default="size"
            Sort criterion.
        
        ladderize_direction: {"left", "right"}, default="right"
            Direction for larger subtrees.
        
        leaf_label_size: float, default=10
            Font size for leaf labels.
        
        linewidth: float, default=1.5
            Branch line width.
        
        height: float, default=0.7
            Figure height per leaf node.
        
        width: float, default=10
            Figure width.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        tree = self.get_tree(metric=metric, outgroup=outgroup, tree_method=tree_method, dtype=dtype)
        plotter = TreePlotter(tree, 
                              ignore_branch_length=ignore_branch_length,
                              innode_label_size=0,
                              ladderize=ladderize,
                              ladderize_by=ladderize_by,
                              ladderize_direction=ladderize_direction,
                              leaf_label_size=leaf_label_size,
                              linewidth=linewidth,
                              width=width,
                              height=height)
        
        #plotter.add_scale_bar(length=0.5, label="Scale", position=(0, 0))
        plotter.plot(figsize=figsize, ax=ax, outfile=outfile)
        return None

class AA_Composition_Multiple_Species_Analysis():
    def __init__(self, data, genetic_code):
        """
        Description
        ----------
        Amino acid composition analysis of multiple species.
        
        Parameters
        ----------
        data: tuple, such as: [("species name1", "cds file1 path"), ("species name2", "cds file2 path"), ...]
            CDS files for a set of species.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        
        Seq1toSeq3 = {v:k for k,v in Seq3toSeq1.items() if k!='*'}
        
        AA_count_dict = {}
        for prefix, file in data:
            AA_count = defaultdict(int)
            for ID, Seq in fastaIO(file):
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
        
    def get_PCA_df(self, dtype="Fraction", outfile=None):
        """
        Description
        ----------
        Return a dataframe of principal component analysis.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
        
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            pca = self.AA_Fraction_pca
        else:
            pca = self.AA_count_pca
        PCA_df = pca.column_correlations
        PCA_df.columns = [f"Dim1 ({pca.eigenvalues_summary.iloc[0, 1]})",
                          f"Dim2 ({pca.eigenvalues_summary.iloc[1, 1]})"]
        PCA_df.index.name = "Species"
        
        if outfile != None:
            PCA_df.to_csv(outfile, sep='\t')
        else:
            return PCA_df
        
    def get_COA_df(self, dtype="Fraction", outfile=None):
        """
        Description
        ----------
        Return a dataframe of correspondence analysis.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
        
        outfile: str, default=None
            A path of outfile.
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
        COA_df.columns = [f"Dim1 ({ca.eigenvalues_summary.iloc[0, 1]})",
                          f"Dim2 ({ca.eigenvalues_summary.iloc[1, 1]})","Type"]
        
        if outfile != None:
            COA_df.to_csv(outfile, sep='\t', index=None)
        else:
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
                     show_xlabels=True,
                     outfile=None):
        """
        Description 
        -----------
        Draw heatmap of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of X-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
             
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (9, df.shape[1]/3)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(df.T, cmap=cmap, ax=ax, 
                    yticklabels=show_ylabels,
                    xticklabels=show_xlabels)
        ax.set_yticklabels(ax.get_yticklabels(),fontsize=ylabels_fontsize, fontstyle=ylabels_fontstyle, rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(),fontsize=xlabels_fontsize, fontstyle=xlabels_fontstyle, rotation=0)
        ax.set_xlabel("")
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                        show_xlabels=True, 
                        outfile=None):
        """
        Description 
        -----------
        Draw clustermap of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
        
        row_cluster: bool, default=True
            Clustering rows.
            
        col_cluster: bool, default=True
            Clustering cols.
        
        outfile: str, default=None
            A path of outfile.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
            
        if figsize == None:
            figsize = (9, df.shape[1]/3)
        
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
        if outfile != None:
            cmp.savefig(outfile, dpi=600)
        return None
    
    def draw_boxplot(self,
                     dtype="Fraction",
                     figsize=None,
                     ax=None,
                     fontstyle='italic',
                     xlabel=None,
                     fontsize=10,
                     outfile=None):
        """
        Description 
        -----------
        Draw boxplot of Count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of labels.
        
        fontsize: float, default=10
            Font size of labels.
        
         ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if xlabel == None and dtype == "Count":
            xlabel = "Count"
        ax.set_xlabel(xlabel)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_pearson_heatmap(self,
                             dtype="Fraction",
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='italic',
                             labels_fontsize=10,
                             ax=None, 
                             outfile=None):
        """
        Description 
        -----------
        
        Draw pearson heatmap of count or Fraction.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        labels_fontstyle: {'normal', 'italic', 'oblique'} default='italic'
            Font style of labels.
            
        labels_fontsize: float, default=10
            Font size of labels.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_COA_plot(self,
                      dtype="Fraction",
                      figsize=(8,8),
                      species_labels_color="black", 
                      species_labels_style="italic", 
                      species_labels_size=8, 
                      species_shapes_color=None, 
                      species_shapes_size=200, 
                      species_labels_ha="left", 
                      species_labels_va="bottom",
                      amino_acid_labels_color="black", 
                      amino_acid_labels_size=8, 
                      amino_acid_shapes_size=100,
                      amino_acid_labels_ha="left",
                      amino_acid_labels_va="bottom",
                      show_species_labels=True,
                      show_species = True,
                      show_amino_acid_labels=True,
                      show_amino_acid = True,
                      title=None,
                      xlabel=None,
                      ylabel=None,
                      title_size=12,
                      xlabel_size=12,
                      ylabel_size=12,
                      ax=None, 
                      outfile=None):
        """
        Description
        ----------
        Draw correspondence analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=(6,4)
            Figure size.
            
        species_labels_color: str, default="black"
            Species labels color.
        
        species_labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Species labels style.
            
        species_labels_size: float, default=8
            Species lables size.
        
        species_shapes_color: dict, default=None,
            Species shaple color, can be highly customized, such as {"species1": "red", "species2":"blue", ... }.
        
        species_shapes_size: float, default=200
            Species shapes size.
        
        species_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
        
        species_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
        
        amino_acid_labels_color: str, default="black"
            Amino acid labels color.
            
        amino_acid_labels_size: float, default=8
            Amino acid labels size.
        
        amino_acid_shapes_size: float, default=100
            Amino acid shapes size.
            
        amino_acid_labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
            
        amino_acid_labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
            
        show_species_labels: bool, or list, default=False
            Show species labels in figture. A list include of species names, or bool value.
        
        show_species: bool, default=True
            Show species in figture.
        
        show_amino_acid_labels: bool, or list, default=True
            Show amino acid labels in figture. A list include of amino acids, or bool value.

        show_amino_acid: bool, default=True
             show amino acids in figture.
                 
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if show_amino_acid:
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
        if show_species:
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
        if show_amino_acid:
            for i in range(0,len(row_coords)):
                Pa = ax.scatter(row_coords[0][i], row_coords[1][i], 
                                c=self.AAColorShapeType[row_coords.index[i]][0], 
                                label=row_coords.index[i],
                                marker=self.AAColorShapeType[row_coords.index[i]][1],
                                s=amino_acid_shapes_size)
                Pa_list.append(Pa)
        
        
        Ps_list = []
        if show_species:
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    
    def draw_PCA_plot(self,
                      dtype="Fraction",
                      figsize=(6,6), 
                      labels_color="black", 
                      labels_style="italic", 
                      labels_size=8, 
                      shapes_color='#E64B35FF', 
                      shapes_type='*',
                      shapes_size=100, 
                      labels_ha="left", 
                      labels_va="bottom",
                      title=None,
                      xlabel=None,
                      ylabel=None,
                      title_size=12,
                      xlabel_size=12,
                      ylabel_size=12,
                      show_labels=True,
                      show_legend=False,
                      ax=None, 
                      outfile=None):
        """
        Description 
        ----------
        Draw principal component analysis of amino acid composition plot.

        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        figsize: tuple, default=(6,6)
            Figure size.
        
        labels_color: str, default="black"
            Species labels color.
        
        labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Species labels style.
            
        labels_size: float, default=8
            Species lables size.
            
        labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
        
        labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
            
        show_labels: bool, or list, default=False
            Show species labels in figture. A list include of species names, or bool value.
            
        shapes_color: dict, default=None
            Species shaple color, can be highly customized, such as {"species1": "red", "species2":"blue", ... }.
        
        shapes_size: float, default=200
            Species shapes size.
        
        shapes_type: dict, default=None
            Species shaple type, can be highly customized, such as {"species1": "*", "species2":"<", ... }.
        
        show_legend: bool, default=True
            Show legend.
        
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
    
    def get_tree(self, metric='euclidean', outgroup="midpoint", tree_method="nj", dtype="Fraction"):
        """
        Description 
        -----------
        Return a tree object from TreeNode class of scikit-bio. 
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
        """
        
        if dtype == "Fraction":
            df = self.AA_Fraction_df
        else:
            df = self.AA_count_df
        
        dist = pdist(np.matrix(df.T), metric=metric)
        dist_matrix = DistanceMatrix(squareform(dist), ids=list(df.T.index))
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
    
    def get_tree_newick_string(self, metric='euclidean', outgroup="midpoint", tree_method="nj", dtype="Fraction"):
        """
        Description 
        -----------
        Return a tree string of newick format. 
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
        """
        tree = self.get_tree(metric, outgroup, tree_method, dtype)
        f = io.StringIO()
        tree.write(f)
        f.seek(0)
        tree_str = f.read()
        f.close()
        return tree_str[:-1]
    
    def draw_tree_plot(self,
                       figsize=(8,8), 
                       ax=None,
                       tree_method="nj",
                       metric='euclidean',
                       outgroup="midpoint",
                       ignore_branch_length=True,
                       innode_label_size=0,
                       ladderize=True,
                       ladderize_by="size",
                       ladderize_direction="right",
                       leaf_label_size=10,
                       linewidth=1.5,
                       width=10,
                       height=0.7, 
                       outfile=None, 
                       dtype="Fraction"):
        """
        Description 
        -----------
        Draw tree plot of amino acid composition by NJ method.
        
        Parameters
        ----------
        dtype: {"Fraction", "Count"} default="Fraction"
            Choose to use amino acid composition Count or Fraction.
            
        tree_method: {nj, upgma, gme, bme} default='nj'
            Method of phylogenetic reconstruction. See also skibio.tree.
        
        metric: {euclidean, cityblock, braycurtis, canberra, chebyshev, 
                 correlation, cosine, dice, hamming, jaccard, jensenshannon, 
                 kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto,
                 russellrao, seuclidean, sokalmichener, sokalsneath,
                 sqeuclidean, yule}, default="euclidean"
            The distance metric to use. 
            Euclidean distance, high sensitivity, suitable for related species; 
            braycurtis distance reduces the impact of extreme values, 
            suitable for distant species or highly variable genes, 
            emphasizing compositional differences, such as ecological data
        
        outgroup: None, midpoint, list, default="midpoint"
            If outgroup == None, return unroot tree.
                    
        figsize: tuple, default=None
            Figure size.
                
        ignore_branch_length: bool, default=True
            Ignore branch lengths for cladogram.
        
        innode_label_size: float, default=0
            Font size for internal node labels.
        
        ladderize: bool, default=True
            Enable ladderize tree sorting.
        
        ladderize_by: {"size", "branch_length"}, default="size"
            Sort criterion.
        
        ladderize_direction: {"left", "right"}, default="right"
            Direction for larger subtrees.
        
        leaf_label_size: float, default=10
            Font size for leaf labels.
        
        linewidth: float, default=1.5
            Branch line width.
        
        height: float, default=0.7
            Figure height per leaf node.
        
        width: float, default=10
            Figure width.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        tree = self.get_tree(metric=metric, outgroup=outgroup, tree_method=tree_method, dtype=dtype)
        plotter = TreePlotter(tree, 
                              ignore_branch_length=ignore_branch_length,
                              innode_label_size=0,
                              ladderize=ladderize,
                              ladderize_by=ladderize_by,
                              ladderize_direction=ladderize_direction,
                              leaf_label_size=leaf_label_size,
                              linewidth=linewidth,
                              width=width,
                              height=height)
        
        #plotter.add_scale_bar(length=0.5, label="Scale", position=(0, 0))
        plotter.plot(figsize=figsize, ax=ax, outfile=outfile)
        return None
    
class Stop_Codon_Analysis():
    def __init__(self, data, genetic_code, incomplete_codon=True):
        """
        Description
        ----------
        Stop codon analysis of multiple species.
        
        Parameters
        ----------
        data: tuple, such as: [("species name1", "cds file1 path"), ("species name2", "cds file2 path"), ...]
            CDS files for a set of species.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        """
        StopCodon_Count_dict = {}
        GeneName2StopCodon_dict = {}
        for prefix, file in data:
            StopCodon_Count = defaultdict(int)
            GeneName2StopCodon = {}
            for ID, Seq in fastaIO(file):
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
        """
        Description
        ----------
        Return a dataframe of Stop codon analysis.
        
        Parameters
        ----------
        outfile: str, default=None
            A path of outfile.
        """
        
        return pd.melt(self.StopCodon_Count_df.reset_index(), id_vars="Stop Codon", value_name="Count", var_name="Group")
        
    def draw_barplot(self, figsize=None, ax=None, palette="Blues", legend_font_style="italic", legend_font_size=10, outfile=None):
        """
        Description 
        ----------
        Draw barplot of Stop codon count.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        palette: palette name, or list, default=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"]
            Colors to use for the different levels of the `hue` variable. Should be something that can be interpreted by :func:`color_palette`.
            Reference value: Set1, Set2, Set3, tab10, tab20, tab20b, tab20c, Dark2.
        
        legend_font_style: {'normal', 'italic', 'oblique'}, default='italic'
            Font stylt of legend.
        
        legend_font_size: 10
            Font size of legend.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        df = self.get_df()
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(x=df["Stop Codon"], y=df["Count"], hue=df["Group"], hue_order=sorted(set(df["Group"])), palette=palette, ax=ax)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, shadow=False, title='', prop={'style':legend_font_style, 'size':legend_font_size})
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                        show_xlabels=True,
                        outfile=None):
        """
        Description 
        -----------
        Draw clustermap of count or Fraction.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        ylabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
        
        xlabels_fontstyle: {'normal', 'italic', 'oblique'}, default='normal'
            Font style of Y-axis labels.
            
        ylabels_fontsize: float, default=10
            Y-axis label font size of figture.
        
        xlabels_fontsize: float, default=10
            X-axis label font size of figture.
        
        show_ylabels: bool, default=True
            Show Y-axis label of figtrue.
            
        show_xlabels: bool, default=True
             Show Y-axis label of figtrue.
        
        row_cluster: bool, default=True
            Clustering rows.
            
        col_cluster: bool, default=True
            Clustering cols.
        
        outfile: str, default=None
            A path of outfile.
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
        if outfile != None:
            cmp.savefig(outfile, dpi=600)
        return None
    
    def draw_boxplot(self,
                     figsize=None,
                     ax=None,
                     fontstyle='italic',
                     fontsize=10, 
                     outfile=None):
        """
        Description 
        -----------
        Draw boxplot of stop codon count.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        fontstyle: {'normal', 'italic', 'oblique'}, default='italic'
            Font style of labels.
        
        fontsize: float, default=10
            Font size of labels.
        
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        sns.boxplot(self.StopCodon_Count_df, orient='h', ax=ax)
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=fontstyle, fontsize=fontsize)
        ax.set_xlabel("Count")
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
                      ax=None, 
                      outfile=None):
        """
        Description 
        ----------
        Draw principal component analysis of stop codon count plot.

        Parameters
        ----------
        figsize: tuple, default=(6,6)
            Figure size.
        
        labels_color: str, default="black"
            Species labels color.
        
        labels_style: {'normal', 'italic', 'oblique'} default='normal'
            Species labels style.
            
        labels_size: float, default=8
            Species lables size.
            
        labels_ha: {"left", "right", "top", "bottom", "center"}, default="left"
            Horizontalalignment of species labels.
        
        labels_va: {"left", "right", "top", "bottom", "center"}, default="bottom"
            Verticalalignment of species labels.
            
        show_labels: bool, or list, default=False
            Show species labels in figture. A list include of species names, or bool value.
            
        shapes_color: dict, default=None
            Species shaple color, can be highly customized, such as {"species1": "red", "species2":"blue", ... }.
        
        shapes_size: float, default=200
            Species shapes size.
        
        shapes_type: dict, default=None
            Species shaple type, can be highly customized, such as {"species1": "*", "species2":"<", ... }.
        
        show_legend: bool, default=True
            Show legend.
        
        title: str, default=None
            Title of figture.
        
        title_size: float, default=12
            Title font size of figture.
        
        xlabel: str, default=None
            X-axis label of figtrue.
            
        xlabel_size: float, default=12
            X-axis label font size of figture.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ylabel_size: float, default=12
            Y-axis label font size of figture.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
            
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
        if outfile != None:
            fig.savefig(outfile, dpi=600)
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
        Description 
        ----------
        Phylogenetic Tree Plotter using scikit-bio's TreeNode.
        
        Parameters
        ----------
        tree: TreeNode
            Input tree (must be rooted).

        height: float, default=0.7
            Figure height per leaf node.
        
        width: float, default=10
            Figure width.
                
        ignore_branch_length: bool, default=True
            Ignore branch lengths for cladogram.
        
        innode_label_size: float, default=0
            Font size for internal node labels.
        
        leaf_label_size: float, default=10
            Font size for leaf labels.
            
        ladderize: bool, default=True
            Enable ladderize tree sorting.
        
        ladderize_by: {"size", "branch_length"}, default="size"
            Sort criterion.
        
        ladderize_direction: {"left", "right"}, default="right"
            Direction for larger subtrees.
        
        linewidth: float, default=1.5
            Branch line width.
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
        """
        Description 
        ----------
        Recursive ladderization of tree structure
        
        Parameters
        ----------
        by: str
            Sorting criterion: "size" (subtree tip count) or 
            "branch_length" (subtree total length)
        direction: str
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
        """        
        Description 
        ----------
        Calculate node coordinates (x, y) via depth-first traversal
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
    
    def plot(self, figsize=(8,8), ax= None, outfile=None):
        """
        Description 
        ----------
        Plot rectangular layout phylogenetic tree
        
        Parameters
        ----------
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        
        Returns
        ----------
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
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def _set_axis_limits(self, ax):
        """Set dynamic axis limits with padding"""
        all_x = [pos[0] for pos in self.node_positions.values()]
        all_y = [pos[1] for pos in self.node_positions.values()]
        
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        
        x_padding = (x_max - x_min) * 0.3
        y_padding = (y_max - y_min) * 0.1
        
        #ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_xlim(x_min - x_padding, x_max + (x_max - x_min) * 0.3)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
    
    def add_node_label(
        self, 
        node_name, 
        label, 
        size= 8,
        **kwargs
    ):
        """
        Description 
        ----------
        Add custom text label to a specific node
        
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
        """Description 
        ----------
        Add evolutionary distance scale bar
        
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
        
    def highlight_clade(
        self, 
        node_names, 
        color, 
        alpha= 0.3,
        **kwargs
    ):
        """
        Description 
        ----------
        Highlight a clade with background rectangle
        
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
        
        all_x = [pos[0] for pos in self.node_positions.values()]
        x_padding = (max(all_x) - min(all_x)) * 1
        width =  max([self.node_positions[k][0] for k in self.node_positions])
        
        # Create background rectangle        
        rect = Rectangle(
            xy=(-0.5, min(y_vals)-0.5),
            #xy=(min(x_vals), min(y_vals)),
            width=width + x_padding,
            #width=max(x_vals) - min(x_vals),
            height=max(y_vals) - min(y_vals) + 1,
            #height=max(y_vals) - min(y_vals),
            color=color,
            alpha=alpha,
            zorder=-10,
            **kwargs
        )
        self._plot_patches.append(rect)

class Sequence_Indices_Analysis():
    def __init__(self, file, genetic_code, optimal_codons, cai_ref_Obs):
        """
        Description
        ----------
        Seqence analysis of various indicators, such as GC3s, ENC, CAI, CBI, Fop, GC, etc.
        
        Parameters
        ----------
        file: str
            A fasta or fasta.gz format file path.
        
        genetic_code: int
            A genetic code id, use `pycubs.CodonTables()` for more details.
        
        cai_ref_Obs: str, dict, default="Escherichia coli"
            get_Obs() function return value.
            Observed number of occurrences of codon in a reference set of genes.
            Or is preset, such as "Escherichia coli", "Bacillus subtilis", "Saccharomyces cerevisiae"
        
        optimal_codons: list, pd.DateFrame, str, default="Escherichia coli"
            It can be the return value from the get_optimal_codons_from_codonw_coafile() function, 
            or return value from get_optimal_codons_from_ENC() function,
            or it can be a custom list containing the best codons,
            or it can be a preset value from the codonw software, such as "Escherichia coli",
            "Bacillus subtilis", "Dictyostelium discoideum", "Aspergillus nidulans", 
            "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans", "Neurospora crassa"
        """
        tmp = []
        index = ['A', 'T', 'G', 'C', 'GC', 'AT', 'GC-skew', 'AT-skew',
                 'A1', 'T1', 'G1', 'C1', 'GC1', 'AT1', 'GC1-skew', 'AT1-skew',
                 'A2', 'T2', 'G2', 'C2', 'GC2', 'AT2', 'GC2-skew', 'AT2-skew',
                 'A3', 'T3', 'G3', 'C3', 'GC3', 'AT3', 'GC3-skew', 'AT3-skew', 'GC12',
                 'A1s', 'T1s', 'G1s', 'C1s', 'GC1s', 'AT1s',
                 'A2s', 'T2s', 'G2s', 'C2s', 'GC2s', 'AT2s',
                 'A3s', 'T3s', 'G3s', 'C3s', 'GC3s', 'AT3s', 'GC12s',
                 'A3s codonW', 'T3s codonW', 'G3s codonW', 'C3s codonW', 'L_sym', 'L_aa',
                 'Aromaticity', 'Hydropathicity', 'ENC', 'CAI', 'CBI', 'Fop']
        for name, seq in fastaIO(file):
            Obs = get_Obs(seq, genetic_code=genetic_code)
            Indices = get_ATGC_Indices(Obs)
            Indices["Aromaticity"] = get_Aromo(Obs)
            Indices["Hydropathicity"] = get_Gravy(Obs)
            Indices["ENC"] = get_ENC(Obs)
            Indices["CAI"] = get_CAI(Obs, ref_Obs=cai_ref_Obs)
            Indices["CBI"] = get_CBI(Obs, optimal_codons=optimal_codons)
            Indices["Fop"] = get_Fop(Obs, optimal_codons=optimal_codons)
            tmp.append(pd.DataFrame({name:Indices}, index=index))
        self.indices_df = pd.concat(tmp, axis=1).T
        
    def draw_pearson_heatmap(self,
                             figsize=None, 
                             cmap="Blues",
                             labels_fontstyle='normal',
                             labels_fontsize=10,
                             ax=None, 
                             indices='all',
                             outfile=None):
        """
        Description 
        -----------
        Draw pearson heatmap.
        
        Parameters
        ----------
        figsize: tuple, default=None
            Figure size.
            
        cmap: matplotlib colormap name or object, or list of colors, default="Blues"
            The mapping from data values to color space, such as 'Greys', 'Purples',
            'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 
            'Spectral', 'coolwarm', 'bwr', 'seismic', 'berlin', 'managua', 'vanimo'
            
        labels_fontstyle: {'normal', 'italic', 'oblique'} default='italic'
            Font style of labels.
            
        labels_fontsize: float, default=10
            Font size of labels.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
            
        indices: str, list, default="all"
            A indices list, such as ["GC1s", "GC2s", "GC3s", "CAI", "CBI", "Fop", "L_aa", 'Aromaticity', 'Hydropathicity'].
        """
        
        if indices == 'all':
            df = self.indices_df
        else:
            df = self.indices_df.loc[:, indices]
        if figsize == None:
            figsize = (df.shape[1]/4+1, df.shape[1]/4)
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
            
        corr = df.corr(method='pearson')
        mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
        sns.heatmap(corr, mask=mask, fmt=".2f", cmap="Blues", ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        ax.set_yticklabels(ax.get_yticklabels(), fontstyle=labels_fontstyle, fontsize=labels_fontsize)
        
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_regression_plot(self, 
                      x_axis,
                      y_axis,
                      figsize=(6,4),
                      show_gene_names=False,
                      gene_names_size=10,
                      gene_names_color="#0A0A0A",
                      point_color="#4F845C",
                      line_color="#C25759", 
                      point_size=20,
                      title=None,
                      xlabel=None,
                      ylabel=None, 
                      ax=None,
                      outfile=None,
                     ):
        """
        Description
        ----------
        Plot data and a linear regression model fit.
        
        Parameters
        ----------
        x_axis: str 
            A sequence indices.
        
        y_axis: str
            A sequence indices.
        
        figsize: tuple, default=None
            Figure size.
        
        show_gene_names: bool, or list, default=False
            Show gene name in plot. A list include of gene name, or bool value.
            
        gene_names_size: float, default=10
            Font size of gene name.
            
        gene_names_color: str, default="#0A0A0A"
            Font color of gene name. 
            
        point_color: str, default="#4F845C",
            Point color.
            
        line_color: str, default="#C25759",
            Strand line color.
            
        point_size: float, default=20
            Point size.

        title: str, default=None
            Title of figture.
            
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        xs = self.indices_df[x_axis].tolist()
        ys = self.indices_df[y_axis].tolist()
        labels = self.indices_df.index.tolist()
        
        PearsonRResult = ss.pearsonr(xs, ys)
        R = PearsonRResult[0]
        P = PearsonRResult[1]
        slope, intercept = np.polyfit(xs, ys, 1)
        
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
        if xlabel==None:
            xlabel = x_axis
        if ylabel==None:
            ylabel = y_axis
        if title == None:
            title = '$y$ = {:.4f}$x$ + {:.4f}; $R$$^2$ = {:.4f}; $P$ = {:.4f}'.format(slope, intercept, pow(R, 2), P)
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
    
    def draw_histogram_plot(self, 
                      x_axis,
                      figsize=(6,4),
                      color=None,
                      title=None,
                      xlabel=None,
                      ylabel="Count", 
                      ax=None,
                      outfile=None,
                     ):
        """
        Description
        ----------
        Draw histogram plot.
        
        Parameters
        ----------
        x_axis: str 
            A sequence indices.
        
        figsize: tuple, default=None
            Figure size.
            
        color: str, default=None
            Bar color.
        
        title: str, default=None
            Title of figture.
            
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        xs = self.indices_df[x_axis].tolist()
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        if xlabel == None:
            xlabel = x_axis
        sns.histplot(xs, ax=ax, color=color)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
        
    def draw_density_plot(self, 
                      x_axis,
                      figsize=(6,4),
                      color=None,
                      title=None,
                      xlabel=None,
                      ylabel="Density", 
                      ax=None,
                      outfile=None,
                     ):
        """
        Description
        ----------
        Draw density plot.
        
        Parameters
        ----------
        x_axis: str 
            A sequence indices.
        
        figsize: tuple, default=None
            Figure size.
            
        color: str, default=None
            Line color.
        
        title: str, default=None
            Title of figture.
            
        xlabel: str, default=None
            X-axis label of figtrue.
            
        ylabel: str, default=None
            Y-axis label of figtrue.
            
        ax: matplotlib Axes, default=None
            Axes object to draw the plot onto, otherwise uses the current Axes.
            
        outfile: str, default=None
            A path of outfile.
        """
        
        xs = self.indices_df[x_axis].tolist()
        if ax == None:
            fig, ax = plt.subplots(figsize=figsize)
        if xlabel == None:
            xlabel = x_axis
        sns.kdeplot(xs, ax=ax, color=color)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if outfile != None:
            fig.savefig(outfile, dpi=600)
        return None
        
def draw_codon_optimization_plot(sequence, ref_Obs, genetic_code, width=30, outfile=None):
    """
    Description
    ----------
    Draw codon optimization plot.
    
    The optimal expression of genes is achieved through the systematic design of the gene itself,
    vector, host system, and culture conditions. Researchers routinely focus on selecting appropriate
    expression vectors and host systems, but often ignore whether the gene itself achieves the best
    match with the vector and host system.
    
    The most basic principle of codon optimization is to replace codons in foreign mRNA sequences 
    with synonymous codons that are frequently used in host cells to ensure that the codons 
    in foreign mRNA sequences are more compatible with the codon usage bias of host cells 
    and avoid rare codons. The current mainstream codon optimization websites basically use 
    optimization CAI as an indicator, for example: https://www.codonbias.cn/heterologous_expression.

    However, we need to know that codons are not the only factors affecting protein expression.
    There are other factors, such as rare codons, GC content, secondary structure (free energy), etc.
    When all the codons in the foreign mRNA sequence are replaced with the best codons of the host cell,
    the protein may not be expressed, because the expression of some proteins requires the presence 
    of rare codons, which delays the progress of the ribosome and provides enough 
    time for the correct folding of the protein.

    However, for a specific amino acid sequence, due to the presence of synonymous codons,
    there may be uncountable candidate mRNA sequences, so we do not provide an optimal solution for codon optimization.
    Instead, we help researchers review the results by visualizing them and choose the 
    combination of the optimal codon and the suboptimal codon.
    
    Parameters
    ----------
    sequence: str 
        A sequence string.
    
    ref_Obs: dict
        get_Obs(), get_Obs_from_emboss_cutfile(), and get_Obs_from_CUBE_file() function return value,
        Observed number of occurrences of codon in a reference set of genes.
    
    genetic_code: int
        A genetic code id, use `pycubs.CodonTables()` for more details.
    
    width: int, default=30
        Amino acid length per row.
    
    outfile: str, default=None
        A path of outfile.
    """
    if len(list(ref_Obs.keys())[0]) == 3:
        tmp = {}
        for aa in ref_Obs:
            tmp.setdefault(Seq3toSeq1[aa], {})
            for c in ref_Obs[aa]:
                tmp[Seq3toSeq1[aa]][c] = ref_Obs[aa][c]
        ref_Obs = tmp
        
    RA = get_Relative_Adaptiveness(ref_Obs)
    start = 0
    sequences = []
    for i in range(int((len(sequence) - (len(sequence) % (width*3)))/(width*3) + 1)):
        if start+width*3 > len(sequence):
            sequences.append(sequence[start: len(sequence)])
        else:
            sequences.append(sequence[start: start + width*3])
        start += width*3
        
    fig, axs = plt.subplots(len(sequences), 1, figsize=(width/2, 2.5*len(sequences)))
    plt.subplots_adjust(hspace=0.3)
    
    max_CAI_sequence = ""
    for index, seq in enumerate(sequences):
        Obs = get_Obs(seq, genetic_code=genetic_code)
        seq_list = []
        for i in range(0, len(seq), 3):
            seq_list.append(seq.upper()[i:i+3])
        x_labels = list(translate(seq, genetic_code=genetic_code))
        
        if len(sequences) == 1:
            axs = [axs]
            
        axs[index].plot(range(0, len(x_labels)), [-1]*len(x_labels))
        axs[index].set_xticks(range(0, len(x_labels)))
        axs[index].set_xticklabels(x_labels)
        axs[index].set_xlim(-1, width)
        axs[index].set_ylim(0, 7)
        axs[index].set_yticks([])

        #pos_dict = defaultdict(dict)
        path = []
        for x, aa in enumerate(x_labels):
            for y, c in enumerate(sorted(ref_Obs[aa].keys(), key=lambda x: ref_Obs[aa][x], reverse=True)):
                #pos_dict[(x, aa)] = (y, c)
                axs[index].text(x, y, s=c, ha='center', va='bottom')
                if y==0:
                    max_CAI_sequence += c
                if seq_list[x] == c:
                    path.append((x, y))
                if seq_list[x] not in ref_Obs[aa]:
                    axs[index].text(x, len(ref_Obs[aa]), s=seq_list[x], ha='center', va='bottom', 
                                    bbox=dict(facecolor='white', edgecolor='black'))
                    path.append((x, len(ref_Obs[aa])))
                    
        paths = []
        for i in range(len(path)):
            if i < len(path)-1:
                paths.append((path[i], path[i+1]))

        for p1, p2 in paths:
             axs[index].annotate('', xy=(p2[0], p2[1]+0.25), xytext=(p1[0], p1[1]+0.25), rotation=0, 
                         arrowprops=dict(color='#00A087FF', arrowstyle='-', connectionstyle='arc', alpha=0.5,
                                         shrinkA=0, shrinkB=0, linewidth=3))

    axs[-1].annotate('', xy=(math.floor(len(seq)/3) + 3, 1), xytext=(len(seq)/3 + 3, 5.5), rotation=0, 
                     arrowprops=dict(shrink=0, edgecolor='#FB8072', facecolor="#FB8072", linewidth=0.5))
    axs[index].text(x=math.floor(len(seq)/3) + 3 + 0.1, y=4.5, s='Low',rotation=90)
    axs[index].text(x=math.floor(len(seq)/3) + 3 + 0.1, y=2.5, s='High',rotation=90)
    axs[index].text(x=math.floor(len(seq)/3) + 3 - 0.5, y=2.5, s='CAI value',rotation=90)

    axs[index].annotate('Raw CDS', xy=(path[-1][0], path[-1][1]+0.25), xytext=(path[-1][0] + 1, path[-1][1]+0.25), rotation=0, 
                arrowprops=dict(color='#00A087FF', arrowstyle='-', connectionstyle='arc', alpha=0.5,
                                shrinkA=0, shrinkB=0, linewidth=3))
    
    axs[index].text(0, -2.5, s="Note:", ha='center', va='bottom')
    axs[index].text(1.5, -2.5, s="Codon", ha='center', va='bottom', bbox=dict(facecolor='white', edgecolor='black'))
    axs[index].text(8.2, -2.5, s="inconsistent with the genetic code table of the target species.", ha='center', va='bottom')
    
    if outfile != None:
            fig.savefig(outfile, dpi=600)
    genetic_code=infer_genetic_code_from_Obs(ref_Obs)
    CAI = get_CAI(Obs=get_Obs(max_CAI_sequence, genetic_code=genetic_code), ref_Obs=ref_Obs)
    print(f">max_CAI_sequence CAI={CAI} Length={len(max_CAI_sequence)}\n{max_CAI_sequence}")
    return None
    
def infer_genetic_code_from_Obs(Obs):
    """
    Description 
    ----------
    Return a genetic code ID from Obs object.
    
    Parameters
    ----------
    Obs: dict
        get_Obs() function return value.
        
    """
    query_table = {c:aa for aa in Obs for c in Obs[aa]}
    Seq1toSeq3 = {v:k for k,v in Seq3toSeq1.items()}
    for i in [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33]:
        table = CodonTables().get(i)
        _status = True
        for codon in table:
            if len(query_table[codon]) == 3:
                if table[codon] != query_table[codon]:
                    _status = False
                    break
            else:
                if table[codon] != Seq1toSeq3[query_table[codon]]:
                    _status = False
                    break
        if _status:
            return i
