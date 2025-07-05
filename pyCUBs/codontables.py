#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2024/05/17"
__all__ = ['CodonTables']
__version__ = "v0.01"

CodonTable1={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "*", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable2={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Met", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "*", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "*", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable3={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Thr", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Thr", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Thr", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Thr", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Met", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"} 

CodonTable4={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable5={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Met", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Ser", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Ser", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable6={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "Gln", 'TGA': "*", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "Gln", 'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable9={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys", 
             'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys", 
             'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp", 
             'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp", 
             'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg", 
             'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg", 
             'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg", 
             'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg", 
             'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser", 
             'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser", 
             'ATA': "Ile", 'ACA': "Thr", 'AAA': "Asn", 'AGA': "Ser", 
             'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Ser", 
             'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly", 
             'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly", 
             'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly", 
             'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable10={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Cys",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable11={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable12={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Ser", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable13={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Met", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Gly",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Gly",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable14={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Tyr", 'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Asn", 'AGA': "Ser",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Ser",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable16={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Leu", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable21={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Met", 'ACA': "Thr", 'AAA': "Asn", 'AGA': "Ser",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Ser",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable22={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "*",   'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Leu", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable23={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "*",   'TCA': "Ser", 'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable24={'TTT': "Phe", 'TCT': "Ser",  'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser",  'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser",  'TAA': "*",   'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser",  'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro",  'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro",  'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro",  'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro",  'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr",  'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr",  'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr",  'AAA': "Lys", 'AGA': "Ser",
              'ATG': "Met", 'ACG': "Thr",  'AAG': "Lys", 'AGG': "Lys",
              'GTT': "Val", 'GCT': "Ala",  'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala",  'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala",  'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala",  'GAG': "Glu", 'GGG': "Gly"}

CodonTable25={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "Gly",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable26={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "*",   'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Ala", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable27={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Gln", 'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Gln", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable28={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Gln", 'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Gln", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable29={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Tyr", 'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Tyr", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable30={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Glu", 'TGA': "*",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Glu", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable31={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Glu", 'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "Glu", 'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Arg",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Arg",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

CodonTable33={'TTT': "Phe", 'TCT': "Ser", 'TAT': "Tyr", 'TGT': "Cys",
              'TTC': "Phe", 'TCC': "Ser", 'TAC': "Tyr", 'TGC': "Cys",
              'TTA': "Leu", 'TCA': "Ser", 'TAA': "Tyr", 'TGA': "Trp",
              'TTG': "Leu", 'TCG': "Ser", 'TAG': "*",   'TGG': "Trp",
              'CTT': "Leu", 'CCT': "Pro", 'CAT': "His", 'CGT': "Arg",
              'CTC': "Leu", 'CCC': "Pro", 'CAC': "His", 'CGC': "Arg",
              'CTA': "Leu", 'CCA': "Pro", 'CAA': "Gln", 'CGA': "Arg",
              'CTG': "Leu", 'CCG': "Pro", 'CAG': "Gln", 'CGG': "Arg",
              'ATT': "Ile", 'ACT': "Thr", 'AAT': "Asn", 'AGT': "Ser",
              'ATC': "Ile", 'ACC': "Thr", 'AAC': "Asn", 'AGC': "Ser",
              'ATA': "Ile", 'ACA': "Thr", 'AAA': "Lys", 'AGA': "Ser",
              'ATG': "Met", 'ACG': "Thr", 'AAG': "Lys", 'AGG': "Lys",
              'GTT': "Val", 'GCT': "Ala", 'GAT': "Asp", 'GGT': "Gly",
              'GTC': "Val", 'GCC': "Ala", 'GAC': "Asp", 'GGC': "Gly",
              'GTA': "Val", 'GCA': "Ala", 'GAA': "Glu", 'GGA': "Gly",
              'GTG': "Val", 'GCG': "Ala", 'GAG': "Glu", 'GGG': "Gly"}

global IDs, Names, Tables, Seq3toSeq1
IDs = [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33]
Names = ["Standard", 
         "Vertebrate Mitochondrial",
         "YeastMitochondrial",
         "Mold Mitochondrial, Protozoan Mitochondrial, Coelenterate Mitochondrial, Mycoplasma, Spiroplasma",
         "Invertebrate Mitochondrial",
         "Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear",
         "Echinoderm Mitochondrial, Flatworm Mitochondrial",
         "Euplotid Nuclear",
         "Bacterial, Archaeal, Plant Plastid",
         "Alternative Yeast Nuclear",
         "Ascidian Mitochondrial",
         "Alternative Flatworm Mitochondrial",
         "Chlorophycean Mitochondrial",
         "Trematode Mitochondrial",
         "Scenedesmus obliquus Mitochondrial",
         "Thraustochytrium Mitochondrial",
         "Rhabdopleuridae Mitochondrial",
         "Candidate Division SR1, Gracilibacteria",
         "Pachysolen tannophilus Nuclear",
         "Karyorelict Nuclear",
         "Condylostoma Nuclear",
         "Mesodinium Nuclear",
         "Peritrich Nuclear",
         "Blastocrithidia Nuclear",
         "Cephalodiscidae Mitochondrial UAA-Tyr",
]

Tables = [CodonTable1,  CodonTable2,  CodonTable3,  CodonTable4,  CodonTable5,  CodonTable6,
          CodonTable9,  CodonTable10, CodonTable11, CodonTable12, CodonTable13, CodonTable14,
          CodonTable16, CodonTable21, CodonTable22, CodonTable23, CodonTable24, CodonTable25,
          CodonTable26, CodonTable27, CodonTable28, CodonTable29, CodonTable30, CodonTable31, CodonTable33,]

Seq3toSeq1 = {"Gly": "G",
              "Ala": "A",
              "Val": "V",
              "Leu": "L",
              "Ile": "I",
              "Glu": "E",
              "Gln": "Q",
              "Asp": "D",
              "Asn": "N",
              "Met": "M",
              "Ser": "S",
              "Thr": "T",
              "Phe": "F",
              "Trp": "W",
              "Tyr": "Y",
              "Arg": "R",
              "His": "H",
              "Cys": "C",
              "Pro": "P",
              "Lys": "K",
              "*":   "*",
}

class CodonTables:
    def __init__(self):
        self.IDs = IDs
        self.Names = Names
        self.Tables = Tables
        self.Seq3toSeq1 = Seq3toSeq1
        
    def get(self, value, aaseq3=True):
        if (value in self.IDs):
            if aaseq3:
                return self.Tables[self.IDs.index(value)]
            else:
                return self._aaSeq3toSeq1(self.Tables[self.IDs.index(value)])
        elif (value in self.Names):
            if aaseq3:
                return self.Tables[self.Names.index(value)]
            else:
                return self._aaSeq3toSeq1(self.Tables[self.Names.index(value)])
        else:
            print('Plase input ID or Name of Codons.')
        
    def __str__(self):
        tmp = "Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes\n\nTranslate Tables/Genetic Codes:\n"
        for ID, Name in zip(self.IDs, self.Names):
            tmp +="{:>2s}".format(str(ID)) +": "+Name+"\n"
        return str(tmp)
    
    def __repr__(self):
        return self.__str__()
    
    def __len__(self):
        return len(self.IDs)
    
    def _aaSeq3toSeq1(self, codontable):
        tmp = {}
        for codon in codontable:
            tmp.setdefault(codon, self.Seq3toSeq1[codontable[codon]])
        return tmp
        
CBI_and_Fop_preset = {"Escherichia coli":
                      ("Ikemura (1985) Mol. Biol. Evol. 2:13-34 (updated by INCBI 1991)", 
                       ['TCT', 'TTC', 'TCC', 'TAC', 'TGC', 'CGT', 'CAC', 'CGC','CTG', 'CCG', 'CAG', 'ACT', 'ATC', 'ACC', 'AAC', 'AGC', 'AAA', 'GTT', 'GCT', 'GGT', 'GAC', 'GGC', 'GAA', 'GCG']
                      ),
                      "Bacillus subtilis":
                      ("Sharp and Devine (1989) Nucl. Acids Res 17:5029-5039)",
                       ['TCT', 'TTC', 'TAC', 'CTT', 'CCT', 'CGT', 'CGC', 'CCA','CAA', 'ACT', 'ATC', 'AAC', 'AAA', 'GTT', 'GCT', 'GGT','GAC', 'GTA', 'GAA']
                      ),
                      "Dictyostelium discoideum":
                      ("Sharp and Devine (1989) Nucl. Acids Res 17:5029-5039)", 
                       ['TTC', 'TAC', 'CGT', 'CTC', 'CAC', 'CCA', 'CAA', 'ATC', 'ACC', 'AAC', 'AAG', 'GGT', 'GTC', 'GCC', 'GAA']
                      ),
                      "Aspergillus nidulans":
                      ("Lloyd and Sharp (1991) Mol. Gen. Genet 230: 288-294", 
                       ['TTC', 'TCC', 'TAC', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CAG', 'ATC', 'ACC', 'AAC', 'AAG', 'GCT', 'GGT', 'GTC', 'GCC', 'GAC', 'GAG']
                      ),
                      "Saccharomyces cerevisiae":
                      ("Sharp and Cowe (1991) Yeast 7:657-678",
                       ['TCT', 'TGT', 'TTC', 'TCC', 'TAC', 'TTG', 'CAC', 'CCA','CAA', 'ATT', 'ACT', 'ATC', 'ACC', 'AAC', 'AGA', 'AAG','GTT', 'GCT', 'GGT', 'GTC', 'GAC', 'GAA']
                      ),
                      "Drosophila melanogaster":
                      ("Shields et al. (1988) Mol Biol Evol 5: 704-716", 
                       ['TTC', 'TCC', 'TAC', 'TGC', 'CGT', 'CCC', 'CAC', 'CGC','CTG', 'CAG', 'ATC', 'ACC', 'AAC', 'AAG', 'GTC', 'GCC','GAC', 'GGC', 'GTG', 'GAG']
                      ),
                      "Caenorhabditis elegans":
                      ("Stenico, Lloyd and Sharp Nuc. Acids Res. 22: 2437-2446(1994)",
                       ['TTC', 'TCC', 'TAC', 'TGC', 'CTT', 'CGT', 'CTC', 'CAC', 'CGC', 'CCA', 'ATC', 'ACC', 'AAC', 'AAG', 'GCT', 'GTC','GCC', 'GAC', 'GGA', 'GAG']
                      ),
                      "Neurospora crassa": 
                      ("Lloyd and Sharp (1993)",
                       ['TCT', 'TTC', 'TCC', 'TAC', 'TGC', 'CGT', 'CTC', 'CCC','CAC', 'CGC', 'CAG', 'ACT', 'ATC', 'ACC', 'AAC', 'AAG','GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GAG']
                      )
                     }

CAI_preset = {"Escherichia coli":
 ("No reference from codonw", {'Phe': {'TTT': 0.296, 'TTC': 1.0},
                               'Ser': {'TCT': 1.0, 'TCC': 0.744, 'TCA': 0.077, 'TCG': 0.017, 'AGT': 0.085, 'AGC': 0.41},
                               'Tyr': {'TAT': 0.239, 'TAC': 1.0},
                               'Cys': {'TGT': 0.5, 'TGC': 1.0},
                               'Leu': {'TTA': 0.02, 'TTG': 0.02, 'CTT': 0.042, 'CTC': 0.037, 'CTA': 0.007, 'CTG': 1.0},
                               '*': {'TAA': 0.0, 'TGA': 0.0, 'TAG': 0.0},
                               'Trp': {'TGG': 1.0},
                               'Pro': {'CCT': 0.07, 'CCC': 0.012, 'CCA': 0.135, 'CCG': 1.0},
                               'His': {'CAT': 0.291, 'CAC': 1.0},
                               'Arg': {'CGT': 1.0, 'CGC': 0.356, 'CGA': 0.004, 'CGG': 0.004, 'AGA': 0.004, 'AGG': 0.002},
                               'Gln': {'CAA': 0.124, 'CAG': 1.0},
                               'Ile': {'ATT': 0.185, 'ATC': 1.0, 'ATA': 0.003},
                               'Thr': {'ACT': 0.965, 'ACC': 1.0, 'ACA': 0.076, 'ACG': 0.099},
                               'Asn': {'AAT': 0.051, 'AAC': 1.0},
                               'Lys': {'AAA': 1.0, 'AAG': 0.253},
                               'Met': {'ATG': 1.0},
                               'Val': {'GTT': 1.0, 'GTC': 0.066, 'GTA': 0.495, 'GTG': 0.221},
                               'Ala': {'GCT': 1.0, 'GCC': 0.122, 'GCA': 0.586, 'GCG': 0.424},
                               'Asp': {'GAT': 0.434, 'GAC': 1.0},
                               'Gly': {'GGT': 1.0, 'GGC': 0.724, 'GGA': 0.01, 'GGG': 0.019},
                               'Glu': {'GAA': 1.0, 'GAG': 0.259}}),
"Bacillus subtilis":
 ('No reference from codonw', {'Phe': {'TTT': 0.571, 'TTC': 1.0},
                               'Ser': {'TCT': 1.0, 'TCC': 0.021, 'TCA': 0.458, 'TCG': 0.021, 'AGT': 0.125, 'AGC': 0.208},
                               'Tyr': {'TAT': 0.5, 'TAC': 1.0},
                               'Cys': {'TGT': 1.0, 'TGC': 1.0},
                               'Leu': {'TTA': 1.0, 'TTG': 0.036, 'CTT': 0.857, 'CTC': 0.143, 'CTA': 0.5, 'CTG': 0.071},
                               '*': {'TAA': 0.0, 'TGA': 0.0, 'TAG': 0.0},
                               'Trp': {'TGG': 1.0},
                               'Pro': {'CCT': 1.0, 'CCC': 0.071, 'CCA': 0.714, 'CCG': 0.143},
                               'His': {'CAT': 1.0, 'CAC': 0.083},
                               'Arg': {'CGT': 1.0, 'CGC': 0.609, 'CGA': 0.022, 'CGG': 0.043, 'AGA': 0.435, 'AGG': 0.022},
                               'Gln': {'CAA': 1.0, 'CAG': 0.214},
                               'Ile': {'ATT': 0.5, 'ATC': 1.0, 'ATA': 0.071},
                               'Thr': {'ACT': 1.0, 'ACC': 0.033, 'ACA': 0.867, 'ACG': 0.2},
                               'Asn': {'AAT': 0.417, 'AAC': 1.0},
                               'Lys': {'AAA': 1.0, 'AAG': 0.097},
                               'Met': {'ATG': 1.0},
                               'Val': {'GTT': 1.0, 'GTC': 0.188, 'GTA': 0.75, 'GTG': 0.438},
                               'Ala': {'GCT': 1.0, 'GCC': 0.025, 'GCA': 0.275, 'GCG': 0.125},
                               'Asp': {'GAT': 0.417, 'GAC': 1.0},
                               'Gly': {'GGT': 0.955, 'GGC': 0.773, 'GGA': 1.0, 'GGG': 0.045},
                               'Glu': {'GAA': 1.0, 'GAG': 0.412}}),
"Saccharomyces cerevisiae":
 ("Sharp and Cowe (1991) Yeast 7:657-678", {'Phe': {'TTT': 0.113, 'TTC': 1.0},
                                            'Ser': {'TCT': 1.0, 'TCC': 0.693, 'TCA': 0.036, 'TCG': 0.005, 'AGT': 0.021, 'AGC': 0.031},
                                            'Tyr': {'TAT': 0.071, 'TAC': 1.0},
                                            'Cys': {'TGT': 1.0, 'TGC': 0.077},
                                            'Leu': {'TTA': 0.117, 'TTG': 1.0, 'CTT': 0.006, 'CTC': 0.003, 'CTA': 0.039, 'CTG': 0.003},
                                            '*': {'TAA': 0.0, 'TGA': 0.0, 'TAG': 0.0},
                                            'Trp': {'TGG': 1.0},
                                            'Pro': {'CCT': 0.047, 'CCC': 0.009, 'CCA': 1.0, 'CCG': 0.002},
                                            'His': {'CAT': 0.245, 'CAC': 1.0},
                                            'Arg': {'CGT': 0.137, 'CGC': 0.002, 'CGA': 0.002, 'CGG': 0.002, 'AGA': 1.0, 'AGG': 0.003},
                                            'Gln': {'CAA': 1.0, 'CAG': 0.007},
                                            'Ile': {'ATT': 0.823, 'ATC': 1.0, 'ATA': 0.003},
                                            'Thr': {'ACT': 0.921, 'ACC': 1.0, 'ACA': 0.012, 'ACG': 0.006},
                                            'Asn': {'AAT': 0.053, 'AAC': 1.0},
                                            'Lys': {'AAA': 0.135, 'AAG': 1.0},
                                            'Met': {'ATG': 1.0},
                                            'Val': {'GTT': 1.0, 'GTC': 0.831, 'GTA': 0.002, 'GTG': 0.018},
                                            'Ala': {'GCT': 1.0, 'GCC': 0.316, 'GCA': 0.015, 'GCG': 0.001},
                                            'Asp': {'GAT': 0.554, 'GAC': 1.0},
                                            'Gly': {'GGT': 1.0, 'GGC': 0.02, 'GGA': 0.002, 'GGG': 0.004},
                                            'Glu': {'GAA': 1.0, 'GAG': 0.016}})
}

if __name__ == '__main__':
    print(CodonTables())
