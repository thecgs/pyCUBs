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
                     
if __name__ == '__main__':
    print(CodonTables())
