#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2024/05/17"
__all__ = ["GetObs", "GetFranction", "GetFrequency", "GetRSCU", "DrawCodonBarplot", "GetCusp", "GetcodonW", 
           "NPA", "DrawNPA", "GetNC", "GetGC3s", "ENC", "DrawENC", "Find4Dtv", "GetPR2", "PR2", "DrawPR2"]
__version__ = "v0.01"

from fastaio import FastaIO
from codontables import CodonTables

def GetObs(Seqs, Genetic_Codes:int, aaSeq3:bool=True):
    """
    Calculate Observed number of occurrences of codon (Obs, 密码子使用频次)
        
    Seqs: input tupe include of list or str.
    Genetic_Codes：
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
    
    Reference website: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
    """
    
    codontable = CodonTables().get(Genetic_Codes, aaSeq3)
    
    if isinstance(Seqs, str):
        Seqs = [Seqs.upper()]

    Codons = {}
    for Seq in Seqs:
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
        Acid = codontable[Codon]
        if Acid in Obs:
            Obs[Acid][Codon] = Codons.get(Codon, 0)
        else:
            Obs.setdefault(Acid, {Codon: Codons.get(Codon, 0)})
            
    return Obs

def GetFranction(Obs):
    """
    Calculate franction of codon (Franction)
    
    Franction represents the proportion of each codon in the codon encoding the amino acid, i.e.
    
    Franction = the number of occurrences of a codon/the number of occurrences of all codons of the amino acid encoded by the codon.
    
    note: Cusp software is consistent with the calculated results
          Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    Franction = {}
    for Acid in Obs:
        Franction.setdefault(Acid, {})
        Total = sum(Obs[Acid].values())
        if Total == 0:
            continue
        for Codon in Obs[Acid]:
            Franction[Acid][Codon] = Obs[Acid][Codon]/Total
    return Franction

def GetFrequency(Obs):
    """
    Calculate frequency of codon (Frequency)
    
    Frequency indicates the Frequency of the codon occurrence in the total gene codon encoding, 
    generally expressed as the number of the codon occurrence in 1000 codons.
    
    Frequency = the number of the codon occurrence *1000/ the total number of all codons of the gene
    
    note: Cusp software is consistent with the calculated results
          Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    Total = sum([sum(Obs[Acid].values()) for Acid in Obs])
    Frequency = {}
    for Acid in Obs:
        Frequency.setdefault(Acid, {})
        for Codon in Obs[Acid]:
            Frequency[Acid][Codon] = Obs[Acid][Codon]/Total*1000
    return Frequency

def GetRSCU(Obs):
    """
    Calculate relative synonymous codon usage (RSCU, 相对同义密码子使用度)
    """
    
    RSCU = {}
    for Acid in Obs:
        RSCU.setdefault(Acid, {})
        Total = sum(Obs[Acid].values())
        if Total == 0:
            continue
        n = len(Obs[Acid])
        for Codon in Obs[Acid]:
            RSCU[Acid][Codon] = Obs[Acid][Codon]*n/Total
    return RSCU

def DrawCodonBarplot(obj, data_type='RSCU', color_preset=['Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c', 'Dark2'][1], width=0.9):
    """
    Draw a codons barplot.
    
    obj: return value of GetObs, GetFranction, GetFrequency, or GetRSCU function 
    data_type: ["Number", "Franction", "Frequency", "RSCU"]
    color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
    width: default=0.9
    
    """
    
    import matplotlib.pyplot as plt
    
    obj = obj.copy()
    del obj['*']
    fig, ax = plt.subplots(1,1, figsize=(8,4), dpi=600)
    cols = plt.colormaps.get_cmap(color_preset).colors
    cex = max([sum(obj[Acid].values()) for Acid in obj])*0.16
    for acid in obj:
        value = 0
        values = []
        colors = []
        codons = []
        for codon, color in zip(obj[acid], cols):
            value += obj[acid][codon]
            values.append(value)
            colors.append(color)
            codons.append(codon)
            
        for y, codon, value, color in zip(reversed(range(len(codons))), reversed(codons), reversed(values), reversed(colors)):
            ax.bar(acid, value, width, label=acid, fc=color)
            plt.text(x=acid, y=(-y*cex-2.2*cex)/3, ha="center", va="center", s=codon, fontdict=dict(fontsize=8, color='black', family='monospace'),
                     bbox={'facecolor': color, 'edgecolor':color, 'pad':1})
    plt.margins(x=0.01)
    plt.ylabel(data_type)
    #plt.savefig(data_type+'.png', bbox_inches='tight')
    #plt.savefig(data_type+'.pdf', bbox_inches='tight')
    return None

def GetCusp(Obs, human_format:bool=False):
    """
    return cusp calculate result.
    
    note: Cusp software is consistent with the calculated results
          Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
    """
    
    ATGC1 = {}
    ATGC2 = {}
    ATGC3 = {}
    ATGC3_2d = {}
    ATGC3_4d = {}
    ATGC3_2d_4d = {}
    for Acid in Obs:
        for Codon in Obs[Acid]:
            if Codon[0] not in ATGC1:
                ATGC1.setdefault(Codon[0], Obs[Acid][Codon])
            else:
                ATGC1[Codon[0]] += Obs[Acid][Codon]
            if Codon[1] not in ATGC2:
                ATGC2.setdefault(Codon[1], Obs[Acid][Codon])
            else:
                ATGC2[Codon[1]] += Obs[Acid][Codon]
            if Codon[2] not in ATGC3:
                ATGC3.setdefault(Codon[2], Obs[Acid][Codon])
            else:
                ATGC3[Codon[2]] += Obs[Acid][Codon]
                    
    Total_base_num = sum(ATGC1.values()) + sum(ATGC2.values())+sum(ATGC3.values())
    GC = (ATGC1.get('C', 0) + ATGC1.get('G', 0) + ATGC2.get('C', 0) + ATGC2.get('G', 0) + ATGC3.get('C', 0) + ATGC3.get('G', 0))/Total_base_num
    GC1 = (ATGC1.get('C', 0) + ATGC1.get('G', 0))/sum(ATGC1.values())
    GC2 = (ATGC2.get('C', 0) + ATGC2.get('G', 0))/sum(ATGC2.values())
    GC3 = (ATGC3.get('C', 0) + ATGC3.get('G', 0))/sum(ATGC3.values())
    
    Franction = GetFranction(Obs)
    Frequency = GetFrequency(Obs)
    
    CupsResult = [{"Coding GC": GC, "1st letter GC": GC1, "2nd letter GC": GC2, "3rd letter GC": GC3},
                  {"Franction": Franction, "Frequency":Frequency, "Number": Obs}
                 ]
    
    if human_format:
        import io
        out = io.StringIO()
        for k in CupsResult[0]:
            print('#'+k, str(round(CupsResult[0][k]*100, 2))+'%', file=out)
        print('\n#Codon AA Fraction Frequency Number', file=out)
        for Acid in CupsResult[1]['Number']:
            for Codon in CupsResult[1]['Franction'][Acid]:
                print(Codon, Acid, 
                      round(CupsResult[1]['Franction'][Acid][Codon], 3),
                      round(CupsResult[1]['Frequency'][Acid][Codon], 3),
                      CupsResult[1]['Number'][Acid][Codon], sep='\t', file=out)
        out.seek(0)
        CupsResult = out.read()
        out.close()
        return CupsResult
    else:
        return CupsResult

def GetcodonW(Obs, human_format:bool=False):
    """
    return codonW software calculate result.
    
    note: codonW  software is consistent with the calculated results
          codonW  website: https://codonw.sourceforge.net/
    """
    
    def X3s(Obs, base=["A", "T", "G", "C"]):
        X3s_Obs = {}
        for Acid in Obs:
            for Codon in Obs[Acid]:
                if Codon[2] == base:
                    X3s_Obs.setdefault(Acid, Obs[Acid])
                    continue
    
        X3s_codons = {}
        for Acid in X3s_Obs:
                if Acid in ['*']:
                    continue
                if len(X3s_Obs[Acid]) < 2:
                    continue
                for Codon in Obs[Acid]:
                    if Codon[2] not in X3s_codons:
                        X3s_codons.setdefault(Codon[2], Obs[Acid][Codon])
                    else:
                        X3s_codons[Codon[2]] += Obs[Acid][Codon]
        return X3s_codons.get(base, 0)/sum(X3s_codons.values())
    
    ATGC = {}
    for Acid in Obs:
        if Acid in ['*']:
            continue
        for Codon in Obs[Acid]:
            if Codon[0] not in ATGC:
                ATGC.setdefault(Codon[0], Obs[Acid][Codon])
            else:
                ATGC[Codon[0]] += Obs[Acid][Codon]
            if Codon[1] not in ATGC:
                ATGC.setdefault(Codon[1], Obs[Acid][Codon])
            else:
                ATGC[Codon[1]] += Obs[Acid][Codon]
            if Codon[2] not in ATGC:
                ATGC.setdefault(Codon[2], Obs[Acid][Codon])
            else:
                ATGC[Codon[2]] += Obs[Acid][Codon]
    
    GC = (ATGC.get('C', 0) + ATGC.get('G', 0))/sum(ATGC.values())
    L_aa = int(sum(ATGC.values())/3)
    
    ATGC3s = {}
    for Acid in Obs:
        if Acid in ['*']:
            continue
        if len(Obs[Acid]) < 2:
            continue
        for Codon in Obs[Acid]:
            if Codon[2] not in ATGC3s:
                ATGC3s.setdefault(Codon[2], Obs[Acid][Codon])
            else:
                ATGC3s[Codon[2]] += Obs[Acid][Codon]
                
    GC3s= ATGC3s.get('G', 0)/sum(ATGC3s.values()) + ATGC3s.get('C', 0)/sum(ATGC3s.values())
    L_sym = sum(ATGC3s.values())
    
    #GCn3 = (GC - G3s -C3s )/(L_aa *3 - L_sym)
    
    A3s = X3s(Obs, base="A")
    T3s = X3s(Obs, base="T")
    G3s = X3s(Obs, base="G")
    C3s = X3s(Obs, base="C")
    
    Nc = GetNC(Obs)
    
    codonWResult = {"A3s":A3s, "T3s":T3s, "G3s":G3s, "C3s":C3s, "GC3s": GC3s, "GC": GC, "Nc": Nc, 'L_sym': L_sym, 'L_aa':L_aa}
    if human_format:
        import io
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
              file=out,
             )
        out.seek(0)
        codonWResult = out.read()
        out.close()
    return codonWResult

def NPA(inputfile, Genetic_Codes, sym=True):
    """
    Neutral plot analysis.
    """
    
    import numpy as np
    import scipy.stats as ss

    def GC123(Obs, sym=True):
        ATGC1 = {}
        ATGC2 = {}
        ATGC3 = {}
        for Acid in Obs:
            if sym:
                if Acid in ['*']:
                    continue
                if len(Obs[Acid]) < 2:
                    continue
            for Codon in Obs[Acid]:
                if Codon[0] not in ATGC1:
                    ATGC1.setdefault(Codon[0], Obs[Acid][Codon])
                else:
                    ATGC1[Codon[0]] += Obs[Acid][Codon]
                if Codon[1] not in ATGC2:
                    ATGC2.setdefault(Codon[1], Obs[Acid][Codon])
                else:
                    ATGC2[Codon[1]] += Obs[Acid][Codon]
                if Codon[2] not in ATGC3:
                    ATGC3.setdefault(Codon[2], Obs[Acid][Codon])
                else:
                    ATGC3[Codon[2]] += Obs[Acid][Codon]

        GC1= ATGC1.get('G', 0)/sum(ATGC1.values()) + ATGC1.get('C', 0)/sum(ATGC1.values())     
        GC2= ATGC2.get('G', 0)/sum(ATGC2.values()) + ATGC2.get('C', 0)/sum(ATGC2.values())     
        GC3= ATGC3.get('G', 0)/sum(ATGC3.values()) + ATGC3.get('C', 0)/sum(ATGC3.values())
        GC12 = (GC1 + GC2)/2
        return {"GC1":GC1, "GC2":GC2, "GC12":GC12, "GC3":GC3}

    GC12 = []
    GC3 = []
    GeneName = []
    for ID, Seq in FastaIO(inputfile):
        Obs = GetObs(Seqs=Seq, Genetic_Codes=Genetic_Codes)
        res = GC123(Obs, sym)
        GC12.append(res["GC12"])
        GC3.append(res["GC3"])
        GeneName.append(ID)

    PearsonRResult = ss.pearsonr(GC3, GC12)
    R = PearsonRResult[0]
    P = PearsonRResult[1]
    slope, intercept = np.polyfit(GC3, GC12, 1)
    return {"sym":sym, "R": R, "P":P,"slope":slope,"intercept":intercept, "GC12":GC12, "GC3":GC3, "GeneName":GeneName}

def DrawNPA(NPAResult, show_label=True):
    """
    Draw NPA plot.
    """
    
    import seaborn as sns
    import matplotlib.pyplot as plt

    xs = NPAResult['GC3']
    ys = NPAResult['GC12']
    labels = NPAResult['GeneName']

    if NPAResult['sym']:
        xlabel = "P$_3$/GC$_3$$_s$"
        ylabel = "P$_1$$_2$/GC$_1$$_,$$_2$$_s$"
    else:
        xlabel = "GC$_3$"
        ylabel = "GC$_1$$_,$$_2$"

    title = "Neutral plot analysis"
    formula = '$y$ = {:.4f}$x$ + {:.4f}; $R$$^2$ = {:.4f}; $P$ = {:.4f}'.format(NPAResult['slope'], NPAResult['intercept'], pow(NPAResult["R"], 2), NPAResult["P"])

    fig, ax = plt.subplots()
    sns.regplot(x=xs, y=ys, 
                fit_reg=True, scatter_kws={"color": "r"}, 
                line_kws={"color":"blue", "linewidth": 2}, 
                label="Regression Line",
                ax=ax)
    if show_label:
        for x,y,l in zip(xs, ys, labels):
            plt.text(x,y,l)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title+"\n"+formula)
    plt.plot(xs, xs, 'b')
    return None

def GetNC(Obs):
    """
    The effective number of codons (NC) (Wright 1990). 
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
        Fxmeans.append(sum(Fx)/len(Fx) if Fx != None else None)
    
    Total = [len(F1_Obs)]
    for FxObs, Fxmean in zip(Fx_Obs, Fxmeans):
        Total.append(len(FxObs)/Fxmean if Fxmean !=None else 0)
        
    ENC = sum(Total)
    return ENC

def GetGC3s(Obs):
    """
    Calculate GC3s.
    """
    ATGC3s = {}
    for Acid in Obs:
        if Acid in ['*']:
            continue
        if len(Obs[Acid]) < 2:
            continue
        for Codon in Obs[Acid]:
            if Codon[2] not in ATGC3s:
                ATGC3s.setdefault(Codon[2], Obs[Acid][Codon])
            else:
                ATGC3s[Codon[2]] += Obs[Acid][Codon]
    GC3s= ATGC3s.get('G', 0)/sum(ATGC3s.values()) + ATGC3s.get('C', 0)/sum(ATGC3s.values())
    return GC3s

def ENC(inputfile, Genetic_Codes):
    """
    Effective number of codons (ENC) analysis
    """
    
    ys = []
    xs = []
    GeneName = []
    for ID, Seq in FastaIO(inputfile):
        Obs = GetObs(Seqs=Seq, Genetic_Codes=Genetic_Codes)
        xs.append(GetGC3s(Obs))
        ys.append(GetNC(Obs))
        GeneName.append(ID)
    return {"GC3s":xs, "ENC":ys, "GeneName":GeneName}

def DrawENC(ENCResult, show_label=False):
    """
    Draw ENC plot.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    def gc2enc(value):
        return 2 + value + 29/(pow(value, 2) + pow(1-value, 2))

    xs = ENCResult["GC3s"]
    ys = ENCResult["ENC"]
    labels = ENCResult["GeneName"]
    
    plt.scatter(x=xs, y=ys, c='r')
    plt.xlim(0,1)
    plt.ylim(0,70)
    plt.xlabel("GC$_3$$_s$")
    plt.ylabel("ENC")
    plt.title("ENC plot analysis")
    
    if show_label:
        for x,y,l in zip(xs,ys,labels):
            plt.text(x, y, l)
    
    plt.plot(np.arange(0, 1, 0.005), gc2enc(np.arange(0, 1, 0.005)), c='b')
    return None

def Find4Dtv(Obs):
    """
    Find the four-codon degenerate amino acids of at the third codon position.
    """
    
    sym4d = {}
    for Acid in Obs:
        if len(Obs[Acid].keys()) > 3:
            m = {}
            for Codon in Obs[Acid].keys():
                if Codon[:2] not in m:
                    m.setdefault(Codon[:2], [Codon])
                else:
                    m[Codon[:2]].append(Codon)
            for k in m:
                if len(m[k]) == 4:
                    sym4d.setdefault(Acid, m[k])
    return sym4d

def GetPR2(Obs):
    """
    Parity rule 2 (PR2) analysis.
    
    reference1: https://pubmed.ncbi.nlm.nih.gov/10570983/
    reference2: https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-015-0456-4 
    """
    
    sym4d = Find4Dtv(Obs)
    sym4d_codon = [i for k in sym4d for i in sym4d[k]]
    ATGC3_4d = {}
    for Acid in Obs:
        for Codon in Obs[Acid]:
            if Codon not in sym4d_codon:
                continue
            if Codon[2] not in ATGC3_4d:
                ATGC3_4d.setdefault(Codon[2], Obs[Acid][Codon])
            else:
                ATGC3_4d[Codon[2]] += Obs[Acid][Codon]
                    
    A3_4d = ATGC3_4d.get('A', 0)/sum(ATGC3_4d.values())
    T3_4d = ATGC3_4d.get('T', 0)/sum(ATGC3_4d.values())
    G3_4d = ATGC3_4d.get('G', 0)/sum(ATGC3_4d.values())
    C3_4d = ATGC3_4d.get('C', 0)/sum(ATGC3_4d.values())
    AT3_bias_4d =A3_4d/(A3_4d+T3_4d)
    GC3_bias_4d =G3_4d/(G3_4d+C3_4d)
    return {"A3/(A3+T3)|4": AT3_bias_4d, "G3/(G3+C3)|4": GC3_bias_4d}
 
def PR2(inputfile, Genetic_Codes):
    """
    Parity rule 2 (PR2) analysis.
    """
    
    ys = []
    xs = []
    GeneName = []
    for ID, Seq in FastaIO(inputfile):
        Obs = GetObs(Seqs=Seq, Genetic_Codes=Genetic_Codes)
        res = GetPR2(Obs)
        ys.append(res["A3/(A3+T3)|4"])
        xs.append(res["G3/(G3+C3)|4"])
        GeneName.append(ID)
    return {"A3/(A3+T3)|4":ys, "G3/(G3+C3)|4":xs, "GeneName":GeneName}

def DrawPR2(PR2Result, show_label=True):
    """
    Draw parity rule 2 (PR2) plot.
    """
    
    import matplotlib.pyplot as plt
    ys = PR2Result["A3/(A3+T3)|4"]
    xs = PR2Result["G3/(G3+C3)|4"]
    labels = PR2Result["GeneName"]

    plt.scatter(xs, ys, c='r')
    plt.axhline(0.5, c='b')
    plt.axvline(0.5, c='b')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("G$_3$/(G$_3$+C$_3$)|4")
    plt.ylabel("A$_3$/(A$_3$+T$_3$)|4")
    plt.title("PR2 plot analysis")

    if show_label:
        for x,y,l in zip(xs, ys, labels):
            plt.text(x,y,l)
    return None
