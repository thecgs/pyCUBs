#!/usr/bin/env python
# coding: utf-8

__author__ = "Author: Guisen Chen; Email: thecgs001@foxmail.com; Date: 2024/05/17"
__all__ = ["GetObs", "GetFranction", "GetFrequency", "GetRSCU", "DrawCodonBarplot", "GetCusp", "GetcodonW", 
           "NPA", "DrawNPA", "GetNC", "GetGC3s", "ENC", "DrawENC", "Find4Dtv", "GetPR2", "PR2", "DrawPR2"]
__version__ = "v0.01"


from fastaio import FastaIO
from codontables import CodonTables

def GetObs(seqences, genetic_codes:int, aaseq3:bool=True):
    """
    Description: 
        
        Calculate Observed number of occurrences of codon.
    
    Optional:
        
        seqences: {str, list} a seqence string, or a seqences list, or a fasta or fasta.gz format file path.
        genetic_codes: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        aaseq3: if value is True, amino acid three letter code.
    """
    
    codontable = CodonTables().get(genetic_codes, aaseq3)
    
    def _FastaIO(inputfile):
        for ID, Seq in FastaIO(inputfile):
            yield Seq
    
    if isinstance(seqences, str):
        import os
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
        Acid = codontable[Codon]
        if Acid in Obs:
            Obs[Acid][Codon] = Codons.get(Codon, 0)
        else:
            Obs.setdefault(Acid, {Codon: Codons.get(Codon, 0)})
            
    return Obs

def GetFranction(Obs):
    """
    Description: 
        
        Calculate franction of codon (Franction).
        Franction represents the proportion of each codon in the codon encoding the amino acid, i.e.
        Franction = the number of occurrences of a codon/the number of occurrences of all codons of the amino acid encoded by the codon.
    
    Optional:
        
        Obs: GetObs function return value.

    Note: 
        
        Cusp software is consistent with the calculated results
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
    Description: 
    
        Calculate frequency of codon (Frequency).
        Frequency indicates the Frequency of the codon occurrence in the total gene codon encoding, 
        generally expressed as the number of the codon occurrence in 1000 codons.
        Frequency = the number of the codon occurrence *1000/ the total number of all codons of the gene
    
    Optional:
        
        Obs: GetObs function return value.
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
    Description: 
    
        Calculate relative synonymous codon usage (RSCU)
    
    Optional:
        
        Obs: GetObs function return value.
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

def DrawCodonBarplot(obj, ylabel='RSCU', title=None, color_preset=["#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"], width=0.9, remove_stop_codon=True):
    """
    Description: 
        
        Draw a codons barplot.
        
    Optional:
        
        obj: return value of GetObs, GetFranction, GetFrequency, or GetRSCU function.
        ylabel: ylabel of plot.
        title: title of plot.
        width: bar spacing width. default=0.9
        color_preset: ["Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "Dark2"]
        remove_stop_codon: remove stop codon.
    """
    
    import matplotlib.pyplot as plt
    
    obj = obj.copy()
    if remove_stop_codon:
        del obj['*']
    fig, ax = plt.subplots(1,1, figsize=(8,4), dpi=600)
    
    if isinstance(color_preset, str):
        cols = plt.colormaps.get_cmap(color_preset).colors
    else:
        cols = color_preset
        
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
    plt.ylabel(ylabel)
    plt.title(title)
    return None

def GetCusp(Obs, human_format:bool=False):
    """
    Description: 
    
        Return cusp calculate result.
        
    Optional:
        
        Obs: GetObs function return value.
        human_format: if value is True, return human readable format.
    
    Note:
        
        Cusp software is consistent with the calculated results
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
    Description: 
    
        Return codonW software calculate result.
        
    Optional:
        
        Obs: GetObs function return value.
        human_format: if value is True, return human readable format.
    
    Note:
        
        Cusp software is consistent with the calculated results
        Cusp website: https://www.bioinformatics.nl/cgi-bin/emboss/cusp
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

def NPA(inputfile, genetic_codes, sym=True):
    """
    Description: 
    
        Neutral plot analysis.
        
    Optional:
        
        inputfile: a fasta or fasta.gz format file.
        genetic_codes: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        sym: only synonymous codons model.
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
        Obs = GetObs(seqences=Seq, genetic_codes=genetic_codes)
        res = GC123(Obs, sym)
        GC12.append(res["GC12"])
        GC3.append(res["GC3"])
        GeneName.append(ID)

    PearsonRResult = ss.pearsonr(GC3, GC12)
    R = PearsonRResult[0]
    P = PearsonRResult[1]
    slope, intercept = np.polyfit(GC3, GC12, 1)
    return {"sym":sym, "R": R, "P":P,"slope":slope,"intercept":intercept, "GC12":GC12, "GC3":GC3, "GeneName":GeneName}

def DrawNPA(NPAResult, show_gene_name=False, gene_name_size=10,
            gene_name_color="#0A0A0A", point_color="#4F845C", line_color="#C25759", 
            title="Neutral plot analysis", xlabel=None, ylabel=None, 
           ):
    """
    Description: 
        
        Draw NPA plot.
    
    Optional:
        
        NPAResult: NPA function return value.
        show_gene_name: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_name_size: font size of gene name. 
        gene_name_color: font color of gene name. 
        point_color: point color.
        line_color: strand line color.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
    
    """
    
    import seaborn as sns
    import matplotlib.pyplot as plt

    xs = NPAResult['GC3']
    ys = NPAResult['GC12']
    labels = NPAResult['GeneName']
    
    if xlabel == None:
        if NPAResult['sym']:
            xlabel = "P$_3$/GC$_3$$_s$"
        else:
            xlabel = "GC$_3$"
            
    if ylabel == None:
        if NPAResult['sym']:
            ylabel = "P$_1$$_2$/GC$_1$$_,$$_2$$_s$"
        else:
            ylabel =  "GC$_1$$_,$$_2$"
    
    formula = '$y$ = {:.4f}$x$ + {:.4f}; $R$$^2$ = {:.4f}; $P$ = {:.4f}'.format(NPAResult['slope'], NPAResult['intercept'], pow(NPAResult["R"], 2), NPAResult["P"])

    fig, ax = plt.subplots()
    sns.regplot(x=xs, y=ys, 
                fit_reg=True, scatter_kws={"color": point_color}, 
                line_kws={"color":line_color, "linewidth": 2}, 
                label="Regression Line",
                ax=ax)
    
    if show_gene_name:
        for x,y,l in zip(xs, ys, labels):
            if isinstance(show_gene_name, bool):
                plt.text(x, y, l, c=gene_name_color, size=gene_name_size)
            else:
                if l in show_gene_name:
                    plt.text(x, y, l, c=gene_name_color, size=gene_name_size)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title+"\n"+formula)
    plt.plot(xs, xs, c=line_color)
    return None

def GetNC(Obs):
    """
    Description: 
    
        The effective number of codons (NC) (Wright 1990). 
        
    Optional:
        
        inputfile: a fasta or fasta.gz format file.
        genetic_codes: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
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
    Description: 
    
        Calculate GC3s.
        
    Optional:
        
        Obs: GetObs function return value.    
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

def ENC(inputfile, genetic_codes):
    """
    Description: 
    
        Effective number of codons (ENC) analysis
        
    Optional:
        
        Obs: GetObs function return value.   
    """
    
    ys = []
    xs = []
    GeneName = []
    for ID, Seq in FastaIO(inputfile):
        Obs = GetObs(seqences=Seq, genetic_codes=genetic_codes)
        xs.append(GetGC3s(Obs))
        ys.append(GetNC(Obs))
        GeneName.append(ID)
    return {"GC3s":xs, "ENC":ys, "GeneName":GeneName}

def DrawENC(ENCResult, show_gene_name=False, gene_name_size=10,
            gene_name_color="#0A0A0A", point_color="#4F845C", line_color="#C25759", 
            title="ENC plot analysis", xlabel="GC$_3$$_s$", ylabel="ENC", 
           ):
    """
    Description: 
        
        Draw ENC plot.
    
    Optional:
        
        ENCResult: ENC function return value.
        show_gene_name: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_name_size: font size of gene name. 
        gene_name_color: font color of gene name. 
        point_color: point color.
        line_color: strand line color.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    def gc2enc(value):
        return 2 + value + 29/(pow(value, 2) + pow(1-value, 2))
    
    xs = ENCResult["GC3s"]
    ys = ENCResult["ENC"]
    labels = ENCResult["GeneName"]
    
    plt.scatter(x=xs, y=ys, c=point_color)
    plt.xlim(0,1)
    plt.ylim(0,70)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if show_gene_name:
        for x,y,l in zip(xs, ys, labels):
            if isinstance(show_gene_name, bool):
                plt.text(x, y, l, c=gene_name_color, size=gene_name_size)
            else:
                if l in show_gene_name:
                    plt.text(x, y, l, c=gene_name_color, size=gene_name_size)
            
    plt.plot(np.arange(0, 1, 0.005), gc2enc(np.arange(0, 1, 0.005)), c=line_color)
    return None

def Find4Dtv(Obs):
    """
    Description: 
    
        Find the four-codon degenerate amino acids of at the third codon position.
        
    Optional:
        
        Obs: GetObs function return value.
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
    Description:
    
        Parity rule 2 (PR2) analysis.
    
    Optional:
        
        Obs: GetObs function return value.
        
    Reference:
        
        https://pubmed.ncbi.nlm.nih.gov/10570983/
        https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-015-0456-4 
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
 
def PR2(inputfile, genetic_codes):
    """
    Description: 
        
        Parity rule 2 (PR2) analysis.
    
    Optional:
        
        inputfile: a fasta or fasta.gz format file.
        genetic_codes: genetic code id, use `import codontables; codontables.CodonTables()` for more details.
        
    Reference: 
        [1] Sueoka N. Intrastrand parity rules of DNA base composition and usage biases of synonymous codons.
            J Mol Evol. 1995 Mar;40(3):318-25. doi: 10.1007/BF00163236. PMID: 7723058.
        [2] Sueoka N. Translation-coupled violation of Parity Rule 2 in human genes is not the cause of heterogeneity of
            the DNA G+C content of third codon position. Gene. 1999 Sep 30;238(1):53-8. doi: 10.1016/s0378-1119(99)00320-0. PMID: 10570983.
    """
    
    ys = []
    xs = []
    GeneName = []
    for ID, Seq in FastaIO(inputfile):
        Obs = GetObs(seqences=Seq, genetic_codes=genetic_codes)
        res = GetPR2(Obs)
        ys.append(res["A3/(A3+T3)|4"])
        xs.append(res["G3/(G3+C3)|4"])
        GeneName.append(ID)
    return {"GeneName":GeneName, "G3/(G3+C3)|4":xs, "A3/(A3+T3)|4":ys}

def DrawPR2(PR2Result, show_gene_name=False, gene_name_size=10,
            gene_name_color="#0A0A0A", point_color="#4F845C", line_color="#C25759", 
            title="PR2 plot analysis", xlabel="G$_3$/(G$_3$+C$_3$)|4", ylabel="A$_3$/(A$_3$+T$_3$)|4", 
           ):
    """
    Description: 
        
        Draw parity rule 2 (PR2) plot.
    
    Optional:
        
        PR2Result: PR2 function return value.
        show_gene_name: {bool, ["gene_name1", "gene_name2", ...]} show gene name in plot.
        gene_name_size: font size of gene name. 
        gene_name_color: font color of gene name. 
        point_color: point color.
        line_color: strand line color.
        title: title of plot.
        xlabel: xlabel of plot.
        ylabel: ylabel of plot.
    
    """
    
    import matplotlib.pyplot as plt
    
    ys = PR2Result["A3/(A3+T3)|4"]
    xs = PR2Result["G3/(G3+C3)|4"]
    labels = PR2Result["GeneName"]
    
    plt.scatter(xs, ys, c=point_color)
    plt.axhline(0.5, c=line_color)
    plt.axvline(0.5, c=line_color)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if show_gene_name:
        for x,y,l in zip(xs, ys, labels):
            if isinstance(show_gene_name, bool):
                plt.text(x, y, l, c=gene_name_color, size=gene_name_size)
            else:
                if l in show_gene_name:
                    plt.text(x, y, l, c=gene_name_color, size=gene_name_size)
    return None
