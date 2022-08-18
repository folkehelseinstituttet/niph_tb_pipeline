#!/usr/bin/env python

''' The script that finishes the TeX report. Called from TB_pipeline.py '''

import sys
import os
import time
from .__init__ import __version__ as VERSION

def CreateFooter(metainfo):
    currentdate = time.strftime("%Y-%b-%d")
    ID = metainfo["ID"].replace("_","\\_")
    footer = '\\rfoot{Isolat ID: %s | Dato: %s }' % (ID, currentdate)
    with open("Latex_template/include/footer.tex","w") as outfile:
        outfile.write(footer)

def CreateInfo(metainfo, covdicsample, samplesindb):
    # DataQual should be assessed from presence of Fastqc / kaiju problem files
    num_variants = covdicsample["VARIANT"]
    bases_lowcov = covdicsample["LOWCOV"]
    bases_het = covdicsample["HET"]
    samples_in_db = str(samplesindb)
    percentage_aligned = covdicsample["ALIGNED"] / (covdicsample["LENGTH"] - covdicsample["MASKED"]) * 100
    percentage_lowqual = covdicsample["LOWCOV"] / (covdicsample["LENGTH"] - covdicsample["MASKED"])

    if os.path.isfile("averagedepth.txt"):
        with open("averagedepth.txt", "rU") as rdfile:
            RD = rdfile.readlines()[0].rstrip("\n")
            try:
                RD = float(RD)
            except NameError:
                RD = 0.0
    else:
        RD = 0.0
    DataQualList = []
    if RD < 30.0:
        DataQualList.append("Lav dybde")
    if os.path.isfile("Mashclassificationproblem"):
        DataQualList.append("Ikke Mycobacterium")
    if os.path.isfile("Mashothermycobacterium"):
        DataQualList.append("Ikke MTB")
    if percentage_aligned < 90.00:
        DataQualList.append("Lav ref. coverage")
    if os.path.isfile("Fastqc_problems"):
        DataQualList.append("Lav kvalitet")
    if len(DataQualList) == 0:
        DataQual = "OK"
    else:
        DataQual = ", ".join(DataQualList)

    ID = metainfo["ID"].replace("_","\\_")
    Barcode = metainfo["Barcode"]
    Location = metainfo["Location"]
    LocatedFrom = metainfo["Source"]
    Isolationdate = metainfo["Isolated"]



    infostring = '''
    ID for pr\\o ve    &  {SampleName}    & Perc. aligned          & {percentage_aligned:.2f}           \\\ \hline
    Variants   & {num_variants}         & Bases low cov       & {bases_lowcov}          \\\ \hline
    Datakvalitet & {DataQual}       & Samples in DB & {samples_in_db}      \\\ \hline
    Read depth  & {Readdepth:.2f}       & Dato  & {currentdate}  \\\ \hline
    '''

    infostringfull = infostring.format(
        SampleName = ID,
        percentage_aligned = percentage_aligned,
        num_variants = num_variants,
        bases_lowcov = bases_lowcov,
        DataQual = DataQual,
        samples_in_db = samples_in_db,
        Readdepth = RD,
        currentdate = time.strftime("%Y-%b-%d"))
    with open("Latex_template/include/info.tex","w") as outfile:
        outfile.write(infostringfull)

def CreateOppsummering(resistens, clusterjanei, species, mashspecies):

    if clusterjanei:
        smitte = 'Pr\\o ven tilh\\o rer et sannsynlig smittecluster, noe som antyder \\textbf{nylig smitte}. '
    else:
        smitte = 'Pr\\o ven er ikke n\\ae rt beslektet med tidligere sekvenserte pr\\o ver. '

    if len(resistens) == 0:
        res = 'Det ble ikke funnet noen resistensmutasjoner. '
    elif len(resistens) == 1:
        res = 'Det ble funnet mutasjoner som indikerer resistens mot %s. ' % list(resistens.keys())[0]
    else:
        reskeys = list(resistens.keys())
        alleres = ", ".join(["\\textbf{%s}" % r for r in reskeys[:-1]]) + " og \\textbf{%s}." % reskeys[-1]
        res = 'Det ble funnet mutasjoner som indikerer resistens mot {alleres} '.format(alleres=alleres)

    if not (os.path.isfile("Mashothermycobacterium") or os.path.isfile("Mashclassificationproblem")):
        oppsummering = 'Pr\\o ven var positiv for \\textbf{Mycobacterium tuberculosis-komplekset (MTC)}. %s %s ' % (res, smitte)
        tophit = species
        oppsummering += 'Beste artstreff fra Mykrobe var \\textbf{%s}. Beste fra MASH var \\textbf{%s}. ' % (species, mashspecies)
    elif os.path.isfile("Mashclassificationproblem"):
        tophit = open("Mashclassificationproblem","rU").read()
        oppsummering = 'Pr\\o ven er ikke i \\textbf{MTBC}. Beste treff fra Mykrobe er \\textbf{%s} og fra MASH \\textbf{%s}. ' % (species, mashspecies)
    else:
        tophit = open("Mashothermycobacterium","rU").read()
        oppsummering = 'Pr\\o ven er ikke \\textbf{Mycobacterium tuberculosis}. Beste treff fra Mykrobe er \\textbf{%s} og fra MASH \\textbf{%s}. ' % (species, mashspecies)
    with open("Latex_template/include/oppsummering.tex","w") as outfile:
        outfile.write(oppsummering)

def CreateTyping(lineage):
    # Now gathered directly from colltyper
    #lineage = lineage.replace("_","-")
    #if lineage is None:
    lineage = lineage.replace("_","\\_")
    mix = EvaluateLineageMix(lineage)
    if lineage == "":
        typing = 'Pr\\o ven ble ikke entydig typet til en bestemt lineage. '
    elif lineage == "4.9" and (os.path.isfile("Mashothermycobacterium") or os.path.isfile("Mashclassificationproblem")):
        typing = 'Pr\\o ven hadde ingen fylogenetiske mutasjoner mot Hr37Rv, indikerer 4.9 eller annen art'
    else:
        typing = 'Pr\\o ven ble funnet\ \\aa\ tilh\\o re lineage %s. ' % (lineage)
    if mix:
        typing += 'Pr\\o ven kan v\\ae re en blanding av ulike samples. '
    with open("Latex_template/include/typing.tex","w") as outfile:
        outfile.write(typing)

def EvaluateLineageMix(lineage):
    '''Takes the lineage text and evaluates whether there is a mix of different lineages (e.g. 4.7 and 3.2)'''
    lintab = lineage.split(" / ")
    # Finding even one non-fit indicates a possible mix
    # Start with the first, iterate over the rest to check all pairs
    for i in range(len(lintab)-1):
        for j in range(i+1, len(lintab)):
            if TwoLineagesMix(lintab[i],lintab[j]):
                return True
    return False

def TwoLineagesMix(i,j):
    '''Evaluate if lineage i and j are compatible or not'''
    # Add catch for BOV and BOV\\_AFRI
    if "BOV" in i and "BOV" in j:
        return False
    if "BOV" in i and ("5" in j or "6" in j):
        return False
    if ("5" in j or "6" in j) and "BOV" in j:
        return False
    isplit = i.split('.')
    jsplit = j.split('.')
    largest = max(len(isplit),len(jsplit))
    if len(isplit) == largest:
        for elem in range(len(jsplit)):
            if not jsplit[elem] == isplit[elem]:
                return True
    else:
        for elem in range(len(isplit)):
            if not isplit[elem] == jsplit[elem]:
                return True
    return False
    
def CreateResistensBokser(resistens):
    # Figure out the appropriate box:
    # DEFAULTS:
    IngenRes = '& $\square$ Ingen resistens predikert \\\ '
    MonoRes = '& $\square$ Mono-resistens \\\ '
    Res = '& $\square$ Resistens \\\ '
    MultiRes = '& $\square$ Multi-resistens (MDR) \\\ '
    XRes = '& $\square$ Omfattende resistens (XDR) \\\ '


    if len(resistens) == 0:
        IngenRes = '& \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Ingen resistens predikert \\\ '

    #elif len(resistens) == 1 and any([drug in resistens for drug in ["Pyrazinamide","Isoniazid","Rifampicin","Ethambutol"]]):
    #    MonoRes = '& \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Mono-resistens \\\ '

    else:
        if 'Isoniazid' in resistens and 'Rifampicin' in resistens and 'Quinolones' in resistens and any([drug in resistens for drug in ["Amikacin","Capreomycin","Kanamycin"]]):
            XRes = '& \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Omfattende resistens (XDR) \\\ '
        elif 'Isoniazid' in resistens and 'Rifampicin' in resistens:
            MultiRes = '&  \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Multi-resistens (MDR) \\\ '
        elif 'Isoniazid' in resistens or 'Rifampicin' in resistens:
            MonoRes = '& \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Mono-resistens \\\ '
        else:
            Res = '& \\cellcolor[HTML]{EFEFEF} $\\text{\\rlap{$\\checkmark$}}\\square$ Resistens \\\ '

    text = '''%s
%s
%s
%s
%s''' % (IngenRes, MonoRes, Res, MultiRes, XRes)
    with open("Latex_template/include/resistensbokser.tex","w") as outfile:
        outfile.write(text)

def CreateResistens(resistensdic):
    # Sort out first line drugs and second line drugs
    # Sort resistant and non-resistant
    # Get specific mutation of resistant

    FirstLine = ["Ethambutol", "Pyrazinamide", "Streptomycin", "Isoniazid", "Rifampicin"] # PASS FORSKJELLER I STAVING NORSK/ENGELSK
    #SecondLine = ["Ofloxacin", "Levofloxacin", "Moxifloxacin", "Amikacin", "Kanamycin", "Capreomycin"] # OBS: Levofloxacin ligger ikke inne i MYKROBE PREDICTOR!! Men lik andre kinoloner
    SecondLine = ["Fluorokinoloner", "Amikacin", "Kanamycin", "Capreomycin"]
    FirstSens = [drug for drug in FirstLine if drug not in resistensdic]
    FirstRes = [drug for drug in FirstLine if drug in resistensdic]
    SecSens = [drug for drug in SecondLine if drug not in resistensdic]
    SecRes = [drug for drug in SecondLine if drug in resistensdic]
    if 'Quinolones' in resistensdic:
        #SecRes += ['Ofloxacin', 'Levofloxacin', 'Moxifloxacin']
        #SecSens.remove('Ofloxacin')
        #SecSens.remove('Levofloxacin')
        #SecSens.remove('Moxifloxacin')
        SecRes += ["Fluorokinoloner"]
        SecSens.remove("Fluorokinoloner")
    Row1SR1 = Row2SR1 = Row3SR1 = Row4SR1 = Row5SR1 = Row1SR2 = Row2SR2 = Row3SR2 = Row4SR2 = '' # = Row5SR2 = Row6SR2 = '' Removed when quinolones collapsed
    
    # Figure out the second column
    if len(FirstSens) == 5:
        Row5SR1 = '\multirow{-5}{*}{Sensitiv}'
    elif len(FirstSens) == 4:
        Row5SR1 = '\cellcolor[HTML]{EFEFEF}Resistent'
        Row4SR1 = '\multirow{-4}{*}{Sensitiv}'
    elif len(FirstSens) == 3:
        Row5SR1 = '\multirow{-2}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row4SR1 = '\cellcolor[HTML]{EFEFEF}'
        Row3SR1 = '\multirow{-3}{*}{Sensitiv}'
    elif len(FirstSens) == 2:
        Row5SR1 = '\multirow{-3}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row4SR1 = Row3SR1 = '\cellcolor[HTML]{EFEFEF}'
        Row2SR1 = '\multirow{-2}{*}{Sensitiv}'
    elif len(FirstSens) == 1:
        Row5SR1 = '\multirow{-4}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row4SR1 = Row3SR1 = Row2SR1 = '\cellcolor[HTML]{EFEFEF}'
        Row1SR1 = 'Sensitiv'
    else:
        Row5SR1 = '\multirow{-5}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row4SR1 = Row3SR1 = Row2SR1 = Row1SR1 = '\cellcolor[HTML]{EFEFEF}'


    #if len(SecSens) == 6:
    #    Row6SR2 = '\multirow{-6}{*}{Sensitiv}'
    #elif len(SecSens) == 5:
    #    Row6SR2 = '\cellcolor[HTML]{EFEFEF}Resistent'
    #    Row5SR2 = '\multirow{-5}{*}{Sensitiv}'
    if len(SecSens) == 4:
        #Row6SR2 = '\multirow{-2}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        #Row5SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row4SR2 = '\multirow{-4}{*}{Sensitiv}'
    elif len(SecSens) == 3:
        #Row6SR2 = '\multirow{-3}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row4SR2 = '\cellcolor[HTML]{EFEFEF}Resistent' # = Row5SR2
        Row3SR2 = '\multirow{-3}{*}{Sensitiv}'
    elif len(SecSens) == 2:
        #Row6SR2 = '\multirow{-4}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        #Row5SR2 = Row4SR2 = Row3SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row4SR2 = '\multirow{-2}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row3SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row2SR2 = '\multirow{-2}{*}{Sensitiv}'
    elif len(SecSens) == 1:
        #Row6SR2 = '\multirow{-5}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        #Row5SR2 = Row4SR2 = Row3SR2 = Row2SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row4SR2 = '\multirow{-3}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row3SR2 = Row2SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row1SR2 = 'Sensitiv'
    else:
        #Row6SR2 = '\multirow{-6}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        #Row5SR2 = Row4SR2 = Row3SR2 = Row2SR2 = Row1SR2 = '\cellcolor[HTML]{EFEFEF}'
        Row4SR2 = '\multirow{-4}{*}{\cellcolor[HTML]{EFEFEF}Resistent}'
        Row3SR2 = Row2SR2 = Row1SR2 = '\cellcolor[HTML]{EFEFEF}'

    # Figure out order of drugs
    FirstTot = FirstSens + ['\cellcolor[HTML]{EFEFEF}' + drug for drug in FirstRes]
    SecTot = SecSens + ['\cellcolor[HTML]{EFEFEF}' + drug for drug in SecRes]
    muts = []
    for drug in FirstSens + FirstRes + SecSens + SecRes:
        if drug == 'Fluorokinoloner':
            lookup = 'Quinolones'
        else:
            lookup = drug
        if lookup in resistensdic:
            muts.append('\cellcolor[HTML]{EFEFEF}' + resistensdic[lookup]) # For example \cellcolor[HTML]{EFEFEF}rpoB (S531L)
        else:
            muts.append('Ingen mutasjon detektert')


    text = '''
Type                    & Tolkning                                           & Antibiotika                              & 
Resistensgen \\scriptsize{{(Aminosyreforandring)}}          \\\ 
\\cline{{1-4}}
                              & {Row1SR1}     & {Fdrug1}      & {muts1}      \\\ 
                              & {Row2SR1}     & {Fdrug2}      & {muts2}      \\\ 
                              & {Row3SR1}     & {Fdrug3}      & {muts3}      \\\ 
                              & {Row4SR1}     & {Fdrug4}      & {muts4}      \\\ 
\\multirow{{-5}}{{*}}{{F\\o rstelinje}}  & {Row5SR1} & {Fdrug5} & {muts5}      \\\ 
\\cline{{1-4}}
                              &  {Row1SR2}    & {Sdrug1}      & {muts6}      \\\ 
                              &  {Row2SR2}    & {Sdrug2}      & {muts7}      \\\ 
                              &  {Row3SR2}    & {Sdrug3}      & {muts8}      \\\ 
\\multirow{{-4}}{{*}}{{Andrelinje}}  & {Row4SR2}   & {Sdrug4}  & {muts9}     \\\ 
\\cline{{1-4}}
%--------------- ^ ADD YOUR TABLE CONTENTS ABOVE ^ ---------------------
'''.format(
    Row1SR1 = Row1SR1,
    Row2SR1 = Row2SR1,
    Row3SR1 = Row3SR1,
    Row4SR1 = Row4SR1,
    Row5SR1 = Row5SR1,
    Row1SR2 = Row1SR2,
    Row2SR2 = Row2SR2,
    Row3SR2 = Row3SR2,
    Row4SR2 = Row4SR2,
    Fdrug1 = FirstTot[0],
    Fdrug2 = FirstTot[1],
    Fdrug3 = FirstTot[2],
    Fdrug4 = FirstTot[3],
    Fdrug5 = FirstTot[4],
    Sdrug1 = SecTot[0],
    Sdrug2 = SecTot[1],
    Sdrug3 = SecTot[2],
    Sdrug4 = SecTot[3],
    muts1 = muts[0],
    muts2 = muts[1],
    muts3 = muts[2],
    muts4 = muts[3],
    muts5 = muts[4],
    muts6 = muts[5],
    muts7 = muts[6],
    muts8 = muts[7],
    muts9 = muts[8])

    with open("Latex_template/include/resistens.tex","w") as outfile:
        outfile.write(text)


def CreateBeslektedeOppsummering(clusterjanei):
    if clusterjanei: # TRUE hvis cluster
        text = 'Pr\\o ven var n\\ae rt beslektet med tidligere sekvenserte isolater, noe som antyder \\textbf{nylig smitte}. '
    else:
        text = 'Pr\\o ven var ikke n\\ae rt beslektet med noen av v\\aa re tidligere sekvenserte isolater. '
    with open("Latex_template/include/beslektede_oppsummering.tex","w") as outfile:
        outfile.write(text)

def CreateBeslektede(relationtoothers):
    # Determine if part of cluster
    # Determine number of closely related
    # Determine number of related

    text = '''
Terskelverdi  & Antall tidligere sekvenserte isolater  \\\ \hline 
N\\ae rt beslektet (0 til 5 mutasjoner forskjell) & \\textbf{%s} isolater \\\ 
Beslektet (0 til 12 mutasjoner forskjell) & \\textbf{%s} isolater \\\ \hline 
''' % (str(relationtoothers["close"]), str(relationtoothers["somewhat"]))

    with open("Latex_template/include/beslektede.tex","w") as outfile:
        outfile.write(text)

def CreatePipelineInfo(metainfo):
    # Get metainfo from master file or configure them here
    if "Sekvensator" in metainfo:
        Sekvensator = metainfo["Sekvensator"]
    else:
        Sekvensator = 'Illumina'
    Metode = 'Helgenom'
    Pipeline = 'NIPH TB pipeline v.%s' % (VERSION)
    Referansegenom = 'H37Rv'

    text = '''
Sekvensator & {maskin} & Metode & {metode} \\\ \hline 
Pipeline & {pipeline} & Referansegenom & {ref} \\\ \hline 
'''.format(
        maskin = Sekvensator,
        metode = Metode,
        pipeline = Pipeline,
        ref = Referansegenom)

    with open("Latex_template/include/pipelinedetaljer.tex","w") as outfile:
        outfile.write(text)

def CreateReport(metainfo, resdic, clusterjanei, lineage, species, speciesmash, relationtoothers, covdicsample, samplesindb):
    CreateFooter(metainfo)
    CreateInfo(metainfo, covdicsample, samplesindb)
    CreateOppsummering(resdic, clusterjanei, species, speciesmash)
    CreateTyping(lineage)
    CreateResistensBokser(resdic)
    CreateResistens(resdic)
    CreateBeslektedeOppsummering(clusterjanei)
    CreateBeslektede(relationtoothers)
    CreatePipelineInfo(metainfo)
