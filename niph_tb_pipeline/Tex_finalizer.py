#!/usr/bin/env python2

''' The script that finishes the TeX report. Called from TB_pipeline.py '''

import sys
import os
import time
from __init__ import __version__ as VERSION

def CreateFooter(metainfo):
    currentdate = time.strftime("%Y-%b-%d")
    ID = metainfo["ID"]
    footer = '\\rfoot{Isolat ID: %s | Dato: %s }' % (ID, currentdate)
    with open("Latex_template/include/footer.tex","w") as outfile:
        outfile.write(footer)

def CreateInfo(metainfo, covdicsample):
    # DataQual should be assessed from presence of Fastqc / kaiju problem files
    if os.path.isfile("Fastqc_problems"):
        DataQual = "Lav output"
    elif os.path.isfile("Kaijuclassificationproblem"):
        DataQual = "Ikke MTB"
    elif os.path.isfile("Kaijucontaminationproblem"):
        DataQual = "Mulig kontaminasjon"
    elif covdicsample < 90.00:
        DataQual = "Lav ref. coverage"
    else:
        DataQual = "OK"

    ID = metainfo["ID"]
    Barcode = metainfo["Barcode"]
    Location = metainfo["Location"]
    LocatedFrom = metainfo["Source"]
    Isolationdate = metainfo["Isolated"]

    if os.path.isfile("averagedepth.txt"):
        with open("averagedepth.txt", "rU") as rdfile:
            RD = rdfile.readlines()[0].rstrip("\n")
    else:
        RD = "N/A"

    infostring = '''
    ID for pr\\o ve    &  {SampleName}    & Barcode          & {Barcode}           \\\ \hline
    Sted   & {Location}         & Isolert fra       & {LocatedFrom}          \\\ \hline
    Datakvalitet & {DataQual}       & Isolasjonsdato & {Isolationdate}      \\\ \hline
    Read depth  & {Readdepth}       & Dato  & {currentdate}  \\\ \hline
    '''

    infostringfull = infostring.format(
        SampleName = ID,
        Barcode = Barcode,
        Location = Location,
        LocatedFrom = LocatedFrom,
        DataQual = DataQual,
        Isolationdate = Isolationdate,
        Readdepth = RD,
        currentdate = time.strftime("%Y-%b-%d"))
    with open("Latex_template/include/info.tex","w") as outfile:
        outfile.write(infostringfull)

def CreateOppsummering(resistens, clusterjanei):

    if clusterjanei:
        smitte = 'Pr\\o ven tilh\\o rer et sannsynlig smittecluster, noe som antyder \\textbf{nylig smitte}.'
    else:
        smitte = 'Pr\\o ven er ikke n\\ae rt beslektet med tidligere sekvenserte pr\\o ver.'

    if len(resistens) == 0:
        res = 'Det ble ikke funnet noen resistensmutasjoner.'
    elif len(resistens) == 1:
        res = 'Det ble funnet mutasjoner som indikerer resistens mot %s.' % resistens.keys()[0]
    else:
        reskeys = resistens.keys()
        alleres = ", ".join(["\\textbf{%s}" % r for r in reskeys[:-1]]) + " og \\textbf{%s}." % reskeys[-1]
        res = 'Det ble funnet mutasjoner som indikerer resistens mot {alleres}'.format(alleres=alleres)

    if not os.path.isfile("Kaijuclassificationproblem"):
        oppsummering = 'Pr\\o ven var positiv for \\textbf{Mycobacterium tuberculosis}. %s %s' % (res, smitte)
    else:
        tophit = open("Kaijuclassificationproblem","rU").read()
        oppsummering = 'Pr\\o ven er ikke \\textbf{Mycobacterium tuberculosis}, men \\textbf{%s}. ' % tophit
    with open("Latex_template/include/oppsummering.tex","w") as outfile:
        outfile.write(oppsummering)

def CreateTyping(lineage):
    # Now gathered directly from colltyper
    #lineage = lineage.replace("_","-")
    #if lineage is None:
    if lineage == "":
        typing = 'Pr\\o ven ble ikke entydig typet til en bestemt lineage. '
    else:
        typing = 'Pr\\o ven ble funnet\ \\aa\ tilh\\o re lineage %s. ' % (lineage)
    with open("Latex_template/include/typing.tex","w") as outfile:
        outfile.write(typing)
    
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
Beslektet (0 til 30 mutasjoner forskjell) & \\textbf{%s} isolater \\\ \hline 
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

def CreateReport(metainfo, resdic, clusterjanei, lineage, relationtoothers, covdicsample):
    CreateFooter(metainfo)
    CreateInfo(metainfo, covdicsample)
    CreateOppsummering(resdic, clusterjanei)
    CreateTyping(lineage)
    CreateResistensBokser(resdic)
    CreateResistens(resdic)
    CreateBeslektedeOppsummering(clusterjanei)
    CreateBeslektede(relationtoothers)
    CreatePipelineInfo(metainfo)
