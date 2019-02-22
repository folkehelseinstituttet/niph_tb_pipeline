#!/usr/bin/env python

'''TB pipeline.
- Make sure files are organized properly
- Get path of previous libraries
- Run FastQC. Flag sample if shit quality
    - Open summary.txt in .fastq.gz. If Basic statistics = PASS -> OK else -> WARN
- Run Kaiju. Flag sample if shit quality
- Run snippy pipeline
- Run mykrobe predictor
- Copy results to global library
- Run snippy-core on all strains
- Use snp-dists and find 20 closest strains (for each strain in current sample)
- Create tree using FastTree
- Create strain report from (1) QC. (2) Mykrobe predictor (3) FastTree

To avoid harmful shell injection, do assert no spaces in any files. NO files are allowed to contain space

CURRENT PROBLEMS:
SWITCH TO MASH INSTEAD OF KAIJU? [X]
ONLY SNIPPY OUTPUT IS WRITTEN TO STDOUT? (SOLUTION - FLUSH EVERY PRINT STATEMENT) [X]
SOME PROBLEMS COLLATING RESULTS FROM MASH/MYKROBE PROBABLY [ ]
NEED TO REMOVE TMP DIRECTORY FROM MYKROBE? (OR NOT - IT TAKES 0KB) [X]
MAKE GLOBAL TREE MUCH SIMPLER TO CALCULATE. GTR / LIKELIHOOD TAKES TOO LONG WITH 1K+ ISOLATES - RAPIDNJ?

'''
import os
import sys
import time
import csv
import re
import Bio
import ete3
from subprocess import Popen, PIPE, call
from Bio import AlignIO
from . import Tex_finalizer


KAIJU_NODES_DMP = "/mnt/kaijudb/nodes.dmp"
KAIJU_DB_DMI = "/mnt/kaijudb/kaiju_db_nr.fmi"
KAIJU_NAMES_DMP = "/mnt/kaijudb/names.dmp"
KAIJU_BIN = "/opt/kaiju/bin"
MASH_REFSEQ_SKETCH = "/media/mash/refseq.genomes.k21s1000.msh"
TB_REF = "/mnt/Reference/M_tuberculosis_H37Rv.gb"
TB_EXCLUDECOLS = "/mnt/Reference/Trimal_excludecolumns.txt"
#MCCORTEX31_PATH = "/opt/Mykrobe-predictor/mccortex/bin/mccortex31"
FIGTREE_EXEC = "java -jar /mnt/FigTree_v1.4.3/lib/figtree.jar"
GLOBAL_COLLECTION = "/mnt/global_collection"
TEX_TEMPLATE_DIR = "/mnt/Latex_template"
NUMBER_NEIGHBORS_IN_TREE = 10

def FindSampleName(sample):
    sampleName = sample
    try:
        assert ' ' not in sampleName
    except AssertionError:
        sys.exit("Could not proceed because some file names contain space. Check sample %s " % sampleName)

    # Further checks to make sure there is no injection in filenames
    try:
        assert re.match("^[\w_\-]+/?$", sampleName)
    except AssertionError:
        sys.exit("Isolate names can only contain letters (a-z, A-Z), numbers (0-9), dash (-) and underscore (_)")

    if sample.endswith("/"):
        sampleName = sampleName[:-1]
    if "_" in sample:
        sampleName = sampleName.split("_")[0]
    return sampleName

def FindReads():
    files = os.listdir(".")
    R1 = [s for s in files if "R1_001.fastq" in s]
    R2 = [s for s in files if "R2_001.fastq" in s]

    if len(R1) > 1 or len(R2) > 1:
        sys.exit("Ambiguity error: More than one file matches 'R1_001.fastq' or 'R2_001.fastq'")
    if len(R1) < 1 or len(R2) < 1:
        R1 = [s for s in files if "R1.fastq" in s]
        R2 = [s for s in files if "R2.fastq" in s]
    if len(R1) < 1 or len(R2) < 1:
        sys.exit("Could not locate reads. Verify correct naming ('R1_001.fastq.gz')")
    
    try:
        assert ' ' not in R1[0]
        assert ' ' not in R2[0]
    except AssertionError:
        sys.exit("Could not proceed because some file names contain space. Check sample %s and %s" % (R1[0], R2[0]))

    # Further checks to make sure there is no injection in filenames
    try:
        assert re.match("^[\w_\-]+(R1_001\.fastq|R1\.fastq)", R1[0])
        assert re.match("^[\w_\-]+(R2_001\.fastq|R2\.fastq)", R2[0])
    except AssertionError:
        print("Found %s and %s - Illegal characters in name" % (R1[0], R2[0]), flush=True)
        sys.exit("Isolate names can only contain letters (a-z, A-Z), numbers (0-9), dash (-) and underscore (_)")

    print("Found R1: %s, \t R2: %s" % (R1, R2), flush=True)
    return {"R1": R1[0], "R2": R2[0]}

def ReadSummary(summary):
    try:
        with open(summary,"rU") as f:
            myf = csv.reader(f,delimiter="\t")
            fastqcdic = {line[1]: line[0] for line in myf}
            try:
                if fastqcdic["Per sequence quality scores"] == "FAIL":
                    return 2
                elif fastqcdic["Per sequence quality scores"] == "WARN":
                    return 1
                else:
                    return 0
            except KeyError:
                return -1

            #if myf[:4] == "PASS":
            #    return 0
            #elif myf[:4] == "WARN":
            #    return 1
            #elif myf[:4] == "FAIL":
            #    return 2
            #else:
            #    return -1
    except FileNotFoundError:
        return -1

def RemoveFastqSuffix(file):
    if file.endswith(".fastq"):
        return file[:-6]
    elif file.endswith(".fastq.gz"):
        return file[:-9]
    elif file.endswith(".fq"):
        return file[:-3]
    elif file.endswith(".fq.gz"):
        return file[:-6]
    else:
        sys.exit("Can't identify suffix of file: %s" % file)

def RunFastQC(R1, R2):

    # Scan output files
    R1nosuf = RemoveFastqSuffix(R1)
    R2nosuf = RemoveFastqSuffix(R2)
    R1zip = R1nosuf + "_fastqc.zip"
    R2zip = R2nosuf + "_fastqc.zip"
    R1sum = R1nosuf + "_fastqc/summary.txt"
    R2sum = R2nosuf + "_fastqc/summary.txt"
    # Check if files exist already:
    if os.path.isfile(R1zip) and os.path.isfile(R2zip):
        print("FastQC results already exists in %s" % (os.getcwd()), flush=True)
        return 0
    
    # Run FastQC
    print("Running cmd: fastqc %s %s in dir %s" % (R1, R2, os.getcwd()), flush=True)
    errorcode = call("fastqc --extract %s %s" % (R1, R2), shell=True)
    if errorcode != 0:
        sys.exit("FastQC did not complete correctly.")

    #call("unzip -o -j %s %s -d '.'" % (R1zip, R1sum), shell=True)
    R1info = ReadSummary(R1sum)
    #call("unzip -o -j %s %s -d '.'" % (R2zip, R2sum), shell=True)
    R2info = ReadSummary(R2sum)
    if R1info != 0 or R2info != 0:
        # The presence of a Fastqc_problems file indicates a problem with the sequence data
        out = open("Fastqc_problems","w")
        if (R1info == 2 or R2info == 2):
            out.write("Error rate > 1%")
        elif (R1info == 1 or R2info == 1):
            out.write("Error rate > .2%")
        else:
            out.write("Problem reading FastQC file")
        out.close()

def splitKaijuReportLine(line):
    linelist = line.split("\t")
    linelist = [s.strip(" \n") for s in linelist]
    return {"Species": linelist[2], "Percentage": float(linelist[0]), "Reads": int(linelist[1])}


def AnalyzeKaijuReport(myfile):
    with open(myfile,"rU") as mf:
        content = mf.readlines()
        mostcontent = splitKaijuReportLine(content[2])
        runnerup = splitKaijuReportLine(content[3])
        # Completely arbitrary cutoff that highest hit has to have 10 times as many as runnerup
        if "Mycobacterium" not in mostcontent["Species"]:
            of = open("Kaijuclassificationproblem","w")
            of.write(mostcontent["Species"] + "\n")
            of.close()
        else:
            if mostcontent["Species"] != "Mycobacterium tuberculosis":
                of = open("Kaijuothermycobacterium","w")
                of.write(mostcontent["Species"] + "\n")
                of.close()
            elif mostcontent["Reads"] < 10 * runnerup["Reads"]:
                if not "Mycobacterium" in runnerup["Species"]:
                    open("Kaijucontaminationproblem","w")
            # Else - Do not write any files. Top species is MTB with at least 10X the n reads of runnerup

def RunKaiju(R1, R2):
    # Run Kaiju
    print("Checking species ID with kaiju")

    # Check if file exists already:
    if os.path.isfile("kaiju.summary"):
        print("Kaiju results already exists in %s" % os.getcwd())
        return 0

    #errorcode = call("%s/kaiju -x -t %s -f %s -i %s -j %s -o kaiju.out -z 8" % (KAIJU_BIN, KAIJU_NODES_DMP, KAIJU_DB_DMI, R1, R2), shell=True)
    errorcode = call("kaiju -x -t %s -f %s -i %s -j %s -o kaiju.out -z 8" % (KAIJU_NODES_DMP, KAIJU_DB_DMI, R1, R2), shell=True)
    errorcode2 = call("kaijuReport -i kaiju.out -o kaiju.summary -t %s -n %s -r species" % (KAIJU_NODES_DMP, KAIJU_NAMES_DMP), shell=True)

    try:
        call("rm kaiju.out",shell=True)
    except:
        sys.exit("Failed to remove kaiju.out file in %s" % os.getcwd())

    # Analyze Kaiju
    AnalyzeKaijuReport("kaiju.summary")

def RunMash(R1, R2):
    '''Convert species identification to MASH instead.'''
    print("Checking species ID with MASH", flush=True)
    if os.path.isfile("mashreport.tab"):
        print("Mash results already exists in %s" % os.getcwd(), flush=True)
        AnalyzeMashTopHit("mashreport.tab")
        return 0

    errorcode1 = call("mash screen -w -p 4 %s %s %s > mashreport.tab" % (MASH_REFSEQ_SKETCH, R1, R2), shell=True)
    errorcode2 = call("sort -gr mashreport.tab > mashreport.sorted.tab", shell=True)
    errorcode3 = call("head mashreport.sorted.tab > mashreport.tab", shell=True)
    errorcode4 = call("rm mashreport.sorted.tab", shell=True)
    AnalyzeMashTopHit("mashreport.tab")

def ReadMashTopHit(mashreport):
    # NOT USED
    with open(mashreport,"rU") as mf:
        data = csv.reader(mf, delimiter="\t")
        tophit = next(data)
        topmatch = CaptureMashHit(tophit[5])

    return topmatch

def SplitMashReportLine(line):
    pattern = "[A-Z][a-z]* [a-z]*( strain)?( BCG)?(-1)?( subsp\. \w+)?( phiX174)?"
    return {"identity": line[0], "sharedhashes":line[1], "medianmultiplicity":line[2], "pvalue":line[3], "queryid":line[4],"querycomment":line[5], "species": re.search(pattern,line[5])[0]}

def AnalyzeMashTopHit(myfile):
    with open(myfile, "rU") as mf:
        content = csv.reader(mf,delimiter="\t")
        mostcontent = SplitMashReportLine(next(content))
        runnerup = SplitMashReportLine(next(content))
        # Switch it up if phiX is top hit
        if "Enterobacteria phage phiX174" in mostcontent["species"]:
            mostcontent, runnerup = runnerup, mostcontent

        if "Mycobacterium" not in mostcontent["species"]:
            of = open("Mashclassificationproblem","w")
            of.write(mostcontent["species"] + "\n")
            of.close()
        else:
            if "Mycobacterium tuberculosis" not in mostcontent["species"]:
                of = open("Mashothermycobacterium","w")
                of.write(mostcontent["species"] + "\n")
                of.close()
        # ELSE: Do not write any files. Top species is MTB. Runner-up info not used (No contamination screen)


def CaptureMashHit(string):
    '''Method for capturing a binomial name from MASH screen results.'''
    pattern = "[A-Z][a-z]* [a-z]*( strain)?( BCG)?(-1)?( subsp\. \w+)?"
    match = re.search(pattern,string)[0]
    return match


def RunSnippy(R1, R2):
    # Check if snippy dir exists already:
    if os.path.isdir("snippy"):
        print("Snippy results already exits in %s" % os.getcwd(), flush=True)
        return 0
    errorcode = call("snippy --outdir ./snippy --ref %s --R1 %s --R2 %s" % (TB_REF, R1, R2), shell=True) # REMOVED --cleanup (Needed for samtools depth - Remove BAM files later in script)

def FindCoverage():
    print("Checking coverage", flush=True)
    if os.path.isfile("averagedepth.txt"):
        print("Average depth already calculated", flush=True)
        return 0
    #errorcode = call("gzip -cd snippy/snps.depth.gz | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR; print sqrt(sumsq/NR - (sum/NR)**2)}' > averagedepth.txt", shell=True)
    errorcode1 = call("samtools depth -aa snippy/snps.bam > snippy/snps.depth", shell=True)
    #errorcode = call("gzip -cd snippy/snps.depth.gz | awk '{sum+=$3} END { print sum/NR}' > averagedepth.txt", shell=True)
    errorcode2 = call("cat snippy/snps.depth | awk '{sum+=$3} END { print sum/NR}' > averagedepth.txt", shell=True)

def CleanupSnippyData():
    print("Cleaning up snippy data", flush=True)
    if not os.path.isdir("snippy/reference/ref/"):
        print("Snippy data already cleaned up", flush=True)
        return 0
    errorcode1 = call("rm -rf snippy/reference/ref/", shell=True)
    errorcode2 = call("rm -rf snippy/reference/genomes/", shell=True)
    errorcode3 = call("rm snippy/reference/ref.fa.*", shell=True)
    #errorcode4 = call("rm snippy/snps.*.tbi", shell=True)
    errorcode4 = call("rm snippy/snps.depth", shell=True)
    errorcode5 = call("rm snippy/snps.*.gz", shell=True)
    errorcode6 = call("rm snippy/snps.consensus*.fa", shell=True)
    errorcode7 = call("rm snippy/snps.bam", shell=True)
    errorcode8 = call("rm snippy/snps.bam.bai", shell=True)
    # Dont remove raw.vcf - Has important heterozygous information
    #errorcode9 = call("rm snippy/snps.raw.vcf", shell=True)

def RunMykrobe(R1, R2, sampleName):
    # Check if mykrobe predictor results already exists:
    if os.path.isfile("mykrobe_output.csv"):
        print("Mykrobe predictor results already exists in %s" % os.getcwd(), flush=True)
        return 0
    #errorcode1 = call("mykrobe predict %s tb --mccortex31_path %s -1 %s %s > mykrobe_output.json" % (sampleName, MCCORTEX31_PATH, R1, R2), shell=True)
    errorcode1 = call("mykrobe predict %s tb --output mykrobe_output.csv --format csv -1 %s %s" % (sampleName, R1, R2), shell=True)
    #errorcode2 = call("json_to_tsv mykrobe_output.json > mykrobe_output.tsv", shell=True)
    #errorcode3 = call("rm -rf atlas", shell=True)
    errorcode4 = call("rm -rf mykrobe", shell=True)
    if os.path.isfile("mykrobe_output.tsv"):
        errorcode5 = call("rm mykrobe_output.json",shell=True)
    if os.path.isdir("tmp"):
        errorcode6 = call("rm -rf tmp", shell=True)

def CollType():
    print("Typing according to Coll (2014) scheme", flush=True)
    if os.path.isfile("colltype.txt"):
        print("Coll type already calculated", flush=True)
        return 0
    #errorcode = call("colltyper -o colltype.txt snippy/snps.vcf", shell=True)
    # USE unfiltered VCF file to find mixed infections
    errorcode = call("colltyper -o colltype.txt snippy/snps.raw.vcf", shell=True)

def sampleAnalysis(sample):
    
    sampleName = FindSampleName(sample)
    print("Sample name: %s" % sampleName, flush=True)
    try:
        os.chdir(sample)
    except:
        sys.exit("Could not access dir: %s" % sample)
    print("Current dir: %s " % os.getcwd(), flush=True)
    myfiles = FindReads()
    try:
        assert ' ' not in myfiles["R1"]
        assert ' ' not in myfiles["R2"]
    except AssertionError:
        sys.exit("Could not proceed because some file names contain space. Check %s and %s" % (myfiles["R1"], myfiles["R2"]))

    RunFastQC(myfiles["R1"], myfiles["R2"])
    #RunKaiju(myfiles["R1"], myfiles["R2"])
    RunMash(myfiles["R1"], myfiles["R2"])
    RunSnippy(myfiles["R1"], myfiles["R2"])
    FindCoverage()
    CleanupSnippyData()
    CollType()
    RunMykrobe(myfiles["R1"], myfiles["R2"], sampleName)
    try:
        os.chdir("..")
    except:
        sys.exit("Failed to go out of directory: %s" % sample)

def CopyToGlobalDir(sample):
    '''Method that copies results to global dir. Needs to be amended because no longer write access'''

    if not os.path.isdir("%s/%s" % (GLOBAL_COLLECTION, sample)):
        call("mkdir %s/%s" % (GLOBAL_COLLECTION, sample), shell=True)
    files = os.listdir(sample)
    for fil in files:
        if fil.endswith(".fastq"):
            continue
        if fil.endswith(".gz"):
            continue
        if fil.endswith(".out"):
            continue
        if fil.endswith("snippy"):
            continue
        if fil.endswith("atlas"):
            continue
        if fil.endswith("tmp"):
            continue
        if fil.endswith("fastqc"):
            continue
        if fil.endswith("template"):
            continue
        if not os.path.isfile("%s/%s/%s" % (GLOBAL_COLLECTION, sample, fil)):
            call("cp -rf %s/%s %s/%s/%s" % (sample, fil, GLOBAL_COLLECTION, sample, fil), shell=True)
    if not os.path.isfile("%s/%s/snps.tab" % (GLOBAL_COLLECTION, sample)):
        call("cp %s/snippy/snps.tab %s/%s/snps.tab" % (sample, GLOBAL_COLLECTION, sample), shell=True)
    if not os.path.isfile("%s/%s/snps.aligned.fa" % (GLOBAL_COLLECTION, sample)):
        call("cp %s/snippy/snps.aligned.fa %s/%s/snps.aligned.fa" % (sample, GLOBAL_COLLECTION, sample), shell=True)
    # call("ln -s %s %s/%s/" % (SNIPPY_REF_DIR, GLOBAL_COLLECTION, sample), shell=True) # NOT NEEDED FOR DOCKER VERSION

def CopyForEasyGlobalDirMove(sample):
    '''Method that copies all relevant GLOBAL_COLLECTION data to one folder, and additionally all relevant reports to another folder'''
    if not os.path.isdir("COPY_TO_TB_PIPELINE_DATABASE/%s" % sample):
        call("mkdir COPY_TO_TB_PIPELINE_DATABASE/%s" % sample, shell=True)
    files = os.listdir(sample)
    for fil in files:
        if fil.endswith(".fastq"):
            continue
        if fil.endswith(".gz"):
            continue
        if fil.endswith(".out"):
            continue
        if fil.endswith("snippy"):
            continue
        if fil.endswith("atlas"):
            continue
        if fil.endswith("tmp"):
            continue
        if fil.endswith("fastqc"):
            continue
        if fil.endswith("template"):
            continue
        if fil.endswith("pdf"):
            call("cp -rf %s/%s COPY_TO_REPORTS/%s" % (sample, fil, fil), shell=True)
        if not os.path.isfile("COPY_TO_TB_PIPELINE_DATABASE/%s/%s" % (sample, fil)):
            call("cp -rf %s/%s COPY_TO_TB_PIPELINE_DATABASE/%s/%s" % (sample, fil, sample, fil), shell=True)

def CopySnippyDataToShallowDir(sample):
    if os.path.isfile("%s/snippy/snps.tab" % sample):
        errorcode1 = call("mv %s/snippy/snps.tab %s/snps.tab" % (sample, sample), shell=True)
    if os.path.isfile("%s/snippy/snps.aligned.fa" % sample):
        errorcode2 = call("mv %s/snippy/snps.aligned.fa %s/snps.aligned.fa" % (sample, sample), shell=True)
    if os.path.isfile("%s/snippy/snps.vcf" % sample):
        errorcode1 = call("mv %s/snippy/snps.vcf %s/snps.vcf" % (sample, sample), shell=True)

def MaskRepetitiveRegions(alnfile):
    '''This method is not yet complete'''
    outfilename = alnfile[:-8] + 'masked.fasta'
    excludefile = TB_EXCLUDECOLS
    with open(excludefile, 'rU') as exfile:
        exclude = exfile.read()
    errorcode = call("trimal -in %s -out %s -selectcols %s " % (alnfile, outfilename, exclude), shell=True)
    return outfilename

def replaceOldSNPalignment(maskedfile, filetoreplace):
    call("snp-sites -o %s %s" % (filetoreplace, maskedfile), shell=True)


def RunSnippyCore(basedir, timestamp):
    # Change so that this analysis is not run in GLOBAL dir but local
    #try:
    #    os.chdir(GLOBAL_COLLECTION)
    #except:
    #    sys.exit("Unable to move to global collection dir: %s" % GLOBAL_COLLECTION)
    print("Running snippy-core", flush=True)

    dirs = next(os.walk(basedir))[1]
    dirs.remove("COPY_TO_REPORTS")
    dirs.remove("COPY_TO_TB_PIPELINE_DATABASE")
    dirs = " ".join(dirs)
    globaldirs = next(os.walk(GLOBAL_COLLECTION))[1]
    if len(globaldirs) == 0:
        globaldirstxt = ""
    else:
        globaldirstxt = " ".join([GLOBAL_COLLECTION + "/" + g for g in globaldirs])

    if not os.path.isfile("snippy-core.log"):
        try:
            #errorcode1 = call("snippy-core --prefix=TB_all_%s --ref=/mnt/Reference/ref.fa --mask=/mnt/Reference/Trimal_excludecolumns.bed %s %s 2> snippy-core.log" % (timestamp, globaldirstxt, dirs), shell=True)
            errorcode1 = call("snippy-core --prefix=TB_all_%s --ref=/mnt/Reference/ref.fa %s %s 2> snippy-core.log" % (timestamp, globaldirstxt, dirs), shell=True)
        except:
            #print("Command failed: snippy-core --prefix=TB_all_%s --ref=/mnt/Reference/ref.fa --mask=/mnt/Reference/Trimal_excludecolumns.bed %s %s 2> snippy-core.log" % (timestamp, globaldirstxt, dirs))
            print("Command failed: snippy-core --prefix=TB_all_%s --ref=/mnt/Reference/ref.fa %s %s 2> snippy-core.log" % (timestamp, globaldirstxt, dirs), flush=True)
            sys.exit("Snippy-core failed for unknown reason. See snippy-core log file")
    # Mask bad regions in full alignment
    # (LEGACY) (No longer needed since snippy-core can do this automatically) NOTE - NOT TRUE - STILL USE
    if os.path.isfile("TB_all_%s.full.aln" % timestamp):
        if not os.path.isfile("TB_all_%s.masked.fasta" % timestamp):
            maskedfile = MaskRepetitiveRegions("TB_all_%s.full.aln" % timestamp)
    
    # Replace basic SNP alignment with masked FULL alignment
    replaceOldSNPalignment("TB_all_%s.masked.fasta" % timestamp, "TB_all_%s.aln" % timestamp)

    # Discard whole-genome alignment
    if os.path.isfile("TB_all_%s.full.aln" % timestamp):
        errorcode2 = call("rm TB_all_%s.full.aln" % timestamp, shell=True)

    # Create continuously evolving global tree
    # REMOVE IN VERSION 0.9 - HAS BECOME TOO SLOW
    #if not os.path.isfile("Global_collection_tree.nwk"):
    #    errorcode3 = call("FastTree -nt -gtr TB_all_%s.aln > Global_collection_tree.nwk" % (timestamp), shell=True)
    
    # Remove masked file
    # (LEGACY) (No longer needed since Snippy-core does this)
    #if os.path.isfile("TB_all_%s.masked.fasta" % timestamp):
    #    errorcode4 = call("rm TB_all_%s.masked.fasta", shell=True)
    
    #errorcode4 = call("mv TB_all* %s" % basedir, shell=True)
    #errorcode5 = call("mv snippy-core.log %s" % basedir, shell=True)
    #try:
    #    os.chdir(basedir)
    #except:
    #    sys.exit("Unable to move back from global collection to %s " % basedir)

def ReadSnippyCoreLog():
    """ LEGACY method. Reads snippy-core log to get % coverage"""
    with open("snippy-core.log","rU") as snippycorelogfile:
        print("Reading snippy-core log", flush=True)
        snippyloglines = snippycorelogfile.readlines()
        covdic = {}
        for l in snippyloglines:
            if "coverage" not in l:
                continue
            else:
                l2 = l.split("\t")[1] # Get information part
                l2s = l2.split(" ")
                covdic[l2s[0]] = float(l2s[4][:-2]) # Remove % and \n
    return covdic

def ReadSnippyCoreCov(timestamp):
    """ NEW method for getting % cov from snippy-core. Uses .txt file"""
    with open("TB_all_%s.txt" % (timestamp), "rU") as txtfile:
        print("Reading coverage information", flush=True)
        data = csv.reader(txtfile, delimiter="\t")
        header = next(data)
        #covdickeys = {key:i for i,key in enumerate(header)}
        covdic = {}
        for line in data:
            covdic[ line[0] ] = { header[i]: int(line[i]) for i in range(1,len(header))  }
        return covdic

def RunSnpDists():
    ''' Should be based on FULL alignment (with masked regions) rather than post snp-dists version'''
    print("Finding distances between all isolates in global collection", flush=True)
    errorcode = call("snp-dists -c TB_all*.masked.fasta > TB_all_dists.csv", shell=True)

def ReadSnpDistsObject(dists):
    header = next(dists)[1:]
    res = {}
    with open("TB_all_dists.phylip","w") as phylipfile:
        writer = csv.writer(phylipfile, delimiter="\t")
        samplesindb = dists.line_num - 1
        writer.writerow([str(samplesindb)])
        for row in dists:
            res[row[0]] = {}
            writer.writerow(row)
            for i in range(len(header)):
                res[row[0]][header[i]] = int(row[i+1])
    return res

def FindNeighbors(sampledists, threshold):
    Neighbors = []
    # Threshold = # closest isolates
    values = sorted(sampledists.values())
    if len(values) <= threshold:
        return [key for key in sampledists]
    else:
        cutoff = values[threshold-1] # e.g. values # 20
        Neighbors = [key for key in sampledists if sampledists[key] <= cutoff]
        return Neighbors[:threshold]

def MakeTree(sample, Neighbors, all_snps):
    os.chdir("./%s" % (sample))
    # Select which Bio.AlignIO seqs should be part of new alignment
    templist = []
    for record in all_snps:
        if record.id == sample or record.id in Neighbors:
            templist.append(record)
    with open("Neighbors.aln","w") as outfile:
        Bio.AlignIO.write(Bio.Align.MultipleSeqAlignment(records=templist,alphabet=Bio.Alphabet.SingleLetterAlphabet()), outfile, "fasta")

    errorcode1 = call("snp-sites -m -o Neighbors_SNPs.aln Neighbors.aln", shell=True)
    errorcode2 = call("FastTree -nt -gtr Neighbors_SNPs.aln > NeighborTree.nwk", shell=True)
    errorcode3 = call("%s -graphic PNG -width 750 -height 400 NeighborTree.nwk NeighborTree.png" % FIGTREE_EXEC, shell=True)
    errorcode4 = call("convert NeighborTree.png -background white -alpha remove -alpha off NeighborTreeWhite.png", shell=True) # Convert is imagemagick...

    try:
        os.chdir("..")
    except:
        sys.exit("Unable to move back from directory %s" % sample)

def CopyTexTemplate():
    errorcode = call("cp -r %s ." % TEX_TEMPLATE_DIR, shell=True)

def HandleMutation(mutation):
    # Promotor mutations can have dash (-) in them that should not be split (minus character)
    mut = mutation.split(':')[0]
    if mut.count('-') > 1:
        # Promotor mutation
        mutation = ''.join(mutation.split('-')[:-1])
        mutation = mutation.replace('_',' (') + ')'
        mutation = mutation.replace('-','')
        return mutation
    else:
        # Protein mutation
        mutation = mutation.split('-')[0]
        mutation = mutation.replace('_',' (') + ')'
        return mutation

def PimpResDic(mykrobetsvfile):
    with open(mykrobetsvfile, "rU") as infile:
        data = csv.reader(infile, delimiter=",")
        header = next(data)
        drugcol = header.index('drug')
        mutcol = [ i for i, word in enumerate(header) if word.startswith('variant') ][0]
        resistensdic = {}
        for row in data:
            # Go from katG_S315T-S315T:66:1:266 to katG (S315T)
            mutation = row[mutcol]
            if mutation == '':
                # Dont enter drug if no resistance found
                continue
            else:
                # Can have multiple mutations split by ;
                muts = mutation.split(';')
                resistensdic[row[drugcol]] = ', '.join([HandleMutation(m) for m in muts])

        return resistensdic

def GetLineage(mykrobetsvfile):
    """
    LEGACY METHOD. NO LONGER USED.
    """
    with open(mykrobetsvfile, "rU") as infile:
        data = csv.reader(infile,delimiter=",")
        header = next(data)
        lineagecol = header.index('lineage')
        return next(data)[lineagecol]

def GetSpeciesMykrobe(mykrobetsvfile):
    '''Get species from Mykrobe predictor'''
    with open(mykrobetsvfile, "rU") as infile:
        data = csv.reader(infile,delimiter=",")
        header = next(data)
        speciescol = header.index('species')
        binomial_name = re.sub("_"," ", next(data)[speciescol])
        return binomial_name

def GetLineageColl(colltyperfile):
    with open(colltyperfile, "rU") as infile:
        data = csv.reader(infile,delimiter="\t")
        header = next(data)
        lineagecol = header.index('Lineage')
        lineage = ""
        for row in data:
            lineage += (row[lineagecol]) + " / "
        if lineage.endswith(" / "):
            return lineage[:-3]
        else:
            return lineage


def IsItInCluster(sample, dists):
    for key in dists:
        if key == sample:
            continue
        if int(dists[key]) <= 30:
            return True
    return False

def NumberRelated(sample, dists):
    related = {"close": 0, "somewhat" : 0}
    for key in dists:
        if key == sample:
            continue
        if int(dists[key]) <= 12:
            if int(dists[key]) > 5:
                related["somewhat"] += 1
            else:
                related["close"] += 1
                related["somewhat"] += 1
    return related



def FinalizeSampleReport(sample, metainfo, resdic, clusterjanei, lineage, species, speciesmash, relationtoothers, covdicsample):
    '''Runs TeX scripts to actually create the report. Needs to be amended due to no write access to GLOBAL_COLLECTION'''
    try:
        os.chdir("./%s" % (sample))
    except:
        sys.exit("Could not move into directory %s" % sample)
        
    CopyTexTemplate()

    # Copy tree file to tex directory
    call("cp NeighborTreeWhite.png Latex_template/imageFiles/tree.png",shell=True)

    # Run TB-finalizer to create tex files
    Tex_finalizer.CreateReport(metainfo, resdic, clusterjanei, lineage, species, speciesmash, relationtoothers, covdicsample)
    # Create the pdf
    try:
        os.chdir("Latex_template")
    except:
        sys.exit("Unable to move in tex directory of %s" % sample)
    call("pdflatex tb-wgs-report.tex > /dev/null 2>&1", shell=True)
    call("cp tb-wgs-report.pdf ../%s-tb-wgs-report.pdf" % sample, shell=True)
    # NOTE: Following needs to be amended due to no write access
    #call("cp tb-wgs-report.pdf %s/%s/%s-tb-wgs-report.pdf" % (GLOBAL_COLLECTION, sample, sample), shell=True)
    try:
        os.chdir("../..")
    except:
        sys.exit("Unable to move back from directory %s" % sample)

def main():
    print("Starting", flush=True)
    # Verify that Metadata.csv exists
    metadataexists = True
    if not os.path.isfile("Metadata.csv"):
        metadataexists = False
        print("WARNING: Unable to locate Metadata.csv in root directory", flush=True)

    basedir = os.getcwd()
    timestamp = time.strftime("%d_%b_%Y")

    dirs = next(os.walk(basedir))[1]
    if "COPY_TO_TB_PIPELINE_DATABASE" in dirs:
        dirs.remove("COPY_TO_TB_PIPELINE_DATABASE")
    if "COPY_TO_REPORTS" in dirs:
        dirs.remove("COPY_TO_REPORTS")
    # Create directory that will hold data to go into global database
    if not os.path.isdir("COPY_TO_TB_PIPELINE_DATABASE"):
        call("mkdir COPY_TO_TB_PIPELINE_DATABASE", shell=True)
    # Create directory that will hold all new reports
    if not os.path.isdir("COPY_TO_REPORTS"):
        call("mkdir COPY_TO_REPORTS", shell=True)

    for sample in dirs:
        if "COPY_TO_TB_PIPELINE_DATABASE" in sample:
            continue
        if "COPY_TO_REPORTS" in sample:
            continue
        print("Running sample %s" % sample, flush=True)

        sampleAnalysis(sample)

        # When done with individual analyses, copy all results to Global dir
        # NOTE: Removed CopyToGlobalDir function due to no write access
        #CopyToGlobalDir(sample)
        CopySnippyDataToShallowDir(sample)
    
    RunSnippyCore(basedir, timestamp)
    RunSnpDists()

    # For each sample, find the 20 closest and draw a FastTree
    with open("TB_all_dists.csv","rU") as distfile:
        print("Reading distances between all isolates in global collection", flush=True)
        dists = csv.reader(distfile, delimiter=",")
        usedists = ReadSnpDistsObject(dists)

    print("Writing NJ tree of all isolates in DB", flush=True)
    try:
        call("rapidnj TB_all_dists.phylip -i pd > Global_collection_tree.nwk", shell=True)
    except:
        print("Unable to write global collection tree", flush=True)

    # Load Biopython alignment of all snps
    with open("TB_all_%s.aln" % timestamp, "rU") as snpfile:
        print("Loading global alignment of SNPs", flush=True)
        all_snps = Bio.AlignIO.read(snpfile, "fasta")

    
    # Load sample metadata
    if metadataexists:
        with open("Metadata.csv","rU") as metainfofile:
            print("Loading metadata", flush=True)
            metainfo = csv.reader(metainfofile,delimiter=",")
            metainfodic = {}
            header = next(metainfo)
            for row in metainfo:
                rowdic = {header[i]: row[i] for i in range(len(row))}
                metainfodic[row[0]] = rowdic
    else:
        metainfodic = {}
        for sample in dirs:
            metainfodic[sample] = {"ID": sample, "Barcode":"","Location":"","Source":"","Isolated":""}

    # with open("snippy-core.log","rU") as snippycorelogfile:
    ## LEGACY - moved to ReadSnippyCoreLog
    #     print("Reading snippy-core log")
    #     snippyloglines = snippycorelogfile.readlines()
    #     covdic = {}
    #     for l in snippyloglines:
    #         if "coverage" not in l:
    #             continue
    #         else:
    #             l2 = l.split("\t")[1] # Get information part
    #             l2s = l2.split(" ")
    #             covdic[l2s[0]] = float(l2s[4][:-2]) # Remove % and \n

    covdic = ReadSnippyCoreCov(timestamp)

    print("Finalizing reports for each strain", flush=True)

    for sample in dirs:
        Neighbors = FindNeighbors(usedists[sample], NUMBER_NEIGHBORS_IN_TREE)
        MakeTree(sample, Neighbors, all_snps)

        # Pimp resdic
        pimpedresdic = PimpResDic('%s/mykrobe_output.csv' % sample)

        # Get lineage
        # Consider implementing COLL scheme instead
        #lin = GetLineage('%s/mykrobe_output.tsv' % sample)
        lin = GetLineageColl("%s/colltype.txt" % sample)
        species = GetSpeciesMykrobe('%s/mykrobe_output.csv' % sample)
        speciesmash = ReadMashTopHit("%s/mashreport.tab" % sample)

        # Find out if sample is part of a cluster.
        # Lowest distance is always to self (0)
        relationtoothers = NumberRelated(sample, usedists[sample])
        clusterjanei = relationtoothers["somewhat"] > 0

        FinalizeSampleReport(sample, metainfodic[sample], pimpedresdic, clusterjanei, lin, species, speciesmash, relationtoothers, covdic[sample])

        # Finally, copy data that goes to global directory/reports
        CopyForEasyGlobalDirMove(sample)
        
    print("Finished")

if __name__ == '__main__':
    main()