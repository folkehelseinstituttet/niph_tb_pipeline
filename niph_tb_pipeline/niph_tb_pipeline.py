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
- Xres box not checked?
- Still need to mask PE/PPE regions
- Save time by not zipping fastqc output? Use --extract flag when running fastq. Can then remove unzip commands [X]
- Colltyper not latest version in image
- Make sure permissions are set so that we have write access on Felles. As of 23 Feb, even olbb in Mint does not have write permission.
- Fixed json_to_tsv [X]

'''
import os
import sys
import time
import csv
import Bio
import ete3
from subprocess import Popen, PIPE, call
from Bio import AlignIO
from . import Tex_finalizer


KAIJU_NODES_DMP = "/mnt/kaijudb/nodes.dmp"
KAIJU_DB_DMI = "/mnt/kaijudb/kaiju_db_nr.fmi"
KAIJU_NAMES_DMP = "/mnt/kaijudb/names.dmp"
KAIJU_BIN = "/opt/kaiju/bin"
TB_REF = "/mnt/Reference/M_tuberculosis_H37Rv.gb"
TB_EXCLUDECOLS = "/mnt/Reference/Trimal_excludecolumns.txt"
MCCORTEX31_PATH = "/opt/Mykrobe-predictor/mccortex/bin/mccortex31"
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
    
    print("Found R1: %s, \t R2: %s" % (R1, R2))
    return {"R1": R1[0], "R2": R2[0]}

def ReadSummary(summary):
    with open(summary,"rU") as f:
        myf = f.read()
        if myf[:4] == "PASS":
            return 0
        elif myf[:4] == "WARN":
            return 1
        elif myf[:4] == "FAIL":
            return 2
        else:
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
        print("FastQC results already exists in %s" % (os.getcwd()))
        return 0
    
    # Run FastQC
    print("Running cmd: fastqc %s %s in dir %s" % (R1, R2, os.getcwd()))
    errorcode = call("fastqc --extract %s %s" % (R1, R2), shell=True)
    if errorcode != 0:
        sys.exit("FastQC did not complete correctly.")

    #call("unzip -o -j %s %s -d '.'" % (R1zip, R1sum), shell=True)
    R1info = ReadSummary(R1sum)
    #call("unzip -o -j %s %s -d '.'" % (R2zip, R2sum), shell=True)
    R2info = ReadSummary(R2sum)
    if R1info != 0 or R2info != 0:
        # The presence of a Fastqc_problems file indicates a problem with the sequence data
        open("Fastqc_problems","w")

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
        if mostcontent["Species"] != "Mycobacterium tuberculosis":
            of = open("Kaijuclassificationproblem","w")
            of.write(mostcontent["Species"] + "\n")
            of.close()
        if mostcontent["Reads"] < 10 * runnerup["Reads"]:
            open("Kaijucontaminationproblem","w")

def RunKaiju(R1, R2):
    # Run Kaiju
    print("Checking species ID with kaiju")

    # Check if file exists already:
    if os.path.isfile("kaiju.summary"):
        print("Kaiju results already exists in %s" % os.getcwd())
        return 0

    errorcode = call("%s/kaiju -x -t %s -f %s -i %s -j %s -o kaiju.out -z 8" % (KAIJU_BIN, KAIJU_NODES_DMP, KAIJU_DB_DMI, R1, R2), shell=True)
    #errorcode = call("kaiju -x -t %s -f %s -i %s -j %s -o kaiju.out -z 8" % (KAIJU_NODES_DMP, KAIJU_DB_DMI, R1, R2), shell=True)
    errorcode2 = call("%s/kaijuReport -i kaiju.out -o kaiju.summary -t %s -n %s -r species" % (KAIJU_BIN, KAIJU_NODES_DMP, KAIJU_NAMES_DMP), shell=True)

    try:
        call("rm kaiju.out",shell=True)
    except:
        sys.exit("Failed to remove kaiju.out file in %s" % os.getcwd())

    # Analyze Kaiju
    AnalyzeKaijuReport("kaiju.summary")


def RunSnippy(R1, R2):
    # Check if snippy dir exists already:
    if os.path.isdir("snippy"):
        print("Snippy results already exits in %s" % os.getcwd())
        return 0
    errorcode = call("snippy --outdir ./snippy --ref %s --pe1 %s --pe2 %s" % (TB_REF, R1, R2), shell=True) # REMOVED --cleanup (Needed for samtools depth - Remove BAM files later in script)

def FindCoverage():
    print("Checking coverage")
    if os.path.isfile("averagedepth.txt"):
        print("Average depth already calculated")
        return 0
    #errorcode = call("gzip -cd snippy/snps.depth.gz | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR; print sqrt(sumsq/NR - (sum/NR)**2)}' > averagedepth.txt", shell=True)
    errorcode = call("gzip -cd snippy/snps.depth.gz | awk '{sum+=$3} END { print sum/NR}' > averagedepth.txt", shell=True)

def CleanupSnippyData():
    print("Cleaning up snippy data")
    errorcode1 = call("rm -rf snippy/reference/ref/", shell=True)
    errorcode2 = call("rm -rf snippy/reference/genomes/", shell=True)
    errorcode3 = call("rm snippy/reference/ref.fa.*", shell=True)
    errorcode4 = call("rm snippy/snps.*.tbi", shell=True)
    errorcode5 = call("rm snippy/snps.*.gz", shell=True)
    errorcode6 = call("rm snippy/snps.consensus*.fa", shell=True)
    errorcode7 = call("rm snippy/snps.bam", shell=True)
    errorcode8 = call("rm snippy/snps.bam.bai", shell=True)
    errorcode9 = call("rm snippy/snps.raw.vcf", shell=True)

def RunMykrobe(R1, R2, sampleName):
    # Check if mykrobe predictor results already exists:
    if os.path.isfile("mykrobe_output.tsv"):
        print("Mykrobe predictor results already exists in %s" % os.getcwd())
        return 0
    errorcode1 = call("mykrobe predict %s tb --mccortex31_path %s -1 %s %s > mykrobe_output.json" % (sampleName, MCCORTEX31_PATH, R1, R2), shell=True)
    #errorcode1 = call("mykrobe predict %s tb -1 %s %s > mykrobe_output.json" % (sampleName, R1, R2), shell=True)
    errorcode2 = call("json_to_tsv mykrobe_output.json > mykrobe_output.tsv", shell=True)
    errorcode3 = call("rm -rf atlas", shell=True)
    errorcode4 = call("rm mykrobe_output.json",shell=True)

def CollType():
    print("Typing according to Coll (2014) scheme")
    if os.path.isfile("colltype.txt"):
        print("Coll type already calculated")
        return 0
    errorcode = call("colltyper -o colltype.txt snippy/snps.vcf", shell=True)

def sampleAnalysis(sample):
    
    sampleName = FindSampleName(sample)
    print("Sample name: %s" % sampleName)
    try:
        os.chdir(sample)
    except:
        sys.exit("Could not access dir: %s" % sample)
    print("Current dir: %s " % os.getcwd())
    myfiles = FindReads()
    try:
        assert ' ' not in myfiles["R1"]
        assert ' ' not in myfiles["R2"]
    except AssertionError:
        sys.exit("Could not proceed because some file names contain space. Check %s and %s" % (myfiles["R1"], myfiles["R2"]))

    RunFastQC(myfiles["R1"], myfiles["R2"])
    RunKaiju(myfiles["R1"], myfiles["R2"])
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
        if not os.path.isfile("%s/%s/%s" % (GLOBAL_COLLECTION, sample, fil)):
            call("cp -rf %s/%s %s/%s/%s" % (sample, fil, GLOBAL_COLLECTION, sample, fil), shell=True)
    if not os.path.isfile("%s/%s/snps.tab" % (GLOBAL_COLLECTION, sample)):
        call("cp %s/snippy/snps.tab %s/%s/snps.tab" % (sample, GLOBAL_COLLECTION, sample), shell=True)
    if not os.path.isfile("%s/%s/snps.aligned.fa" % (GLOBAL_COLLECTION, sample)):
        call("cp %s/snippy/snps.aligned.fa %s/%s/snps.aligned.fa" % (sample, GLOBAL_COLLECTION, sample), shell=True)
    # call("ln -s %s %s/%s/" % (SNIPPY_REF_DIR, GLOBAL_COLLECTION, sample), shell=True) # NOT NEEDED FOR DOCKER VERSION

def CopySnippyDataToShallowDir(sample):
    errorcode1 = call("mv %s/snippy/snps.tab %s/snps.tab" % (sample, sample), shell=True)
    errorcode2 = call("mv %s/snippy/snps.aligned.fa %s/snps.aligned.fa" % (sample, sample), shell=True)

def MaskRepetitiveRegions(alnfile):
    '''This method is not yet complete'''
    outfilename = alnfile[:-8] + 'masked.fasta'
    excludefile = TB_EXCLUDECOLS
    with open(excludefile, 'rU') as exfile:
        exclude = exfile.read()
    errorcode = call("trimal -in %s -selectcols %s -out %s" % (alnfile, exclude, outfilename), shell=True)
    return outfilename

def replaceOldSNPalignment(maskedfile, filetoreplace):
    call("snp-sites -o %s %s" % (filetoreplace, maskedfile), shell=True)


def RunSnippyCore(basedir, timestamp):
    # Change so that this analysis is not run in GLOBAL dir but local
    #try:
    #    os.chdir(GLOBAL_COLLECTION)
    #except:
    #    sys.exit("Unable to move to global collection dir: %s" % GLOBAL_COLLECTION)
    if not os.path.isfile("snippy-core.log"):
        errorcode1 = call("snippy-core --prefix=TB_all_%s --ref=/mnt/Reference/ref.fa %s/* 2> snippy-core.log" % (timestamp, GLOBAL_COLLECTION), shell=True)
    
    # Mask bad regions in full alignment
    maskedfile = MaskRepetitiveRegions("TB_all_%s.full.aln" % timestamp)
    replaceOldSNPalignment(maskedfile, "TB_all_%s.aln" % timestamp)

    # Discard whole-genome alignment
    if os.path.isfile("TB_all_%s.full.aln" % timestamp):
        errorcode2 = call("rm TB_all_%s.full.aln" % timestamp, shell=True)

    # Create continuously evolving global tree
    if not os.path.isfile("Global_collection_tree.nwk"):
        errorcode3 = call("FastTree -nt -gtr TB_all_%s.aln > Global_collection_tree.nwk" % (timestamp), shell=True)
    
    # Remove masked file
    if os.path.isfile("TB_all_%s.masked.fasta" % timestamp):
        errorcode4 = call("rm TB_all_%s.masked.fasta", shell=True)
    #errorcode4 = call("mv TB_all* %s" % basedir, shell=True)
    #errorcode5 = call("mv snippy-core.log %s" % basedir, shell=True)
    #try:
    #    os.chdir(basedir)
    #except:
    #    sys.exit("Unable to move back from global collection to %s " % basedir)

def RunSnpDists():
    print("Finding distances between all isolates in global collection")
    errorcode = call("snp-dists -c TB_all*.aln > TB_all_dists.csv", shell=True)

def ReadSnpDistsObject(dists):
    header = next(dists)[1:]
    res = {}
    for row in dists:
        res[row[0]] = {}
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
        data = csv.reader(infile, delimiter="\t")
        header = data.next()
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
        data = csv.reader(infile,delimiter="\t")
        header = data.next()
        lineagecol = header.index('lineage')
        return data.next()[lineagecol]

def GetLineageColl(colltyperfile):
    with open(colltyperfile, "rU") as infile:
        data = csv.reader(infile,delimiter="\t")
        header = data.next()
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
        if int(dists[key]) <= 30:
            if int(dists[key]) > 5:
                related["somewhat"] += 1
            else:
                related["close"] += 1
                related["somewhat"] += 1
    return related



def FinalizeSampleReport(sample, metainfo, resdic, clusterjanei, lineage, relationtoothers, covdicsample):
    try:
        os.chdir("./%s" % (sample))
    except:
        sys.exit("Could not move into directory %s" % sample)
        
    CopyTexTemplate()

    # Copy tree file to tex directory
    call("cp NeighborTreeWhite.png Latex_template/imageFiles/tree.png",shell=True)

    # Run TB-finalizer to create tex files
    Tex_finalizer.CreateReport(metainfo, resdic, clusterjanei, lineage, relationtoothers, covdicsample)
    # Create the pdf
    try:
        os.chdir("Latex_template")
    except:
        sys.exit("Unable to move in tex directory of %s" % sample)
    call("pdflatex tb-wgs-report.tex", shell=True)
    call("cp tb-wgs-report.pdf ../%s-tb-wgs-report.pdf" % sample, shell=True)
    call("cp tb-wgs-report.pdf %s/%s/%s-tb-wgs-report.pdf" % (GLOBAL_COLLECTION, sample, sample), shell=True)
    try:
        os.chdir("../..")
    except:
        sys.exit("Unable to move back from directory %s" % sample)

def main():
    print("Starting")
    # Verify that Metadata.csv exists
    metadataexists = True
    if not os.path.isfile("Metadata.csv"):
        metadataexists = False
        print("WARNING: Unable to locate Metadata.csv in root directory")

    basedir = os.getcwd()
    timestamp = time.strftime("%d_%b_%Y")

    dirs = next(os.walk(basedir))[1]

    for sample in dirs:
        print("Running sample %s" % sample)

        sampleAnalysis(sample)

        # When done with individual analyses, copy all results to Global dir
        CopyToGlobalDir(sample)
        CopySnippyDataToShallowDir(sample)
    
    RunSnippyCore(basedir, timestamp)
    RunSnpDists()
    # For each sample, find the 20 closest and draw a FastTree
    with open("TB_all_dists.csv","rbU") as distfile:
        print("Reading distances between all isolates in global collection")
        dists = csv.reader(distfile, delimiter=",")
        usedists = ReadSnpDistsObject(dists)

    # Load Biopython alignment of all snps
    with open("TB_all_%s.aln" % timestamp, "rbU") as snpfile:
        print("Loading global alignment of SNPs")
        all_snps = Bio.AlignIO.read(snpfile, "fasta")

    
    # Load sample metadata
    if metadataexists:
        with open("Metadata.csv","rbU") as metainfofile:
            print("Loading metadata")
            metainfo = csv.reader(metainfofile,delimiter=",")
            metainfodic = {}
            header = metainfo.next()
            for row in metainfo:
                rowdic = {header[i]: row[i] for i in range(len(row))}
                metainfodic[row[0]] = rowdic
    else:
        metainfodic = {}
        for sample in dirs:
            metainfodic[sample] = {"ID": sample, "Barcode":"","Location":"","Source":"","Isolated":""}

    with open("snippy-core.log","rbU") as snippycorelogfile:
        print("Reading snippy-core log")
        snippyloglines = snippycorelogfile.readlines()
        covdic = {}
        for l in snippyloglines:
            if "coverage" not in l:
                continue
            else:
                l2 = l.split("\t")[1] # Get information part
                l2s = l2.split(" ")
                covdic[l2s[0]] = float(l2s[4][:-2]) # Remove % and \n

    print("Finalizing reports for each strain")

    for sample in dirs:
        Neighbors = FindNeighbors(usedists[sample], NUMBER_NEIGHBORS_IN_TREE)
        MakeTree(sample, Neighbors, all_snps)

        # Pimp resdic
        pimpedresdic = PimpResDic('%s/mykrobe_output.tsv' % sample)

        # Get lineage
        # Consider implementing COLL scheme instead
        #lin = GetLineage('%s/mykrobe_output.tsv' % sample)
        lin = GetLineageColl("%s/colltype.txt" % sample)

        # Find out if sample is part of a cluster.
        # Lowest distance is always to self (0)
        relationtoothers = NumberRelated(sample, usedists[sample])
        clusterjanei = relationtoothers["somewhat"] > 0

        FinalizeSampleReport(sample, metainfodic[sample], pimpedresdic, clusterjanei, lin, relationtoothers, covdic[sample])
        
    print("Finished")

if __name__ == '__main__':
    main()