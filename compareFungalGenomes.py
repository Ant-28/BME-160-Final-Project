#!/usr/bin/env python3
# Group Members: Ananthajit ("Ananth") Srikanth (asrikan1), Layla Myers (laamyers), Teresa Joseph (tkjoseph)

# Required modules: Bio

import Bio
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import Align

# TODO: Download two fungal genomes from NCBI
# todo: explain comments

import sys
import pdfplumber
import os
import statistics


class pdfWriter:
    
    def __init__(self, inFile = ""):
        self.file = inFile

    def filterToText(self,outputFile, temp, linesToIgnore = 5):
        with pdfplumber.open(self.file) as pdf:
            with open(temp,'w') as tempFile:
                
                for page in range(len(pdf.pages)):
                    
                    tempFile.write(pdf.pages[page].extract_text())

        with open(temp,'r') as tempFile:
            with open(outputFile, 'w') as output:
                for i in range(linesToIgnore):
                    line = tempFile.readline()
                while True:
                    line = tempFile.readline()
                    if line == "":
                        break
                    elif not line.isspace():
                        output.write(line)
        if os.path.exists(temp):
            os.remove(temp)

class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN and defaulting to Unicode.'''
        if self.fname == '':
            # default to unicode for special character support
            sys.stdin.reconfigure(encoding='utf-8')
            sys.stdout.reconfigure(encoding='utf-8')

            return sys.stdin
        else:
            return open(self.fname, encoding='utf-8')

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence



# Complete this asap
class genomeDownloader:
    '''Add genome from NCBI'''
    def __init__(self, genome):
        self.genome = genome
    def genomeDownloader(self, filename):
        with open(filename, "w") as genomeFile:
            pass


class fungalGenome(FastAreader):
    
    ITS1sForward = ['TCCGTAGGTGAACCTGCGG','CTTGGTCATTTAGAGGAAGTAA', 'CTCGGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGAACTAA',
                'CCCGGTCATTTAGAGGAAGTAA', 'CTAGGCTATTTAGAGGAAGTAA', 'CTTAGTTATTTAGAGGAAGTAA',
                'CTACGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGTCGTAA','GCATCGATGAAGAACGC', 'ATCGATGAAGAACGCAG']
    ITS1sReverse = [ 'GCTGCGTTCTTCATCGATGC','GCTGCGTTCTTCATCGATGG', 'GCTACGTTCTTCATCGATGC',
                     'GCTGCGTTCTTCATCGATGT', 'ACTGTGTTCTTCATCGATGT', 'GCTGCGTTCTTCATCGTTGC',
                     'GCGTTCTTCATCGATGC', 'GGGCGCAATGTGCGTTCAAA', 'AAACTCTGTCGTGCTGGGGATA','GATTGAATGGCTTAGTGAGG'
                     ]
    
    #overhangs = ['TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG']
    
    def __init__(self,  primerFile = "primers.pdf", linesToIgnore = 5, fname = ""):
        super().__init__(fname = fname)
        
        self.primersGen(primerFile, linesToIgnore)
        
     
    # ITS2s = [] # wat
    
    def primersGen(self, primerFile, linesToIgnore):
        primerPDF = pdfWriter(primerFile)
    
        primerPDF.filterToText("primers.txt", "temp.txt", linesToIgnore)
        
        with open("primers.txt") as primers:
            while True:
                line = primers.readline()
                if line == "":
                    break
                else:
                    lineList = line.strip().split()
                    if len(lineList) == 6:
                        self.ITS1sForward.append(lineList[2])
                        self.ITS1sReverse.append(lineList[5])
                    else:
                        self.ITS1sForward.append(lineList[2])
                        
        
    
    def returnChromosome(self, chromosome):
        count = 0
        for header, sequence in self.readFasta():
            count += 1
            if count == chromosome:
                storeSeq = sequence
            
        return storeSeq.upper()
            
    def sequenceMatch(self, inSeq, otherSeq):
        match = 0 
        assert len(inSeq) == len(otherSeq)
        for ind in range(len(inSeq)):
            
            if inSeq[ind] == otherSeq[ind]:
                match += 1
            
        return match/len(inSeq)
    
    def findITS(self):
        chromsList = []
        chrom = 0
        for header, sequence in self.readFasta():
            
            chrom += 1
            for primer in self.ITS1sForward:
                checkSeq = Seq(primer).reverse_complement()
                
            
                
                sequence = sequence.upper()
                
                
                ITSPos = SeqUtils.nt_search(sequence, str(checkSeq))
                ITSPos = ITSPos[1:]
                if ITSPos:
                    chromsList.append((chrom, ITSPos))
                        
                        
                        
            for primer in self.ITS1sReverse:
                
                checkSeq = str(Seq(primer).reverse_complement())
                sequence = sequence.upper()
                newSeq = str(Seq(sequence).reverse_complement())
                
                    
                
                
                
                temp = SeqUtils.nt_search(newSeq, checkSeq)
                temp = ITSPos[1:]
                ITSPos = []
                if temp:
                    for ITS in temp:
                        ITSPos.append(len(sequence) - ITS)
                    chromsList.append((-chrom, ITSPos))        
            
                
        
        return chromsList
    
    
    def compareGenome(self, other, comparisonLength = 1000):
        chromsList = self.findITS()
        if not chromsList:
            print("Error: No Primers Bound to genome.")
            return
        
        tempChroms = []
        for potentialITS in chromsList:
            tempChroms.append(potentialITS[0])
            
        chromToUse = statistics.mode(tempChroms)
        print("ITS Chromosome: ", chromToUse )
        tempITSIndices = []
        for potentialITS in chromsList:
            tempITSIndices.append(potentialITS[1][0])
        ITSIndex = round(statistics.mean(tempITSIndices))
        
        barcode = self.returnChromosome(chromToUse)
        
        barcode = barcode[ITSIndex:(ITSIndex + comparisonLength)]
        try:
            otherChrom = other.returnChromosome(chromToUse)
            aligner = Align.PairwiseAligner()
            aligner.mode = 'local'
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(otherChrom, barcode)
            print(f'Overall match = {alignments.score/comparisonLength * 100.:3f}%')
            
            return 
        except UnboundLocalError:
            print("Error: Chromosome out of Bounds, genomes cannot be compared")
            print("Species unlikely to be the same")
            return
        
        
        
my_old_file = fungalGenome(fname = 'flavus.fna')
my_file = fungalGenome(fname = 'asper.fna')

my_file.compareGenome(my_old_file)