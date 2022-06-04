#!/usr/bin/env python3
# Group Members: Ananthajit ("Ananth") Srikanth (asrikan1), Layla Myers (laamyers), Teresa Joseph (tkjoseph)

# Required modules: Bio, sys, pdfplumber, os, statistics, datetime


from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import Align
import sys
import pdfplumber
import os
import statistics
from datetime import datetime



class PdfWriter:
    '''
    pdfWriter: write a PDF to text file
    initialize: infile (does not use stdin)
    
    filterToText: convert file input to text file
    '''
    def __init__(self, inFile = ""):
        '''Initialize pdf file'''
        self.file = inFile

    def filterToText(self,outputFile, temp, linesToIgnore = 5):
        '''Convert pdf file to text, ignore number of lines in 
        linesToIgnore'''
        
        # write temporary file
        with pdfplumber.open(self.file) as pdf:
            with open(temp,'w') as tempFile:
                
                for page in range(len(pdf.pages)):
                    
                    tempFile.write(pdf.pages[page].extract_text())
        # write text file without extra spaces

            
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
        
        # remove temporary file
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





class FungalGenome(FastAreader):
    '''
    FungalGenome: read in an entire fungal genome in fastA format and compare 
    to another genome
    Initialization: primer PDF file, number of lines to ignore in PDF (see PdfWriter)
    and genome fasta name
    
    Methods: primersGen - generate primers
    
    '''
    
    
    # ITS1 forward and reverse sequences, mostly from Illumina
    ITS1sForward = ['TCCGTAGGTGAACCTGCGG','CTTGGTCATTTAGAGGAAGTAA', 'CTCGGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGAACTAA',
                'CCCGGTCATTTAGAGGAAGTAA', 'CTAGGCTATTTAGAGGAAGTAA', 'CTTAGTTATTTAGAGGAAGTAA',
                'CTACGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGTCGTAA','GCATCGATGAAGAACGC', 'ATCGATGAAGAACGCAG']
    ITS1sReverse = [ 'GCTGCGTTCTTCATCGATGC','GCTGCGTTCTTCATCGATGG', 'GCTACGTTCTTCATCGATGC',
                     'GCTGCGTTCTTCATCGATGT', 'ACTGTGTTCTTCATCGATGT', 'GCTGCGTTCTTCATCGTTGC',
                     'GCGTTCTTCATCGATGC', 'GGGCGCAATGTGCGTTCAAA', 'AAACTCTGTCGTGCTGGGGATA','GATTGAATGGCTTAGTGAGG'
                     ]
    
    #overhangs = ['TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG']
    
    def __init__(self,  primerFile = "primers.pdf", linesToIgnore = 5, fname = ""):
        '''Initialize fasta, and generate primer text file'''
        super().__init__(fname = fname) # initialize using fastAreader init
        
        self.primersGen(primerFile, linesToIgnore) # generate primers
        
     
    # ITS2s = [] 
    
    def primersGen(self, primerFile, linesToIgnore = 5):
        '''Generate primers from primer PDF file'''
        
        primerPDF = PdfWriter(primerFile)
        
        primerPDF.filterToText("primers.txt", "temp.txt", linesToIgnore) # removes first 5 lines as there are no primers
        
        with open("primers.txt") as primers:
            while True:
                line = primers.readline()
                if line == "":
                    break
                else:
                    lineList = line.strip().split()
                    if len(lineList) == 6:
                        # based on text file, forward primer is
                        # third item when split
                        # and reverse is sixth
                        self.ITS1sForward.append(lineList[2])
                        self.ITS1sReverse.append(lineList[5])
                    else:
                        # there are more forward primers than reverse primers
                        self.ITS1sForward.append(lineList[2])
                        
        
    
    def returnChromosome(self, chromosome):
        '''Return nth chromosome (as chromosomes are sorted in order in genomic fasta)'''
        
        count = 0
        for header, sequence in self.readFasta():
            count += 1
            if count == chromosome:
                storeSeq = sequence
            
        return storeSeq.upper()
            
    
    def findITS(self):
        '''Find all ITS positions that each primer can bind to'''
        chromsList = []
        chrom = 0
        for header, sequence in self.readFasta():
            
            chrom += 1 # chromosome number
            # check forward sequences
            for primer in self.ITS1sForward:
                # reverse complement of primer binds to ITS1 region
                checkSeq = str(Seq(primer).reverse_complement())
                
            
                
                sequence = sequence.upper()
                
                #nt_search searches for indices of checkSeq in sequence
                ITSPos = SeqUtils.nt_search(sequence, checkSeq)
                
                ITSPos = ITSPos[1:]
                if ITSPos: # check if ITSPos is not empty
                    chromsList.append((chrom, ITSPos))
                       
                        
            # check reverse sequence            
            for primer in self.ITS1sReverse:
                
                checkSeq = str(Seq(primer).reverse_complement())
                sequence = sequence.upper()
                newSeq = str(Seq(sequence).reverse_complement())
                
                    
                
                
                #nt_search searches for indices of checkSeq in newSeq
                temp = SeqUtils.nt_search(newSeq, checkSeq)
                temp = ITSPos[1:]
                ITSPos = []
                # reverse position is length of sequence - position on reverse 
                # complement
                if temp: # check if temp is emptu
                    for ITS in temp:
                        ITSPos.append(len(sequence) - ITS)
                        
                    # ITS is approximately 200 bases long, so find approximate
                    # start of ITS
                    chromsList.append((-chrom, ITSPos - 200))    
                
                    
               
                
        
        return chromsList
    
    
    def compareGenome(self, other, comparisonLength = 1000):
        '''Compare self (expected to be reference genome)
            to another genome to detect if they are the same species
            comparisonLength estimates the ITS1-5.8S-ITS2 region to be 
            about 1000 bases long'''
        chromsList = self.findITS() # return chromosomes/position
         
        
        if not chromsList: # no primers bound
            print("Error: No Primers Bound to genome.")
            return
        
        tempChroms = []
        
        for potentialITS in chromsList:
            tempChroms.append(abs(potentialITS[0]))
            
        chromToUse = statistics.mode(tempChroms) # find chromosome on which
        # its exists, mode is used in case a primer binds to another chromosome
        # so it checks the modal frequency of the chromosome 
        
        
        print("ITS Chromosome: ", chromToUse )
        tempITSIndices = []
        
        # take the first index of each primer position and average it
        # (some primers bind to more than one region)
        for potentialITS in chromsList:
            if abs(potentialITS[0]) == chromToUse:
                tempITSIndices.append(potentialITS[1][0])
        # round to get integer index                
        ITSIndex = round(statistics.mean(tempITSIndices))
        
        
        # store a barcode sequence (ITS1 - 5.8S - ITS2 region)
        barcode = self.returnChromosome(chromToUse)
        
        barcode = barcode[ITSIndex:(ITSIndex + comparisonLength)]
        # take about 1000 base long region
        
        try:
            
            # pairwise alignment algorithm, penalize gaps used to align sequences
            
            otherChrom = other.returnChromosome(chromToUse)
            aligner = Align.PairwiseAligner()
            aligner.mode = 'local'
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -0.5
            # aligner is case sensitive
            alignments = aligner.align(otherChrom.upper(), barcode.upper())
           
            match = alignments.score/comparisonLength
            
            match2 = match*100
            print(f'Overall match = {match2:.1f}%')
            
            if match > 0.97:
                print("Species are the same")
            else:
                print("Species are different")
            
            return 
        except UnboundLocalError:
            # chromosome indexed out of bounds
            print("Error: Chromosome out of Bounds or does not exist, genomes cannot be compared")
            print("Species unlikely to be the same")
            return
        
        
def main(referenceGenome, genomeToCheck, compLength = 1000):
    '''main: compares a genome to a reference genome to check if the two
    are the same
    input: reference and genome to check, comparison length 
    output: chromosome on which the ITS sequence exists
        % match of region between genomes
        whether the species are the same
        time taken to generate comparison
        '''
    # measure times to show computational efficiency
    time1 = datetime.now()
    
    # path of this script
    filepath = os.path.dirname(os.path.realpath(__name__))
    
    # path to genome folder (assumes there is a genomes folder in the same
    # location of script, containing genomic fasta files)
    filepathRef = filepath + os.sep + 'genomes' + os.sep + referenceGenome 
    
    filepathCheck = filepath + os.sep + 'genomes' + os.sep + genomeToCheck 
    
    check = FungalGenome(fname = filepathCheck)
    reference = FungalGenome(fname = filepathRef)
    
    # genome comparison
    reference.compareGenome(check, compLength)
    
    
    time2 = datetime.now()
    time3 = time2 - time1
    # calculate minutes and seconds and print time taken to run genome comparison
    time3Mins = time3.seconds // 60
    time3Secs = time3.seconds % 60
    
    
    print(f"Time taken: {time3Mins} minutes and {time3Secs} seconds")
    
    
if __name__ == "__main__":
    
      # cases tested
      print("Saccharomyces cerevisiae (2014 Reference) vs Saccharomyces cerevisiae (c.1999)")
      main("SaccharomycesCerevisiae2014.fna", "SaccharomycesCerevisiae1999.fsa")
     
      print()
     
      print("Saccharomyces cerevisiae (2014 Reference) vs Saccharomyces paradoxus")
      main("SaccharomycesCerevisiae2014.fna", "SaccharomycesParadoxus.fna")
     
      print()
     
      print("Saccharomyces cerevisiae (2014 Reference) vs Agaricus bisporus")
      main("SaccharomycesCerevisiae2014.fna", "AgarBisporus.fna")
     
      print()
     
      print("Saccharomyces cerevisiae (2014 Reference) vs Saccharomyces cerevisiae (c.1999)")
     
      # note to print
      print("Note: the 2011 Genome also failed to bind to ITS primers")
      main("SaccharomycesCerevisiae2014.fna", "SaccharomycesCerevisiae2011.fna")
     
      print()
     
      print("Aspergillus Fumigatus vs Aspergillus Flavus")
      main("AsperFumigatus.fna", "AsperFlavus.fna")
     
    