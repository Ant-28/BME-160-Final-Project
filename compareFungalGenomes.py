#!/usr/bin/env python3
# Group Members: Ananthajit ("Ananth") Srikanth (asrikan1), Layla Myers (laamyers), Teresa Joseph (tkjoseph)

# Required modules: Bio

import Bio
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqUtils

# TODO: Download two fungal genomes from NCBI
import sys


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


class genomeDownloader:
    '''Add genome from NCBI'''
    def __init__(self, genome):
        self.genome = genome
    def genomeDownloader(self, filename):
        with open(filename, "w") as genomeFile:
            pass

class fungalGenome(FastAreader):
    
    
    ITS1s_forward = ['CTTGGTCATTTAGAGGAAGTAA', 'CTCGGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGAACTAA',
                'CCCGGTCATTTAGAGGAAGTAA', 'CTAGGCTATTTAGAGGAAGTAA', 'CTTAGTTATTTAGAGGAAGTAA',
                'CTACGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGTCGTAA', 'TCCGTAGGTGAACCTGCGG', 'CTTGGTCATTTAGAGGAAGTAA'
                ,'GGAAGTAAAAGTCGTAACAAGG', 'TACGTCCCTGCCCTTTGTAC', 'GTCCCTGCCCTTTGTACACA', 'CTGCCCTTTGTACACACCGC', 'ACACACCGCCCGTCGCTACT']
    ITS1s_reverse = [ 'GCTGCGTTCTTCATCGATGC','GCTGCGTTCTTCATCGATGG', 'GCTACGTTCTTCATCGATGC',
                     'GCTGCGTTCTTCATCGATGT', 'ACTGTGTTCTTCATCGATGT', 'GCTGCGTTCTTCATCGTTGC',
                     'GCGTTCTTCATCGATGC', 'GGGCGCAATGTGCGTTCAAA']
    ITS2s = [] # wat
    overhangs = ['TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG']
    
    def returnChromosome(self, chromosome):
        count = 0
        for header, sequence in self.readFasta():
            count += 1
            if count == chromosome:
                storeSeq = sequence
            
        return storeSeq
            
    def sequenceMatch(self, inSeq, otherSeq):
        match = 0 
        assert len(inSeq) == len(otherSeq)
        for ind in range(len(inSeq)):
            
            if inSeq[ind] == otherSeq[ind]:
                match += 1
            
        return match/len(inSeq)
    
    def findChromosome(self):
        chromsList = []
        for primer in self.ITS1s_forward:
            checkSeq = Seq(primer).reverse_complement()
            chrom = 0
        
            for header, sequence in self.readFasta():
                sequence = sequence.upper()
                chrom += 1
                if str(checkSeq) in sequence:
                    chromsList.append((chrom, sequence.index(str(checkSeq))))
        for primer in self.ITS1s_reverse:
            checkSeq = Seq(primer).reverse_complement()
            chrom = 0
            newSeq = str(Seq(sequence).reverse_complement())
            for header, sequence in self.readFasta():
                
                sequence = sequence.upper()
                chrom += 1
                if str(checkSeq) in newSeq:
                    chromsList.append((chrom, len(newSeq) - newSeq.index(str(checkSeq))))        
                
        return chromsList
    
my_file = fungalGenome('GCA_022626425.1_ASM2262642v1_genomic.fna')

print(my_file.findChromosome())

