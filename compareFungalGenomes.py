#!/usr/bin/env python3
# Group Members: Ananthajit ("Ananth") Srikanth (asrikan1), Layla Myers (laamyers), Teresa Joseph (tkjoseph)

# Required modules: Bio


from Bio import Entrez

# TODO: Download two fungal genomes from NCBI

class fungalGenome:
    ITS1s = ['CTTGGTCATTTAGAGGAAGTAA', 'CTCGGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGAACTAA',
                'CCCGGTCATTTAGAGGAAGTAA', 'CTAGGCTATTTAGAGGAAGTAA', 'CTTAGTTATTTAGAGGAAGTAA',
                'CTACGTCATTTAGAGGAAGTAA', 'CTTGGTCATTTAGAGGTCGTAA', 'GCTGCGTTCTTCATCGATGC',
                'GCTGCGTTCTTCATCGATGG', 'GCTACGTTCTTCATCGATGC', 'GCTGCGTTCTTCATCGATGT',
                'ACTGTGTTCTTCATCGATGT', 'GCTGCGTTCTTCATCGTTGC', 'GCGTTCTTCATCGATGC']
    ITS2s = []
    overhangs = [‘TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG’, ‘GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG’]
    pass #hi there :)



