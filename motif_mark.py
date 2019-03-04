#!/usr/bin/env python
'''Motif mapping module for BI 625 assignment'''

import argparse
import re
import cairo

def get_arguments():
    parser = argparse.ArgumentParser(
        description="reads in fasta file and splicing regulation motifs, to visualize mapping")
    parser.add_argument("-f", "--filename", help="name of fasta file",
                        required=True, type=str)
    parser.add_argument("-m", "--motifs", help="file containing motifs, 1 per line",
                        required=True, type=str)
    return parser.parse_args()

def build_exon(line, base):
    '''This function is designed to build an exon dictionary entry
    from an uppercase segment in a fasta file.'''
    exon = ''
    start = base
    while line[base].isupper():
        exon += line[base]
        base += 1
    fin = base - 1
    return base, start, exon, fin

def build_intron(line, base):
    '''This function is designed to build an intron dictionary entry
    from an uppercase segment in a fasta file.'''
    intron = ''
    start = base
    while line[base].islower():
        intron += line[base]
        base += 1
    fin = base - 1
    return base, start, intron, fin

def parse_fa(fa_file):
    '''(file) -> dict, dict, dict
    This function is going to parse a fasta file where introns appear
    in lowercase and exons appear in uppercase. It will return an
    intron and an exon dictionary where the key represents the order
    of the segment.'''
    ln, base = 0, 0 # init counters
    exons={} #key will be start coordinate, value will be seq, end coordinate
    introns={} #key will be start coordinate, value will be seq, end coordinate
    with open(fa_file) as fa:
        for line in fa:
            ln+=1
            line = line.strip('\n')
            if ln%2 == 1: # header
                header = line.split( )
                gene = header[0]
                loc = header[1].split(:)
                stat = header[2] + header[3] # (reverse complement)
                chr = loc[0]
                pos = loc[1].split(-)
                base = pos[0]
            if ln%2 == 0: # sequence line
                while base < len(line):
                    if line[base].isupper():
                        base, start, exon, fin = build_exon(line, base)
                        exons[start] = [exon, fin]
                    if line[base].islower():
                        base, start, intron, fin = build_intron(line, base)
                        introns[start] = [intron, fin]
    return segments, introns, exons

def id_motif(m_file, introns, exons):
    '''(file, dict, dict) -> dict
    This function takes in a sequence as a string and uses regex to find motifs.
    Motifs are identified and returned with coordinates (start position?).'''
    motif_dict = {} #define acceptable versions of motifs
    motif_coords = {} #dict for motif mapping
    # parse motif file
    with open(m_file) as motifs:
        for motif in motifs:
            motif_dict[motif] = motif
            for base in range(len(motif)): # Populate motif variants of 1 Y or N
                motif_dict[motif[0:base]+'W'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'S'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'M'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'K'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'R'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'Y'+motif[base+1:len(motif)]] = motif
                motif_dict[motif[0:base]+'N'+motif[base+1:len(motif)]] = motif
    # find motifs in introns
    for intron in introns:
        for key in motif_dict:
            # use re.findall to identify motifs
            if motif_dict[key] not in motif_coords:
                motif_coords[motif_dict[key]]= []
            # motif_coords[motif_dict[key]].append(start_coordinates)
    # find motifs in exons
    for exon in exons:
        for key in motif_dict:
            # use re.findall to identify motifs
            if motif_dict[key] not in motif_coords:
                motif_coords[motif_dict[key]]= []
            # motif_coords[motif_dict[key]].append(start_coordinates)
    return motif_coords

def draw_motifs(s_dict, m_dict, i_dict, e_dict):
    '''(dict,dict,dict) -> svg
    This function uses dictionaries generated from parse_fa() and id_motif. Dictionaries
    are: s_dict (segment dictionary, key = order found, value = true start position),
    m_dict (motif dictionary, key = motif, value = list of start positions),
    i_dict (intron dictionary, key = true start pos, value = sequence, end pos),
    e_dict (exon dictionary, key = true start pos, value = sequence, end pos).
    Function generates an SVG image of the gene including introns,
    exons and motif mapping, using pycairo to draw.'''
    #add code
    return None

def main():
    '''documentation'''
    seg_dict, intron_dict, exon_dict = parse_fa(args.filename)
    motif_coords = id_motif(args.motifs, intron_dict, exon_dict)
    draw_motifs(seg_dict, motif_coords, intron_dict, exon_dict)
    return None
