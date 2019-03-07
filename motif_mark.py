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

def parse_header(line):
    '''Parses fasta file header to obtain read info'''
    header = line.split(' ')
    print(header)
    gene = header[0]
    gene = gene[1:]
    loc = header[1].split(':')
    stat = header[2] + header[3] # (reverse complement)
    chrm = loc[0]
    pos = loc[1].split('-')
    base = int(pos[0])
    end = int(pos[1])
    return gene, chrm, base, end, stat

def build_exon(line, base, chr):
    '''This function is designed to build an exon dictionary entry
    from an uppercase segment in a fasta file.'''
    exon = ''
    start = base
    while line[chr].isupper():
        exon += line[chr]
        base += 1
        chr += 1
    fin = base - 1
    return base, start, exon, fin

def build_intron(line, base, chr):
    '''This function is designed to build an intron dictionary entry
    from an uppercase segment in a fasta file.'''
    intron = ''
    start = base
    while line[chr].islower(): #broken, index out of range?
        intron += line[chr]
        base += 1
        chr += 1
    fin = base - 1
    intron = intron.upper()
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
    genes = {} #key will be gene abbr, value will be chr, start coor, end coor
    with open(fa_file) as fa:
        for line in fa:
            ln+=1
            chr = 0
            print(ln)
            line = line.strip('\n')
            if line.startswith('>'): # header
                gene, chrm, base, end, stat = parse_header(line)
                #have reverse complement info under stat if wanted
                if gene not in genes:
                    genes[gene] = [chr, base, end]
                new = True
            else: # sequence line
                while chr < len(line):
                    if line[chr].isupper():
                        base, start, exon, fin, chr = build_exon(line, base, chr)
                        exons[start] = [exon, fin]
                    elif line[chr].islower():
                        base, start, intron, fin, chr = build_intron(line, base, chr)
                        introns[start] = [intron, fin]
                new = False
    return introns, exons

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

def draw_motifs(m_dict, i_dict, e_dict):
    '''(dict,dict,dict) -> svg
    This function uses dictionaries generated from parse_fa() and id_motif. Dictionaries
    are: m_dict (motif dictionary, key = motif, value = list of start positions),
    i_dict (intron dictionary, key = true start pos, value = sequence, end pos),
    e_dict (exon dictionary, key = true start pos, value = sequence, end pos).
    Function generates an SVG image of the gene including introns,
    exons and motif mapping, using pycairo to draw.'''
    #add code
    return None

def main():
    '''documentation'''
    args = get_arguments()
    intron_dict, exon_dict = parse_fa(args.filename)
    print(intron_dict, exon_dict)
    #motif_coords = id_motif(args.motifs, intron_dict, exon_dict)
    #draw_motifs(motif_coords, intron_dict, exon_dict)
    return None


if __name__ == '__main__':
    main()
