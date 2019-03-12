#!/usr/bin/env python
'''Motif mapping module for BI 625 assignment'''

import argparse
import re
import numpy as np
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
    gene = header[0]
    gene = gene[1:]
    loc = header[1].split(':')
    if len(header) == 4:
        stat = header[2] + header[3] # (reverse complement)
    else:
        stat = ''
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
    return base, start, exon, fin, chr

def build_intron(line, base, chr):
    '''This function is designed to build an intron dictionary entry
    from an uppercase segment in a fasta file.'''
    intron = ''
    start = base
    while chr < len(line) and line[chr].islower():
        intron += line[chr]
        base += 1
        chr += 1
    fin = base - 1
    intron = intron.upper()
    return base, start, intron, fin, chr

def parse_fa(fa_file):
    '''(file) -> dict, dict, dict
    This function is going to parse a fasta file where introns appear
    in lowercase and exons appear in uppercase. It will return an
    intron and an exon dictionary where the key represents the order
    of the segment.'''
    ln, base = 0, 0 # init counters
    new = True
    exons={} #key will be start coordinate, value will be seq, end coordinate
    introns={} #key will be start coordinate, value will be seq, end coordinate
    genes = {} #key will be gene abbr, value will be chr, start coor, end coor
    with open(fa_file) as fa:
        for line in fa:
            ln+=1
            line = line.strip('\n')
            if line.startswith('>'): # header
                if new == False:
                    new = True
                    char = 0
                    while char < len(prev_line):
                        if prev_line[char].isupper():
                            base, start, exon, fin, char = build_exon(prev_line, base, char)
                            exons[start] = [exon, fin]
                        elif prev_line[char].islower():
                            base, start, intron, fin, char = build_intron(prev_line, base, char)
                            introns[start] = [intron, fin]
                gene, chrm, base, end, stat = parse_header(line)
                if gene not in genes:
                    genes[gene] = [chrm, base, end]
                prev_line = ''
            else: # sequence line
                prev_line += line
                new = False
    char = 0
    while char < len(prev_line):
        if prev_line[char].isupper():
            base, start, exon, fin, char = build_exon(prev_line, base, char)
            exons[start] = [exon, fin]
        elif prev_line[char].islower():
            base, start, intron, fin, char = build_intron(prev_line, base, char)
            introns[start] = [intron, fin]
    return introns, exons, genes

def make_motif_pattern(motif):
    '''(string) -> string
    This function reads in a motif and builds a regex expression that will
    recognize all equivalent motifs.'''
    motif_pattern = ''
    trans_dict = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "U": "T",
        "M": "[AC]",
        "R": "[AG]",
        "W": "[AT]",
        "S": "[CG]",
        "Y": "[CT]",
        "K": "[GT]",
        "V": "[ACG]",
        "H": "[ACT]",
        "D": "[AGT]",
        "B": "[CGT]",
        "N": "[GATC]"}
    for char in motif:
        motif_pattern += trans_dict[char]
    return motif_pattern

def id_motif(m_file, introns, exons):
    '''(file, dict, dict) -> dict
    This function takes in a sequence as a string and uses regex to find motifs.
    Motifs are identified and returned with coordinates (start position?).'''
    motif_dict={}
    motif_coords = {} #dict for motif mapping
    # parse motif file
    with open(m_file) as motifs:
        for motif in motifs:
            motif = motif.strip('\n')
            motif_dict[motif] = make_motif_pattern(motif.upper())
    for motif in motif_dict:
        pattern = re.compile(motif_dict[motif])
        if motif not in motif_coords:
                motif_coords[motif]= []
        for intron in introns:
            coords = []
            for match in re.finditer(motif_dict[motif], introns[intron][0]):
                real_c = tuple(np.add(match.span(), (intron,intron)))
                coords.append(real_c) # adjusted coordinates of each match as np array
            for coord in coords:
                motif_coords[motif].append(coord)
        for exon in exons:
            coords = []
            for match in re.finditer(motif_dict[motif], exons[exon][0]):
                real_c = tuple(np.add(match.span(), (intron,intron)))
                coords.append(real_c) #coordinates of each match in tuple form
            for coord in coords:
                motif_coords[motif].append(coord)
    return motif_coords

def id_exons(exon_dict, gene_info):
    '''(dict, list) -> list
    This function is meant to be used by draw_motifs() to id exons for drawing.
    exon_dict (key = true start pos, value = sequence, end pos)
    gene_info [chromosome (str), start pos (int), end pos (int)]'''
    exons = []
    start = gene_info[1]
    end = gene_info[2]
    chrom = gene_info[0]
    for entry in exon_dict:
        if entry > start and entry < end:
            exons.append((entry,exon_dict[entry][1])) #append tuple of exon coordinates to list
    return exons

def pick_motifs(m_list, gene_info):
    '''This function is meant to be used by draw_motifs() to id motifs
    for drawing on a particular gene.
    m_list (list of motif coordinates in tuple form)
    gene_info [chromosome (str), start pos (int), end pos (int)]'''
    motifs = []
    start = gene_info[1]
    end = gene_info[2]
    chrom = gene_info[2]
    for coordinates in m_list:
        if coordinates[0] > start and coordinates[0] < end:
            motifs.append(coordinates)
    return motifs

def draw_motifs(m_dict, i_dict, e_dict, g_dict):
    '''(dict,dict,dict,dict) -> svg
    This function uses dictionaries generated from parse_fa() and id_motif. Dictionaries
    are: m_dict (motif dictionary, key = motif, value = list of start positions),
    i_dict (intron dictionary, key = true start pos, value = sequence, end pos),
    e_dict (exon dictionary, key = true start pos, value = sequence, end pos),
    g_dict (gene dictionary, key = gene name, value = chromosome, start pos, end pos).
    Function generates an SVG image of the gene including introns,
    exons and motif mapping, using pycairo to draw.'''
    g, m = 0, 0 # init gene & motif counter
    w = 1500 # set pixel width
    h = 1000 * len(g_dict) # pixel height = frame x each gene to draw
    colors = [[0.9,0.1,0.1], #red
              [0.4, 0.9, 0.4], #green
              [0.2, 0.23, 0.9], #blue
              [0.7, 0.7, 0.2], #yellow
              [0.2,0.7,0.5], #teal
              [0.57, 0.2, 0.7], #purple
              [0.7, 0.45, 0.2], #orange
              [0.1, 0.9, 0.1], #vv green
              [0.2, 0.7, 0.7] #light blue
             ]    # Up to 9 motif colors available, will repeat colors after 9
    for motif in m_dict:
        m_dict[motif].append(colors[m%9])
        m+=1
    with cairo.SVGSurface('motif_map.svg', w, h) as surface:
        context = cairo.Context(surface)
        for gene in g_dict:
            scale = 1300 / (g_dict[gene][2] - g_dict[gene][1]) #scale gene drawing to fit frame
            #draw gene outline
            context.set_source_rgb(0,0,0) #drawing a black line
            context.set_line_width(10)
            context.move_to(100, 500 + 1000*g) # start position
            context.line_to(1400, 500 + 1000*g) # draw line to end pos
            context.stroke()
            #draw exons
            exon_list = id_exons(e_dict, g_dict[gene])
            for exon in exon_list:
                start = exon[0] - g_dict[gene][1]
                end = exon[1] - exon[0]
                context.rectangle((100+start*scale), (450+1000*g), ((end)*scale), 100)
                context.fill()
            #draw motifs
            for motif in m_dict:
                context.set_source_rgb(m_dict[motif][-1][0],m_dict[motif][-1][1],m_dict[motif][-1][2])
                motif_list = pick_motifs(m_dict[motif][:-1], g_dict[gene])
                for site in motif_list:
                    start = site[0] - g_dict[gene][1]
                    end = site[1] - site[0]
                    context.rectangle((100+start*scale), (450+1000*g), (end*scale), 100)
                    context.fill()
            #final things
            g += 1 # update gene count
    return None

def main():
    '''documentation'''
    args = get_arguments()
    intron_dict, exon_dict, gene_dict = parse_fa(args.filename)
    motif_coords = id_motif(args.motifs, intron_dict, exon_dict)
    draw_motifs(motif_coords, intron_dict, exon_dict, gene_dict)
    print('made .svg drawing')
    return None


if __name__ == '__main__':
    main()
