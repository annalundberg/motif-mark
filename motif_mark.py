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
                coords.append(match.span()) #coordinates of each match in tuple form
            for coord in coords:
                motif_coords[motif].append(coord)
        for exon in exons:
            coords = []
            for match in re.finditer(motif_dict[motif], exons[exon][0]):
                coords.append(match.span()) #coordinates of each match in tuple form
            for coord in coords:
                motif_coords[motif].append(coord)
    return motif_coords

def draw_motifs(m_dict, i_dict, e_dict, g_dict):
    '''(dict,dict,dict,dict) -> svg
    This function uses dictionaries generated from parse_fa() and id_motif. Dictionaries
    are: m_dict (motif dictionary, key = motif, value = list of start positions),
    i_dict (intron dictionary, key = true start pos, value = sequence, end pos),
    e_dict (exon dictionary, key = true start pos, value = sequence, end pos),
    g_dict (gene dictionary, key = gene name, value = chromosome, start pos, end pos).
    Function generates an SVG image of the gene including introns,
    exons and motif mapping, using pycairo to draw.'''
    g = 0 # init gene counter
    w = 1500 # set pixel width
    h = 1000 * len(g_dict) # pixel height = frame x each gene to draw
    colors = [[0.9,0.1,0.1], #red
              [0.4, 0.9, 0.4], #green
              [0.2, 0.23, 0.9], #blue
              [0.7, 0.7, 0.2], #yellow
              [0.2,0.7,0.5], #teal
              [0.57, 0.2, 0.7], #purple
              [0.7, 0.45, 0.2], #orange
              [0.1, 0.9, 0.1], #very green
              [0.2, 0.7, 0.7] #light blue
             ]    # Up to 9 motif colors available, will repeat colors after 9
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

            g += 1 # update gene count
    return None

def main():
    '''documentation'''
    args = get_arguments()
    intron_dict, exon_dict, gene_dict = parse_fa(args.filename)
    print(intron_dict, '\n', exon_dict, '\n', gene_dict)
    motif_coords = id_motif(args.motifs, intron_dict, exon_dict)
    print(motif_coords)
    #draw_motifs(motif_coords, intron_dict, exon_dict, gene_dict)
    return None


if __name__ == '__main__':
    main()
