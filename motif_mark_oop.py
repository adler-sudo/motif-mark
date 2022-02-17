#!/usr/bin/env python

# import modules
import cairo
import math
import re
import itertools
import random

import argparse
import os

from Bioinfo import DNAbases, IUPACbases


# argparse
def parse_args():
    """
    argparse
    """
    parser = argparse.ArgumentParser(description='Generate image representing locations of motifs within gene.')
    parser.add_argument('-f','--input_file',help='Input fasta file.',default='test.fa')
    parser.add_argument('-m','--motif_file',help='File containing one motif per line.',default='test_motif.txt')
    args = parser.parse_args()
    return args

# classes
class Motif:
    """
    the motif class
    """
    def __init__(self, motif):
        self.motif = motif.upper()
        self.combos = self.generate_combos(self.motif)
        
    # methods
    def print_motif(self):
        print(self.motif)

    def generate_combos(self,motif):
        pattern = ''.join([IUPACbases[c] for c in motif])
        # add lowercase
        # TODO: do we need to take into account combos of upper AND lower case (binding at start or stop)
        lower = pattern.lower()
        patterns = [pattern, lower]
        return patterns

class Gene:
    """
    the gene class
    """
    def __init__(self, id_line, sequence):
        self.gene = self.extract_gene(id_line)
        self.location = self.extract_location(id_line)
        self.sequence = sequence
        self.exons = self.identify_exons(sequence)
        self.matches = {}
        self.length = len(self.sequence)

    def extract_gene(self,id_line):
        gene = id_line.split(' ')[0][1:]
        return gene

    def extract_location(self,id_line):
        chr = id_line.split(' ')[1]
        return chr

    def identify_exons(self,sequence):
        exons = [(m.start(0),m.end(0)) for m in re.finditer(pattern='[A-Z]+',string=sequence)]
        return exons

    def identify_matches(self,motif,sequence):
        """
        identify motif matches in sequence

        Parameters:
        -----------
        motif : class Motif
            object of class type Motif
        """
        
        # TODO: should store motif objects within gene object (like tokens)
        # TODO: ^this will allow us to access their attributes later
        for pattern in motif.combos:
            matches = []
            # TODO: need to account for overlapping (for example 'aaaaaa' would only return one match for 'aaaa')
            pattern_matches = [(m.start(0),m.end(0)) for m in re.finditer(pattern,sequence)]
            if len(pattern_matches) > 0:
                for pattern_match in pattern_matches:
                    matches.append(pattern_match)

            # include motif and matches in master dict
            if len(matches) > 0:
                self.matches[pattern] = matches


# motif reader
def read_in_motifs(motif_file:str):
    """
    reads in all motifs
    """
    f = open(motif_file,'r')
    motifs = []
    while True:
        motif = f.readline().rstrip()
        # break at end of file
        if motif == '':
            break
        motifs.append(motif)

    f.close()
    
    return motifs

# pycairo
def generate_pycairo_legend(context,motif_color_dict,x,y):
    """
    generate pycairo legend
    """
    
    for k in motif_color_dict:
        context.set_line_width(2)
        context.set_source_rgb(
            motif_color_dict[k][0],
            motif_color_dict[k][1],
            motif_color_dict[k][2])
        context.move_to(x-10,y-5)
        context.line_to(x-10,y+5)
        context.stroke()
        context.move_to(x,y)
        context.set_source_rgb(0.5,0.5,0.5)
        context.set_font_size(12)
        context.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.show_text(k)
        # shift location down
        y += 20


def generate_pycairo(master_list,output_file:str):
    """
    generate pycairo image

    Parameters:
    -----------
    # TODO: need to adjust this once we decide which way we're going
    # master_dict : dict
    #     Contains genes and corresponding positions of matches. Matches are contained in tuples representing
    #     beginning and end location. See structure below for details.

    #     Structure:

    #         {
    #                 <motif_pattern>:
    #                 [(<begin1>,<end1>),(<begin2>,<end2>),...,(<beginn>,<endn>)]
    #         }
    """

    # TODO: go back and functionalize this. may make sense to have pycairo class

    # motif color dict for storing colors
    motif_color_dict = {}

    # set surface
    surface_width = 800
    surface_height = len(master_list) * 100
    surface = cairo.SVGSurface('plot.svg',surface_width,surface_height) # TODO: come back and fix surface later
    context = cairo.Context(surface)
    
    # define starting height
    height = 0
    start = 800 * .05 # TODO: can use this to center our images a bit
    # TODO: ^ need to incorpoarte into our motif location points

    # loop through master dict
    for gene_class_object in master_list:
        # generate gene title
        context.set_source_rgb(0.5,0.5,0.5)
        context.set_font_size(12)
        context.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.move_to(surface_width*.25,height+10)
        context.show_text(gene_class_object.gene)

        # generate gene location subtitle
        context.set_source_rgb(0.5,0.5,0.5)
        context.set_font_size(10)
        context.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        context.move_to(surface_width*.25,height+20)
        context.show_text(gene_class_object.location)

        # generate gene representation
        context.set_line_width(1)
        context.set_source_rgb(0.2,0.2,0.2)
        context.move_to(0,height+50) 
        context.line_to(gene_class_object.length, height+50)
        context.stroke()

        # generate exon representation
        # TODO: put length of exon in gene class (not sure best way to handle)
        context.rectangle(
            gene_class_object.exons[0][0],
            height+25,
            gene_class_object.exons[0][1]-gene_class_object.exons[0][0],
            50)
        context.fill()

        # generate marks
        for pattern in gene_class_object.matches:
            matches = gene_class_object.matches[pattern]

            if pattern in motif_color_dict:
                context.set_source_rgb(
                    motif_color_dict[pattern][0],
                    motif_color_dict[pattern][1],
                    motif_color_dict[pattern][2]
                )
            else:
                # set color for specific motif match
                red = random.random()
                blue = random.random()
                green = random.random()
                
                context.set_source_rgb(red,blue,green)
                motif_color_dict[pattern] = [red,blue,green]
                
            # TODO: put length of motif in motif class (not sure best way to handle)
            context.set_line_width(len(pattern) / 4)

            for match in matches:
                context.move_to(match[0],height+25)
                context.line_to(match[0],height+75)
                context.stroke()

        height += 100

    # generate legend
    generate_pycairo_legend(context,motif_color_dict,surface_width*0.8,surface_height*0.1)
    # save as png
    surface.write_to_png(output_file)
    surface.finish()

# output name
def generate_output_filename(input_file:str):
    """
    generate output file
    """
    output_file = os.path.splitext(input_file)[0]+'.png'
    return output_file

# program
def main():
    """
    main program
    """
    # set random seed
    random.seed(0)

    # read in args
    args = parse_args()
    output_file = generate_output_filename(args.input_file)

    # read in motifs
    motifs = read_in_motifs(args.motif_file)

    # define master list
    master_list = []

    # open fasta
    open_fasta = open(args.input_file)

    while True:
        id_line = open_fasta.readline().rstrip()
        # evalaute if end of file
        if id_line == '':
            break
        sequence = open_fasta.readline().rstrip()
        # instantiate new gene
        gene = Gene(id_line,sequence)
        # process each motif through gene
        for motif in motifs:
            motif_object = Motif(motif)
            gene.identify_matches(motif_object, sequence)
        
        master_list.append(gene)

        # generate pycairo image
        generate_pycairo(
            master_list,
            output_file=output_file)

        
    
    open_fasta.close()
    print(master_list)
    
    

# run program
if __name__ == '__main__':
    main()
