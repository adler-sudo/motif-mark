#!/usr/bin/env python

# import modules
from turtle import position
import cairo
import math
import re
import itertools
import random

import argparse
import os

from Bioinfo import DNAbases, IUPACbases, oneline_fasta


# argparse
def parse_args():
    """
    argparse
    """
    parser = argparse.ArgumentParser(description='Generate image representing locations of motifs within gene.')
    parser.add_argument('-f','--input_file',help='Input fasta file.',default='test.fa')
    parser.add_argument('-m','--motif_file',help='File containing one motif per line.',default='test_motif.txt')
    parser.add_argument('-d','--output_dir',help='Specify alternative output directory.')
    parser.add_argument('-n','--new_name',help='Specify alternative name for output png file.')
    args = parser.parse_args()
    return args

# classes
class Motif:
    """
    The motif class
    """
    def __init__(self, motif):
        """
        motif : str
        """
        self.motif = motif
        self.upper_motif = motif.upper()
        self.combos = self.generate_combos(self.upper_motif)
        self.length = len(motif)
        
    # methods
    def print_motif(self,motif:str) -> str:
        return self.motif

    def generate_combos(self,motif:str) -> list:
        """
        Uses IUPACbases to generate the unique motif pattern.
        
        Parameters:
        -----------
        motif : str
            The motif as it appears in the input motifs file.
        
        Returns:
        --------
        pattern : list
            List of all patterns of this motif."""
        pattern = [''.join([IUPACbases[c] for c in motif])]
        return pattern

class Gene:
    """
    The Gene class
    """
    def __init__(self, id_line, sequence:str):
        self.gene = self.extract_gene(id_line)
        self.location = self.extract_location(id_line)
        self.sequence = sequence
        self.exons = self.identify_exons(sequence)
        self.matches = {}
        self.length = len(self.sequence)

    def extract_gene(self,id_line:str) -> str:
        """
        Extracts the gene from the ID line.
        
        Parameters:
        -----------
        id_line : str
            The ID line for the gene as seen in the fasta file.
        
        Returns:
        --------
        gene : str
            The gene name as it appears in the fasta file.
        """
        gene = id_line.split(' ')[0][1:]
        return gene

    def extract_location(self,id_line:str) -> str:
        """
        Extracts the gene location from the ID line.
        
        Parameters:
        -----------
        id_line : str
            The ID line for the gene as seen in the fasta file.
            
        Returns:
        --------
        chr : str
            The location of the gene extracted from the ID line of the fasta file.
        """
        chr = id_line.split(' ')[1]
        return chr

    def identify_exons(self,sequence:str) -> list:
        """
        Identifies the exon locations in the fasta sequence for the gene - indicated by capital bases.
        
        Parameters:
        -----------
        sequence : str
            The entire sequence of the gene as it appears in the fasta file
        
        Returns:
        --------
        exons : list
            The start and end position of each exon in the fasta sequence. Each unit is composed of a tuple: (start position, end position)
        """
        exons = [(m.start(0),m.end(0)) for m in re.finditer(pattern='[A-Z]+',string=sequence)]
        return exons

    def store_motif(self,motif:Motif) -> None:
        """
        Stores the motif object in the gene class attribute matches. This is a dictionary with motif class as key and empty list as value.

        Parameters:
        -----------
        motif : Motif
            Motif class object
        """
        self.matches[motif] = []
        return None
    
    def identify_matches(self,motif:Motif,sequence:str) -> None:
        """
        Identify the matches of a motif in the current sequence.

        Parameters:
        -----------
        motif : Motif
            Motif class object being evaluated.
        sequence : str
            String of the gene glass object
        """
        # store the motif
        self.store_motif(motif)
        for pattern in motif.combos:
            matches = []
            # TODO: need to account for overlapping (for example 'aaaaaa' would only return one match for 'aaaa')
            pattern = '(?={0})'.format(pattern)
            pattern_matches = [m.start(0) for m in re.finditer(pattern,sequence,re.IGNORECASE)]
            if len(pattern_matches) > 0:
                for pattern_match in pattern_matches:
                    self.matches[motif].append(pattern_match)
        return None

# motif reader
def read_in_motifs(motif_file:str) -> list:
    """
    Generates a list with all of the motifs from the motif input file.

    Parameters:
    -----------
    motif_file : str
        The input motif file to extract motif sequences from

    Returns:
    --------
    motifs : list
        List of each of the motifs in the motif input file. Motifs exist as strings in the list.
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

class GeneCollection():
    """
    Collection of Genes to be drawn in pycairo.
    """
    def __init__(self):
        self.gene_list = []
        self.longest_gene = 0
        self.num_genes = 0

    def add_new_gene(self,gene:Gene) -> None:
        """
        Execute series of steps to process new Gene class object
        
        Parameters:
        -----------
        gene : Gene
            Gene class object to be added.
        """
        self.gene_list, self.num_genes = self.append_gene(gene,self.gene_list)
        self.longest_gene = self.evaluate_gene_length(self.longest_gene,gene.length)

    def evaluate_gene_length(self,longest_gene:int,gene_length:int) -> int:
        """
        Evaluate new gene length and change longest_gene if longer than current.
        
        Parameters:
        -----------
        gene_length : int
            Length of gene being added to collection
        
        Returns:
        --------
        gene_length : int
            Length of longest gene in collection
        """
        new_gene_length = max(longest_gene,gene_length)
        return new_gene_length

    def append_gene(self,gene:Gene,gene_list:list) -> tuple:
        """
        Adds new gene to gene_list
        
        Parameters:
        -----------
        gene : Gene
            Gene to be added to list
        gene_list : list
            List of genes 
            
        Returns:
        --------
        gene_list : list
            List of genes with new gene appended
        """
        gene_list.append(gene)
        num_genes = len(gene_list)
        return gene_list, num_genes


#### PYCAIRO ####

class MotifCairo(cairo.Context):
    """
    Class extension of cairo.Context
    """
    def __init__(self,surface):
        super().__init__()
        self.motif_color_dict = {}
        self.surface = surface

    def generate_pycairo_legend(self,motif_color_dict:dict,x:int,y:int) -> None:
        """
        Generate the pycairo legend for the output figure.

        Parameters:
        -----------
        motif_color_dict : dict
            Dicitonary holding the colors associated with each motif. Motif is they key and color is the value.
        x : int
            Beginning x coordinate for pycairo image.
        y : int
            Beginning y coordinate for pycairo image.
        """
        exons_per_row = 4
        exon_counter = 0

        self.move_to(x,y-20)
        self.set_source_rgb(0.5,0.5,0.5)
        self.set_font_size(14)
        self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        self.show_text('Legend')

        for k in self.motif_color_dict:
            self.set_line_width(2)
            self.set_source_rgb(
                self.motif_color_dict[k][0],
                self.motif_color_dict[k][1],
                self.motif_color_dict[k][2])
            self.move_to(x,y-5)
            self.line_to(x,y)
            self.stroke()
            self.move_to(x+10,y)
            self.set_source_rgb(0.5,0.5,0.5)
            self.set_font_size(12)
            self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
            self.show_text(k)
            # shift location right
            x += 120
        return None

    def generate_motif_color(self) -> tuple[float,float,float]:
        """
        Generates and sets random color to be used in pycairo image
        
        Parameters:
        -----------
        context : cairo.Context
            Context object specifying where to draw the image

        Returns:
        --------
        red,blue,green : tuple
            The red blue and green values generated
        """
        red = random.uniform(0,1)
        blue = random.uniform(0,1)
        green = random.uniform(0,1)
        self.set_source_rgb(red,blue,green)
        return red,blue,green

    def set_mark_characteristics(self,motif_object:Motif) -> None:
        """
        Pull motif color from the dictionary and set to source rgb
        """
        if motif_object.motif not in self.motif_color_dict:
            self.motif_color_dict[motif_object.motif] = self.generate_motif_color()
        self.set_source_rgb(
            self.motif_color_dict[motif_object.motif][0],
            self.motif_color_dict[motif_object.motif][1],
            self.motif_color_dict[motif_object.motif][2]
        )
        self.set_line_width(motif_object.length)

    def generate_gene_representation(self,start,gene_class_object,gene_center):
        self.set_line_width(1)
        self.set_source_rgb(0.2,0.2,0.2)
        self.move_to(start,gene_center) 
        self.line_to(gene_class_object.length+start, gene_center)
        self.stroke()

    def generate_exons(self,gene_class_object,start,init_height):
        for exon in gene_class_object.exons:
            self.set_source_rgb(0.2,0.2,0.2)
            self.rectangle(
                exon[0]+start,
                init_height+25,
                exon[1]-exon[0],
                50)
            self.fill()
            self.rectangle(
                exon[0]+start,
                init_height+25,
                exon[1]-exon[0],
                50)
            self.set_source_rgb(1,0,1)
            self.set_line_width(1)
            self.stroke()
        return None

    def generate_new_gene_title(self,surface_width:int,height:int,gene:str,chr_location:str) -> None:
        """
        Generates title for new gene
        """
        # generate gene title
        self.set_source_rgb(0.5,0.5,0.5)
        self.set_font_size(12)
        self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        self.move_to(surface_width*.25,height+10)
        self.show_text(gene)

        # generate gene location subtitle
        self.set_source_rgb(0.5,0.5,0.5)
        self.set_font_size(10)
        self.select_font_face('Arial',cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
        self.move_to(surface_width*.25,height+20)
        self.show_text(chr_location)

    def generate_pycairo_image(
        self,
        gene_collection:GeneCollection,
        output_file:str,
        surface_width:int,
        surface_height:int,
        init_height:int=100) -> None:
        """
        Generate the pycairo image

        Parameters:
        -----------
        master_dict : dict
            Dictionary containing the list of gene objects 'master_list' and 'longest_genes'
        output_file : str
            The output file name to write the image to
        """
        # define start location for gene representation
        start = surface_width * .1
        
        # define mark height
        motif_height = 10

        # loop through master dict
        for gene_class_object in gene_collection.gene_list:
            
            # use to avoid overlaps, adjusting y axis of motif marks
            motif_position_dict = {}

            gene_center = init_height + 50

            self.generate_new_gene_title(surface_width,init_height,gene_class_object.gene,gene_class_object.location)
            self.generate_gene_representation(start,gene_class_object,gene_center)
            self.generate_exons(gene_class_object,start,init_height)

            # generate marks
            for motif_object in gene_class_object.matches:

                self.set_mark_characteristics(motif_object)

                # shift y axis location if overlap
                for match_loc in gene_class_object.matches[motif_object]:
                    stagger_height_adjustment = 0
                    while True:
                        positions_covered = set(range(match_loc,match_loc+len(motif_object.motif)+10)) # plus ten gives padding to marks
                        if stagger_height_adjustment in motif_position_dict.keys():
                            if positions_covered.isdisjoint(motif_position_dict[stagger_height_adjustment]):
                                motif_position_dict[stagger_height_adjustment].update(positions_covered)   
                                break                          
                            else:
                                stagger_height_adjustment += 3
                        else:
                            motif_position_dict[stagger_height_adjustment] = positions_covered
                            break
                        
                    self.move_to(match_loc+start,gene_center-stagger_height_adjustment)
                    self.line_to(match_loc+start,gene_center-stagger_height_adjustment+2)
                    self.stroke()

            # adjust height for next gene
            init_height += 100

        # generate legend
        self.generate_pycairo_legend(self.motif_color_dict,surface_width*0.1,50)
        # save as png
        self.surface.write_to_png(output_file)
        self.surface.finish()
        return None

def generate_image_dimensions(gene_collection:GeneCollection):
    """
    Defines width and height of the surface object
    """
    # define starting height
    height = 100

    # set surface
    surface_width = gene_collection.longest_gene + gene_collection.longest_gene * .2
    surface_height = len(gene_collection.gene_list) * 100 + height

    return surface_width, surface_height

# output name
def generate_output_filename(input_file:str,name:str=None,output_dir:str=None) -> str:
    """
    Generate the output filename given the input file. The output file name should have the same name but with .png extension.

    input_file : str
        The input fasta file to be evaluated.
    name : str
        If the user has specified a new name for their file then this will be a value other than None
    output_dir : str
        If the user has specified a new directory for their file other than default this will be a value other than None.

    Returns:
    --------
    output_file : str
        The output file to write the image to
    """
    extension = '.png'
    
    # define output directory
    if output_dir is not None:
        output_dir = os.path.dirname(output_dir)
        if not os.path.exists(output_dir):
            raise FileNotFoundError(f'{output_dir} directory does not exist. Please choose valid directory.')
    else:
        output_dir = os.path.join(os.path.dirname(input_file),'../output')
    
    # define file
    if name is not None:
        name = os.path.splitext(os.path.basename(name))[0]
        name = name + extension
        output_file = os.path.join(output_dir,name)
    else:
        name = os.path.splitext(os.path.basename(input_file))[0]
        name = name + extension
        output_file = os.path.join(output_dir,name)
    return output_file

# program
def main():
    """
    main program
    """
    # random seed for consistent colors
    random.seed(100)

    args = parse_args()
    
    oneline_file = oneline_fasta(args.input_file)
    output_file = generate_output_filename(args.input_file,args.new_name,args.output_dir)

    motifs = read_in_motifs(args.motif_file)

    open_fasta = open(oneline_file)

    gene_collection = GeneCollection()
    
    while True:
        id_line = open_fasta.readline().rstrip()
        # evalaute if end of file
        if id_line == '':
            break
        sequence = open_fasta.readline().rstrip()
        # instantiate new gene
        gene = Gene(id_line,sequence)
        gene_collection.add_new_gene(gene)
        # process each motif through gene
        for motif in motifs:
            motif_object = Motif(motif)
            gene.identify_matches(motif_object, sequence)

    surface_width,surface_height = generate_image_dimensions(gene_collection)
    
    surface = cairo.SVGSurface('plot.svg',surface_width,surface_height)
    context = MotifCairo(surface)
    
    context.generate_pycairo_image(
        gene_collection,
        output_file=output_file,
        surface_width=surface_width,
        surface_height=surface_height)

    open_fasta.close()

    # remove temporary fasta and plot
    os.remove(oneline_file)
    os.remove('plot.svg')
    print('Your image is ready my majesty <bows>')
    

# run program
if __name__ == '__main__':
    main()
