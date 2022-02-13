#!/usr/bin/env python

# import modules
import cairo
import math
import re

import argparse
import os


# argparse
def parse_args():
    """
    argparse
    """
    parser = argparse.ArgumentParser(description='Generate image representing locations of motifs within gene.')
    parser.add_argument('-f','--input_file',help='Input fasta file.')
    parser.add_argument('-m','--motif_file',help='File containing one motif per line.')
    args = parser.parse_args()
    return args

# pycairo
def generate_pycairo(output_file:str):
    """
    generate pycairo image
    """
    surface = cairo.SVGSurface('plot.svg',100,100)

    context = cairo.Context(surface)
    context.set_line_width(1)
    context.move_to(0,50)
    context.line_to(25, 50)
    context.stroke()
    context.move_to(50,50)
    context.line_to(100,50)
    context.stroke()

    context.rectangle(25,25,50,50)
    context.fill()
    surface.write_to_png(output_file)
    surface.finish()

# output name
def generate_output_filename(input_file:str):
    """
    generate output file
    """
    output_file = os.path.splitext(input_file)[0]+'.png'
    return output_file

# run program
if __name__ == '__main__':

    # read in args
    args = parse_args()
    output_file = generate_output_filename(args.input_file)

    # generate image
    generate_pycairo(output_file)
