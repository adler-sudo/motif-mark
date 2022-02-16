#!/usr/bin/env python

# import modules
from motif_mark_oop import *


## TEST CLASSES ##
def test_motif_class():
    """
    tests the motif class
    """
    motif = Motif(motif='YYYY')
    assert motif.motif == 'YYYY', 'motif attribute not working'
    assert motif.evaluate_base('Y') == ['C','T','U'], 'error in evaluate_base() method of Motif class'
    assert motif.evaluate_base('T') == ['T','U'], 'error in evaluate_base() method of Motif class'
    
    assert motif.generate_combos('YT') == ['CT', 'CU', 'TT', 'TU', 'UT', 'UU', 'ct', 'cu', 'tt', 'tu', 'ut', 'uu'], 'error in generate_combos() method of Motif class'

def test_gene_class(id_line='>THOMAS something./something\.',sequence='aaaaAAAAaaaa'):
    motif = Motif(motif='aaaa')
    gene = Gene(id_line,sequence)

    assert gene.extract_gene('>THOMAS something') == 'THOMAS', 'error in extract_gene() method of Gene class'
    assert gene.identify_exons('aaaaAAAAaaaa') == [(4,8)], 'error in identify_exons() method of Gene class'

    gene.identify_matches(motif,sequence='aaAAAAaaaa')
    assert gene.matches == {'AAAA': [(2, 6)], 'aaaa': [(6, 10)]}, 'error in identify_matches() method of Gene class'


## TEST FUNCTIONS ##
def test_read_in_motifs(motif_file):
    assert read_in_motifs(motif_file) == ['ygcy', 'GCAUG', 'catag', 'YYYYYYYYYY'], 'error in read_in_motifs() function'

def test_generate_output_filename(input_file):
    # TODO: need to make this more robust
    assert generate_output_filename(input_file) == 'data/input/Figure_1.png', 'error in generate_output_filename() function'

def test_generate_pycairo(output_file):
    pass

## DEFINE MAIN PROGRAM ##
def main():
    test_motif_file = 'data/input/Fig_1_motifs.txt'
    test_input_file = 'data/input/Figure_1.fasta'

    test_motif_class()
    test_gene_class()

    test_read_in_motifs(test_motif_file)
    test_generate_output_filename(test_input_file)

    print('All tests passed successfully you wizard')


## RUN PROGRAM ##
if __name__ == '__main__':
    main()