#!/usr/bin/env python

# import modules
import unittest
import os

from motif_mark_oop import * # TODO: change this


# test files
TEST_MOTIF_FILE = os.path.join(os.path.dirname('../data/input/'),'Fig_1_motifs.txt')

# build test class
class TestMotif(unittest.TestCase):

    def setUp(self) -> None:
        self.m = 'yY'
        self.motif = Motif(self.m)
        return super().setUp()

    def tearDown(self) -> None:
        """
        def tearDown(self):
            Clean up temp files
        
            # Remove temp files after tests
            os.unlink(self.temp_lib1.name)
            os.unlink(self.temp_lib2.name)
            os.unlink(self.temp_lib3.name)
            # remove temp dir
            shutil.rmtree(self.temp_dir)
        """
        return super().tearDown()

    def test_print_motif(self):
        # TODO: build motif class object in each test case
        # TODO: standard test case and then extremes (ie empty string)
        # TODO: what if motif contains character not in dictionary
        # TODO: implimenting exceptions use self.assertRaises
        self.assertEqual(self.motif.print_motif('yY'),'yY')
        
    def test_generate_combos(self):
        self.assertEqual(self.motif.generate_combos('YY'),['[CTU][CTU]','[ctu][ctu]'])

class TestGene(unittest.TestCase):
    
    def setUp(self) -> None:
        self.id_line = '>THOMAS chr1:1-2'
        self.sequence = 'aaaaAAAAaaaa'
        self.gene = Gene(self.id_line,self.sequence)
        # generate motif class
        self.m = 'aaaa'        
        self.motif = Motif(self.m)
        return super().setUp()
    
    def test_extract_gene(self):
        self.assertEqual(self.gene.extract_gene(self.id_line),'THOMAS')

    def test_extract_location(self):
        self.assertEqual(self.gene.extract_location(self.id_line), 'chr1:1-2')

    def test_identify_exons(self):
        self.assertListEqual(self.gene.identify_exons(self.sequence),[(4,8)])

    def test_store_motif(self):
        self.gene.store_motif(self.motif)
        self.assertEqual(self.gene.matches,{self.motif:[]})

    def test_identify_matches(self):
        self.gene.identify_matches(self.motif,self.sequence)
        self.assertDictEqual(self.gene.matches,{self.motif:[4, 0, 8]})

class TestAdditionalMethods(unittest.TestCase):

    def setUp(self) -> None:
        self.id_line = '>THOMAS chr1:1-2'
        self.sequence = 'aaaaAAAAaaaa'
        self.gene = Gene(self.id_line,self.sequence)
        self.m = 'aaaa'        
        self.motif = Motif(self.m)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_read_in_motifs(self):
        self.assertListEqual(read_in_motifs(TEST_MOTIF_FILE),['ygcy', 'GCAUG', 'catag', 'YYYYYYYYYY'])

    def test_update_longest_gene(self):
        master_dict = {'master_list':[self.gene],'longest_gene':3}
        new_master_dict = update_longest_gene(master_dict,self.gene)
        self.assertEqual(master_dict['longest_gene'],len(self.sequence))

if __name__ == '__main__':
    unittest.main()