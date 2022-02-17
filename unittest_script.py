#!/usr/bin/env python

# import modules
import unittest
from motif_mark_oop import Motif


# build test class
class TestMotif(unittest.TestCase):

    def test_print_motif(self):
        motif = Motif('yY')
        self.assertEqual(motif.motif,'YY')

    def test_generate_combos(self):
        pass

    def test_evaluate_base(self):
        pass

if __name__ == '__main__':
    unittest.main()