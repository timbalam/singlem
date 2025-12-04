#!/usr/bin/env python3

#=======================================================================
# Authors: Rossen Zhao, Tim Lamberton
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import os.path
import sys
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.otu_table import OtuTable
from singlem.condense  import Condenser
from singlem.archive_otu_table import ArchiveOtuTable
from singlem.pipe import QUERY_BASED_ASSIGNMENT_METHOD

path_to_script = 'singlem'
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','condense')

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_apply_expectation_maximization_core_trivial(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        # str.split('gene    sample    sequence    num_hits    coverage    taxonomy    read_names    nucleotides_aligned  taxonomy_by_known? read_unaligned_sequences equal_best_hit_taxonomies taxonomy_assignment_method')
        otus.data = [
            ['g1', 'sample1', 'seq1', 1, 1.05,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.05},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split1(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            ['g1', 'sample1', 'seq1', 11,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD],
            ['g1', 'sample1', 'seq2', 11,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 1, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split2(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            ['g1', 'sample1', 'seq1', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            ['g1', 'sample1', 'seq2', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 2.2},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split3(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            ['g1', 'sample1', 'seq1', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            ['g1', 'sample1', 'seq2', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2','Root; d__Bacteria; p;c;o;f;g; tax3'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax2': 1.1},
            species_to_coverage
        )

    def test_apply_expectation_maximization_core_split4(self):
        otus = ArchiveOtuTable()
        otus.fields = ArchiveOtuTable.FIELDS_VERSION4
        otus.data = [
            ['g1', 'sample1', 'seq1', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2'],QUERY_BASED_ASSIGNMENT_METHOD],
            ['g1', 'sample1', 'seq1', 1,1.1,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax1','Root; d__Bacteria; p;c;o;f;g; tax2','Root; d__Bacteria; p;c;o;f;g; tax3'],QUERY_BASED_ASSIGNMENT_METHOD],
            ['g1', 'sample1', 'seq1', 1,1.2,'','','','','',['Root; d__Bacteria; p;c;o;f;g; tax5','Root; d__Bacteria; p;c;o;f;g; tax4'],QUERY_BASED_ASSIGNMENT_METHOD]
        ]
        species_to_coverage, best_hit_taxonomy_sets = Condenser()._apply_species_expectation_maximization_core(otus, 0, {'Bacteria': ['g1']}, min_genes_for_whitelist=0, proximity_cutoff=0)
        self.assertEqual(
            {'Root; d__Bacteria; p;c;o;f;g; tax1': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax2': 1.1, 'Root; d__Bacteria; p;c;o;f;g; tax4': 0.6, 'Root; d__Bacteria; p;c;o;f;g; tax5': 0.6},
            species_to_coverage
        )
        self.assertEqual(sorted([
            ['Root; d__Bacteria; p;c;o;f;g; tax1', 'Root; d__Bacteria; p;c;o;f;g; tax2'],
            ['Root; d__Bacteria; p;c;o;f;g; tax1', 'Root; d__Bacteria; p;c;o;f;g; tax2', 'Root; d__Bacteria; p;c;o;f;g; tax3'],
            ['Root; d__Bacteria; p;c;o;f;g; tax5', 'Root; d__Bacteria; p;c;o;f;g; tax4']
        ]), sorted(best_hit_taxonomy_sets))

if __name__ == "__main__":
    import logging
    # logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    unittest.main()
