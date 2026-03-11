#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2025 Tim Lamberton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Tim Lamberton"
__copyright__ = "Copyright 2025"
__credits__ = ["Tim Lamberton"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import os
import sys
import argparse
import logging
import itertools
from collections import defaultdict

#import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.otu_table_collection import StreamingOtuTableCollection
from singlem.metapackage import Metapackage
from singlem.condense import Condenser
from singlem.taxonomy import *

def debug_write_archive(archive_otu_table, metapackage_path, hits_table):

    otus = StreamingOtuTableCollection()
    otus.add_archive_otu_table_file(archive_otu_table.strip())
    
    if metapackage_path:
        logging.info("Using the metapackage at {}".format(metapackage_path))
        metapackage = Metapackage.acquire(metapackage_path)
    elif not metapackage:
        # Neither were specified, so use the default set of packages
        logging.info("Using default SingleM metapackage")
        metapackage = Metapackage.acquire_default()
    
    if metapackage.version < 3:
        raise Exception("Condense function now only works with version 3+ metapackages.")

    markers = {} # set of markers used to the domains they target
    
    for spkg in metapackage.singlem_packages:
        # ensure v3 packages
        if not spkg.version in [3,4]:
            raise Exception("Only works with v3 or v4 singlem packages.")
        marker_name = spkg.graftm_package_basename()
        markers[marker_name] = spkg.target_domains()
    
    with(open(hits_table, "w")) as f:
        otu_hits_table = {"sample": [],
                          "marker": [],
                          "domain": [],
                          "coverage": []}
        taxon_hits_table = {}
        num_otus = 0
        for sample, sample_otus in otus.each_sample_otus(generate_archive_otu_table=True):
            sample_otus = Condenser()._remove_off_target_otus(sample_otus, markers)
            Condenser()._convert_diamond_best_hit_ids_to_taxonomies(metapackage, sample_otus)

            anc_to_child_to_count = defaultdict(lambda: defaultdict(lambda: 0))
            otu_to_best_hits = []
            for otu in sample_otus:
                best_hit_taxonomies = otu.equal_best_hit_taxonomies()
                good_taxonomies = otu.good_taxonomies()
                if (
                        best_hit_taxonomies is not None or good_taxonomies is not None
                        and otu.taxonomy_assignment_method()
                        in (QUERY_BASED_ASSIGNMENT_METHOD, DIAMOND_ASSIGNMENT_METHOD)
                        ):
                    
                    child_to_anc = {}
                    clean_best_hits = []
                    for taxonomies in (best_hit_taxonomies, good_taxonomies):
                        clean_hits = []
                        if taxonomies is not None:
                            for tax in taxonomies:
                                clean_tax = TaxonomyUtils.clean_taxonomy_string(tax)
                                clean_hits.append(clean_tax)

                                # track ancestors of best hit taxa.
                                child_tax = None
                                for anc_tax in TaxonomyUtils.ancestor_taxonomies(clean_tax):
                                    if child_tax is not None:
                                        child_to_anc[child_tax] = anc_tax
                                    child_tax = anc_tax
                        
                        clean_best_hits.append(clean_hits)

                    # count unique parent-child pairs.
                    for child_tax, anc_tax in child_to_anc.items():
                        anc_to_child_to_count[anc_tax][child_tax] += 1
                    
                    otu_to_best_hits.append((child_to_anc, clean_best_hits, otu.marker, otu.coverage))

            new_num_otus = len(otu_to_best_hits)
            new_markers = [""] * new_num_otus
            new_domains = [""] * new_num_otus
            new_coverages = ["0"] * new_num_otus
            new_taxons = {key: ["0"] * new_num_otus for key in taxon_hits_table.keys()}
            for i, (child_to_anc, (best_hits, good_hits), marker, coverage) in enumerate(otu_to_best_hits):
                new_markers[i] = marker
                new_domains[i] = ";".join(markers[marker])
                new_coverages[i] = f"{coverage:.3}"

                anc_hits = []
                for child_tax, anc_tax in child_to_anc.items():
                    if child_tax not in anc_to_child_to_count:
                        # child is a leaf
                        continue
                    
                    max_child_child_count = max(anc_to_child_to_count[child_tax].values())

                    if anc_to_child_to_count[anc_tax][child_tax] > max_child_child_count:
                        anc_hits.append(child_tax)

                for best_hit_tax in best_hits:
                    if best_hit_tax not in taxon_hits_table:
                        taxon_hits_table[best_hit_tax] = ["0"] * num_otus
                        new_taxons[best_hit_tax] = ["0"] * new_num_otus
                    new_taxons[best_hit_tax][i] = "1"
                for good_tax in good_hits:
                    if good_tax not in taxon_hits_table:
                        taxon_hits_table[good_tax] = ["0"] * num_otus
                        new_taxons[good_tax] = ["0"] * new_num_otus
                    if new_taxons[good_tax][i] == "1":
                        new_taxons[good_tax][i] = "3"
                    else:
                        new_taxons[good_tax][i] = "2"
                for anc_tax in anc_hits:
                    if anc_tax not in taxon_hits_table:
                        taxon_hits_table[anc_tax] = ["0"] * num_otus
                        new_taxons[anc_tax] = ["0"] * new_num_otus
                    if new_taxons[anc_tax][i] == "0":
                        new_taxons[anc_tax][i] = "4"
            
            otu_hits_table["sample"] += [sample] * new_num_otus
            otu_hits_table["marker"] += new_markers
            otu_hits_table["domain"] += new_domains
            otu_hits_table["coverage"] += new_coverages
            for taxon, hits in new_taxons.items():
                taxon_hits_table[taxon] += hits
    
        f.write("\t".join(itertools.chain(otu_hits_table.keys(), taxon_hits_table.keys())))
        f.write("\n")
        for row in zip(*otu_hits_table.values(), *taxon_hits_table.values()):
            f.write("\t".join(row))
            f.write("\n")

    
def debug_write_props(sample_otus, otu_to_taxon_to_props, genes_to_domains, f):
    #[{taxon -> prop}]
    num_otus = len(otu_to_taxon_to_props)
    taxon_to_otu_to_prop = {"marker": [""] * num_otus,
                            "domain": [""] * num_otus,
                            #"sequence": [""] * num_otus,
                            "coverage": ["0"] * num_otus}
    for i, (otu, taxon_to_props) in enumerate(zip(sample_otus, otu_to_taxon_to_props)):
        taxon_to_otu_to_prop["marker"][i] = otu.marker
        taxon_to_otu_to_prop["domain"][i] = ";".join(genes_to_domains[otu.marker])
        #taxon_to_otu_to_prop["sequence"][i] = otu.sequence
        taxon_to_otu_to_prop["coverage"][i] = f"{otu.coverage:.3}"
        for tax, prop in taxon_to_props.items():
            if tax not in taxon_to_otu_to_prop:
                taxon_to_otu_to_prop[tax] = ["0"] * num_otus
            taxon_to_otu_to_prop[tax][i] = f"{float(prop):.3}"
    
    f.write("\t".join(taxon_to_otu_to_prop.keys()))
    f.write("\n")
    for row in zip(*taxon_to_otu_to_prop.values()):
        f.write("\t".join(row))
        f.write("\n")


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--input-archive-otu-table', help="Output hits from this table", required=True)
    parent_parser.add_argument('--metapackage', help = 'Set of SingleM packages to use [default: use the default set]')
    parent_parser.add_argument('--hits-table', help = "TSV output file path", required = True)
    
    args = parent_parser.parse_args()
    debug_write_archive(args.input_archive_otu_table, args.metapackage, args.hits_table)
