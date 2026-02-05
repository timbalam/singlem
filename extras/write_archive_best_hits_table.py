#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2020 Ben Woodcroft
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

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2022"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import os
import sys
import argparse
import logging
import itertools

#import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.otu_table_collection import StreamingOtuTableCollection
from singlem.metapackage import Metapackage
from singlem.condense import Condenser


def debug_write_archive(archive_otu_tables, metapackage_path, output_dir):

    otus = StreamingOtuTableCollection()
    if archive_otu_tables:
        for o in archive_otu_tables:
            otus.add_archive_otu_table_file(o.strip())
    
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
    
    os.makedirs(output_dir, exist_ok = True)
    for sample, sample_otus in otus.each_sample_otus(generate_archive_otu_table=True):
        sample_otus = Condenser()._remove_off_target_otus(sample_otus, markers)
        Condenser()._convert_diamond_best_hit_ids_to_taxonomies(metapackage, sample_otus)

        with(open(os.path.join(output_dir, f"{sample}_hits.tsv"), "w")) as f:
            debug_write_best_hits(sample_otus, markers, f)

    
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

def debug_write_best_hits(sample_otus, genes_to_domains, f):
    sample_otus = list(sample_otus)
    num_otus = len(sample_otus)
    taxon_to_otu_to_hit = {"marker": [""] * num_otus,
                        "domain": [""] * num_otus,
                        #"sequence": [""] * num_otus,
                        "coverage": ["0"] * num_otus}
    for i, otu in enumerate(sample_otus):
        taxon_to_otu_to_hit["marker"][i] = otu.marker
        taxon_to_otu_to_hit["domain"][i] = ";".join(genes_to_domains[otu.marker])
        #taxon_to_otu_to_hit["sequence"][i] = otu.sequence
        taxon_to_otu_to_hit["coverage"][i] = f"{otu.coverage:.3}"
        for best_hit_tax in itertools.chain(otu.equal_best_hit_taxonomies(),
                                            otu.good_taxonomies()):
            if best_hit_tax not in taxon_to_otu_to_hit:
                taxon_to_otu_to_hit[best_hit_tax] = ["0"] * num_otus
            taxon_to_otu_to_hit[best_hit_tax][i] = "1"
    
    f.write("\t".join(taxon_to_otu_to_hit.keys()))
    f.write("\n")
    for row in zip(*taxon_to_otu_to_hit.values()):
        f.write("\t".join(row))
        f.write("\n")


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--input-archive-otu-tables', '--input-archive-otu-table', nargs = '+', help = "Condense from these archive tables", required = True)
    parent_parser.add_argument('--metapackage', help = 'Set of SingleM packages to use [default: use the default set]')
    parent_parser.add_argument('--output-dir', help = "output directory", required = True)
    
    args = parent_parser.parse_args()

    debug_write_archive(args.input_archive_otu_tables, args.metapackage, args.output_dir)
