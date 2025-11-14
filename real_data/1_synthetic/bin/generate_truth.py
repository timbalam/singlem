#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft, 2025 Tim Lamberton
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
__copyright__ = "Copyright 2022, 2025"
__credits__ = ["Ben Woodcroft", "Tim Lamberton"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os

import polars as pl
import pyarrow
import tempfile
import extern

#sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--coverage-file', required=True, help='Path to coverage file')
    parent_parser.add_argument('--gtdb-bac-tax', required=True, help='Path to GTDB taxonomy file')
    parent_parser.add_argument('--gtdb-ar-tax', required=True, help='Path to GTDB taxonomy file')
    parent_parser.add_argument('--output-condensed', required=True, help='Path to output file in singlem condensed format')

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    output1 = os.path.abspath(args.read1)
    output2 = os.path.abspath(args.read2)
    output_condensed = os.path.abspath(args.output_condensed)

    # Read coverages
    bot_coverages = pl.read_excel(
        args.coverage_file,
        engine = 'xlsx2csv',
        read_options = dict(
            skip_lines = 20
        )
    )
    top_coverages = pl.read_excel(
        args.coverage_file,
        engine = 'xlsx2csv',
        schema_overrides = bot_coverages.schema,
        read_options = dict(
            skip_lines = 2,
            n_rows = 16,
            null_values = "x",
            ignore_errors = True
        )
    )
    coverages = pl.concat([
        top_coverages.select(
            species = 'Organism Name',
            phylum = 'Phylum',
            abundance = 'Genome molecules /   ng gDNA'
        ),
        bot_coverages.select(
            species = 'Organism Name',
            phylum = 'Phylum',
            abundance = 'Genome molecules /   ng gDNA'
        )
    ])
    coverages = coverages.filter(pl.col('abundance') > 0]
    logging.info(f"Read {len(coverages)} coverages > 0.")

    bac = pl.read_csv(args.gtdb_bac_tax, separator = '\t',
                      has_header = False,
                      names = ['accession', 'gtdb_taxonomy'])
    ar = pl.read_csv(args.gtdb_ar_tax, separator = '\t',
                     has_header = False,
                     names = ('accession', 'gtdb_taxonomy'))
    taxonomies = pl.concat([
        bac.select('gtdb_taxonomy'),
        ar.select('gtdb_taxonomy')
    ]).with_columns(
        (pl.col('gtdb_taxonomy')
            .str.extract_groups('p__(<phylum>[^;]+).*;s__(<species>.+)')
          .alias('fields')
          .to_frame()
          .unnest('fields')
        )

    coverages_taxonomies = coverages.join(
        taxonomies,
        on = ["species", "phylum"]
        how = "inner"
    ).select(
        sample = pl.lit("SRR606249"),
        coverage = "abundance",
        taxonomy = "gtdb_taxonomy"
    )

    coverages_taxonomies.write_csv(args.output_condensed,
                                   separator = '\t')

