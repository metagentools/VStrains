#!/usr/bin/env python3

import argparse

__author__ = "Runpeng Luo"
__copyright__ = ""
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = ""
__status__ = ""

PARSER = argparse.ArgumentParser(
    prog="ns", 
    description="""Construct full-length viral strains under deno vo approach from contigs and assembly graph, 
    currently supports SPAdes and Flye""")

PARSER.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph, .gfa format')
PARSER.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, .paths format')
PARSER.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, help=("minimum strain coverage (optional)"))
PARSER.add_argument('-minlen' '--minimum_contig_length', dest='min_len', default=250, type=int, help=("minimum initial contig length for strains (default 250)"))
PARSER.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, DEBUG_MODE only')
PARSER.add_argument('-o', '--output_dir', dest='output_dir', default='acc/', type=str, help='output directory (default: acc/)')