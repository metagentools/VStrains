#!/usr/bin/env python3

import argparse
import sys
import os
import time
import numpy
from datetime import date

from utils import (
    ns_SPAdes,
    ns_Flye,
)
__author__ = "Runpeng Luo"
__copyright__ = ""
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = ""
__status__ = ""

def run(args):
    numpy.seterr(all='raise')
    RUNNER = {
        "spades": ns_SPAdes.run,
        "flye": ns_Flye.run,
    }
    RUNNER[args.assembler](args)


def main():
    parser = argparse.ArgumentParser(prog="ns", 
        description="""Construct full-length viral strains under deno vo approach 
        from contigs and assembly graph, currently supports SPAdes and Flye""")

    parser.add_argument(
        '-a', '--assembler', dest='assembler', type=str, required=True, choices=['spades', 'flye'],
        help="name of the assembler used. [spades, flye]")

    parser.add_argument(
        '-g', '--graph', dest='gfa_file', type=str, required=True, 
        help='path to the assembly graph, (.gfa format)')

    parser.add_argument(
        '-p', '--path', dest='path_file', type=str, required=False,
        help='contig file from SPAdes (.paths format), only required for SPAdes. e.g., contigs.paths')

    parser.add_argument(
        '-i', '--info', dest='info_file', type=str, required=False,
        help='contig information file from Flye (.txt format), only required for Flye. e.g., assembly_info.txt')

    parser.add_argument(
        '-mc', '--minimum_coverage', dest='min_cov', type=int, default=0,
        help=("minimum node coverage cutoff, [default: auto]"))

    parser.add_argument(
        '-ml', '--minimum_contig_length', dest='min_len', default=250, type=int, 
        help=("minimum initial contig length for strains [default: 250]"))

    parser.add_argument(
        '-r', "--reference_fa", dest='ref_file', type=str, 
        help='path to the reference strain, .fasta format, DEBUG_MODE only')

    parser.add_argument(
        '-o', '--output_dir', dest='output_dir', default='acc/', type=str, 
        help='path to the output directory [default: acc/]')
    
    args = parser.parse_args()

    # parsing arguments, sanity check
    if (not args.gfa_file) or (not os.path.exists(args.gfa_file)):
        print("\nPath to the assembly graph is required, (.gfa format)")
        print("Please ensure the path is correct")
        print("\nExiting...\n")
        sys.exit(1)

    if args.assembler.lower() == 'spades':
        if (not args.path_file) or (not os.path.exists(args.path_file)):
            print("\nPath to Contig file from SPAdes (.paths format) is required for SPAdes assmbler option. e.g., contigs.paths")
            print("\nExiting...\n")
            sys.exit(1)   
        
    elif args.assembler.lower() == 'flye':
        if (not args.info_file) or (not os.path.exists(args.info_file)):
            print("\nPath to Contig information file from Flye (.txt format) is required for Flye. e.g., assembly_info.txt")
            print("\nExiting...\n")
            sys.exit(1)
    else:
        print("\nPlease make sure to provide the correct assembler type (SPAdes or Flye).")
        print("\nExiting...\n")
        sys.exit(1)
    
    if args.min_len < 0 or args.min_cov < 0:
        print("\nPlease make sure to provide the correct option (invalid value for min_len or min_cov).")
        print("\nExiting...\n")
        sys.exit(1)

    if args.output_dir[-1] == "/":
        args.output_dir = args.output_dir[:-1]

    # initialize output directory
    os.makedirs(args.output_dir, exist_ok=True)
    try:
        os.makedirs(args.output_dir+"/gfa/")
        os.makedirs(args.output_dir+"/tmp/")
        os.makedirs(args.output_dir+"/paf/")
    except OSError as _:
        print("\nCurrent output directory is not empty")
        print("Please empty/re-create the output directory")
        print("\nExiting...\n")
        sys.exit(1)

    # all good
    start = time.time()
    run(args)
    elapsed = time.time() - start
    print("\nresult is stored in {0}/strain.fasta".format(args.output_dir))
    print("Finished: ", date.today().strftime("%B %d, %Y"))
    print("Elapsed time: ", elapsed)
    print("\nExiting...\n")
    return 0

if __name__ == "__main__":
    main()