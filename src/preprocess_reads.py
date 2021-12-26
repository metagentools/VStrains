import sys, os
import json
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort
import subprocess
import argparse
import time


usage = 'Build variation graph using SPAdes metaviral mode, with input pair-end reads and store the graph.'

def main():
    parser = argparse.ArgumentParser(prog='preprocess_reads.py', description=usage)
    parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    parser.add_argument('-spades', '--spades_path', dest='spades', type=str, required=True, help='absolute path to spades executable')
    parser.add_argument('-t', '--threads', default=8, help="Set number of threads used for spades.")
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, required=True)
    args = parser.parse_args()
    
    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()
    
    filepath = os.path.dirname(os.path.abspath(__file__))
    spades = args.spades


    if spades:
        # subprocess.check_call("rm -rf ../{0} && mkdir ../{0}".format(args.output_dir))
        print(filepath)
        subprocess.check_call("rm -rf {0} && mkdir {0}".format(
            args.output_dir), shell=True)

        subprocess.check_call(spades + " --metaviral -1 {0} -2 {1} -o {2}".format(
            args.forward, args.reverse, args.output_dir), shell=True)
    else:
        print("SPAdes executable path haven't specified.")

    t1_stop = time.perf_counter()
    t2_stop = time.process_time()

    print("\preprocess reads completed")
    print("Elapsed time: {:.1f} seconds".format(t1_stop-global_t1_start))
    print("CPU process time: {:.1f} seconds".format(t2_stop-global_t2_start))
    return


if __name__ == "__main__":
    main()