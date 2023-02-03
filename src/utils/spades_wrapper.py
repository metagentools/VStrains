import os
import subprocess
import argparse
import time

__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, vsAware Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="spades_wrapper.py",
        description="""Build assembly graph&contig using SPAdes --careful mode, 
        with input pair-end reads and store the graph.""",
    )
    parser.add_argument(
        "-f",
        "--forward",
        dest="forward",
        type=str,
        required=True,
        help="Forward reads, fastq format",
    )
    parser.add_argument(
        "-r",
        "--reverse",
        dest="reverse",
        type=str,
        required=True,
        help="Reverse reads, fastq format",
    )
    parser.add_argument(
        "-spades",
        "--spades_path",
        dest="spades",
        type=str,
        required=True,
        help="absolute path to spades executable",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="thread_count",
        default=8,
        help="Set number of threads used for SPAdes.",
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", type=str, required=True
    )
    args = parser.parse_args()

    global_t1_start = time.perf_counter()
    global_t2_start = time.process_time()

    filepath = os.path.dirname(os.path.abspath(__file__))
    spades = args.spades

    if spades:
        print(filepath)
        subprocess.check_call(
            "rm -rf {0} && mkdir {0}".format(args.output_dir), shell=True
        )

        subprocess.check_call(
            spades
            + " -1 {0} -2 {1} --careful -t {3} -o {4}".format(
                args.forward, args.reverse, args.thread_count, args.output_dir
            ),
            shell=True,
        )
    else:
        print("SPAdes executable path haven't specified.")

    t1_stop = time.perf_counter()
    t2_stop = time.process_time()

    print("\SPAdes assembly completed")
    print("Elapsed time: {:.1f} seconds".format(t1_stop - global_t1_start))
    print("CPU process time: {:.1f} seconds".format(t2_stop - global_t2_start))
