#!/usr/bin/env python3
import argparse
import subprocess
import sys
import random

def main():
    parser = argparse.ArgumentParser(
            prog="sampling",
            description="""Sampling the pairend fastq file""",
        )

    parser.add_argument(
        "-s",
        "--sampling_ratio",
        dest="sratio",
        type=int,
        required=True,
        help="sampling ratio, 2 for sampling half the dataset, etc.,"
    )
    parser.add_argument(
        "-f",
        "--forward",
        dest="fwd",
        type=str,
        required=True,
        help="forward .fastq file",
    )

    parser.add_argument(
        "-r",
        "--reverse",
        dest="rve",
        type=str,
        required=True,
        help="reverse .fastq file",
    )

    parser.add_argument(
        "-of",
        "--out_forward",
        dest="ofwd",
        type=str,
        required=True,
        help="output forward .fastq file",
    )

    parser.add_argument(
        "-or",
        "--out_reverse",
        dest="orve",
        type=str,
        required=True,
        help="output reverse .fastq file",
    )

    args = parser.parse_args()

    if 1/args.sratio <= 0 or 1/args.sratio >= 1:
        print("error ratio, please input a valid ratio")
        sys.exit(1)

    subprocess.check_call("echo "" > {0}".format(args.ofwd), shell=True)
    subprocess.check_call("echo "" > {0}".format(args.orve), shell=True)

    with open(args.fwd, 'r') as fwd:
        with open(args.rve, 'r') as rve:
            with open(args.ofwd, 'w') as ofwd:
                with open(args.orve, 'w') as orve:
                    flines = fwd.readlines()
                    rlines = rve.readlines()
                    n = len(flines) // 4
                    k = 0
                    print("total number of reads: ", n)
                    for i in range(n):
                        if random.random() > 1/args.sratio:
                            continue
                        k += 1
                        for fcurr in flines[i*4:i*4+4]:
                            ofwd.write(fcurr)
                        for rcurr in rlines[i*4:i*4+4]:
                            orve.write(rcurr)
                    print("sample {0} reads given ratio {1}".format(k, args.sratio))
                    orve.close()
                ofwd.close()
            rve.close()
        fwd.close()
    
    return
    


if __name__ == '__main__':
    sys.exit(main())
