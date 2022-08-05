#!/usr/bin/env python3
import subprocess
import argparse
from quast_evaluation import quast_eval

import logging
import threading

usage = "Run multiple benchmark case in parallel"
author = "Runpeng Luo"


def para_eval(config, t_rank):
    logging.info("Thread %s: starting reconstruction", t_rank)

    name, gfa_addr, contig_addr, ref_addr, odir_addr, log_addr, overlap = config
    comm = "time python src/hap_construction.py -gfa {0} -c {1} -overlap {2} -ref {3} -o {4} > {5}".format(
        gfa_addr, contig_addr, overlap, ref_addr, odir_addr, log_addr
    )
    subprocess.check_call(comm, shell=True)

    logging.info("Thread %s: finishing reconstruction, start quast evaluation", t_rank)

    quast_eval(
        "{0}strain.fasta".format(odir_addr), ref_addr, "quast_{0}/".format(name), t_rank
    )

    logging.info("Thread %s: finishing quast evaluation", t_rank)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="benchmark_eval.py", description=usage)
    parser.add_argument(
        "-c",
        "--config_file",
        dest="config_file",
        type=str,
        required=True,
        help="benchmark configuration file",
    )
    args = parser.parse_args()

    assert args.config_file

    # remove all the acc directory
    subprocess.check_call("rm -rf eval_result && mkdir eval_result/", shell=True)

    configs = []
    with open(args.config_file, "r") as config_file:
        count = int(config_file.readline()[:-1])
        config_file.readline()  # delimiter
        for i in range(count):
            name = config_file.readline()[:-1]
            gfa_addr = config_file.readline()[:-1]
            contig_addr = config_file.readline()[:-1]
            ref_addr = config_file.readline()[:-1]
            odir_addr = config_file.readline()[:-1]
            log_addr = config_file.readline()[:-1]
            overlap = config_file.readline()[:-1]
            config_file.readline()

            configs.append(
                (name, gfa_addr, contig_addr, ref_addr, odir_addr, log_addr, overlap)
            )
        config_file.close()

    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")

    logging.info("Main    : before creating thread")
    threads = []
    for i, config in enumerate(configs):
        logging.info("Main    : create and start thread %d.", i)
        thread = threading.Thread(
            target=para_eval,
            args=(
                config,
                i,
            ),
        )
        threads.append(thread)
        thread.start()

    for i, thread in enumerate(threads):
        logging.info("Main    : before joining thread %d.", i)
        thread.join()
        logging.info("Main    : thread %d done", i)

    print("All benchmark is finished, start cleanup")
    # relocate all the log file
    subprocess.check_call(
        "mv hap_*.log eval_result/ && mv acc_* eval_result/", shell=True
    )
    subprocess.check_call("mv quast_* eval_result/", shell=True)
    print("Finished")
