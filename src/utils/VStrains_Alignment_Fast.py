#!/usr/bin/env python3
import argparse
import os
import time
import subprocess
import numpy
import sys

__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, VStrains Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"

rev_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}

def reverse_seq(seq: str):
    return "".join(rev_dict[x] for x in reversed(seq))

def single_end_read_mapping(seq: str, kmer_htable: dict, index2seqlen: list, split_len: int, len_index2id: int):
    nodes = numpy.zeros(len_index2id, dtype=int)
    coords = [sys.maxsize for _ in range(len_index2id)]
    kindices = [sys.maxsize for _ in range(len_index2id)]

    rlen = len(seq)
    for i in range(rlen - split_len + 1):
        kmer = seq[i : i + split_len]
        if kmer in kmer_htable:
            # found a collide node
            for (rid, rcord) in kmer_htable[kmer]:
                nodes[rid] += 1
                coords[rid] = min(coords[rid], rcord)
                kindices[rid] = min(kindices[rid], i)

    
    saturates = []
    L = 0
    R = 0
    for i, v in enumerate(nodes):
        if coords[i] == sys.maxsize or kindices[i] == sys.maxsize:
            continue
        L = max(coords[i], coords[i] - kindices[i])
        R = min(coords[i] + index2seqlen[i] - 1, coords[i] - kindices[i] + rlen - 1)
        saturate = R - L - (split_len - 1) + 1
        expected = (min(rlen, index2seqlen[i]) - split_len + 1) * (rlen - split_len) / rlen
        if v >= max(min(saturate, expected), 1):
            # print(i,v,"passed")
            saturates.append(i)
    return saturates     


def main():
    print(
        "----------------------Paired-End Information Alignment----------------------"
    )
    parser = argparse.ArgumentParser(
        prog="pe_info",
        description="""Align Paired-End reads to nodes in graph to obtain strong links""",
    )

    parser.add_argument(
        "-g", "--gfa,", dest="gfa", type=str, required=True, help="graph, .gfa format"
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        dest="dir",
        type=str,
        required=True,
        help="output directory",
    )

    parser.add_argument(
        "-f", "--forward", dest="fwd", required=True, help="forward read, .fastq"
    )

    parser.add_argument(
        "-r", "--reverse", dest="rve", required=True, help="reverse read, .fastq"
    )

    parser.add_argument(
        "-k",
        "--kmer_size",
        dest="kmer_size",
        type=int,
        default=128,
        help="unique kmer size",
    )

    args = parser.parse_args()

    # initialize output directory
    if args.dir[-1] == "/":
        args.dir = args.dir[:-1]
    subprocess.check_call("rm -rf {0}".format(args.dir), shell=True)
    os.makedirs(args.dir, exist_ok=True)

    glb_start = time.time()


    # get gfa node informations
    index2id = []
    index2seq = []
    index2seqlen = []
    
    with open(args.gfa, "r") as gfa:
        for Line in gfa:
            splited = (Line[:-1]).split("\t")
            if splited[0] == "S":
                index2id.append(splited[1])
                index2seq.append(splited[2])
                index2seqlen.append(len(splited[2]))
        gfa.close()

    split_len = args.kmer_size + 1

    # construct hash table for gfa nodes with chunck kmer
    kmer_htable = {}
    for i, seq in enumerate(index2seq):
        seqlen = index2seqlen[i]
        for sub_i in range(seqlen - split_len + 1):
            kmer = seq[sub_i : sub_i + split_len]
            rev_kmer = reverse_seq(kmer)
            if kmer in kmer_htable:
                # not unique
                kmer_htable[kmer].append((i, sub_i))
            else:
                # unique
                kmer_htable[kmer] = [(i, sub_i)]
            
            if rev_kmer in kmer_htable:
                # not unique
                kmer_htable[rev_kmer].append((i, sub_i))
            else:
                # unique
                kmer_htable[rev_kmer] = [(i, sub_i)]

    # init nodes pairwise relationship
    len_index2id = len(index2id)
    node_mat = numpy.zeros((len_index2id, len_index2id), dtype=int)
    short_mat = numpy.zeros((len_index2id, len_index2id), dtype=int)

    n_reads = 0
    short_reads = 0
    used_reads = 0

    print("Start aligning reads to gfa nodes")
    fwd_fd = open(args.fwd, "r")
    rve_fd = open(args.rve, "r")
    fwd_reads = fwd_fd.readlines()
    rve_reads = rve_fd.readlines()
    fwd_fd.close()
    rve_fd.close()

    total_size = min(len(fwd_reads) // 4, len(rve_reads) // 4)
    for read_idx in range(total_size):
        if read_idx % 100000 == 0:
            print("Number of processed reads: ", read_idx)
        [_, fseq, _, _] = [s[:-1] for s in fwd_reads[read_idx * 4 : (read_idx + 1) * 4]]
        [_, rseq, _, _] = [s[:-1] for s in rve_reads[read_idx * 4 : (read_idx + 1) * 4]]
        if fseq.count("N") or rseq.count("N"):
            n_reads += 1
        elif len(fseq) < split_len or len(rseq) < split_len:
            short_reads += 1
        else:
            used_reads += 1
            # valid read pair
            lefts = single_end_read_mapping(fseq, kmer_htable, index2seqlen, split_len, len_index2id)
            rights = single_end_read_mapping(rseq, kmer_htable, index2seqlen, split_len, len_index2id)
            
            k = 0
            for i in lefts:
                for i2 in lefts[k:]:
                    short_mat[i][i2] += 1
                k += 1

            k = 0
            for j in rights:
                for j2 in rights[k:]:
                    short_mat[j][j2] += 1
                k += 1

            for i in lefts:
                for j in rights:
                    node_mat[i][j] += 1



    out_file = "{0}/pe_info".format(args.dir)
    out_file2 = "{0}/st_info".format(args.dir)
    subprocess.check_call("touch {0}; echo " " > {0}".format(out_file), shell=True)
    subprocess.check_call("touch {0}; echo " " > {0}".format(out_file2), shell=True)
    with open(out_file, "w") as outfile:
        with open(out_file2, "w") as outfile2:
            for i in range(len_index2id):
                for j in range(len_index2id):
                    outfile.write(
                        "{0}:{1}:{2}\n".format(
                            index2id[i], index2id[j], node_mat[i][j]
                        )
                    )
                    outfile2.write(
                        "{0}:{1}:{2}\n".format(
                            index2id[i], index2id[j], short_mat[i][j]
                        )
                    )
            outfile2.close()
        outfile.close()

    glb_elapsed = time.time() - glb_start
    print("Global time elapsed: ", glb_elapsed)
    print("result stored in: ", out_file)


if __name__ == "__main__":
    main()
    sys.exit(0)
