#!/usr/bin/env python3
import argparse
import os
import time
import subprocess
import numpy
import sys
from datetime import date

import multiprocessing

__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, vsAware Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"

def process_paf_file(index2id, index2reflen, len_index2id, read_ids, fwd_paf_file, rve_paf_file, tid):
    print("Thread {0} start".format(tid))
    print("current pid: {0}".format(os.getpid()))
    start = time.time()

    node_mat = numpy.zeros((len_index2id, len_index2id), dtype=int)
    short_mat = numpy.zeros((len_index2id, len_index2id), dtype=int)

    id2index = {}
    for i in range(len_index2id):
        id2index[index2id[i]] = i

    read2index = {}
    index2read = numpy.array([(k,fwdlen,revlen) for (k, _, _, fwdlen,revlen) in read_ids], dtype=int)

    conf_alns_f = [None for _ in index2read]
    # numpy.array([None for _ in index2read], dtype=object)
    conf_cords_f = [None for _ in index2read]
    # numpy.array([None for _ in index2read], dtype=object)

    conf_alns_r = [None for _ in index2read]
    # numpy.array([None for _ in index2read], dtype=object)
    conf_cords_r = [None for _ in index2read]
    # numpy.array([None for _ in index2read], dtype=object)
    
    for i, (glb_index, f_local_inds, r_local_inds, _, _) in enumerate(read_ids):
        read2index[glb_index] = i
        conf_alns_f[i] = [[] for _ in range(f_local_inds)]
        # numpy.array([[] for _ in range(f_local_inds)], dtype=object)
        conf_cords_f[i] = [[] for _ in range(f_local_inds)]
        # numpy.array([[] for _ in range(f_local_inds)], dtype=object)
        conf_alns_r[i] = [[] for _ in range(r_local_inds)]
        # numpy.array([[] for _ in range(r_local_inds)], dtype=object)
        conf_cords_r[i] = [[] for _ in range(r_local_inds)]
        # numpy.array([[] for _ in range(r_local_inds)], dtype=object)
    
    for file in [fwd_paf_file, rve_paf_file]:
        with open(file, "r") as fwd_paf:
            file_count = 0
            for line in fwd_paf:
                if round((time.time() - start) % 60) == 1:
                    print("Thread {0}: Processed {1} mapping up to now.".format(tid, file_count))
                if line == "\n":
                    break
                splited = (line[:-1]).split("\t")
                seg_no = splited[0]
                [glb_seg_no, sub_no] = seg_no.split("_")
                ref_no = str(splited[5])
                ref_start_coord = int(splited[7]) # 0-based
                nm = int(splited[10]) - int(splited[9])
                if nm == 0 and int(splited[10]) == 128:
                    if file == fwd_paf_file:
                        conf_alns_f[read2index[int(glb_seg_no)]][int(sub_no)].append(id2index[ref_no])
                        conf_cords_f[read2index[int(glb_seg_no)]][int(sub_no)].append(ref_start_coord)
                    else:
                        conf_alns_r[read2index[int(glb_seg_no)]][int(sub_no)].append(id2index[ref_no])
                        conf_cords_r[read2index[int(glb_seg_no)]][int(sub_no)].append(ref_start_coord)

                # cur_status = None

                # if file == fwd_paf_file:
                #     cur_status = conf_alns_f[read2index[int(glb_seg_no)]][int(sub_no)]
                # else:
                #     cur_status = conf_alns_r[read2index[int(glb_seg_no)]][int(sub_no)]
                # if nm == 0 and int(splited[10]) == 128:
                #     if cur_status != -1:
                #         # another perfect match, ambiguous, mark as empty string
                #         cur_status = -2
                #     else:
                #         # first perfect match ever found
                #         cur_status = id2index[ref_no]
                #         if file == fwd_paf_file:
                #             conf_cords_f[read2index[int(glb_seg_no)]][int(sub_no)] = ref_start_coord
                #         else:
                #             conf_cords_r[read2index[int(glb_seg_no)]][int(sub_no)] = ref_start_coord

                # if file == fwd_paf_file:
                #     conf_alns_f[read2index[int(glb_seg_no)]][int(sub_no)] = cur_status
                # else:
                #     conf_alns_r[read2index[int(glb_seg_no)]][int(sub_no)] = cur_status
                file_count += 1
            fwd_paf.close()

    subprocess.check_call(
        "rm {0}".format(fwd_paf_file), shell=True
    )
    subprocess.check_call(
        "rm {0}".format(rve_paf_file), shell=True
    )
    nonunique_counter = 0
    def retrieve_single_end_saturation(glb_index, conf_alns, conf_cords, rlen=250, ks=128):
        nodes = numpy.zeros(len_index2id, dtype=int)
        coords = [None for _ in range(len_index2id)]
        kindices = [None for _ in range(len_index2id)]
        for i, sub_aln_statuses in enumerate(conf_alns[glb_index]):
            if len(sub_aln_statuses) > 1:
                nonunique_counter += 1
            for j, sub_aln_status in enumerate(sub_aln_statuses):
                nodes[sub_aln_status] += 1
                if coords[sub_aln_status] == None:
                    coords[sub_aln_status] = conf_cords[glb_index][i][j]
                else:
                    coords[sub_aln_status] = min(coords[sub_aln_status], conf_cords[glb_index][i][j])
                if kindices[sub_aln_status] == None:
                    kindices[sub_aln_status] = i
                else:
                    kindices[sub_aln_status] = min(kindices[sub_aln_status], i)
        saturates = []
        L = 0
        R = 0
        for i, v in enumerate(nodes):
            if coords[i] == None or kindices[i] == None:
                continue
            # if coords[i] > 0 and kindices[i] > 0:
            #     continue
            L = max(coords[i], coords[i] - kindices[i])
            R = min(coords[i]+index2reflen[i] - 1, coords[i] - kindices[i] + rlen - 1)
            saturate = R - L - 127 + 1
            expected = (min(rlen, index2reflen[i]) - ks + 1) * (rlen - ks)/rlen
            # print("current node: ",index2id[i], "kmer-count:", v)
            # print("saturate: ", saturate)
            # print("L: ", L, "R: ", R)
            # print("coords on ref: ", coords[i])
            # print("kindex on read: ", kindices[i])
            # at most 10 error in first 5 or last 5 kindexes
            if v >= max(min(saturate, expected), 1):
                # print(i,v,"passed")
                saturates.append(i)
        return saturates

    for (glb_id, fwdlen, revlen) in index2read:
        glb_index = read2index[glb_id]
        if round((time.time() - start) % 60) == 0:
            print("{0}, thread: {1} current read id: {2}".format(date.today().strftime("%B %d, %Y"), tid, glb_id))
        
        
        lefts = retrieve_single_end_saturation(glb_index, conf_alns_f, conf_cords_f, fwdlen)
        rights = retrieve_single_end_saturation(glb_index, conf_alns_r, conf_cords_r, revlen)
        
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

        # free up memory
        conf_alns_f[glb_index] = None
        conf_alns_r[glb_index] = None

    elapsed = time.time() - start
    print("thread: {0} found non unique kmer count: {1}".format(tid, nonunique_counter))
    print("thread: {0} time spent for processing paf file: {1}".format(tid, elapsed))
    return node_mat, short_mat

def batch_split(fwd_file: str, rve_file: str, temp_dir: str, batch_size: int, do_split: bool):
    """split the read file into several 
    Args:
        fwd_file (str): _description_
        rve_file (str): _description_
        batch_size (int): _description_
    Returns:
        list: list of batch files
    """
    split_len = 128
    n_reads = 0
    short_reads = 0
    used_reads = 0
    fkmer = 0
    rkmer = 0

    temp_file_fwd = None
    temp_file_rve = None
    local_reads = 0
    local_list = []
    batch_count = 0
    read_summary = []
    sub_files = []
    # forward reverse read processing
    with open(fwd_file, "r") as fwd:
        with open(rve_file, "r") as rve:
            fwd_reads = fwd.readlines()
            rev_reads = rve.readlines()
            total_size = min(len(fwd_reads) // 4, len(rev_reads) // 4)
            # marker_test = 1
            # total_size = min(marker_test, total_size) 
            for i in range( total_size):
                if i % batch_size == 0:
                    print("Processed {0} reads up to now.".format(i))
                [_, fseq, _, feval] = [
                    s[:-1] for s in fwd_reads[i * 4 : (i + 1) * 4]
                ]
                [_, rseq, _, reval] = [
                    s[:-1] for s in rev_reads[i * 4 : (i + 1) * 4]
                ]
                if fseq.count("N") or rseq.count("N"):
                    n_reads += 1
                    continue
                if len(fseq) < 128 or len(rseq) < 128:
                    short_reads += 1
                    continue
                used_reads += 1
                local_reads += 1
                local_list.append((fseq, feval, rseq, reval))
                if local_reads == batch_size or (local_reads > 0 and i == total_size - 1):
                    # file creation
                    sub_fwd_filename = "{0}/temp_forward_{1}.fastq".format(temp_dir, batch_count)
                    sub_rve_filename = "{0}/temp_reverse_{1}.fastq".format(temp_dir, batch_count)
                    subprocess.check_call("touch {0}; echo " " > {0}".format(sub_fwd_filename), shell=True)
                    subprocess.check_call("touch {0}; echo " " > {0}".format(sub_rve_filename), shell=True)
                    temp_file_fwd = open(sub_fwd_filename, "w")
                    temp_file_rve = open(sub_rve_filename, "w")
                    
                    read_ids = []
                    if do_split:
                        for j, (fseq, feval, rseq, reval) in enumerate(local_list):
                            fread_id_subs = len(fseq) - split_len + 1
                            rread_id_subs = len(rseq) - split_len + 1
                            prefix_name = "@{0}_".format(j)
                            # forward
                            for sub_i in range(len(fseq) - split_len + 1):
                                subfread = fseq[sub_i : sub_i + split_len]
                                subfeval = feval[sub_i : sub_i + split_len]
                                temp_file_fwd.write(prefix_name + "{0} /1\n".format(sub_i))
                                temp_file_fwd.write(subfread + "\n")
                                temp_file_fwd.write("+\n")
                                temp_file_fwd.write(subfeval + "\n")
                            fkmer += len(fseq) - split_len + 1
                            # reverse
                            for sub_i in range(len(rseq) - split_len + 1):
                                subrread = rseq[sub_i : sub_i + split_len]
                                subreval = reval[sub_i : sub_i + split_len]
                                temp_file_rve.write(prefix_name + "{0} /2\n".format(sub_i))
                                temp_file_rve.write(subrread + "\n")
                                temp_file_rve.write("+\n")
                                temp_file_rve.write(subreval + "\n")
                            rkmer += len(rseq) - split_len + 1
                            read_ids.append((j, fread_id_subs, rread_id_subs, len(fseq), len(rseq)))
                    else:
                        for j, (fseq, feval, rseq, reval) in enumerate(local_list):
                            prefix_name = "@{0}_".format(j)
                            temp_file_fwd.write(prefix_name + "{0} /1\n".format(0))
                            temp_file_fwd.write(fseq + "\n")
                            temp_file_fwd.write("+\n")
                            temp_file_fwd.write(feval + "\n")
                            
                            temp_file_rve.write(prefix_name + "{0} /2\n".format(0))
                            temp_file_rve.write(rseq + "\n")
                            temp_file_rve.write("+\n")
                            temp_file_rve.write(reval + "\n")
                            read_ids.append((j, 1, 1, len(fseq), len(rseq)))
                    temp_file_fwd.close()
                    temp_file_rve.close()
                    read_summary.append(read_ids)
                    sub_files.append((sub_fwd_filename, sub_rve_filename))
                    local_reads = 0
                    local_list = []
                    batch_count += 1 
        fwd.close()
        rve.close()

    print("total number of reads (before): ", total_size)
    print("total reads containing N: ", n_reads)
    print("total reads too short [<128]: ", short_reads)
    print("total number of reads (used): ", used_reads)
    print("total number of forward reads kmer: ", fkmer)
    print("total number of reverse reads kmer: ", rkmer)
    return read_summary, sub_files

def minimap_alignment(fasta_file, sub_files, temp_dir):
    paf_files = []
    for i, (sub_fwd_filename, sub_rve_filename) in enumerate(sub_files):
        print("minimap reads {0},{1} to graph..".format(sub_fwd_filename, sub_rve_filename))
        start = time.time()
        sub_fwd_paf = "{0}/temp_fwd_aln_{1}.paf".format(temp_dir, i)
        subprocess.check_call(
            "minimap2 -c -t 16 {0} {1} > {2}".format(
                fasta_file, sub_fwd_filename, sub_fwd_paf
            ),
            shell=True,
        )
        # -B 40 -O 20,50 -E 30,10 -z 1,1 -k 27 -w 18 -s 256 
        subprocess.check_call(
            "rm {0}".format(sub_fwd_filename), shell=True
        )

        sub_rve_paf = "{0}/temp_rve_aln_{1}.paf".format(temp_dir, i)
        subprocess.check_call(
            "minimap2 -c -t 16 {0} {1} > {2}".format(
                fasta_file, sub_rve_filename, sub_rve_paf
            ),
            shell=True,
        )
        subprocess.check_call(
            "rm {0}".format(sub_rve_filename), shell=True
        )

        paf_files.append((sub_fwd_paf, sub_rve_paf))
        elapsed = time.time() - start
        print("Time spent for minimap2: ", elapsed)
    return paf_files

def main():
    print("----------------------Paired-End Information Alignment----------------------")
    parser = argparse.ArgumentParser(
        prog="pe_info",
        description="""Align Paired-End reads to nodes in graph to obtain strong links""",
    )

    parser.add_argument(
        "-g", "--gfa,", dest="gfa", type=str, required=True, help="graph, .gfa format"
    )

    parser.add_argument(
        "-o", "--output_dir", dest="dir", type=str, required=True, help="output directory",
    )

    parser.add_argument(
        "-f", "--forward", dest="fwd", required=True, help="forward read, .fastq"
    )

    parser.add_argument(
        "-r", "--reverse", dest="rve", required=True, help="reverse read, .fastq"
    )

    args = parser.parse_args()

    # initialize output directory
    if args.dir[-1] == "/":
        args.dir = args.dir[:-1]
    subprocess.check_call("rm -rf {0}".format(args.dir), shell=True)
    os.makedirs(args.dir, exist_ok=True)

    glb_start = time.time()
    tmp_g2s_file = "{0}/temp_graph_seq.fasta".format(args.dir)

    # convert gfa to fasta file
    index2id = []
    index2reflen = []
    with open(args.gfa, "r") as gfa:
        with open(tmp_g2s_file, "w") as fasta:
            for Line in gfa:
                splited = (Line[:-1]).split("\t")
                if splited[0] == "S":
                    fasta.write(">{0}\n{1}\n".format(splited[1], splited[2]))
                    index2id.append(splited[1])
                    index2reflen.append(len(splited[2]))
            fasta.close()
        gfa.close()

    # split reads to several batches
    read_summary, sub_files = batch_split(args.fwd, args.rve, args.dir, 40000, True)
    # minimap2 reads to fasta file
    paf_files = minimap_alignment(tmp_g2s_file, sub_files, args.dir)

    len_index2id = len(index2id)
    # process paf files
    local_mats = None
    args_proc = [(index2id, index2reflen, len_index2id, read_summary[i], paf_files[i][0], paf_files[i][1], i) for i in range(len(paf_files))]
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        local_mats = pool.starmap(process_paf_file, args_proc)
        pool.close()
        pool.join()

    print("All processes have finished their job, combine the result.")
    # combine all the outputs
    glb_node_mat = numpy.sum(numpy.array([mat for (mat,_) in local_mats]), axis=0)
    glb_strand_mat = numpy.sum(numpy.array([mat for (_,mat) in local_mats]), axis=0)
    out_file = "{0}/pe_info".format(args.dir)
    out_file2 = "{0}/st_info".format(args.dir)
    subprocess.check_call("touch {0}; echo " " > {0}".format(out_file), shell=True)
    with open(out_file, "w") as outfile:
        with open(out_file2, "w") as outfile2:
            for i in range(len_index2id):
                for j in range(len_index2id):
                    outfile.write("{0}:{1}:{2}\n".format(index2id[i], index2id[j], glb_node_mat[i][j]))
                    outfile2.write("{0}:{1}:{2}\n".format(index2id[i], index2id[j], glb_strand_mat[i][j]))
            outfile2.close()
        outfile.close()

    glb_elapsed = time.time() - glb_start
    print("Global time elapsed: ", glb_elapsed)
    print("result stored in: ", out_file)


if __name__ == "__main__":
    main()
    sys.exit(0)