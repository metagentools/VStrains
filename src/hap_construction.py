#!/usr/bin/env python3

import sys, os
import subprocess
import argparse
# import re

# import graph_tool
from graph_tool.topology import all_circuits, all_shortest_paths, min_spanning_tree
from graph_tool.topology import transitive_closure, is_DAG
from graph_tool.draw import graph_draw
from graph_tool.all import Graph

from math import ceil
import numpy
import heapq
import matplotlib.pyplot as plt
import seaborn
import pandas

from collections import deque


from graph_converter import *
from search_algos import *

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = False
TEMP_DIR = "acc/"

def main():
    """
    --------------------------------------------OVERALL FLOW----------------------------------------
    Input: Graph, contig
    operation:
    -> START
    -> flip graph [DONE]
    Output --> graph_L0.gfa

    -> tip removal based on minimap2 [DONE] only apply if graph is not dag
    Output --> t_graph_L1.gfa

    -> cycle detection and node partition
    Output --> nc_graph_L2p.gfa

    -> node depth rebalance + assign edge flow [DONE] 
    Output --> dt_graph_L2.gfa

    -> removed all the node that less than threshold.
    Output --> sdt_graph_L3.gfa

    -> contig coverage rebalance to minimum edge flow along the contig

    -> split graph branch if supported by contig and cov difference < threshold
    Output --> bsdt_graph_L4.gfa

    -> graph compactification (without the contig involved path)
    Output --> cbsdt_graph_L5.gfa

    -> construct contig clique graph
    Output --> cliq_graph.png

    -> based on the contig clique grah topology, concat the contig via the following priority
    PRIORITY
    HIGH
        --> self-cycle contig if no out-in edges exist
        --> self-cycle contig, with in-out from same alter contig
        --> 
    LOW
    -> END
    ------------------------------------------------------------------------------------------------
    """
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, default=100, help=("minimum coverage for strains"))
    parser.add_argument('-minlen', '--minimum_strain_length', dest='min_len', default=4000, type=int, help=("minimum strain length"))
    parser.add_argument('-maxlen', '--maximum_strain_length', dest='max_len', default=12000, type=int, help=("maximum strain length"))
    parser.add_argument('-overlap', '--vertex_overlap', dest='overlap', default=127, type=int, help=("adjacent vertex overlap in the graph"))
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, DEBUG_MODE only')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, help='output directory')
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')
    args = parser.parse_args()
    if not args.gfa_file:
        print("gfa file is not imported")
        return 1

    if not args.output_dir:
        print("output directory unspecified, use default dir: acc/")
    else:
        global TEMP_DIR
        TEMP_DIR = args.output_dir
    
    subprocess.check_call("rm -rf {0} && mkdir {0}".format(TEMP_DIR), shell=True)

    print("----------------------------------INPUT---------------------------------------")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    contig_dict, node_to_contig_dict, edge_to_contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_len, args.min_cov)
    
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}graph_L0.gfa".format(TEMP_DIR))
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict, "{0}graph_L0.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
        contig_dict_to_fasta(graph, contig_dict, simp_node_dict, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
        minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    
    print("----------------------------------TIP REMOVAL---------------------------------------")
    # apply tip removal for cyclic graph.
    if not graph_is_DAG(graph0, simp_node_dict0):
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(graph0, simp_node_dict0, simp_edge_dict0, TEMP_DIR, args.overlap)

        contig_dict_fix(graph0, simp_node_dict0, simp_edge_dict0, contig_dict, args.overlap)
    else:
        print("Graph is DAG, tip removal skipped.")
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1.gfa".format(TEMP_DIR))

    print("------------------------------NODE PARTITION-----------------------------------")
    # store the no-cycle nodes in nc_graph_L3p.gfa
    noncyc_nodes = None
    simple_paths = None
    # is dag may not work if the graph is not connected as several parts
    if not graph_is_DAG(graph1, simp_node_dict1):
        noncyc_nodes, simple_paths = node_partition(graph1, simp_node_dict1, simp_edge_dict1, TEMP_DIR)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth has been stored in the dict
    # ratio: normalised balanced node depth / previous node depth
    _, _, _, ratio = coverage_rebalance(graph1, simp_node_dict1, simp_edge_dict1)
    
    if noncyc_nodes != None and simple_paths != None:
        print("rebalance linear subgraph now..")
        graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}nc_graph_L2p.gfa".format(TEMP_DIR))
        coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc)
        print("Done, start coverage merge")

        for no, node in simp_node_dictnc.items():
            cnode = simp_node_dict1[no]
            merge_dp = graphnc.vp.dp[node] + graph1.vp.dp[cnode]
            graph1.vp.dp[cnode] = merge_dp   
    else:
        print("no linear subgraph available..")


    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}dt_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}dt_graph_L2.gfa".format(TEMP_DIR))
    assign_edge_flow(graph2, simp_node_dict2, simp_edge_dict2)

    print("-------------------------------GRAPH SIMPLIFICATION AND REBALANCE-----------------------------------")
    mediandp = numpy.median([graph2.vp.dp[node] for node in simp_node_dict2.values()])
    print("MEDIAN NODE DEPTH: ", mediandp)
    # TODO more
    THRESHOLD= mediandp/20
    ids = []
    for no, node in simp_node_dict2.items():
        if graph2.vp.dp[node] < THRESHOLD:
            ids.append(no)
    print("nodes that less than THRESHOLD: ", THRESHOLD, " cov be removed: ", list_to_string(ids))
    graph_simplification(graph2, simp_node_dict2, simp_edge_dict2, node_to_contig_dict, edge_to_contig_dict, THRESHOLD)

    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    coverage_rebalance(graph3, simp_node_dict3, simp_edge_dict3)
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)

    print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    # re-evaluate the contig coverage
    contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict)

    print("-----------------------GRAPH BRANCH SPLIT & COMPACTIFICATION-------------------------------")
    total_removed_branch = 0
    num_split = -1
    while num_split != 0:
        num_split = graph_splitting(graph3, simp_node_dict3, simp_edge_dict3, contig_dict, THRESHOLD, strict_mode=False)
    
        graph_to_gfa(graph3, simp_node_dict3, simp_edge_dict3, "{0}bsdt_graph_L4.gfa".format(TEMP_DIR))
        graph4, simp_node_dict4, simp_edge_dict4 = flipped_gfa_to_graph("{0}bsdt_graph_L4.gfa".format(TEMP_DIR))

        coverage_rebalance(graph4, simp_node_dict4, simp_edge_dict4)

        contig_cov_fix(graph4, simp_node_dict4, simp_edge_dict4, contig_dict)

        simp_path_dict = simple_paths_to_dict(graph4, simp_node_dict4, simp_edge_dict4, contig_dict, args.overlap)
        simp_path_compactification(graph4, simp_node_dict4, simp_edge_dict4, simp_path_dict, args.overlap)

        graph_to_gfa(graph4, simp_node_dict4, simp_edge_dict4, "{0}cbsdt_graph_L5.gfa".format(TEMP_DIR))
        graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}cbsdt_graph_L5.gfa".format(TEMP_DIR))
        assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
        if num_split != 0:
            total_removed_branch += num_split
            graph3 = graph5
            simp_node_dict3 = simp_node_dict5
            simp_edge_dict3 = simp_edge_dict5
    print("Total branch removed: ", total_removed_branch)
    # print("-----------------------------> contig variation plot")
    # # TODO
    # contig_variation_span(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.overlap, args.max_len, TEMP_DIR)
    
    print("-------------------------------CONTIG CLIQUE GRAPH BUILD-----------------------------------")
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}cbsdt_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    
    cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict = contig_clique_graph_build(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.max_len, THRESHOLD, args.overlap)
    draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), TEMP_DIR, "cliq_graphL1.png")
    print("-------------------------------CONTIG PAIRWISE CONCATENATION-----------------------------------")
    concat_contig_dict = contig_pairwise_concatenation(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, 
    cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict, 
    args.min_cov, args.min_len, args.max_len, args.overlap, THRESHOLD, TEMP_DIR)

    # partial result
    contig_dict_to_fasta(graph5, concat_contig_dict, simp_node_dict5, args.overlap, "{0}concat_contig.fasta".format(TEMP_DIR))
    contig_dict_to_path(concat_contig_dict, "{0}concat_contig.paths".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}concat_contig.fasta".format(TEMP_DIR), "{0}concat_contig_to_strain.paf".format(TEMP_DIR))

    print("-------------------------------------STRAIN EXTENSION--------------------------------------")
    extended_contig_dict = strain_extension(graph5, simp_node_dict5, simp_edge_dict5, concat_contig_dict, args.min_cov, args.max_len, THRESHOLD, args.overlap)
    # graph_to_gfa(graph5, simp_node_dict5, simp_edge_dict5, "{0}stecbsdt_graph_L7.gfa".format(TEMP_DIR))

    # print("----------------------------------FINAL STRAIN EXTRACTION----------------------------------")
    # final_strain_extraction(graph5, simp_node_dict5, simp_edge_dict5, extended_contig_dict, THRESHOLD, args.overlap)
    # print("-------------------------------LOCAL SEARCH OPTIMISATION-----------------------------------")
    # local_search_optimisation(graph5, simp_node_dict5, simp_edge_dict5, concat_contig_dict, args.min_cov, args.min_len, args.max_len, args.overlap, THRESHOLD)
    # stat
    contig_dict_to_fasta(graph5, extended_contig_dict, simp_node_dict5, args.overlap, "{0}extended_contig.fasta".format(TEMP_DIR))
    contig_dict_to_path(extended_contig_dict, "{0}extended_contig.paths".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}extended_contig.fasta".format(TEMP_DIR), "{0}extended_contig_to_strain.paf".format(TEMP_DIR))  
    return 0

curr_path = []
def node_partition(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_dir):
    """
    partition the graph into cyclic node and linear intersect cyclic + linear node

    store the cyc node into one graph object, where noncyc node into another object
    """
    def noncyc_path(graph: Graph, src, tgt):
        all_noncyc_paths = []
        visited = {}
        for v in graph.vertices():
            visited[v] = False

        def dfs(u, v):
            global curr_path
            visited[u] = True
            curr_path.append(u)
            # print(path_to_id_string(graph, curr_path, graph.vp.id[u]))
            if u == v:
                print(path_to_id_string(graph, curr_path, "path found"))
                all_noncyc_paths.append(curr_path[:])
            else:
                for next in u.out_neighbors():
                    if not visited[next]:
                        dfs(next, v)
            curr_path = curr_path[:-1]
            visited[u] = False
            return
        dfs(src, tgt)
        return all_noncyc_paths
            
    # for any source and target node pair in the graph, the tranversing non-cycle/simple path
    # would only include the linear&cycle intersection nodes.

    noncyc_nodes = set()
    simple_paths = []

    srcs = []
    tgts = []
    for node in simp_node_dict.values():
        if node.in_degree() == 0 and node.out_degree() == 0:
            print_vertex(graph, node, "isolated node")
        elif node.in_degree() == 0:
            srcs.append(node)
        elif node.out_degree() == 0:
            tgts.append(node)
        else:
            None
    if srcs == [] or tgts == []:
        print("no src or tgt be found on the graph")
        return None, None

    for src in srcs:
        for tgt in tgts:
            paths = noncyc_path(graph, src, tgt)
            for p in paths:
                [noncyc_nodes.add(graph.vp.id[n]) for n in p]
                simple_paths.append([graph.vp.id[n] for n in p])

    print(list_to_string(list(noncyc_nodes), "non-cyclic+intersection ids"))
    [print(p) for p in simple_paths]

    global curr_path
    del curr_path
    
    # partitate into linear graph
    noncyc_graph = graph.copy()
    simp_node_dict_noncyc, simp_edge_dict_noncyc = graph_to_dict(noncyc_graph)
    graph_color_other_to_gray(noncyc_graph, simp_node_dict_noncyc, noncyc_nodes)
    graph_to_gfa(noncyc_graph, simp_node_dict_noncyc, simp_edge_dict_noncyc, "{0}nc_graph_L2p.gfa".format(temp_dir))

    return noncyc_nodes, simple_paths

def tip_removal(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_dir, overlap):
    """
    retrieve all the source/tail simple path, and merge them into adjacent neighbor path if possible
    
    the collapse step can be done before node depeth rebalance, since it only regards to
    matching score within node seq len

    if is the case, then spades contig may also be modified.
    """
    def remove_tip(graph: Graph, simp_node_dict: dict, from_node, to_path):
        """
        collapse the node with the given path, increment given path depth, remove related information
        about the node.
        """
        graph.vp.color[from_node] = 'gray'
        pending_dp = graph.vp.dp[from_node]
        for node in to_path:
            graph.vp.dp[node] += pending_dp
        simp_node_dict.pop(graph.vp.id[from_node])
        for e in from_node.all_edges():
            graph.ep.color[e] = 'gray'
        print(path_to_id_string(graph, to_path, "Tip Node {0} collapsed to path".format(graph.vp.id[from_node])))
        return

    def cand_collapse_path(graph: Graph, from_node, to_paths, temp_dir, overlap, accept_rate = 0.95):
        """
        use minimap2 -c to evaluation the node-path similarity, sort based on matching score in DESC order
        
        return: the most similar path if there exist a path with score >= accept rate, else return None
        """
        ref_loc = "{0}ref.fa".format(temp_dir)
        query_loc = "{0}query.fa".format(temp_dir)
        overlap_loc = "{0}overlap.paf".format(temp_dir)
        subprocess.check_call('touch {0}; touch {1}'.format(ref_loc, query_loc), shell=True)
        
        id_path_dict = {}
        for id, path in list(enumerate(to_paths)):
            id_path_dict[id] = path

        # retrieve all the path information and save into ref.fa
        with open(ref_loc, 'w') as ref_file:
            for id, path in id_path_dict.items():
                name = ">" + str(id) + "\n"
                seq = path_to_seq(graph, path, id, overlap) + "\n"
                ref_file.write(name)
                ref_file.write(seq)
            ref_file.close()

        # save from node info to query.fa
        with open(query_loc, 'w') as query_file:
            name = ">" + graph.vp.id[from_node] + "\n"
            seq = path_to_seq(graph, [from_node], name, overlap) + "\n"
            query_file.write(name)
            query_file.write(seq)
            query_file.close()

        # minimap to obtain matching score for all node-path
        id_evalscore = {}
        minimap_api(ref_loc, query_loc, overlap_loc)
        with open(overlap_loc, 'r') as overlap_file:
            for Line in overlap_file:
                splited = Line.split('\t')
                path_no = int(splited[5])
                nmatch = int(splited[9])
                nblock = int(splited[10])
                if path_no not in id_evalscore:
                    id_evalscore[path_no] = [nmatch/nblock]
                else:
                    id_evalscore[path_no].append(nmatch/nblock)
            overlap_file.close()
        
        # remove temp file
        subprocess.check_call('rm {0}; rm {1}; rm {2}'.format(ref_loc, query_loc, overlap_loc), shell=True)
        
        id_evalscore_sum = []
        for id, scores in id_evalscore.items():
            mean_score = numpy.mean(scores) if len(scores) != 0 else 0
            id_evalscore_sum.append((id, mean_score))
        
        best_match = sorted(id_evalscore_sum, key=lambda t: t[1], reverse=True)
        print("Tip Node: ", graph.vp.id[from_node], best_match)
        if len(best_match) == 0:
            return None
        elif best_match[0][1] >= accept_rate:
            return id_path_dict[best_match[0][0]]
        else:
            return None
    
    is_removed = True
    # get all the source simple path
    src_nodes = []
    tgt_nodes = []
    isolated_node = []
    for node in simp_node_dict.values():
        if node.in_degree() + node.out_degree() == 0:
            print_vertex(graph, node, "isolated node")
            isolated_node.append(node)
        elif node.in_degree() == 0:
            src_nodes.append(node)
        elif node.out_degree() == 0:
            tgt_nodes.append(node) 
        else:
            None
    
    # src node collapse
    src_nodes = sorted(src_nodes, key=lambda x: graph.vp.dp[x])
    for src in src_nodes:
        print("--------------------------src: {0} --------------".format(graph.vp.id[src]))
        src_len = path_len(graph, [src], overlap)
        potential_paths = []
        # path retrieve
        for out_branch in src.out_neighbors():
            if graph.vp.id[out_branch] not in simp_node_dict:
                continue
            print("current out branch: ", graph.vp.id[out_branch])
            for in_tgt in out_branch.in_neighbors():
                if graph.vp.id[in_tgt] == graph.vp.id[src]:
                    # coincidence path
                    continue
                if graph.vp.id[in_tgt] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                print("current in tgt: ", graph.vp.id[in_tgt])
                potential_paths.extend(paths_to_tgt(graph, simp_node_dict, src, in_tgt, overlap, src_len))
        cand_path = cand_collapse_path(graph, src, potential_paths, temp_dir, overlap)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, src, cand_path)
            is_removed = False
        else:
            print_vertex(graph, src, "Tip cannot be removed, no matching path")

    # target node collapse
    tgt_nodes = sorted(tgt_nodes, key=lambda x: graph.vp.dp[x])
    for tgt in tgt_nodes:
        print("--------------------------tgt: {0} --------------".format(graph.vp.id[tgt]))
        tgt_len = path_len(graph, [tgt], overlap)
        potential_paths = []
        # path retrieve
        for in_branch in tgt.in_neighbors():
            if graph.vp.id[in_branch] not in simp_node_dict:
                continue
            print("current in branch: ", graph.vp.id[in_branch])
            for out_src in in_branch.out_neighbors():
                if graph.vp.id[out_src] == graph.vp.id[tgt]:
                    # coincidence path
                    continue
                if graph.vp.id[out_src] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                print("current out src: ", graph.vp.id[out_src])
                potential_paths.extend(paths_from_src(graph, simp_node_dict, tgt, out_src, overlap, tgt_len))
        cand_path = cand_collapse_path(graph, tgt, potential_paths, temp_dir, overlap)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, tgt, cand_path)
            is_removed = False
        else:
            print_vertex(graph, tgt, "Tip cannot be removed, no matching path")

    return is_removed

def paths_to_tgt(graph: Graph, simp_node_dict: dict, self_node, tgt, overlap, maxlen):
    """
    retrieve all the path from any node to tgt node
    within maxlen restriction, in reverse direction
    """
    def dfs_rev(graph: Graph, v, curr_path: list, maxlen, visited, all_path):
        visited[v] = True
        curr_path.insert(0, v)
        curr_len = path_len(graph, curr_path, overlap)
        if curr_len >= maxlen:
            all_path.append(list(curr_path))
        else:
            for u in v.in_neighbors():
                if not visited[u]:
                    dfs_rev(graph, u, curr_path, maxlen, visited, all_path)
        curr_path.pop(0)
        visited[v] = False
        return
    visited = {}
    for u in graph.vertices():
        if graph.vp.id[u] not in simp_node_dict:
            visited[u] = True
        else:
            visited[u] = False
    visited[self_node] = True
    all_path = []
    dfs_rev(graph, tgt, [], maxlen, visited, all_path)   
    return all_path

def paths_from_src(graph: Graph, simp_node_dict: dict, self_node, src, overlap, maxlen):
    """
    retrieve all the path from src node to any node 
    within maxlen restriction, in straight direction
    """
    def dfs_rev(graph: Graph, u, curr_path: list, maxlen, visited, all_path):
        visited[u] = True
        curr_path.append(u)
        curr_len = path_len(graph, curr_path, overlap)
        if curr_len >= maxlen:
            all_path.append(list(curr_path))
        else:
            for v in u.out_neighbors():
                if not visited[v]:
                    dfs_rev(graph, v, curr_path, maxlen, visited, all_path)
        curr_path.pop(-1)
        visited[u] = False
        return
    visited = {}
    for u in graph.vertices():
        if graph.vp.id[u] not in simp_node_dict:
            visited[u] = True
        else:
            visited[u] = False
    visited[self_node] = True
    all_path = []
    dfs_rev(graph, src, [], maxlen, visited, all_path)
    return all_path

def simp_path_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
    """
    reduce all the contig to a single node, and keep all the potential src/tgt edge.

    1. reduce the coverage for each involving node by the amount of contig cov
    2. reconnect end-to-end nodes to the contig node
    """
    graph_backup = graph.copy()
    simp_node_dict_backup = simp_node_dict.copy()

    contig_node_dict = {}
    contig_info = []
    # reduce all the contig to a single node from the graph
    for cno, (contig, clen, ccov) in list(contig_dict.items()):
        src = contig[0]
        tgt = contig[-1]
        id = src + "00" + cno + "00" + tgt
        cseq = path_to_seq(graph_backup, [simp_node_dict_backup[n] for n in contig], cno, overlap)
        kc = numpy.median([graph_backup.vp.kc[simp_node_dict_backup[u]] for u in contig])
        in_edges = list((graph_backup.vp.id[e.source()], src) for e in simp_node_dict_backup[src].in_edges())
        out_edges = list((tgt, graph_backup.vp.id[e.target()],) for e in simp_node_dict_backup[tgt].out_edges())
        

        for i in range(len(contig)):
            no = contig[i]
            popnode = simp_node_dict.pop(no)
            graph.vp.color[popnode] = 'gray'
            if i != len(contig) - 1:
                e = simp_edge_dict.pop((contig[i], contig[i+1]))
                graph.ep.color[e] = 'gray'

        cv = graph.add_vertex()
        graph.vp.seq[cv] = cseq
        graph.vp.dp[cv] = ccov
        graph.vp.kc[cv] = kc
        graph.vp.id[cv] = id
        graph.vp.color[cv] = 'black'
        simp_node_dict[id] = cv
        contig_node_dict[cno] = id

        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])
    
    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u,v) in in_edges:
            # print("Previous concat: ", (u,v))
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                ue = graph.add_edge(simp_node_dict[u], node)
                graph.ep.overlap[ue] = overlap
                graph.ep.color[ue] = 'black'
                simp_edge_dict[(u, graph.vp.id[node])] = ue
                # print_edge(graph, ue, "reappend edge")
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict:
                    ue = graph.add_edge(in_node, node)
                    graph.ep.overlap[ue] = overlap
                    graph.ep.color[ue] = 'black'
                    simp_edge_dict[(graph.vp.id[in_node], graph.vp.id[node])] = ue
                    # print_edge(graph, ue, "reappend edge")             

        for (u,v) in out_edges:
            # print("Previous concat: ", (u,v))
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                ve = graph.add_edge(node, simp_node_dict[v])
                graph.ep.overlap[ve] = overlap
                graph.ep.color[ve] = 'black'
                simp_edge_dict[(graph.vp.id[node], v)] = ve
                # print_edge(graph, ve, "reappend edge")
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict:
                    ve = graph.add_edge(node, out_node)
                    graph.ep.overlap[ve] = overlap
                    graph.ep.color[ve] = 'black'
                    simp_edge_dict[(graph.vp.id[node], graph.vp.id[out_node])] = ve
                    # print_edge(graph, ve, "reappend edge")
    return contig_node_dict

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strict=False):
    def expectation_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
        for (u,v),e in simp_edge_dict.items():

            if graph.ep.color[e] != 'black':
                #forbidden edge
                continue

            u_node = simp_node_dict[u]
            u_node_out_degree = len([e for e in u_node.out_edges() if graph.ep.color[e] == 'black'])
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])

            v_node = simp_node_dict[v]
            v_node_in_degree = len([e for e in v_node.in_edges() if graph.ep.color[e] == 'black'])
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])

            flow = 0
            if u_node_out_degree == 1 and v_node_in_degree == 1:
                flow = max(graph.vp.dp[u_node], graph.vp.dp[v_node])
            elif u_node_out_degree > 1 and v_node_in_degree == 1:
                flow = max(graph.vp.dp[v_node], (graph.vp.dp[v_node]/u_out_sum)*graph.vp.dp[u_node])
            elif u_node_out_degree == 1 and v_node_in_degree > 1:
                flow = max(graph.vp.dp[u_node], (graph.vp.dp[u_node]/v_in_sum)*graph.vp.dp[v_node])
            else:
                flow = max((graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node])
            graph.ep.flow[e] = max(graph.ep.flow[e], flow)
            # print_edge(graph, e, "{0}".format(flow))
        return

    def maximization_node_depth(graph: Graph, simp_node_dict: dict):
        """
        return True if any updates happened
        """
        is_update = False
        for no, node in simp_node_dict.items():
            us = [e.source() for e in node.in_edges() if graph.ep.color[e] == 'black']
            u_out_degrees = numpy.sum([(len([e for e in u_node.out_edges() if graph.ep.color[e] == 'black'])) for u_node in us]) 
            node_in_degree = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
            in_neighbor_dp_sum = -1
            if u_out_degrees == node_in_degree:
                in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])
            
            vs = [e.target() for e in node.out_edges() if graph.ep.color[e] == 'black']
            v_in_degrees = numpy.sum([(len([e for e in v_node.in_edges() if graph.ep.color[e] == 'black'])) for v_node in vs]) 
            node_out_degree = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
            out_neighbor_dp_sum = -1
            if v_in_degrees == node_out_degree:
                out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])

            curr_dp = graph.vp.dp[node]
            inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
            outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
            graph.vp.dp[node] = numpy.max([curr_dp, inflow, outflow, in_neighbor_dp_sum, out_neighbor_dp_sum])
            # print_vertex(graph, node, "prev dp: {0}".format(curr_dp))
            # print("inflow: ", inflow, "outflow: ", outflow, "in n sum: ", in_neighbor_dp_sum, " out n sum: ", out_neighbor_dp_sum)
            if curr_dp != graph.vp.dp[node]:
                is_update = True
        return is_update
    
    # set cutoff delta
    if graph_is_DAG(graph, simp_node_dict):
        cutoff = 0.0001 if strict else 0.0001 * len(simp_node_dict)
    else:
        cutoff = 0.0001 * len(simp_node_dict)
    
    # store previous node depth
    prev_dp_dict = {}
    for no, v in simp_node_dict.items():
        prev_dp_dict[no] = graph.vp.dp[v]
    print("cutoff delta: ", cutoff)

    # EM optimization
    sum_delta = 0
    is_update = True
    while is_update:
        # E Step
        expectation_edge_flow(graph, simp_node_dict, simp_edge_dict)
        # Validation
        sum_delta = 0.0
        for no, node in simp_node_dict.items():
            inflow = round(numpy.sum([graph.ep.flow[e] for e in node.in_edges() if graph.ep.color[e] == 'black']), 4)
            node_in_degree = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
            outflow = round(numpy.sum([graph.ep.flow[e] for e in node.out_edges() if graph.ep.color[e] == 'black']), 4)
            node_out_degree = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
            if node_in_degree == 0 or node_out_degree == 0:
                continue
            else:
                dominator = (inflow + outflow)/ 2
                if dominator != 0.0:
                    sum_delta += round((abs(inflow - outflow))/dominator, 4)
        # print(sum_delta)
        if sum_delta < cutoff:
            break
        # M Step
        is_update = maximization_node_depth(graph, simp_node_dict)

    # final evaluation
    ratios = [(graph.vp.dp[u] / prev_dp_dict[no]) for no, u in simp_node_dict.items()]
    sum_depth_before = numpy.sum([dp for dp in prev_dp_dict.values()])
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Sum Ratio: ", sum_ratio, "Ave Ratio: ", numpy.mean(ratios), "Max Ratio: ", numpy.max(ratios), "Min Ratio: ", numpy.min(ratios), "Delta: ", sum_delta)

    ratio_div = sum_ratio if numpy.min(ratios) != 1.0 else 1.0
    print("selected ratio: ", ratio_div)
    for node in simp_node_dict.values():
        graph.vp.dp[node] = graph.vp.dp[node] / ratio_div

    node_ratio_dict = {}
    for no in prev_dp_dict.keys():
        if prev_dp_dict[no] != 0:
            node_ratio_dict[no] = (graph.vp.dp[simp_node_dict[no]] / prev_dp_dict[no]) if prev_dp_dict[no] != 0 else 0
    curr_dp_dict = {}
    for no, v in simp_node_dict.items():
        curr_dp_dict[no] = graph.vp.dp[v]

    return prev_dp_dict, curr_dp_dict, node_ratio_dict, ratio_div

def allowed_concat_init(graph: Graph, contig_dict: dict, simp_node_dict: dict, max_len, threshold, overlap):
    """
    Decide whether any contig pair should be connected
    """
    self_concat_off = graph_is_DAG(graph, simp_node_dict)
    impossible_concat_dict = {}
    graph.vp.prev = graph.new_vertex_property("string", val="")

    sp_path_dict = {}

    for no in contig_dict.keys():
        impossible_concat_dict[no] = set()
    for tail_cno, [tail_contig, tail_clen, tail_ccov] in contig_dict.items():
        for head_cno, [head_contig, head_clen, head_ccov] in contig_dict.items():
            print("---------> tail cno: ", tail_cno, " vs head cno: ", head_cno)
            print("tail cov: {0} - head cov: {1}".format(tail_ccov, head_ccov))

            tail_node = simp_node_dict[contig_dict[tail_cno][0][-1]]
            head_node = simp_node_dict[contig_dict[head_cno][0][0]]

            # in-out degree reachable filter
            if tail_node.out_degree() == 0 or head_node.in_degree() == 0:
                impossible_concat_dict[tail_cno].add(head_cno)
            
            # contig intersection filter
            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len or (set(tail_contig)).intersection(set(head_contig)):
                    impossible_concat_dict[tail_cno].add(head_cno)
            else:
                if self_concat_off:
                    impossible_concat_dict[tail_cno].add(head_cno)
            # Final check
            if head_cno not in impossible_concat_dict[tail_cno]:
                sp, plen = dijkstra_sp(graph, simp_node_dict, tail_node, head_node, min(tail_ccov, head_ccov), threshold, overlap)
                if sp != None:
                    total_len = 0
                    if head_cno == tail_cno:
                        if plen == 0:
                            total_len = head_clen - overlap
                        else:
                            total_len = head_clen + plen - 2*overlap
                        print("total cyclic shortest length: ", total_len)
                    else:
                        if plen == 0:
                            total_len = head_clen + tail_clen - overlap
                        else:
                            total_len = head_clen + tail_clen + plen - 2*overlap
                        print("total linear shortest length: ", total_len)
                    if total_len >= max_len:
                        print("even sp exceed the maxlen: ", max_len)
                        impossible_concat_dict[tail_cno].add(head_cno)
                    else:
                        print("SP length within upper bound max len")
                        sp_path_dict[(tail_cno, head_cno)] = (sp, plen)
                else:
                    print("SP not found, impossible")
                    impossible_concat_dict[tail_cno].add(head_cno)
            else:
                print("impossible, removed")

    all_contig_ids = contig_dict.keys()
    contig_concat_plans = {}
    for key, item in impossible_concat_dict.items():
        contig_concat_plans[key] = all_contig_ids - item
    return contig_concat_plans, sp_path_dict

def contig_clique_graph_build(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, max_len, threshold, overlap):

    def similarity(cov1, cov2):
        return 1 - (abs(cov1-cov2)/(cov1+cov2))
    def similarity_e(e, cliq_graph):
        return similarity(cliq_graph.vp.ccov[e.source()], cliq_graph.vp.ccov[e.target()])
    cliq_graph = Graph(directed=True)
    cliq_graph.vp.cno = cliq_graph.new_vertex_property("string", val="")
    cliq_graph.vp.clen = cliq_graph.new_vertex_property("int32_t")
    cliq_graph.vp.ccov = cliq_graph.new_vertex_property("double")
    cliq_graph.vp.text = cliq_graph.new_vertex_property("string")
    cliq_graph.vp.color = cliq_graph.new_vertex_property("string")
    
    cliq_graph.ep.color = cliq_graph.new_edge_property("string")
    cliq_graph.ep.slen = cliq_graph.new_edge_property("int32_t")
    cliq_graph.ep.text = cliq_graph.new_edge_property("string")

    cliq_node_dict = {}
    cliq_edge_dict = {}
    
    for cno, (contig, clen, ccov) in contig_dict.items():
        contig_node = cliq_graph.add_vertex()
        cliq_graph.vp.cno[contig_node] = cno
        cliq_graph.vp.clen[contig_node] = clen
        cliq_graph.vp.ccov[contig_node] = ccov
        cliq_graph.vp.color[contig_node] = 'black'
        cliq_graph.vp.text[contig_node] = cno + ":" + str(clen) + ":" + str(round(ccov, 2))
        cliq_node_dict[cno] = contig_node

    concat_plan, sp_path_dict = allowed_concat_init(graph, contig_dict, simp_node_dict, max_len, threshold, overlap)
    
    for tail_cno, head_cnos in concat_plan.items():
        print("------> tail cno: ", tail_cno, " can concat with following head cnos: ", head_cnos)
        src_node = cliq_node_dict[tail_cno]
        for head_cno in head_cnos:
            tgt_node = cliq_node_dict[head_cno]
            contig_edge = cliq_graph.add_edge(src_node, tgt_node)
            cliq_graph.ep.slen[contig_edge] = int(sp_path_dict[(tail_cno, head_cno)][1])
            cliq_graph.ep.color[contig_edge] = 'black'
            cliq_graph.ep.text[contig_edge] = str(cliq_graph.ep.slen[contig_edge])
            cliq_edge_dict[(tail_cno, head_cno)] = contig_edge

    # edge reduction
    # ensure every contig have at most 1 in connection and 1 out connection
    if graph_is_DAG(graph, simp_node_dict):
        max_edges = []
        for no, node in cliq_node_dict.items():
            print("current contig evaluating: ", no)
            # in connection
            ine = [e for e in node.in_edges() if cliq_graph.ep.color[e] == 'black']
            if len(ine) > 1:
                ine = sorted(ine, 
                key=lambda e: similarity_e(e, cliq_graph), reverse=True)
                max_sim = similarity_e(ine[0], cliq_graph)
                print("MAX IN EDGE: {0}->{1}, plen: {2}".
                format(cliq_graph.vp.cno[ine[0].source()], cliq_graph.vp.cno[ine[0].target()], cliq_graph.ep.slen[ine[0]]))
                print("MAX SIM: ", max_sim)
                for e in ine:
                    cs = similarity_e(e, cliq_graph)
                    if cs < max_sim:
                        cliq_graph.ep.color[e] = 'gray'
                        print("drop edge {0} -> {1}, sim: {2}, slen: {3}".format(cliq_graph.vp.cno[e.source()], cliq_graph.vp.cno[e.target()], cs, cliq_graph.ep.slen[e]))
                    else:
                        max_edges.append(e)
            # out connection
            oute = [e for e in node.out_edges() if cliq_graph.ep.color[e] == 'black']
            if len(oute) > 1:
                oute = sorted(oute, 
                key=lambda e: similarity_e(e, cliq_graph), reverse=True)
                max_sim = similarity_e(oute[0], cliq_graph)
                print("MAX OUT EDGE: {0}->{1}, plen: {2}".format(cliq_graph.vp.cno[oute[0].source()], cliq_graph.vp.cno[oute[0].target()], cliq_graph.ep.slen[oute[0]]))
                print("MAX SIM: ", max_sim)
                for e in oute:
                    cs = similarity_e(e, cliq_graph)
                    if cs < max_sim:
                        cliq_graph.ep.color[e] = 'gray'    
                        print("drop edge {0} -> {1}, sim: {2}, slen: {3}".format(cliq_graph.vp.cno[e.source()], cliq_graph.vp.cno[e.target()], cs, cliq_graph.ep.slen[e]))     
                    else:
                        max_edges.append(e)

        # recover node in/out connection if edges be removed lead to no connection
        sorted_sim_edges_f = sorted([e for e in cliq_graph.edges() if cliq_graph.ep.color[e] != 'black'], key=lambda e: similarity_e(e, cliq_graph), reverse=True)
        for fe in sorted_sim_edges_f:
            u = fe.source()
            v = fe.target()
            print("current check edge {0} -> {1}".format(cliq_graph.vp.cno[u], cliq_graph.vp.cno[v]))
            uoute = [e for e in u.out_edges() if cliq_graph.ep.color[e] == 'black']
            vine = [e for e in v.in_edges() if cliq_graph.ep.color[e] == 'black']
            if len(uoute) == 0 and len(vine) == 0:
                cliq_graph.ep.color[fe] = 'black'
                print("reappend edge {0} -> {1}, {2}".format(cliq_graph.vp.cno[u], cliq_graph.vp.cno[v], similarity_e(fe, cliq_graph)))
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    
    # graph transitive reduction
    print("total edges: ", len(cliq_edge_dict))
    transitive_graph_reduction(cliq_graph, cliq_node_dict, cliq_edge_dict)
    print("total edges after reduction: ", len([e for e in cliq_edge_dict.values() if cliq_graph.ep.color[e] == 'black']))
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    return cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict

def contig_pairwise_concatenation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
cliq_graph: Graph, cliq_node_dict: dict, cliq_edge_dict: dict, sp_path_dict: dict, 
min_cov, min_len, max_len, overlap, threshold, tempdir):
    """
    1. concat the self cycle and reduce some self cycles if self len < minlen
    2. concat the most confident contig pair in order
    3. if cycle exist within the graph, break the cycle by removing the most weak edge.
    4. gradually concat from source nodes until all the contig be concated.
    """
    # helper functions
    def contig_pair_reduction(cno1, cno2, cand_cov, cand_path, cand_len, cno_mapping: dict):
        """
        1. reduce the graph and cliq graph via founded path and cand_cov, 
        for cliq graph, then duplicate/remove
        the cno1/2 node and merge into a single node with cand_cov, keep 
        all the connections other than 1-2
        """
        # cliq graph reduction
        cnode1 = cliq_node_dict[cno1]
        cliq_graph.vp.ccov[cnode1] -= cand_cov
        cnode2 = cliq_node_dict[cno2]
        cliq_graph.vp.ccov[cnode2] -= cand_cov

        print("L1 node cov after deduction: ", cliq_graph.vp.ccov[cnode1])
        print("L2 node cov after deduction: ", cliq_graph.vp.ccov[cnode2])

        print("merge L1 to L2, keep the L1 in edges and L2 out edges only")
        cno_merged = cno1 + "->" + cno2

        if cno1 in cno_mapping:
            if cno1 in cno_mapping[cno1]:
                cno_mapping[cno1].remove(cno1)
            cno_mapping[cno1].add(cno_merged)
        else:
            cno_mapping[cno1] = {cno_merged}

        if cno2 in cno_mapping:
            if cno2 in cno_mapping[cno2]:
                cno_mapping[cno2].remove(cno2)
            cno_mapping[cno2].add(cno_merged)
        else:
            cno_mapping[cno2] = {cno_merged}

        if cno_merged not in cno_mapping:
            cno_mapping[cno_merged] = {cno_merged}

        prev1 = contig_dict.pop(cno1)
        prev2 = contig_dict.pop(cno2)
        contig_dict[cno_merged] = [prev1[0]+[graph.vp.id[v] for v in cand_path]+prev2[0], cand_len, cand_cov]
        
        cnode_merged = cliq_graph_add_node(cliq_graph, cliq_node_dict, 
        cno_merged, cand_len, cand_cov, 
        cno_merged + ":" + str(cand_len) + ":" + str(cand_cov))

        # removed the used up contig node
        cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno1, cnode1)
        cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno2, cnode2)

        # removed the related L1 edges
        for edge1out in cnode1.out_edges():
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict,
            cliq_graph.vp.cno[edge1out.source()], cliq_graph.vp.cno[edge1out.target()], edge1out)

        # replace the related L1 in edges
        for edge1in in cnode1.in_edges():
            if cliq_graph.ep.color[edge1in] != 'black':
                continue
            src1 = edge1in.source()
            src1cno = cliq_graph.vp.cno[src1]
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            src1cno, cliq_graph.vp.cno[edge1in.target()], edge1in)
            
            # do not add self cycle edges
            if src1cno != cno_merged:
                cliq_graph_add_edge(cliq_graph, cliq_edge_dict, src1cno, src1, 
                cno_merged, cnode_merged, cliq_graph.ep.slen[edge1in], cliq_graph.ep.text[edge1in])
        
        # removed the related L2 in edges
        for edge2in in cnode2.in_edges():
            if cliq_graph.ep.color[edge2in] != 'black':
                continue
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            cliq_graph.vp.cno[edge2in.source()], cliq_graph.vp.cno[edge2in.target()], edge2in)
        
        # replace the related L2 out edges
        for edge2out in cnode2.out_edges():
            if cliq_graph.ep.color[edge2out] != 'black':
                continue
            tgt2 = edge2out.target()
            tgt2cno = cliq_graph.vp.cno[tgt2]
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            cliq_graph.vp.cno[edge2out.source()], tgt2cno, edge2out)
            if cno_merged != tgt2cno:
                cliq_graph_add_edge(cliq_graph, cliq_edge_dict, cno_merged, cnode_merged,
                tgt2cno, tgt2, cliq_graph.ep.slen[edge2out], cliq_graph.ep.text[edge2out])
        return
    def buffer_concatenation(concat_buffer: list, cno_mapping: dict):
        for (cno1, cno2, cov, delta) in concat_buffer:
            print("------------------------------------------------------------------------")
            print("-->Before mapping: {0}: {1} - {2}: {3}".format(cno1, cno_mapping[cno1], cno2, cno_mapping[cno2]))
            for pcno in list(cno_mapping[cno1]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno1].remove(pcno)
                    print("pcno removed: ", pcno)
            for pcno in list(cno_mapping[cno2]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno2].remove(pcno)  
                    print("pcno removed: ", pcno)
            pairs = []
            for cno1m in cno_mapping[cno1]: 
                for cno2m in cno_mapping[cno2]:
                    if (cno1m, cno2m) in cliq_edge_dict:
                        pairs.append((cno1m, cno2m))
            if not pairs:
                print("contig has been used from previous step: {0} {1}".format(cno1, cno2))
                continue

            if cno1 in contig_usage:
                contig_usage[cno1] += 1
            if cno2 in contig_usage:
                contig_usage[cno2] += 1

            cno1m, cno2m = min(pairs, key=lambda p: 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[0]]] - cov, 2) + 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[1]]] - cov, 2))
            
            print("-->PAIR UP {0} - {1}, cov: {2}, diff: {3}".format(cno1m, cno2m, cov, delta))
            src = simp_node_dict[contig_dict[cno1m][0][-1]]
            tgt = simp_node_dict[contig_dict[cno2m][0][0]]
            cand_path, plen = dijkstra_sp(graph, simp_node_dict, src, tgt, cov, threshold, overlap)
            cand_len = get_concat_len(cno1m, contig_dict[cno1m][1], cno2m, contig_dict[cno2m][1], plen, overlap)

            contig_pair_reduction(cno1m, cno2m, cov, cand_path, cand_len, cno_mapping)

        return
    ###################################################################################################        
    contig_usage = {}
    for cno in contig_dict.keys():
        contig_usage[cno] = 0

    # retrieve all the self cycle first
    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                if cliq_graph.vp.clen[contig_node] < min_len:
                    # remove self cycle edge with self cycle + outer connection feature when clen < minlen
                    cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_graph.edge(contig_node, contig_node))
                    print("remove self edge+outer connection {0} -> {0} with node length < minlen".format(cno))
            else:
                print("Circular PATH: ", cno)
                cov = cliq_graph.vp.ccov[contig_node]
                src = simp_node_dict[contig_dict[cno][0][-1]]
                tgt = simp_node_dict[contig_dict[cno][0][0]]
                cand_path, plen = dijkstra_sp(graph, simp_node_dict, src, tgt, cov, threshold, overlap)
                cand_len = get_concat_len(cno, contig_dict[cno][1], cno, contig_dict[cno][1], plen, overlap)

                if cand_path != None and cand_len != None:
                    contig_dict[cno][0].extend([graph.vp.id[n] for n in cand_path])
                    contig_dict[cno][1] = cand_len
                    contig_dict[cno][2] = cov
                    cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_edge_dict[(cno, cno)])
                    print(path_to_id_string(graph, [simp_node_dict[cno] for cno in contig_dict[cno][0]] + cand_path, "cov: {0}".format(cov)))
                else:
                    print("Path not found, error")

    cno_mapping = {}
    for id in cliq_node_dict.keys():
        cno_mapping[id] = {id}

    concat_buffer = []
    for (cno1, cno2) in cliq_edge_dict.keys():
        print("---------------------------------------------------------------")
        delta = abs(cliq_graph.vp.ccov[cliq_node_dict[cno1]] - cliq_graph.vp.ccov[cliq_node_dict[cno2]])
        cov = min(cliq_graph.vp.ccov[cliq_node_dict[cno1]], cliq_graph.vp.ccov[cliq_node_dict[cno2]])
        print("{0} - {1}, cov: {2}, diff: {3}".format(cno1, cno2, cov, delta))
        if delta < threshold:
            concat_buffer.append((cno1, cno2, cov, delta))
    
    concat_buffer = sorted(concat_buffer, key=lambda tuple: tuple[3])
    print("all most confident sorted concats are: ", concat_buffer)
    buffer_concatenation(concat_buffer, cno_mapping)

    for cno, contig in cliq_node_dict.items():
        contig, clen, ccov = contig_dict[cno]
        print_contig(cno, clen, ccov, contig)
    for (u, v), e in cliq_edge_dict.items():
        print("EDGE: ", u, v)

    if not graph_is_DAG(cliq_graph, cliq_node_dict):
        print("------>graph contains cycles, cyclic concatenation processing")
        cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
        cycles = all_circuits(cliq_graph, True)
        for cyc in cycles:
            print("cycle: ", [cliq_graph.vp.cno[n] for n in cyc])
            skip_cycle = False
            all_edges = {}
            for i in range(len(cyc)):
                u = cyc[i]
                v = cyc[(i+1) % len(cyc)]
                e = cliq_graph.edge(u, v)
                if cliq_graph.ep.color[e] != 'black':
                    skip_cycle = True
                    break
                all_edges[(cliq_graph.vp.cno[u],cliq_graph.vp.cno[v])] = cliq_graph.edge(u, v)
            if skip_cycle:
                continue
            removed_edge_tuple = min(all_edges.items(), key=lambda e_tuple: cliq_graph.ep.slen[e_tuple[1]])
            print("Removed edge: ", removed_edge_tuple[0][0], "->", removed_edge_tuple[0][1])
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, removed_edge_tuple[0][0], removed_edge_tuple[0][1], removed_edge_tuple[1])

    print("--------------Start graduate concatentation------------------")
    for cno, cnos in cno_mapping.items():
        print("cno: ", cno, "maps to: ", cnos)
    # concat 
    # process all the non self-cycle contig until no isolated node.
    # direct pair two adjacent contig only if coverage difference is within pairwise threshold
    while True:
        # clean up the cno_mapping
        cno_mapping = {}
        for id in cliq_node_dict.keys():
            cno_mapping[id] = {id}
        # all the src contig
        L1_contigs = {}
        for no, contig_node in list(cliq_node_dict.items()):
            if cliq_graph.vp.color[contig_node] != 'black':
                continue
            ind = len([e for e in contig_node.in_edges() if cliq_graph.ep.color[e] == 'black'])
            outd = len([e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black'])
            if ind == 0 and outd != 0:
                # src node
                L1_contigs[no] = contig_node

        L2_contigs = {}
        st_pairs = {}
        for no, contig_node in list(L1_contigs.items()):
            out_edges = [e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black']
            for out_e in out_edges:
                out_contig_node = out_e.target()
                if cliq_graph.vp.color[out_contig_node] != 'black':
                    continue
                outid = cliq_graph.vp.cno[out_contig_node]
                if not outid in L2_contigs:
                    L2_contigs[outid] = out_contig_node
                st_pairs[(no, outid)] = (contig_node, out_contig_node)
        
        if not L1_contigs or not L2_contigs:
            print("no more linear concatentation left, break")
            break

        print("L1: ", [no for no in L1_contigs.keys()])
        print("L2: ", [no for no in L2_contigs.keys()])
        
        # detect the most confident contig pair first
        concat_buffer = []
        for (cno1, cno2), (node1, node2) in st_pairs.items():
            delta = abs(cliq_graph.vp.ccov[node1] - cliq_graph.vp.ccov[node2])
            concat_buffer.append((cno1, cno2, min(cliq_graph.vp.ccov[node1], cliq_graph.vp.ccov[node2]), delta))
        
        concat_buffer = sorted(concat_buffer, key=lambda tuple: tuple[3])
        print("all most confident sorted concats are: ", concat_buffer)
        buffer_concatenation(concat_buffer, cno_mapping)

    print("--------------> contig after concatentations...")
    for cno, contig_node in list(cliq_node_dict.items()):
        print("--------------------------------------------------------")
        if cno in contig_usage:
            count = contig_usage[cno]
            print("cno: ", cno, "original usage count: ", count)
            if count > 0:
                print("original contig has been used, remove it to reduce the potential duplication ratio")
                cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno, contig_node)
                contig_dict.pop(cno)
                continue
        contig, clen, ccov = contig_dict[cno]
        cflows = contig_flow(graph, simp_edge_dict, contig)
        print_contig(cno, clen, ccov, contig)
        if len(cflows) > 0:
            print("mean: ", numpy.mean(cflows), 
            " median: ", numpy.median(cflows), 
            " min: ", numpy.min(cflows))

    for (u, v), e in cliq_edge_dict.items():
        print("potential uncaught error: EDGE: ", u, v)
    
    # simplify the graph
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), tempdir, "cliq_graphL2.png")

    return contig_dict

def strain_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
min_cov, max_len, threshold, overlap):
    """
    Extend the strain length
    1. map all the strain back into the graph
    """
    # find all the srcs-targets
    src_nodes = [node for node in graph.vertices() if node.in_degree() == 0 and node.out_degree() != 0]
    tgt_nodes = [node for node in graph.vertices() if node.in_degree() != 0 and node.out_degree() == 0]
    
    # create a global source node & sink node, and concat all the curr src to the global source node
    global_src = graph.add_vertex()
    graph.vp.id[global_src] = 'global_src'
    graph.vp.dp[global_src] = 0
    graph.vp.color[global_src] = 'black'
    simp_node_dict[graph.vp.id[global_src]] = global_src
    for src in src_nodes:
        e = graph.add_edge(global_src, src)
        graph.ep.flow[e] = graph.vp.dp[src]
        graph.ep.color[e] = 'black'
        graph.ep.overlap[e] = overlap
        graph.vp.dp[global_src] += graph.ep.flow[e]
    print("srcs: ", list_to_string([graph.vp.id[src] for src in src_nodes]))

    global_sink = graph.add_vertex()
    graph.vp.id[global_sink] = 'global_sink'
    graph.vp.dp[global_sink] = 0
    graph.vp.color[global_sink] = 'black'
    simp_node_dict[graph.vp.id[global_sink]] = global_sink
    for tgt in tgt_nodes:
        e = graph.add_edge(tgt, global_sink)
        graph.ep.flow[e] = graph.vp.dp[tgt]
        graph.ep.color[e] = 'black'
        graph.ep.overlap[e] = overlap
        graph.vp.dp[global_sink] += graph.ep.flow[e]
    print("sinks: ", list_to_string([graph.vp.id[tgt] for tgt in tgt_nodes]))

    # extend all the contig from both end
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        print_contig(cno, clen, ccov, contig, "-----> current extending contig: ")
        if simp_node_dict[contig[0]].in_degree() == 0 and simp_node_dict[contig[-1]].out_degree() == 0:
            print("not extensible")
        elif simp_node_dict[contig[0]] in list(simp_node_dict[contig[-1]].out_neighbors()):
            print("self cycle, not entensible")
        else:
            print("extensible")
            contig_head = simp_node_dict[contig[0]]
            if contig_head.in_degree() != 0 and global_src not in contig_head.in_neighbors():
                print("---> {0} ~~> contig head: {1}".format(graph.vp.id[global_src], graph.vp.id[contig_head]))
                sp, _ = dijkstra_sp(graph, simp_node_dict, global_src, contig_head, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [[graph.vp.id[n] for n in sp] + contig_dict[cno][0], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")

            contig_tail = simp_node_dict[contig[-1]]  
            if contig_tail.out_degree() != 0 and global_sink not in contig_tail.out_neighbors():
                print("---> contig tail: {0} ~~> {1}".format(graph.vp.id[contig_tail], graph.vp.id[global_sink]))
                sp, _ = dijkstra_sp(graph, simp_node_dict, contig_tail, global_sink, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [contig_dict[cno][0] + [graph.vp.id[n] for n in sp], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")
    
    # udp, node depth with related contig cov be deducted
    graph.vp.udp = graph.new_vertex_property("double")
    for no, node in simp_node_dict.items():
        graph.vp.udp[node] = graph.vp.dp[node]
    print("------------------> all the extended contigs: ")
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        print_contig(cno, clen, ccov, contig, "-----> extended contig ")
        call = [simp_node_dict[n] for n in contig]
        if global_src in simp_node_dict[contig[0]].in_neighbors():
            call.insert(0, global_src)
        if global_sink in simp_node_dict[contig[-1]].out_neighbors():
            call.append(global_sink)
        graph_reduction_c(graph, call, ccov, threshold)

    for v in graph.vertices():
        if graph.vp.udp[v] < threshold:
            graph.vp.color[v] = 'gray'
    for e in graph.edges():
        if graph.ep.flow[e] < threshold:
            graph.ep.color[e] = 'gray'

    partial_used_nodes = []
    free_nodes = []
    overused_nodes = []
    pre_usages = {}
    for no, node in simp_node_dict.items():
        print("----------------------------------------------------")
        if no == 'global_src' or no == 'global_sink':
            continue
        ratio = round(((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2)
        print("Node: {0}, full: {1}, left: {2}, usage: {3}, {4}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), ratio, graph.vp.color[node]))
        if ratio < 100 and graph.vp.udp[node] > threshold:
            partial_used_nodes.append(no)
        if ratio <= 0:
            free_nodes.append(no)
        if ratio > 100 and graph.vp.udp[node] < -threshold:
            overused_nodes.append(no)
        pre_usages[no] = ratio
    overall_pre_usage = numpy.mean(list(pre_usages.values()))
    print("Free nodes: ", list_to_string(free_nodes))
    print("Partial used nodes: ", list_to_string(partial_used_nodes))
    print("Usage: ", overall_pre_usage)

    # analytics part
    global TEMP_DIR
    seaborn.set_theme(style="darkgrid")
    plt.figure(figsize=(128,64))

    graph.vp.color[global_sink] = 'gray'
    graph.vp.color[global_src] = 'gray'
    dps_store = {}
    for node in graph.vertices():
        dps_store[graph.vp.id[node]] = graph.vp.dp[node]
        graph.vp.dp[node] = graph.vp.udp[node]
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}ecbsdt_graph_L6.gfa".format(TEMP_DIR))
    graph.vp.color[global_sink] = 'black'
    graph.vp.color[global_src] = 'black'
    for node in graph.vertices():
        graph.vp.dp[node] = dps_store[graph.vp.id[node]]
    df = pandas.DataFrame(
        {'Id': pre_usages.keys(),
        'pre': pre_usages.values(),
        })
    ax = seaborn.barplot(x='Id', y='pre', data=df)
    for container in ax.containers:
        ax.bar_label(container)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    plt.title('Bar plot of Node Usage')
    plt.savefig("{0}barplot_usage.png".format(TEMP_DIR))

    return contig_dict

def final_strain_extraction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold, overlap):
    # st path FIXME later
    global_src = simp_node_dict['global_src']
    global_sink = simp_node_dict['global_sink']
    print("Overall usage below bound, start st path finding")
    for src in global_src.out_neighbors():
        print("global connect src: ", graph.vp.id[src], "udp: ", graph.vp.udp[src], graph.vp.color[src])
    for tgt in global_sink.in_neighbors():
        print("global connect tgt: ", graph.vp.id[tgt], "udp: ", graph.vp.udp[tgt], graph.vp.color[tgt])
    path = []
    curr_id = 0
    while path != None:
        path, plen = dijkstra_sp_v2(graph, simp_node_dict, global_src, global_sink, overlap)
        if path != None:
            # FIXME
            cflows = contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in path[1:-1]])
            # stat
            print("mincov: ", numpy.min(cflows), 
            "meancov: ", numpy.mean(cflows), 
            "mediancov: ", numpy.median(cflows),
            "maxcov: ", numpy.max(cflows))
            for i in range(len(cflows)):
                print("u: ", graph.vp.id[path[i+1]], "v: ", graph.vp.id[path[i+2]], "flow: ", round(cflows[i],2))
            
            d = pandas.Series({"coverage": cflows})
            ax = seaborn.countplot(x="coverage", data=d)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
            plt.title('Count plot of Path {0}'.format("st" + str(curr_id)))
            plt.savefig("{0}countplot_{1}.png".format(TEMP_DIR, "st" + str(curr_id)))

            # end stat
            redcov = numpy.min(cflows)
            graph_reduction_c(graph, path, redcov, threshold)
            for v in graph.vertices():
                if graph.vp.udp[v] < threshold:
                    graph.vp.color[v] = 'gray'
            for e in graph.edges():
                if graph.ep.flow[e] < threshold:
                    graph.ep.color[e] = 'gray'
            contig_dict["st" + str(curr_id)] = [[graph.vp.id[n] for n in path[1:-1]], plen, redcov]
            curr_id += 1
        else:
            print("Path not found")

    # double check usage
    post_usages = {}
    post_free_nodes = []
    post_partial_used_nodes = []
    for no, node in simp_node_dict.items():
        print("----------------------------------------------------")
        if no == 'global_src' or no == 'global_sink':
            continue
        ratio = round(((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2)
        print("Node: {0}, full: {1}, left: {2}, usage: {3}, {4}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), ratio, graph.vp.color[node]))
        if ratio < 100 and graph.vp.udp[node] > threshold:
            post_partial_used_nodes.append(no)
        if ratio <= 0:
            post_free_nodes.append(no)
        post_usages[no] = ratio
    overall_post_usage = numpy.mean(list(post_usages.values()))
    print("Free nodes: ", list_to_string(post_free_nodes))
    print("Partial used nodes: ", list_to_string(post_partial_used_nodes))
    print("Overall Usage Post: {0}".format(overall_post_usage))
    # do some plotting
    # df = pandas.DataFrame(
    #     {'Id': [i for i in range(len(pre_usages.keys()))],
    #     'pre': pre_usages.values(),
    #     'post': post_usages.values()
    #     })
    # tidy = df.melt(id_vars='Id').rename(columns=str.title)
    # ax = seaborn.barplot(x='Id', y='Value', hue='Variable', data=tidy)
    # for container in ax.containers:
    #     ax.bar_label(container)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    # plt.title('Bar plot of Node Usage')
    # plt.savefig("{0}barplot_usage.png".format(TEMP_DIR))
    # post process
    for node in graph.vertices():
        graph.vp.dp[node] = graph.vp.udp[node]
    simp_node_dict.pop(graph.vp.id[global_src])
    graph.vp.color[global_src] = 'gray'
    simp_node_dict.pop(graph.vp.id[global_sink])
    graph.vp.color[global_sink] = 'gray'

def local_search_optimisation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, min_cov, min_len, max_len, overlap, threshold):
    return

if __name__ == "__main__":
    main()