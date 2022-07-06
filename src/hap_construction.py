#!/usr/bin/env python3

import sys, os
import subprocess
import argparse

from graph_tool.topology import all_circuits
from graph_tool.draw import graph_draw
from graph_tool.all import Graph

import numpy
import heapq
import matplotlib.pyplot as plt
import seaborn
import pandas

from collections import deque
from itertools import product

from graph_converter import *
from search_algos import *

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = False
TEMP_DIR = "acc/"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, default=100, help=("minimum coverage for strains"))
    parser.add_argument('-minlen', '--minimum_strain_length', dest='min_len', default=2500, type=int, help=("minimum strain length"))
    parser.add_argument('-maxlen', '--maximum_strain_length', dest='max_len', default=12000, type=int, help=("maximum strain length"))
    parser.add_argument('-overlap', '--vertex_overlap', dest='overlap', default=127, type=int, help=("adjacent vertex overlap in the graph"))
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, DEBUG_MODE only')
    parser.add_argument('-o', '--output_dir', dest='output_dir', default='acc/', type=str, help='output directory (default: acc/)')
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')
    args = parser.parse_args()
    global TEMP_DIR
    TEMP_DIR = args.output_dir
    
    subprocess.check_call("rm -rf {0} && mkdir {0}".format(TEMP_DIR), shell=True)

    print("----------------------------------INPUT---------------------------------------")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    contig_dict, _, _ = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_len)
    
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}graph_L0.gfa".format(TEMP_DIR))
    print("----------------------------------TIP REMOVAL---------------------------------------")
    # apply tip removal for cyclic graph.
    if not graph_is_DAG(graph0, simp_node_dict0):
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(graph0, simp_node_dict0, simp_edge_dict0, TEMP_DIR, args.overlap)

        contig_dict_fix(graph0, simp_node_dict0, contig_dict, args.overlap)
    else:
        print("Graph is DAG, tip removal skipped.")
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1.gfa".format(TEMP_DIR))

    print("-------------------------------ACYCLIC TRANSLATION-----------------------------------")
    #FIXME
    removed_edges = cyclic_to_dag(graph1, simp_node_dict1, contig_dict, args.overlap)
    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}t_graph_L1_dag.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1_dag.gfa".format(TEMP_DIR))

    print("-------------------------------GRAPH SIMPLIFICATION AND REBALANCE-----------------------------------")
    #FIXME
    mediandp = numpy.median([graph1.vp.dp[node] for node in simp_node_dict1.values()])
    THRESHOLD = mediandp/20
    # 20000 * 0.001 
    # numpy.quantile([graph1.vp.dp[node] for node in graph1.vertices()], 0.05)
    # mediandp/20
    print("MEDIAN NODE DEPTH: ", mediandp, "threshold: ", THRESHOLD)
    graph_simplification(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, THRESHOLD)

    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}dt_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}dt_graph_L2.gfa".format(TEMP_DIR))
    coverage_rebalance(graph2, simp_node_dict2, simp_edge_dict2, True)

    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)
    
    print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict)
    draw_edgeflow(graph3, simp_edge_dict3, TEMP_DIR, 'Bar plot of Edge flow', 'barplot_edge_flow.png')
    
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict3, "{0}graph_L0.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph3, contig_dict, simp_node_dict3, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))

    print("-----------------------GRAPH BRANCH SPLIT & COMPACTIFICATION-------------------------------")
    grapha = graph3
    simp_node_dicta = simp_node_dict3
    simp_edge_dicta = simp_edge_dict3
    total_removed_branch_nt = 0
    total_removed_branch_t = 0
    iterCount = 'A'
    num_split = 0
    trivial_split_count = 0
    while True:
        # trivial branch split
        prev_ids = list(simp_node_dicta.keys())
        trivial_split_count, id_mapping = graph_split_trivial(grapha, simp_node_dicta, simp_edge_dicta)

        graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)

        contig_dict_resol(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, id_mapping, prev_ids, args.overlap)

        # non-trivial branch split
        num_split = graph_splitting(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, THRESHOLD)
        graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))

        simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, args.overlap)

        graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        grapha, simp_node_dicta, simp_edge_dicta = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(grapha, simp_node_dicta, simp_edge_dicta)

        contig_dict_simp(grapha, simp_edge_dicta, contig_dict)
        trim_contig_dict(grapha, simp_node_dicta, contig_dict, args.overlap)
        contig_cov_fix(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            iterCount = chr(ord(iterCount) + 1)
        else:
            coverage_rebalance(grapha, simp_node_dicta, simp_edge_dicta, True)
            graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
            break
    print("Total non-trivial branches removed: ", total_removed_branch_nt, " total trivial branches removed: ", total_removed_branch_t)

    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}post_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph5, contig_dict, simp_node_dict5, args.overlap, "{0}post_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}post_contigs.fasta".format(TEMP_DIR), "{0}post_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}post_contigs_to_strain.paf".format(TEMP_DIR))

    strain_dict = extract_cand_path2(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.overlap, THRESHOLD)
    
    contig_dict_simp(graph5, simp_edge_dict5, strain_dict)    
    contig_dict_to_fasta(graph5, strain_dict, simp_node_dict5, args.overlap, "{0}final_contig.fasta".format(TEMP_DIR))
    contig_dict_to_path(strain_dict, "{0}final_contig.paths".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}final_contig.fasta".format(TEMP_DIR), "{0}final_contig_to_strain.paf".format(TEMP_DIR))
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
    cutoff = 0.00001 * len(simp_node_dict) if strict else 0.0001 * len(simp_node_dict)
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
        sum_dp = numpy.sum([graph.vp.dp[node] for node in graph.vertices()])
        dom_delta = 0.0
        deltas = []
        for no, node in simp_node_dict.items():
            node_in_degree = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
            node_out_degree = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
            if node_in_degree == 0:
                inflow = graph.vp.dp[node]
            else:
                inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges() if graph.ep.color[e] == 'black'])
            
            if node_out_degree == 0:
                outflow = graph.vp.dp[node]
            else:
                outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges() if graph.ep.color[e] == 'black'])

            dominator = (inflow + outflow)/ 2
            if dominator != 0.0:
                sum_delta += abs(inflow - outflow)
                dom_delta += dominator
                deltas.append((no, (graph.vp.dp[node] * abs(inflow - outflow))/dominator))

        if strict:
            sum_delta = sum([k[1] for k in deltas]) / sum_dp
        else:
            sum_delta = sum_delta / dom_delta
        print("sum delta: ", sum_delta, "worst delta: ", sorted(deltas, key=lambda p: p[1], reverse=True)[:10])
        if sum_delta < cutoff:
            break
        # M Step
        is_update = maximization_node_depth(graph, simp_node_dict)

    # final evaluation
    ratios = [(graph.vp.dp[u] / prev_dp_dict[no]) for no, u in simp_node_dict.items()]
    sum_depth_before = numpy.sum([dp for dp in prev_dp_dict.values()])
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Sum Ratio: ", sum_ratio, "Ave Ratio: ", numpy.mean(ratios), "Max Ratio: ", numpy.max(ratios), "Min Ratio: ", numpy.min(ratios), "Delta: ", sum_delta)

    ratio_div = sum_ratio
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
    1. two contig cannot concat if:
        1. two contig have intermediate parallel sharing nodes
        2. no path exist between two contig
    """
    self_concat_off = graph_is_DAG(graph, simp_node_dict)
    impossible_concat_dict = {}
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
            intersects = None
            cend = None
            status = None
            if tail_cno != head_cno:
                isParallel, intersects, cend, status = check_contig_intersection(tail_cno, tail_contig, head_cno, head_contig)
                if isParallel:
                    impossible_concat_dict[tail_cno].add(head_cno)
                if status != None:                 
                    if status != 'n' and status != 'o':
                        if status != 'b':
                            # end-to-end overlap case, forward direction
                            sp_path_dict[(tail_cno, head_cno)] = None
                            continue
                        else:
                            impossible_concat_dict[tail_cno].add(head_cno)
            else:
                if self_concat_off:
                    impossible_concat_dict[tail_cno].add(head_cno)
            if head_cno not in impossible_concat_dict[tail_cno]:
                total_len = 0
                if intersects != None and cend != None:
                    # approx overlap length
                    intersect_len = sum([len(graph.vp.seq[simp_node_dict[id]]) for id in intersects]) - overlap * (len(intersects) + cend)
                    total_len = -intersect_len
                sp, plen, pmark = dijkstra_sp(graph, tail_node, head_node, min(tail_ccov, head_ccov), threshold, overlap)
                if sp != None:
                    if head_cno == tail_cno:
                        if plen == 0:
                            total_len += head_clen - overlap
                        else:
                            total_len += head_clen + plen - 2*overlap
                        print("total cyclic shortest length: ", total_len)
                    else:
                        if plen == 0:
                            total_len += head_clen + tail_clen - overlap
                        else:
                            total_len += head_clen + tail_clen + plen - 2*overlap
                        print("total linear shortest length: ", total_len)
                    if total_len >= max_len:
                        print("even sp exceed the maxlen: ", max_len)
                        impossible_concat_dict[tail_cno].add(head_cno)
                    else:
                        print("SP length within upper bound max len")
                        sp_path_dict[(tail_cno, head_cno)] = (sp, plen, pmark)
                else:
                    print("SP not found, impossible")
                    impossible_concat_dict[tail_cno].add(head_cno)

    all_contig_ids = contig_dict.keys()
    contig_concat_plans = {}
    for key, item in impossible_concat_dict.items():
        contig_concat_plans[key] = all_contig_ids - item
    return contig_concat_plans, sp_path_dict

def contig_clique_graph_build(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, max_len, threshold, overlap):
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
    
    cno2cno_adjMtx_csv = []
    cno2cno_adjMtx_csv.append([cno+"/"+str(round(contig_dict[cno][2])) for cno in concat_plan.keys()])
    cno2cno_adjMtx_csv[0].insert(0, "X")
    cno2cno_adjMtx = []

    cno_to_index = {}
    cno_to_index_csv = {}
    for i, cno in enumerate(cno2cno_adjMtx_csv[0]):
        if i > 0:
            cno_to_index_csv[cno.split("/")[0]] = i
            cno_to_index[cno.split("/")[0]] = i - 1

    for tail_cno, head_cnos in concat_plan.items():
        print("------> tail cno: ", tail_cno, " can concat with following head cnos: ", head_cnos)
        curr_row_csv = [" " for _ in concat_plan.keys()]
        curr_row_csv.insert(0, tail_cno+"/"+str(round(contig_dict[tail_cno][2])))
        curr_row = [[sys.maxsize,'X'] for _ in concat_plan.keys()]

        src_node = cliq_node_dict[tail_cno]
        for head_cno in head_cnos:
            plen = int(sp_path_dict[(tail_cno, head_cno)][1]) if sp_path_dict[(tail_cno, head_cno)] != None else -1
            abs_cdif = abs(contig_dict[tail_cno][2]-contig_dict[head_cno][2])
            curr_row_csv[cno_to_index_csv[head_cno]] = str(round(abs_cdif)) + ":" + str(plen)
            curr_row[cno_to_index[head_cno]] = [abs_cdif, 'W']
            tgt_node = cliq_node_dict[head_cno]
            contig_edge = cliq_graph.add_edge(src_node, tgt_node)
            cliq_graph.ep.slen[contig_edge] = plen
            cliq_graph.ep.color[contig_edge] = 'black'
            cliq_graph.ep.text[contig_edge] = str(cliq_graph.ep.slen[contig_edge])
            cliq_edge_dict[(tail_cno, head_cno)] = contig_edge
        cno2cno_adjMtx.append(curr_row)
        cno2cno_adjMtx_csv.append(curr_row_csv)
    return cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict, cno2cno_adjMtx, cno_to_index

def clique_graph_clean(cliq_graph: Graph, cliq_node_dict: dict, cliq_edge_dict: dict, 
adj_matrix, cno_to_index: dict, threshold):
    """
    adj matrix, elem in 5 colors: X, gray(invalid), white(has connection), red(candidate connection), blue(fixed connection)
    FIXME not all contig should be concated, check the diff < threshold, overlap contig pair have highest priority
    """

    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                remove_self_cycle = False
                if not remove_self_cycle:
                    for inn in contig_node.in_neighbors():
                        if inn != contig_node and contig_node in (inn.in_neighbors()):
                            print("remove sc set 01, ")
                            remove_self_cycle = True
                            break
                if not remove_self_cycle:
                    for onn in contig_node.out_neighbors():
                        if onn != contig_node and contig_node in (onn.out_neighbors()):
                            print("remove sc set 02, ")
                            remove_self_cycle = True
                            break    
                if remove_self_cycle:
                    # remove self cycle edge with self cycle + outer connection feature
                    adj_matrix[cno_to_index[cno]][cno_to_index[cno]][1] = 'G'
                    print("remove self edge+outer connection {0} -> {0}".format(cno))
                else:
                    print("Keep self cycle since no feedback from neighbor ", cno)
    index_to_cno = {}
    for cno, i in cno_to_index.items():
        index_to_cno[i] = cno

    for rowId, row in enumerate(adj_matrix):
        for colId, [flow, color] in enumerate(row):
            if color == 'W':
                if flow >= threshold:
                    print("manual forbidden concat: ", index_to_cno[rowId], index_to_cno[colId], flow)
                    adj_matrix[rowId][colId][1] = 'G'
                else:
                    print("Potential concat: ", index_to_cno[rowId], index_to_cno[colId], flow)
    has_changes = True
    dim = len(adj_matrix)
    # iterate the adjacent matrix in rowwise and columnwise
    while has_changes:
        print("looping")
        has_changes = False
        # row wise
        for rowId in range(dim):
            row = get_row(adj_matrix, rowId)
            hasFixed = False
            colId = None
            minFlow = sys.maxsize
            for cid, [flow, color] in enumerate(row):
                if color == 'B':
                    # current row fixed
                    hasFixed = True
                    break
                elif color == 'W':
                    if flow < minFlow:
                        colId = cid
                        minFlow = flow
                elif color == 'X' or color == 'G':
                    # invalid, skip
                    None
                else:
                    print("Error: ", rowId, index_to_cno[rowId], cid, index_to_cno[cid])
                    assert color != 'R'
            if not hasFixed:
                if colId != None:
                    # found a minimum block to assign
                    adj_matrix[rowId][colId][1] = 'R'
                    has_changes = True
                else:
                    # no more spot for the row
                    None
        for colId in range(dim):
            col = get_col(adj_matrix, colId)
            # hasFixed
            rowId = None
            minFlow = sys.maxsize
            # if only one red+blue among the column, then assign it to blue/nochange, otherwise select the cand flow red -> blue, 
            # (also challenge the blue) and recolor the other red to False
            cands = []
            for rid, [flow, color] in enumerate(col):
                if color == 'R' or color == 'B':
                    cands.append((rid, [flow, color]))
            if len(cands) == 0:
                # relax
                None
            elif len(cands) == 1:
                rid, [flow, color] = cands[0]
                if color == 'R':
                    adj_matrix[rid][colId] = [flow, 'B']
                    has_changes = True
            else:
                mrid, [mflow, _] = min(cands, key=lambda p: p[1][0])
                for rid, [flow, color] in cands:
                    adj_matrix[rid][colId] = [flow, 'G']
                adj_matrix[mrid][colId] = [mflow, 'B']
                has_changes = True
    
    for rowId, row in enumerate(adj_matrix):
        print([c for _, c in row])
        for colId, [abs_diff, color] in enumerate(row):
            if color != 'X':
                e = cliq_edge_dict[(index_to_cno[rowId], index_to_cno[colId])]
                if color != 'B':
                    cliq_graph.ep.color[e] = 'gray'
                else:
                    if abs_diff < threshold:
                        cliq_graph.ep.color[e] = 'black'
                    else:
                        cliq_graph.ep.color[e] = 'gray'
                    print("Decided Edge: {0} -> {1}, asb_diff: {2}".format(index_to_cno[rowId], index_to_cno[colId], abs_diff))
    
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    return cliq_graph, cliq_node_dict, cliq_edge_dict

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

            cno1m, cno2m = min(pairs, key=lambda p: 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[0]]] - cov, 2) + 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[1]]] - cov, 2))
            
            print("-->PAIR UP {0} - {1}, cov: {2}, diff: {3}".format(cno1m, cno2m, cov, delta))
            if (cno1m.split('->')[-1], cno2m.split('->')[0]) in sp_path_dict:
                cno1m_rel = cno1m.split('->')[-1]
                cno2m_rel = cno2m.split('->')[0]
                if sp_path_dict[(cno1m_rel, cno2m_rel)] != None:
                    cand_path, plen, pmark = sp_path_dict[(cno1m_rel, cno2m_rel)]
                else:
                    # cno1m -> end overlap with cno2m
                    intersect = set(contig_dict[cno1m][0]).intersection(set(contig_dict[cno2m][0]))
                    print("{0} -> {1} end to end intersects: {2}".format(cno1m, cno2m, intersect))
                    intermediate_nodes_index = [False for _ in contig_dict[cno1m][0]]
                    for i in [contig_dict[cno1m][0].index(e) for e in intersect]:
                        intermediate_nodes_index[i] = True
                    # in order
                    cand_path = [simp_node_dict[node_id] for i, node_id in enumerate(contig_dict[cno1m][0]) if intermediate_nodes_index[i]]
                    plen = path_len(graph, cand_path, overlap)
                    contig1m, clen1m, ccov1m = contig_dict[cno1m]
                    contig2m, clen2m, ccov2m = contig_dict[cno2m]
                    contig1m = contig1m[:-len(cand_path)]
                    clen1m = path_len(graph, [simp_node_dict[node_id] for node_id in contig1m], overlap)
                    contig_dict[cno1m] = [contig1m, clen1m, ccov1m]

                    contig2m = contig2m[len(cand_path):]
                    clen2m = path_len(graph, [simp_node_dict[node_id] for node_id in contig2m], overlap)
                    contig_dict[cno2m] = [contig2m, clen2m, ccov2m]
                    print("overlap contig concat: ", cno1m, contig1m, cno2m, contig2m, intersect)
            else:
                src = simp_node_dict[contig_dict[cno1m][0][-1]]
                tgt = simp_node_dict[contig_dict[cno2m][0][0]]
                cand_path, plen, pmark = dijkstra_sp(graph, src, tgt, cov, threshold, overlap)
            cand_len = get_concat_len(cno1m, contig_dict[cno1m][1], cno2m, contig_dict[cno2m][1], plen, overlap)

            contig_pair_reduction(cno1m, cno2m, cov, cand_path, cand_len, cno_mapping)

        return
    ###################################################################################################        

    # retrieve all the self cycle first
    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                # remove self cycle edge with self cycle + outer connection feature
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_graph.edge(contig_node, contig_node))
                print("should be removed already ... remove self edge+outer connection {0} -> {0}, legacy".format(cno))
            else:
                print("Circular PATH: ", cno)
                cov = cliq_graph.vp.ccov[contig_node]
                if (cno, cno) in sp_path_dict:
                    cand_path, plen, pmark = sp_path_dict[(cno, cno)]
                else:
                    src = simp_node_dict[contig_dict[cno][0][-1]]
                    tgt = simp_node_dict[contig_dict[cno][0][0]]
                    cand_path, plen, pmark = dijkstra_sp(graph, src, tgt, cov, threshold, overlap)
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


    for (u, v), e in cliq_edge_dict.items():
        print("potential uncaught error: EDGE: ", u, v)
    
    # simplify the graph
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    # draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), tempdir, "cliq_graphL3.png")

    return contig_dict

def strain_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
min_cov, max_len, threshold, overlap):
    """
    Extend the strain length
    1. map all the strain back into the graph
    """
    print("--------------Start strain extension------------------")
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap)

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
                sp, _, _ = dijkstra_sp(graph, global_src, contig_head, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [[graph.vp.id[n] for n in sp] + contig_dict[cno][0], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")

            contig_tail = simp_node_dict[contig[-1]]  
            if contig_tail.out_degree() != 0 and global_sink not in contig_tail.out_neighbors():
                print("---> contig tail: {0} ~~> {1}".format(graph.vp.id[contig_tail], graph.vp.id[global_sink]))
                sp, _, _ = dijkstra_sp(graph, contig_tail, global_sink, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [contig_dict[cno][0] + [graph.vp.id[n] for n in sp], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")
            
            # # re-assign the strain cov to min flow among the contig
            # assert len(contig_dict[cno][0]) >= 1
            # if len(contig_dict[cno][0]) == 1:
            #     redcov = graph.vp.dp[contig_dict[cno][0][0]]
            # else:
            #     redcov = numpy.min(contig_flow(graph, simp_edge_dict, contig_dict[cno][0]))
            contig_dict[cno][2] = ccov
    remove_global_source_sink(graph, global_src, global_sink)
    return contig_dict

def final_strain_extraction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, usage_dict: dict, threshold, overlap):
    """
    extract rest of partial used node path from global source to sink if any
    """
    print("--------------Start strain final extraction------------------")
    global_src = simp_node_dict['global_src']
    global_sink = simp_node_dict['global_sink']

    for no, [_, _, state] in usage_dict.items():
        if state == 'full' or state == 'over':
            node = simp_node_dict[no]
            graph.vp.color[node] = 'gray'
            for edge in node.all_edges():
                graph.ep.color[edge] = 'gray'
    
    path = []
    curr_id = 0
    while path != None:
        path, plen = dijkstra_sp_v2(graph, simp_node_dict, global_src, global_sink, overlap)
        if path != None:
            # FIXME
            cflows = contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in path[1:-1]])
            # # stat
            # print("mincov: ", numpy.min(cflows), 
            # "meancov: ", numpy.mean(cflows), 
            # "mediancov: ", numpy.median(cflows),
            # "maxcov: ", numpy.max(cflows))
            
            # d = pandas.Series({"coverage": cflows})
            # ax = seaborn.countplot(x="coverage", data=d)

            # ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
            # plt.title('Count plot of Path {0}'.format("st" + str(curr_id)))
            # plt.savefig("{0}countplot_{1}.png".format(TEMP_DIR, "st" + str(curr_id)))
            # # end stat
            redcov = numpy.min(cflows)
            for node in path[1:-1]:
                node_id = graph.vp.id[node]
                usage_dict[node_id][0] += redcov
                usage_dict[node_id][2] = node_status(usage_dict[node_id][0], 
                    usage_dict[node_id][1], threshold)
                if usage_dict[node_id][2] == 'full' or usage_dict[node_id][2] == 'over':
                    graph.vp.color[node] = 'gray'
                    for edge in node.all_edges():
                        graph.ep.color[edge] = 'gray'
            contig_dict["st" + str(curr_id)] = [[graph.vp.id[n] for n in path[1:-1]], plen, redcov]
            curr_id += 1
        else:
            print("Path not found")
    for no, [used, cap, state] in usage_dict.items():
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(used), round(cap), state))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    # # double check usage
    # post_usages = {}
    # post_free_nodes = []
    # post_partial_used_nodes = []
    # for no, node in simp_node_dict.items():
    #     print("----------------------------------------------------")
    #     if no == 'global_src' or no == 'global_sink':
    #         continue
    #     ratio = round(((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2)
    #     print("Node: {0}, full: {1}, left: {2}, usage: {3}, {4}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), ratio, graph.vp.color[node]))
    #     if ratio < 100 and graph.vp.udp[node] > threshold:
    #         post_partial_used_nodes.append(no)
    #     if ratio <= 0:
    #         post_free_nodes.append(no)
    #     post_usages[no] = ratio
    # overall_post_usage = numpy.mean(list(post_usages.values()))
    # print("Free nodes: ", list_to_string(post_free_nodes))
    # print("Partial used nodes: ", list_to_string(post_partial_used_nodes))
    # print("Overall Usage Post: {0}".format(overall_post_usage))
    # # do some plotting
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
    # plt.savefig("{0}barplot_usage_post.png".format(TEMP_DIR))
    # # post process
    # for node in graph.vertices():
    #     graph.vp.dp[node] = graph.vp.udp[node]
    # simp_node_dict.pop(graph.vp.id[global_src])
    # graph.vp.color[global_src] = 'gray'
    # simp_node_dict.pop(graph.vp.id[global_sink])
    # graph.vp.color[global_sink] = 'gray'
    remove_global_source_sink(graph, global_src, global_sink)
    return

def local_search_optimisation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
overlap, threshold):   
    """
    need to check whether current path flow below maximum flow 
    """
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, store_dict=True)
    
    usage_dict = {}
    # init status to all the node
    node_to_contig_dict, _ = contig_map_node(contig_dict)

    full_count = 0
    partial_count = 0
    over_count = 0
    for no, node in simp_node_dict.items():
        if no not in node_to_contig_dict:
            usage_dict[no] = [0, graph.vp.dp[node], 'partial']
            partial_count += 1
        else:
            sum_used = 0
            for cno in node_to_contig_dict[no]:
                sum_used += contig_dict[cno][2]
            if sum_used - graph.vp.dp[node] > threshold:
                # overused
                usage_dict[no] = [sum_used, graph.vp.dp[node], 'over']
                over_count += 1
            else:
                if abs(sum_used - graph.vp.dp[node]) < threshold:
                    usage_dict[no] = [sum_used, graph.vp.dp[node], 'full']
                    full_count += 1
                else:
                    usage_dict[no] = [sum_used, graph.vp.dp[node], 'partial']
                    partial_count += 1
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(usage_dict[no][0]), round(usage_dict[no][1]), usage_dict[no][2]))

    print("overflow nodes: ", over_count)
    print("full nodes: ", full_count)
    print("partial used nodes: ", partial_count)
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    swap = True
    itercount = 0
    while swap:
        swap = False
        itercount += 1
        print("iteration: ", itercount)
        for sno, snode in simp_node_dict.items():
            for eno, enode in simp_node_dict.items():
                # check if sno & eno adjacent
                if sno == eno or snode in enode.all_neighbors():
                    continue
                # not adjacent nodes, check any bounded over-path and partial-path exist
                over_paths, partial_paths = path_replacement_account(graph, simp_node_dict, simp_edge_dict, usage_dict, snode, enode)
                if len(over_paths) <= 0 or len(partial_paths) <= 0:
                    continue
                print("-----*curr s: {0} t: {1}".format(sno, eno))
                all_cnos = set()
                all_paths = {}
                prev_cno_mapping = {}
                for i, op in enumerate(over_paths):
                    min_bound = min([graph.vp.dp[n] for n in op])
                    print(path_to_id_string(graph, op, "over cap: {0}".format(min_bound)))

                    cno_set = None
                    for n in op:
                        if graph.vp.id[n] in node_to_contig_dict:
                            acc_set = set(node_to_contig_dict[graph.vp.id[n]])
                            cno_set = acc_set if cno_set == None else cno_set.intersection(acc_set)
                        else:
                            cno_set = set()
                            break
                    print("involved strain: ", [(cno, contig_dict[cno][2]) for cno in cno_set])
                    all_cnos = all_cnos.union(cno_set)
                    all_paths[('o', i)] = (op, min_bound)
                    for cno in cno_set:
                        prev_cno_mapping[cno] = ('o', i)

                for i, pp in enumerate(partial_paths):
                    min_bound = min([graph.vp.dp[n] for n in pp])
                    print(path_to_id_string(graph, pp, "partial cap: {0}".format(min_bound)))
                    cno_set = None
                    for n in pp:
                        if graph.vp.id[n] in node_to_contig_dict:
                            acc_set = set(node_to_contig_dict[graph.vp.id[n]])
                            cno_set = acc_set if cno_set == None else cno_set.intersection(acc_set)
                        else:
                            cno_set = set()
                            break
                    print("involved strain: ", [(cno, contig_dict[cno][2]) for cno in cno_set])
                    all_cnos = all_cnos.union(cno_set)
                    all_paths[('p', i)] = (pp, min_bound)
                    for cno in cno_set:
                        prev_cno_mapping[cno] = ('p', i)
                
                all_comb = list(product(list(all_paths.keys()), repeat=len(all_cnos)))
                all_cnos = list(all_cnos)

                prev_state_dict = {}
                for [path, _] in all_paths.values():
                    for n in path:
                        if graph.vp.id[n] not in prev_state_dict:
                            prev_state_dict[graph.vp.id[n]] = copy.deepcopy(usage_dict[graph.vp.id[n]])
                prev_overcount = len([k for k, [_, _, s] in prev_state_dict.items() if s == 'over'])
                print([(cno, contig_dict[cno][2]) for cno in all_cnos])
                optim_score = None
                optim_overcount = None
                optim_comb = []
                for i, comb in enumerate(all_comb):
                    # select the best combination that provides optimal capacity usage
                    comb_subscores = {}
                    for key, [_, cap] in all_paths.items():
                        comb_subscores[key] = [0, cap]
                    for j in range(len(all_cnos)):
                        flow = contig_dict[all_cnos[j]][2]
                        comb_subscores[comb[j]][0] += flow
                    #FIXME how to assign the score for one comb, include the check on over-node reduction
                    score_flow = pow(sum([flow/cap for [flow, cap] in comb_subscores.values()]) - len(comb_subscores), 2)
                    print("comb: ", [(path_to_id_string(graph, all_paths[key][0], str(all_paths[key][1])), "cno: " + all_cnos[i], contig_dict[all_cnos[i]][2]) for i, key in enumerate(comb)])
                    new_state_dict = test_comb(graph, comb, usage_dict, all_cnos, all_paths, prev_cno_mapping, contig_dict, threshold)
                    score_overcount = len([k for k, [_, _, s] in new_state_dict.items() if s == 'over'])
                    print("prev states: ", prev_state_dict)
                    print("new states: ", new_state_dict)
                    print("previous over node count: ", len([k for k, [_, _, s] in prev_state_dict.items() if s == 'over']))
                    print("new over node count: ", score_overcount)

                    if optim_score == None:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount
                    elif score_overcount < optim_overcount:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount 
                    elif score_overcount == optim_overcount and score_flow < optim_score:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount                     

                print("Optim comb: ", [(path_to_id_string(graph, all_paths[key][0], str(all_paths[key][1])), "cno: " + all_cnos[i], contig_dict[all_cnos[i]][2]) for i, key in enumerate(optim_comb)])
                print("Optim score: ", optim_score)
                print("Optim overnode: ", optim_overcount)
                if prev_overcount <= optim_overcount:
                    print("no more over usage node be reduced")
                    continue
                swap = True
                # replacement
                for i, key in enumerate(optim_comb):
                    # update contig dict
                    cno = all_cnos[i]
                    prev_p = [graph.vp.id[n] for n in all_paths[prev_cno_mapping[cno]][0]]
                    new_p = [graph.vp.id[n] for n in all_paths[key][0]]
                    print("prev: {0}, {1}, {2}".format(cno, list_to_string(contig_dict[cno][0]), list_to_string(prev_p)))
                    contig = contig_replacement_c(contig_dict[cno][0], prev_p, new_p)
                    print("after: {0}, {1}, {2}".format(cno, list_to_string(contig), list_to_string(new_p)))
                    clen = path_len(graph, contig, overlap)
                    contig_dict[cno] = [contig, clen, contig_dict[cno][2]]
                    # update usage dict
                    #TODO
                    for prev_id in prev_p:
                        # if prev_id in new_p:
                        #     continue
                        usage_dict[prev_id][0] -= contig_dict[cno][2]
                        prev_state = usage_dict[prev_id][2]
                        usage_dict[prev_id][2] = node_status(usage_dict[prev_id][0], usage_dict[prev_id][1], threshold)
                        print("prev-state changed: {0} {1}->{2}".format(prev_id, prev_state, usage_dict[prev_id][2]))
                    for new_id in new_p:
                        # if new_id in prev_p:
                        #     continue
                        usage_dict[new_id][0] += contig_dict[cno][2]
                        prev_state = usage_dict[new_id][2]
                        usage_dict[new_id][2] = node_status(usage_dict[new_id][0], usage_dict[new_id][1], threshold)
                        print("new-state changed: {0} {1}->{2}".format(new_id, prev_state, usage_dict[new_id][2]))
                # update the node to contig mapping
                node_to_contig_dict, _ = contig_map_node(contig_dict)

    over_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'over'])
    full_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'full'])
    partial_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'partial'])
    print("overflow nodes: ", over_count)
    print("full nodes: ", full_count)
    print("partial used nodes: ", partial_count)
    for no, [used, cap, state] in usage_dict.items():
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(used), round(cap), state))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    return usage_dict

def split_contig(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold):
    """
    split contig out
    """
    dup_record = {}
    reduced_edges = []
    reduced_vertices = []
    for no in simp_node_dict.keys():
        dup_record[no] = 0
    for cno in contig_dict.keys():
        contig, _, ccov = contig_dict[cno]
        contig = list(contig)
        print("--->split contig: ", cno)
        for i in range(1,len(contig) - 1):
            no = contig_dict[cno][0][i]
            prev_no = contig_dict[cno][0][i-1]
            prev_old_no = contig[i-1]
            ie = simp_edge_dict[(prev_old_no, no)]
            graph.ep.flow[ie] -= ccov
            reduced_edges.append(ie)           
            if graph.ep.flow[ie] <= threshold:
                graph_remove_edge(graph, simp_edge_dict, graph.vp.id[ie.source()], graph.vp.id[ie.target()])

            old_vertex = simp_node_dict[no]

            graph.vp.dp[old_vertex] -= ccov
            new_vertex = graph_add_vertex(graph, simp_node_dict, str(no) + "X" + str(dup_record[no]), 
                ccov, graph.vp.seq[old_vertex], graph.vp.kc[old_vertex])
            graph_add_edge(graph, simp_edge_dict, simp_node_dict[prev_no], 
                prev_no, new_vertex, graph.vp.id[new_vertex], graph.ep.overlap[ie], ccov)
            if i == len(contig) - 2:
                next_no = contig_dict[cno][0][i+1]
                oe = simp_edge_dict[(no, next_no)]
                graph.ep.flow[oe] -= ccov
                reduced_edges.append(oe)
                if graph.ep.flow[oe] <= threshold:
                    graph_remove_edge(graph, simp_edge_dict, graph.vp.id[oe.source()], graph.vp.id[oe.target()])
                graph_add_edge(graph, simp_edge_dict, new_vertex, graph.vp.id[new_vertex], 
                    simp_node_dict[next_no], next_no, graph.ep.overlap[ie], ccov) 
            
            reduced_vertices.append(old_vertex)
            if graph.vp.dp[old_vertex] <= threshold:
                graph_remove_vertex(graph, simp_node_dict, no)
            contig_dict[cno][0][i] = graph.vp.id[new_vertex]
            dup_record[no] += 1
    # print("DELTA: ", threshold)
    # for e in sorted(graph.edges(), key=lambda ee: graph.ep.flow[ee]):
    #     if e in reduced_edges:
    #         print("contig split")
    #     else:
    #         print("original")
    #     print_edge(graph, e, "lower than delta" if graph.ep.flow[e] < threshold 
    #         else "greater than delta")
    
    # for v in sorted(graph.vertices(), key=lambda vv: graph.vp.dp[vv]):
    #     if v in reduced_vertices:
    #         print("contig split")
    #     else:
    #         print("original")
    #     print_vertex(graph, v, "lower than delta" if graph.vp.dp[v] < threshold 
    #         else "greater than delta")
    return


def extract_cand_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, strain_dict: dict, overlap, threshold, itercount='A'):
    global TEMP_DIR
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, True)
    for e in sorted(graph.edges(), key=lambda a: graph.ep.flow[a]):
        print_edge(graph, e)
    lb = min([graph.ep.flow[e] for e in graph.edges()])
    ub = lb*2
    print("[{0},{1})".format(lb, ub))
    threshold = ub/20
    print("local threshold: ", threshold)

    min_eset = {}
    for e in graph.edges():
        if graph.ep.flow[e] > lb and graph.ep.flow[e] <= ub:
            min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e

    x = [graph.ep.flow[e] for e in min_eset.values()]
    regions, bins = numpy.histogram(x)
    # find peak region
    peak_region_index = list(regions).index(max(regions))
    peak_lb = bins[peak_region_index+1] - threshold
    # peak_lb = bins[peak_region_index]
    # peak_ub = bins[peak_region_index+1]
    peak_ub = bins[peak_region_index] + threshold
    #################################STAT###################################
    print(regions)
    print(bins)
    print(peak_lb, bins[peak_region_index], bins[peak_region_index+1], peak_ub)
    plt.figure(figsize=(128,64))
    plt.hist(x, bins)
    plt.axvline(x=peak_lb, color='r')
    plt.axvline(x=bins[peak_region_index], color='b')
    plt.axvline(x=bins[peak_region_index+1], color='b')
    plt.axvline(x=peak_ub, color='r')
    plt.title("min_eset")
    plt.savefig("{0}{1}".format(TEMP_DIR, "min_eset.png"))
    ########################################################################
    peak_eset = dict()
    peak_vset = set()
    connect_set = dict()
    for k, e in min_eset.items():
        if graph.ep.flow[e] >= peak_lb and graph.ep.flow[e] < peak_ub:
            peak_eset[k] = e
    for (u,v), e in peak_eset.items():
        connect_set[e] = []
        peak_vset.add(u)
        peak_vset.add(v)    
    print(list_to_string(list(peak_vset), "Peak related nodes: "))

    # connectivity
    for (au, av), e_tail in peak_eset.items():
        for (bu, bv), e_head in peak_eset.items():
            if au == bu and av == bv:
                continue
            if reachable(graph, simp_node_dict, simp_node_dict[av], simp_node_dict[bu]):
                connect_set[e_tail].append(e_head)
    
    # transitive reduction
    for k in connect_set:
        for i in connect_set:
            for j in connect_set:
                if k != i and i != j and k != j:
                    if k in connect_set[i] and j in connect_set[k] and j in connect_set[i]:
                        connect_set[i].remove(j)
    
    graphm = Graph(directed=True)
    graphm.vp.id = graphm.new_vertex_property("string")
    graphm.vp.flow = graphm.new_vertex_property("int32_t")
    graphm.vp.text = graphm.new_vertex_property("string")
    graphm.ep.color = graphm.new_edge_property("string")
    graphm_vdict = dict()
    graphm_edict = dict()
    for u in connect_set.keys():
        vertex = graphm.add_vertex()
        graphm.vp.id[vertex] = str(graph.vp.id[u.source()]) + "->" + str(graph.vp.id[u.target()])
        graphm.vp.flow[vertex] = graph.ep.flow[u]
        graphm.vp.text[vertex] = graphm.vp.id[vertex] + ":" + str(graphm.vp.flow[vertex])
        graphm_vdict[u] = vertex
    for u, vs in connect_set.items():
        for v in vs:
            e = graphm.add_edge(graphm_vdict[u], graphm_vdict[v])
            graphm_edict[(u,v)] = e
    
    #################################STAT###################################
    output_size = 120 * (len(list(graphm.vertices()) )+ len(list(graphm.edges())))
    vsize= 30
    graph_draw(g=graphm, output="{0}{1}".format(TEMP_DIR, "graphm.png"), bg_color="white", 
    vertex_text=graphm.vp.text, vertex_size=vsize, vertex_font_size=int(vsize * 0.8), 
    output_size=(output_size, output_size))
    print("min peak graph has been stored in: {0}{1}".format(TEMP_DIR, "graphm.png"))
    ########################################################################
    adjMtx = [[] for _ in connect_set.keys()]
    e2i = {}
    i2e = {}
    for i, e in enumerate(connect_set.keys()):
        e2i[e] = i
        i2e[i] = e
    for u, vs in connect_set.items():
        row = [[sys.maxsize,'X'] for _ in connect_set.keys()]
        for v in vs:
            abs_dif = abs(graph.ep.flow[u] - graph.ep.flow[v])
            row[e2i[v]] = [abs_dif, 'W']
        adjMtx[e2i[u]] = row
    
    has_changes = True
    dim = len(adjMtx)
    # iterate the adjacent matrix in rowwise and columnwise
    while has_changes:
        print("looping")
        has_changes = False
        # row wise
        for rowId in range(dim):
            row = get_row(adjMtx, rowId)
            hasFixed = False
            colId = None
            minFlow = sys.maxsize
            for cid, [flow, color] in enumerate(row):
                if color == 'B':
                    # current row fixed
                    hasFixed = True
                    break
                elif color == 'W':
                    if flow < minFlow:
                        colId = cid
                        minFlow = flow
                elif color == 'X' or color == 'G':
                    # invalid, skip
                    None
                else:
                    print("Error: ", rowId, i2e[rowId], cid, i2e[cid])
                    assert color != 'R'
            if not hasFixed:
                if colId != None:
                    # found a minimum block to assign
                    adjMtx[rowId][colId][1] = 'R'
                    has_changes = True
                else:
                    # no more spot for the row
                    None
        for colId in range(dim):
            col = get_col(adjMtx, colId)
            # hasFixed
            rowId = None
            minFlow = sys.maxsize
            # if only one red+blue among the column, then assign it to blue/nochange, otherwise select the cand flow red -> blue, 
            # (also challenge the blue) and recolor the other red to False
            cands = []
            for rid, [flow, color] in enumerate(col):
                if color == 'R' or color == 'B':
                    cands.append((rid, [flow, color]))
            if len(cands) == 0:
                # relax
                None
            elif len(cands) == 1:
                rid, [flow, color] = cands[0]
                if color == 'R':
                    adjMtx[rid][colId] = [flow, 'B']
                    has_changes = True
            else:
                mrid, [mflow, _] = min(cands, key=lambda p: p[1][0])
                for rid, [flow, color] in cands:
                    adjMtx[rid][colId] = [flow, 'G']
                adjMtx[mrid][colId] = [mflow, 'B']
                has_changes = True
    adj_list = {}
    for e in peak_eset.values():
        adj_list[e] = None
    for rowId, row in enumerate(adjMtx):
        for colId, [_, color] in enumerate(row):
            if color != 'X':
                e = graphm_edict[(i2e[rowId], i2e[colId])]
                if color != 'B':
                    graphm.ep.color[e] = 'gray'
                else:
                    graphm.ep.color[e] = 'black'
                    if adj_list[i2e[rowId]] != None:
                        print("Error, branch detected")
                    else:
                        adj_list[i2e[rowId]] = i2e[colId]
    
    for e in sorted(graphm.edges()):
        if graphm.ep.color[e] != 'black':
            graphm.remove_edge(e)
    #################################STAT###################################
    output_size = 120 * (len(list(graphm.vertices()) )+ len(list(graphm.edges())))
    vsize = 30
    graph_draw(g=graphm, output="{0}{1}".format(TEMP_DIR, "graphm_v2.png"), bg_color="white", 
    vertex_text=graphm.vp.text, vertex_size=vsize, vertex_font_size=int(vsize * 0.8), 
    output_size=(output_size, output_size))
    print("min peak graph has been stored in: {0}{1}".format(TEMP_DIR, "graphm_v2.png"))
    ########################################################################
    epaths = []
    for u in adj_list.keys():
        if graphm_vdict[u].in_degree() != 0:
            # intermediate node
            continue
        epath = []
        acc = u
        while acc != None:
            epath.append(acc)
            acc = adj_list[acc]
        epaths.append(epath)
    for i, epath in enumerate(sorted(epaths, key=len, reverse=True)):
        strain = []
        print("--->{0}{1} Strain finding.. length: {2}".format(itercount, i, len(epath)))
        #FIXME
        ccov = max([graph.ep.flow[e] for e in epath])
        if epath[0].source() != global_src:
            # generate path from global src to epath 1st node
            s = global_src
            t = epath[0].source()
            sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
            strain.extend(sp)
            strain.append(t)
        for i in range(len(epath) - 1):
            s = epath[i].target()
            t = epath[i+1].source()
            if s == t:
                strain.append(s)
            else:
                # find intermediate path between s_e to t_e
                sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
                strain.append(s)
                strain.extend(sp)
                strain.append(t)
        if epath[-1].target() != global_sink:
            # generate path from epath last node to global sink
            s = epath[-1].target()
            t = global_sink
            sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
            strain.append(s)
            strain.extend(sp)
        print("Strain - ccov: {0}, {1}".format(ccov, path_to_id_string(graph, strain)))
        # strain_dict[itercount + str(i)] = [[graph.vp.id[n] for n in strain], path_len(graph, strain, overlap), ccov]
        for node in strain:
            eval_score(graph.vp.id[node], graph.vp.dp[node] - ccov, threshold, lb, True)
        # graph_reduction_c(graph, strain, ccov)
        for e in epath:
            strain_dict[graph.vp.id[e.source()]] = [[graph.vp.id[e.source()]], path_len(graph, [e.source()], overlap), graph.vp.dp[e.source()]]
            strain_dict[graph.vp.id[e.target()]] = [[graph.vp.id[e.target()]], path_len(graph, [e.target()], overlap), graph.vp.dp[e.target()]]
        break


    return

def extract_cand_path2(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, threshold, itercount='A'):
    def dict_to_hist(graph: Graph, contig_dict: dict):
        contig_hist = {}
        for (cno, [contig, clen, ccov]) in contig_dict.items():
            lb = ccov - threshold
            ub = ccov + threshold
            min_eset = {}
            min_vset = set()
            for e in graph.edges():
                if graph.ep.flow[e] > lb and graph.ep.flow[e] <= ub:
                    # if (reachable(graph, simp_node_dict, e.target(), simp_node_dict[contig[0]]) or
                    #     reachable(graph, simp_node_dict, simp_node_dict[contig[-1]], e.source())):
                        min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
                        min_vset.add(e.source())
                        min_vset.add(e.target())
            x = [graph.ep.flow[e] for e in min_eset.values()]
            regions, bins = numpy.histogram(x)
            contig_hist[cno] = [clen, ccov, regions, bins]
        return contig_hist
    strain_dict = {}
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, True)
    contig_hist = {}
    sorted_contigdict = sorted(contig_dict.items(), key=lambda v: v[1][2])
    x = [cov for [_, _, cov] in contig_dict.values()]
    regions, bins = numpy.histogram(x)
    print(regions)
    print(bins)
    plt.figure(figsize=(64,32))
    plt.hist(x)

    plt.title("contig_cov")
    plt.savefig("{0}{1}".format(TEMP_DIR, "contig_cov.png"))
    for (cno, [contig, clen, ccov]) in sorted_contigdict:
        print("-----------------------------------*****")
        print_contig(cno, clen, ccov, contig)
        lb = ccov - threshold
        ub = ccov + threshold
        min_eset = {}
        min_vset = set()
        for e in graph.edges():
            if graph.ep.flow[e] > lb and graph.ep.flow[e] <= ub:
                # if (reachable(graph, simp_node_dict, e.target(), simp_node_dict[contig[0]]) or
                #     reachable(graph, simp_node_dict, simp_node_dict[contig[-1]], e.source())):
                    min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
                    min_vset.add(e.source())
                    min_vset.add(e.target())
        x = [graph.ep.flow[e] for e in min_eset.values()]
        regions, bins = numpy.histogram(x)
        print(regions)
        print(bins)
        print("most related vertices: ", path_to_id_string(graph, list(min_vset)))

        contig_hist[cno] = [clen, ccov, regions, bins]

    while len(contig_dict.keys()) > 0:
        print("----------------------------------------------------------------------------")
        cno = max(contig_dict.keys(), key=lambda k: sum(contig_hist[k][2]))
        # cno = max(contig_dict.keys(), key=lambda k: path_len(graph, [simp_node_dict[id] for id in contig_dict[k][0]], overlap))
        # cno = min(contig_dict.keys(), key=lambda k: contig_dict[k][2])
        print("REGION: ", contig_hist[cno][2])
        print("BIN: ", contig_hist[cno][3])
        contig, clen, ccov = contig_dict.pop(cno)
        lb = ccov - threshold
        ub = ccov + threshold
        print(cno, "->bound: ", lb, ccov, ub, "***" ,list_to_string(contig))

        if ccov < threshold:
            print("current contig {0} is used previously {1}".format(ccov, threshold))
            continue

        sphead, _, _ = dijkstra_sp_v3(graph, global_src, simp_node_dict[contig[0]], ccov, threshold, overlap, lb)
        sptail, _, _ = dijkstra_sp_v3(graph, simp_node_dict[contig[-1]], global_sink, ccov, threshold, overlap, lb)
        strain = []
        if sphead != None:
            strain.extend(sphead)
        strain.extend([simp_node_dict[n] for n in contig])
        if sptail != None:
            strain.extend(sptail)

        score = []
        for flow in contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in strain]):
            s = eval_score(flow, ccov, threshold, lb)
            score.append(s)
        ccov = path_cov(graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in strain])
        print(path_to_id_string(graph, strain, "strain, len{0}ccov{1}".format(len(strain), ccov)))
        print("related edges score: ", score)

        print("cand strain found")
        graph_reduction_c(graph, strain, ccov)
        strain_dict[itercount + cno] = [[graph.vp.id[n] for n in strain], path_len(graph, strain, overlap), ccov]
        
        contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict)
        contig_hist = dict_to_hist(graph, contig_dict)

    strains = sorted(strain_dict.items(), key=lambda k: k[1][2], reverse=True)
    strain_dict = {}
    for [cno, [contig, clen, ccov]] in strains:
        if ccov >= threshold:
            strain_dict[cno] = [contig, clen, ccov]
    return strain_dict

if __name__ == "__main__":
    main()