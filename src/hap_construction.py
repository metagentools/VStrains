#!/usr/bin/env python3

import sys, os
import subprocess
import argparse

from graph_tool.topology import all_circuits
from graph_tool.draw import graph_draw
from graph_tool.all import Graph

import numpy
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
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, required=True, help='contig file from SPAdes, paths format')
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
    contig_dict = get_contig(args.contig_file, simp_node_dict, simp_edge_dict, args.min_len)
    
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}graph_L0.gfa".format(TEMP_DIR))
    print("----------------------------------TIP REMOVAL---------------------------------------")
    # apply tip removal for cyclic graph.
    if not graph_is_DAG(graph0, simp_node_dict0):
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(graph0, simp_node_dict0, args.overlap)

        contig_dict_fix(graph0, simp_node_dict0, contig_dict, args.overlap)
    else:
        print("Graph is DAG, tip removal skipped.")
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1.gfa".format(TEMP_DIR))

    print("-----------------------------DELTA ESTIMATION-----------------------------")
    delta_estimation(graph1, simp_node_dict1, simp_edge_dict1)
    print("-------------------------------GRAPH SIMPLIFICATION & REBALANCE-----------------------------------")
    #FIXME
    mediandp = numpy.median([graph1.vp.dp[node] for node in simp_node_dict1.values()])
    THRESHOLD = mediandp/20
    # 20000 * 0.001 
    # numpy.quantile([graph1.vp.dp[node] for node in graph1.vertices()], 0.05)
    # mediandp/20
    print("MEDIAN NODE DEPTH: ", mediandp, "threshold: ", THRESHOLD)
    graph_simplification(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, THRESHOLD)

    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}st_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}st_graph_L2.gfa".format(TEMP_DIR))

    coverage_rebalance_s(graph2, simp_node_dict2, simp_edge_dict2, True)

    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)
    
    print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict)
    draw_edgeflow(graph3, simp_edge_dict3, TEMP_DIR, 'Bar plot of Edge flow', 'barplot_edge_flow.png')
    
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict3, "{0}graph_L0.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), 
            "{0}temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph3, simp_node_dict3, contig_dict, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    # end stat

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

        # fix trivial splitted contig
        contig_dict_remapping(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, id_mapping, prev_ids, args.overlap)
        # non-trivial branch split
        num_split = graph_splitting(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, THRESHOLD, args.overlap)
    
        graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))

        simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, args.overlap)

        graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        grapha, simp_node_dicta, simp_edge_dicta = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(grapha, simp_node_dicta, simp_edge_dicta)

        contig_dup_removed(grapha, simp_edge_dicta, contig_dict)
        trim_contig_dict(grapha, simp_node_dicta, contig_dict, args.overlap)
        contig_cov_fix(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            iterCount = chr(ord(iterCount) + 1)
        else:
            coverage_rebalance_s(grapha, simp_node_dicta, simp_edge_dicta, True)
            graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
            break
    print("Total non-trivial branches removed: ", total_removed_branch_nt, " total trivial branches removed: ", total_removed_branch_t)

    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}post_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph5, simp_node_dict5, contig_dict, args.overlap, "{0}post_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}post_contigs.fasta".format(TEMP_DIR), "{0}post_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}post_contigs_to_strain.paf".format(TEMP_DIR))
    # end stat

    print("-----------------------CONTIG PATH EXTENSION-------------------------------")
    strain_dict = extract_cand_path2(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.overlap, THRESHOLD)
    
    print("-----------------------FINAL CLEAN UP-------------------------------")
    contig_dup_removed(graph5, simp_edge_dict5, strain_dict)  
    trim_contig_dict(graph5, simp_node_dict5, strain_dict, args.overlap)  
    contig_dict_to_fasta(graph5, simp_node_dict5, strain_dict, args.overlap, "{0}final_contigs.fasta".format(TEMP_DIR))
    contig_dict_to_path(strain_dict, "{0}final_contigs.paths".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}final_contigs.fasta".format(TEMP_DIR), "{0}final_contigs_to_strain.paf".format(TEMP_DIR))
    return 0

def delta_estimation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    global TEMP_DIR
    print("test delta")
    xs = []
    ys = []
    cxs = []
    cys = []
    for node in graph.vertices():
        if sum([x.out_degree() for x in node.in_neighbors()]) == node.in_degree() and sum([y.in_degree() for y in node.out_neighbors()]) == node.out_degree():
            lv = sum([graph.vp.dp[n] for n in node.in_neighbors()])
            rv = sum([graph.vp.dp[n] for n in node.out_neighbors()])
            m = graph.vp.dp[node]
            xs.append(lv)
            ys.append(m)
            xs.append(rv)
            ys.append(m)
            cxs.append(m)
            cys.append(m)
    
    fig = plt.figure(figsize=(48,36))
    plt.scatter(xs, ys, s=40)
    plt.plot(cxs, cys, linewidth=3.0)
    plt.title("delta_scatter_plot")
    plt.savefig("{0}{1}".format(TEMP_DIR, "delta_scatter_plot.png"))
    
    fig2 = plt.figure(figsize=(48,36))
    plt.hist([b-a for b, a in zip(ys,xs)], bins=len(ys))
    plt.title("delta_hist_plot")
    plt.savefig("{0}{1}".format(TEMP_DIR, "delta_hist_plot.png"))
    return

curr_path = []
def node_partition(graph: Graph, simp_node_dict: dict):
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

    # init path acc variable
    global curr_path, TEMP_DIR
    curr_path = []
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

    # partitate into linear graph
    noncyc_graph = graph.copy()
    simp_node_dict_noncyc, simp_edge_dict_noncyc = graph_to_dict(noncyc_graph)
    graph_color_other_to_gray(noncyc_graph, simp_node_dict_noncyc, noncyc_nodes)
    graph_to_gfa(noncyc_graph, simp_node_dict_noncyc, simp_edge_dict_noncyc, "{0}nc_graph_L2p.gfa".format(TEMP_DIR))

    return noncyc_nodes, simple_paths

def tip_removal(graph: Graph, simp_node_dict: dict, overlap):
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

    def cand_collapse_path(graph: Graph, from_node, to_paths, temp_dir, overlap, accept_rate = 0.98):
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
    global TEMP_DIR
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
        cand_path = cand_collapse_path(graph, src, potential_paths, TEMP_DIR, overlap)
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
        cand_path = cand_collapse_path(graph, tgt, potential_paths, TEMP_DIR, overlap)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, tgt, cand_path)
            is_removed = False
        else:
            print_vertex(graph, tgt, "Tip cannot be removed, no matching path")

    return is_removed

def coverage_rebalance_s(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strict=False):
    global TEMP_DIR
    print("------------------------------NODE PARTITION-----------------------------------")
    # store the no-cycle nodes in nc_graph_L3p.gfa
    noncyc_nodes = None
    simple_paths = None
    # graphtool is_DAG() may not work if the graph is not connected as several parts
    if not graph_is_DAG(graph, simp_node_dict):
        noncyc_nodes, simple_paths = node_partition(graph, simp_node_dict, TEMP_DIR)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth has been stored in the dict
    # ratio: normalised balanced node depth / previous node depth
    _, _, _, ratio = coverage_rebalance(graph, simp_node_dict, simp_edge_dict, strict)
    
    if noncyc_nodes != None and simple_paths != None:
        print("rebalance linear subgraph now..")
        graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}nc_graph_L2p.gfa".format(TEMP_DIR))
        coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc, strict)
        print("Done, start coverage merge")

        for no, node in simp_node_dictnc.items():
            cnode = simp_node_dict[no]
            merge_dp = graphnc.vp.dp[node] + graph.vp.dp[cnode]
            graph.vp.dp[cnode] = merge_dp   
    else:
        print("no linear subgraph available..")

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strict=False):
    def expectation_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
        for (u,v),e in simp_edge_dict.items():

            if graph.ep.color[e] != 'black':
                #forbidden edge
                continue

            u_node = simp_node_dict[u]
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])

            v_node = simp_node_dict[v]
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])

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
    cutoff = 0
    if graph_is_DAG(graph, simp_node_dict):
        cutoff = 0.00001 * len(simp_node_dict) if strict else 0.0001 * len(simp_node_dict)
    else:
        cutoff = 0.0001 * len(simp_node_dict) if strict else 0.01 * len(simp_node_dict)
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

def extract_cand_path2(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, threshold, itercount='A'):
    def dict_to_hist(graph: Graph, contig_dict: dict):
        contig_hist = {}
        for (cno, [_, clen, ccov]) in contig_dict.items():
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
        # check if contig can be concated self-t to self-s, if so then cycle
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
            if s == 'P4':
                print("P4: ", flow)
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