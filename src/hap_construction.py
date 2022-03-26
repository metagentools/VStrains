#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
# from graph_tool import GraphView, _in_degree
# import graph_tool
# from graph_tool.search import dfs_iterator
from graph_tool.topology import all_circuits, all_shortest_paths
from graph_tool.draw import graph_draw

import subprocess
from graph_tool.all import Graph
# from graph_tool.topology import is_DAG
import argparse
import numpy
import heapq

from graph_converter import *

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
    -> tip removal based on minimap2 [DONE] 
    Output --> t_graph_L1.gfa
    -> split graph branch based on contig [DONE]
    Output --> st_graph_L2.gfa
    * bubble removal based on minimap2 [FORBIDDEN], since bubble removal may not suitable for metagenomics
    -> cycle detection and node partition
    Output --> c_graph_L3p.gfa, nc_graph_L3p.gfa
    -> node depth rebalance + assign edge flow [DONE] 
    Output --> dst_graph_L3.gfa
    -> compact the reassigned/splitted contig into a single node. Contig coverage be re-assigned to max node dp [TODO] 
    Output --> cdst_graph_L4.gfa
    -> compact all the simple path if exist. [DONE] 
    Output --> scdst_graph_L5.gfa
    -> construct contig clique graph
    -> remove all the cycle as the cand strains from the graph, also reduce 
    the stored dp(A) by amount of cycle strain coverage if A is related 
    to the cycle strain. [TODO]
    Output --> cyc_strain.paths, cyc_strain.fasta, cR_scdbt_graph_L6.gfa
    -> reassgin A's depth to dp(A), and perform node depth rebalance [TODO]
    -> s-t path to enumerate all the linear strain [TODO]
    Output --> TODO
    -> find the best s-t path option that used up all the node coverage. [TODO]
    Output --> TODO
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
    contig_dict, node_to_contig_dict, edge_to_contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_len)
    
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}graph_L0.gfa".format(TEMP_DIR))
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict, "{0}graph_L0.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
        contig_dict_to_fasta(graph, contig_dict, simp_node_dict, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
        minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    
    print("----------------------------------TIP REMOVAL---------------------------------------")
    tip_removed = False
    while not tip_removed:
        tip_removed = tip_removal(graph0, simp_node_dict0, simp_edge_dict0, TEMP_DIR, args.overlap)
    
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1.gfa".format(TEMP_DIR))

    print("--------------------------CONTIG RELATED NODE SPLIT-----------------------------")
    contig_dict_fix(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, args.overlap)
    # new_contig_dict = graph_splitting(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, args.overlap, TEMP_DIR, args.min_cov)
    
    # del contig_dict
    for cno, [contig, clen, ccov] in contig_dict.items():
        print("---------------------------------------------------------------")
        print_contig(cno, clen, ccov, contig)

    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}st_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}st_graph_L2.gfa".format(TEMP_DIR))
    print("------------------------------NODE PARTITION-----------------------------------")
    # store the no-cycle nodes in nc_graph_L3p.gfa
    noncyc_nodes = None
    simple_paths = None
    # is dag may not work if the graph is not connected as several parts
    if not graph_is_DAG(graph2, simp_node_dict2):
        noncyc_nodes, simple_paths = node_partition(graph2, simp_node_dict2, simp_edge_dict2, TEMP_DIR)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth use been store in the dict
    # ratio: normalised balanced node depth / previous node depth
    _, _, _, ratio = coverage_rebalance(graph2, simp_node_dict2, simp_edge_dict2, contig_dict)
    
    if noncyc_nodes != None and simple_paths != None:
        print("rebalance linear subgraph now..")
        graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}nc_graph_L3p.gfa".format(TEMP_DIR))
        coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc, contig_dict)
        print("Done, start coverage merge")

        for no, node in simp_node_dictnc.items():
            cnode = simp_node_dict2[no]
            merge_dp = graphnc.vp.dp[node] + graph2.vp.dp[cnode]
            graph2.vp.dp[cnode] = merge_dp   
    else:
        print("no linear subgraph available..")


    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}dst_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}dst_graph_L3.gfa".format(TEMP_DIR))
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)

    # print("-------------------------------GRAPH SIMPLIFICATION AND REBALANCE-----------------------------------")
    # print("MAX NODE DEPTH: ", numpy.max([graph3.vp.dp[node] for node in simp_node_dict3.values()]))
    # ids = []
    # for no, node in simp_node_dict3.items():
    #     if graph3.vp.dp[node] < 50:
    #         ids.append(no)
    # print("nodes that less than 50 cov be removed: ", list_to_string(ids))
    # graph_simplification(graph3, simp_node_dict3, simp_edge_dict3, node_to_contig_dict, edge_to_contig_dict, 50)

    # graph_to_gfa(graph3, simp_node_dict3, simp_edge_dict3, "{0}cdst_graph_L4.gfa".format(TEMP_DIR))
    # graph4, simp_node_dict4, simp_edge_dict4 = flipped_gfa_to_graph("{0}cdst_graph_L4.gfa".format(TEMP_DIR))

    # coverage_rebalance(graph4, simp_node_dict4, simp_edge_dict4, contig_dict)
    # assign_edge_flow(graph4, simp_node_dict4, simp_edge_dict4)

    # print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    # # re-evaluate the contig coverage
    # for cno, [contig, clen, ccov] in contig_dict.items():
    #     print("---------------------------------------------------------------")
    #     if len(contig) > 1:
    #         contig_dict[cno][2] = numpy.min(contig_flow(graph3, simp_edge_dict3, contig))
    #     else:
    #         contig_dict[cno][2] = graph3.vp.dp[simp_node_dict3[contig[0]]]
    #     print_contig(cno, clen, contig_dict[cno][2], contig)

    # print("-------------------------------GRAPH COMPACTIFICATION-----------------------------------")
    # # Compact all the contig into single node and saved to contig_node_dict, Store once, then compact
    # # the rest of the simple paths.

    # contig_nodes = []
    # for [contig, _, _] in contig_dict.values():
    #     contig_nodes.extend(contig)

    # simp_path_dict = simple_paths_to_dict(graph4, simp_node_dict4, simp_edge_dict4, contig_nodes, args.overlap)
    # simp_path_compactification(graph4, simp_node_dict4, simp_edge_dict4, simp_path_dict, args.overlap)

    # graph_to_gfa(graph4, simp_node_dict4, simp_edge_dict4, "{0}scdst_graph_L5.gfa".format(TEMP_DIR))
    # graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}scdst_graph_L5.gfa".format(TEMP_DIR))
    # assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)

    # print("-------------------------------GRAPH BRANCH SPLIT------------------------------------------")

    # graph_splitting(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.overlap, TEMP_DIR, args.min_cov)
    
    # print("-------------------------------CONTIG CLIQUE GRAPH BUILD-----------------------------------")
    # if args.ref_file:
    #     map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}scdst_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    
    # cliq_graph, cliq_node_dict, cliq_edge_dict, rand_path_dict = contig_clique_graph_build(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.min_cov, args.min_len, args.max_len, args.overlap, TEMP_DIR)
    # print("-------------------------------CONTIG PAIRWISE CONCATENATION-----------------------------------")
    # contig_pairwise_concatenation(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.min_cov, args.min_len, args.max_len, args.overlap)

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
    graph_to_gfa(noncyc_graph, simp_node_dict_noncyc, simp_edge_dict_noncyc, "{0}nc_graph_L3p.gfa".format(temp_dir))

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

def bubble_removal():
    return

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

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, strict=False):
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

    def contig_branch_split(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold):
        """
        for any N-N branch with k contig be involved, split into (N-k)-(N-k) subranches, 1-1, 1-1, ..., 1-1 simple path if coverage
        is matched.
        """
        return
    
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
        print(sum_delta)
        if sum_delta < cutoff:
            break
        # M Step
        is_update = maximization_node_depth(graph, simp_node_dict)

    # final evaluation
    ratios = [(graph.vp.dp[u] / prev_dp_dict[no]) for no, u in simp_node_dict.items()]
    sum_depth_before = numpy.sum([dp for dp in prev_dp_dict.values()])
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Sum Ratio: ", sum_ratio, "Ave Ratio: ", numpy.mean(ratios), "Max Ratio: ", numpy.max(ratios), "Min Ratio: ", numpy.min(ratios), "Delta: ", sum_delta)

    for node in simp_node_dict.values():
        graph.vp.dp[node] = graph.vp.dp[node] / sum_ratio

    node_ratio_dict = {}
    for no in prev_dp_dict.keys():
        if prev_dp_dict[no] != 0:
            node_ratio_dict[no] = (graph.vp.dp[simp_node_dict[no]] / prev_dp_dict[no]) if prev_dp_dict[no] != 0 else 0
    curr_dp_dict = {}
    for no, v in simp_node_dict.items():
        curr_dp_dict[no] = graph.vp.dp[v]

    return prev_dp_dict, curr_dp_dict, node_ratio_dict, sum_ratio

def allowed_concat_init(graph: Graph, contig_dict: dict, simp_node_dict: dict, max_len, overlap):
    
    self_concat_off = graph_is_DAG(graph, simp_node_dict)
    impossible_concat_dict = {}
    graph.vp.prev = graph.new_vertex_property("string", val="")

    rand_path_dict = {}

    for no in contig_dict.keys():
        impossible_concat_dict[no] = set()
    for tail_cno, [tail_contig, tail_clen, _] in contig_dict.items():
        for head_cno, [head_contig, head_clen, _] in contig_dict.items():
            print("tail cno: ", tail_cno, " vs ", head_cno)

            tail_node = simp_node_dict[contig_dict[tail_cno][0][-1]]
            head_node = simp_node_dict[contig_dict[head_cno][0][0]]

            if tail_node.out_degree() == 0 or head_node.in_degree() == 0:
                impossible_concat_dict[tail_cno].add(head_cno)

            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len or (set(tail_contig)).intersection(set(head_contig)):
                    impossible_concat_dict[tail_cno].add(head_cno)
            else:
                if self_concat_off:
                    impossible_concat_dict[tail_cno].add(head_cno)

            if head_cno not in impossible_concat_dict[tail_cno]:
                # final check
                src = contig_dict[tail_cno][0][-1]
                tgt = contig_dict[head_cno][0][0]
                reached, rec_path = reachable(graph, simp_node_dict, simp_node_dict[src], tail_cno, simp_node_dict[tgt], head_cno, contig_dict)
                if not reached:
                    impossible_concat_dict[tail_cno].add(head_cno)
                else:
                    rand_path_dict[(tail_cno, head_cno)] = rec_path

    all_contig_ids = contig_dict.keys()
    contig_concat_plans = {}
    for key, item in impossible_concat_dict.items():
        contig_concat_plans[key] = all_contig_ids - item
    return contig_concat_plans, rand_path_dict

def contig_clique_graph_build(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, min_cov, min_len, max_len, overlap, tempdir):
    cliq_graph = Graph(directed=True)
    cliq_graph.vp.cno = cliq_graph.new_vertex_property("string", val="")
    cliq_graph.vp.clen = cliq_graph.new_vertex_property("int32_t")
    cliq_graph.vp.ccov = cliq_graph.new_vertex_property("double")
    cliq_graph.vp.text = cliq_graph.new_vertex_property("string")

    cliq_node_dict = {}
    cliq_edge_dict = {}
    
    for cno, (contig, clen, ccov) in contig_dict.items():
        contig_node = cliq_graph.add_vertex()
        cliq_graph.vp.cno[contig_node] = cno
        cliq_graph.vp.clen[contig_node] = clen
        cliq_graph.vp.ccov[contig_node] = ccov
        cliq_graph.vp.text[contig_node] = cno + ":" + str(clen) + ":" + str(round(ccov))
        cliq_node_dict[cno] = contig_node

    concat_plan, rand_path_dict = allowed_concat_init(graph, contig_dict, simp_node_dict, max_len, overlap)
    
    for tail_cno, head_cnos in concat_plan.items():
        print("tail cno: ", tail_cno, " can concat with following head cnos: ", head_cnos)
        src_node = cliq_node_dict[tail_cno]
        for head_cno in head_cnos:
            tgt_node = cliq_node_dict[head_cno]
            contig_edge = cliq_graph.add_edge(src_node, tgt_node)
            cliq_edge_dict[(tail_cno, head_cno)] = contig_edge
    
    # further filter out the case when shortest path concat still overbound
    # for (src_cno, tgt_cno), e in cliq_edge_dict.items():
    #     return

    graph_draw(g=cliq_graph, output="{0}cliq_graph.png".format(tempdir), bg_color="white", vertex_text=cliq_graph.vp.text, vertex_size=20, vertex_font_size=15, edge_font_size=15, output_size=(900, 900))


    return cliq_graph, cliq_node_dict, cliq_edge_dict, rand_path_dict


# def contig_pairwise_concatenation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, min_cov, min_len, max_len, overlap):

#     concat_plan = allowed_concat_init(contig_dict, simp_node_dict, max_len)
#     for tail_cno, pairs in concat_plan.items():
#         print("tail cno: ", tail_cno, " can concat with following head cnos: ", pairs)
#         for head_cno in pairs:
#             print("head cno: ", head_cno)
#             src = contig_dict[tail_cno][0][-1]
#             tgt = contig_dict[head_cno][0][0]
#             p = distance_search(graph, simp_node_dict, contig_dict[tail_cno][0], src, contig_dict[head_cno][0], tgt, overlap)
    
#     return

def reachable(graph: Graph, simp_node_dict: dict, src, src_cno, tgt, tgt_cno, contig_dict: dict):
    """
    determine whether src can possibly reach the tgt
    """
    print("reachable check: {0} - {1}".format(graph.vp.id[src], graph.vp.id[tgt]))
    visited = {}
    for no in simp_node_dict.keys():
        visited[no] = False
    for c in contig_dict[src_cno][0]:
        visited[c] = True
    for c in contig_dict[tgt_cno][0]:
        visited[c] = True

    visited[graph.vp.id[src]] = False
    visited[graph.vp.id[tgt]] = False

    queue = [src]
    reached = False
    while queue:
        curr = queue.pop()
        visited[graph.vp.id[curr]] = True
        if curr == tgt:
            reached = True
            break
        for out in curr.out_neighbors():
            if not visited[graph.vp.id[out]]:
                graph.vp.prev[out] = graph.vp.id[curr]
                queue.append(out)
    
    rec_path = None

    if not reached:
        print("not reachable")
    else:
        rec_path = [graph.vp.id[tgt]]
        node = tgt
        while graph.vp.prev[node] != graph.vp.id[src]:
            prev_id = graph.vp.prev[node]
            rec_path.insert(0, prev_id)
            node = simp_node_dict[prev_id]
        rec_path.insert(0, graph.vp.id[src])
        print("PATH: ", list_to_string(rec_path))
        for cno, (contig, _, _) in contig_dict.items():
            if cno in [src_cno, tgt_cno]:
                # self involved, tolerant
                continue
            if all(c in rec_path for c in contig):
                reached = False
                print("cno: {0} is contained in the path".format(cno))

    # clean prev
    for node in simp_node_dict.values():
        graph.vp.prev[node] = ""

    if reached:
        print("reachable")
    return reached, rec_path

def distance_search(graph: Graph, simp_node_dict: dict, contig_covered_node_ids: set, source, sink, max_len: int, overlap: int):
    """
    Compute minimal distance and its path between source node and sink node
    optimise the function with contig overlap check TODO
    FIXME re-define shortest path
    """
    def dfs_helper(graph: graph, u, sink, visited):
        """
        Return None if no shortest path is founded
        """
        if graph.vp.id[u] == sink:
            return [u]
        elif u.out_degree() == 0:
            return None
        else:
            ss_path = None
            for v in u.out_neighbors():
                if not visited[v] and graph.vp.id[v] in simp_node_dict:
                    visited[v] = True
                    cmp_path = dfs_helper(graph, v, sink, visited)
                    if cmp_path != None:
                        # path lead to sink
                        cmp_path.insert(0, u)
                        cmp_len = path_len(graph, cmp_path, overlap)
                        if cmp_len < max_len:
                            if ss_path == None:
                                ss_path = cmp_path
                            else:
                                ss_len = path_len(graph, ss_path, overlap)
                                ss_path = ss_path if ss_len < cmp_len else cmp_path
                    visited[v] = False
        return ss_path

    print("source: ", source, "sink: ", sink)
    s_path = None
    s_len = 0
    print("start ssp")
    visited = {}
    for u in graph.vertices():
        visited[u] = False
    # avoid double contig path
    for s in contig_covered_node_ids:
        visited[simp_node_dict[s]] = True
    # avoid cycle
    visited[simp_node_dict[source]] = True
    visited[simp_node_dict[sink]] = False

    u = simp_node_dict[source]

    s_path = dfs_helper(graph, u, sink, visited)
    print(path_to_id_string(graph, s_path, "path: ") if s_path != None else "path: ")
    if s_path == None:
        print("Path not found")
    elif len(s_path) >= 2:
        s_path = s_path[1:-1]
        # compute directed path between ith tail to tth head (distance = len * #nodes in the path - (#nodes in the path - 1) * overlap)
        s_len = path_len(graph, s_path, overlap)
        print("shortest path found")
    else:
        s_path = None
        print("error path found")
    return s_path, s_len

def strain_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strain_dict: dict, min_cov, max_len, overlap):
    """
    Extend the strain length
    1. map all the strain back into the graph
    """

    return

if __name__ == "__main__":
    main()