#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
# from graph_tool import GraphView, _in_degree
# import graph_tool
# from graph_tool.search import dfs_iterator
from graph_tool.topology import all_circuits, all_shortest_paths


import subprocess
from graph_tool.all import Graph
from graph_tool.topology import is_DAG
import argparse
import numpy

from graph_converter import *

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = True
TEMP_DIR = "acc/"

def main():
    """
    --------------------------------------------OVERALL FLOW----------------------------------------
    Input: Graph, contig
    operation:
    -> START
    -> flip graph [DONE]
    -> increment node dp based on contig [DONE] 
    Output --> graph_L0.gfa
    -> tip removal based on minimap2 [DONE] 
    Output --> t_graph_L1.gfa
    -> bubble removal based on minimap2 [TODO] 
    Output --> bt_graph_L2.gfa
    -> cycle detection and node partition
    Output --> c_graph_L3p.gfa, nc_graph_L3p.gfa
    -> node depth rebalance + assign edge flow [DONE] 
    Output --> dbt_graph_L3.gfa
    -> compact the reassigned/splitted contig into a single node. [TODO] 
    Output --> cdbt_graph_L4.gfa
    -> compact all the simple path if exist. [DONE] 
    Output --> scdbt_graph_L5.gfa
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
        TEMP_DIR = args.output_dir
    
    subprocess.check_call("rm -rf {0} && mkdir {0}".format(TEMP_DIR), shell=True)

    print("----------------------------------INPUT---------------------------------------")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    contig_dict, node_to_contig_dict, edge_to_contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_len)

    print("--------------------------CONTIG RELATED NODE FIX-----------------------------")
    # contig_node_cov_rise(graph, simp_node_dict, contig_dict, node_to_contig_dict)

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
    
    print("--------------------------------BUBBLE REMOVAL--------------------------------------")
    bubble_removal()
    
    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}bt_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}bt_graph_L2.gfa".format(TEMP_DIR))

    print("------------------------------NODE PARTITION-----------------------------------")
    # store the no-cycle nodes in nc_graph_L3p.gfa
    noncyc_nodes = None
    simple_paths = None
    if not is_DAG:
        noncyc_nodes, simple_paths = node_partition(graph2, simp_node_dict2, simp_edge_dict2, TEMP_DIR)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth use been store in the dict
    # ratio: normalised balanced node depth / previous node depth
    _, _, _, ratio = coverage_rebalance(graph2, simp_node_dict2, simp_edge_dict2)
    
    if noncyc_nodes != None and simple_paths != None:
        print("rebalance linear subgraph now..")
        graphnc, simp_node_dictnc, simp_edge_dictnc = gfa_to_graph("{0}nc_graph_L3p.gfa".format(TEMP_DIR))
        coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc)
        print("Done, start coverage merge")

        for no, node in simp_node_dictnc.items():
            cnode = simp_node_dict2[no]
            merge_dp = graphnc.vp.dp[node] * ((ratio - 1) / ratio) + graph2.vp.dp[cnode]
            graph2.vp.dp[cnode] = merge_dp   
    else:
        print("no linear subgraph available..")


    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}dbt_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}dbt_graph_L3.gfa".format(TEMP_DIR))
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)

    print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    # re-evaluate the contig coverage
    contig_dict_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict, args.overlap)
    if DEBUG_MODE:
        for cno, (contig, clen, ccov) in contig_dict.items():
            print("---------------------------------------------------------------")
            print_contig(cno, clen, ccov, contig)

    print("-------------------------------GRAPH COMPACTIFICATION-----------------------------------")
    ## Compact all the contig into single node and saved to contig_node_dict, Store once, then compact
    ## the rest of the simple paths.

    contig_node_dict = contig_compactification(graph3, simp_node_dict3, simp_edge_dict3, contig_dict, args.overlap)
    graph_to_gfa(graph3, simp_node_dict3, simp_edge_dict3, "{0}cdbt_graph_L4.gfa".format(TEMP_DIR))
    graph4, simp_node_dict4, simp_edge_dict4 = flipped_gfa_to_graph("{0}cdbt_graph_L4.gfa".format(TEMP_DIR))
    assign_edge_flow(graph4, simp_node_dict4, simp_edge_dict4)

    simp_path_dict = simple_paths_to_dict(graph4, simp_node_dict4, simp_edge_dict4, contig_node_dict, args.overlap)
    contig_compactification(graph4, simp_node_dict4, simp_edge_dict4, simp_path_dict, args.overlap)
    
    if DEBUG_MODE:
        for cno, node_id in contig_node_dict.items():
            print_vertex(graph4, simp_node_dict4[node_id], cno)
    
    graph_to_gfa(graph4, simp_node_dict4, simp_edge_dict4, "{0}scdbt_graph_L5.gfa".format(TEMP_DIR))
    graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}scdbt_graph_L5.gfa".format(TEMP_DIR))
    assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
    
    print("-------------------------------CONTIG PAIRWISE CONCATENATION-----------------------------------")



def contig_dict_fix(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
    """
    fix the contig dict, reassign the contig coverage to minimum used **edge flow** coverage,
    we not use the minimum used **node** since the min node may still be shared by multiple contigs.
    however, min edge flow may still not be the true contig coverage

    """
    for cno, (contig, _, _) in list(contig_dict.items()):
        print("------------------------cno: {0}".format(cno))
        if all([no in simp_node_dict for no in contig]):
            contig_dict[cno][2] = round(numpy.min(contig_flow(graph, simp_edge_dict, contig)), 2)
        else:
            subcontigs = []
            curr_contig = []
            addLast = False
            for no in contig:
                if no in simp_node_dict:
                    curr_contig.append(no)
                    addLast = True
                else:
                    addLast = False
                    if curr_contig != []:
                        subcontigs.append(curr_contig[:])
                    curr_contig = []
            if addLast:
                subcontigs.append(curr_contig[:])
            contig_dict.pop(cno)

            for i, subc in enumerate(subcontigs):
                sublen = path_len(graph, [simp_node_dict[c] for c in subc], overlap)
                subcov = round(numpy.min(contig_flow(graph, simp_edge_dict, subc)), 2)
                contig_dict[cno + "_" + str(i)] = (subc, sublen, subcov)
    return

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

def contig_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
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
        
        contig_reduction(graph, contig, cno, clen, ccov, simp_node_dict, simp_edge_dict, 10)

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
        for (u,_) in in_edges:
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                ue = graph.add_edge(simp_node_dict[u], node)
                graph.ep.overlap[ue] = overlap
                graph.ep.color[ue] = 'black'
                simp_edge_dict[(u, graph.vp.id[node])] = ue
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict:
                    ue = graph.add_edge(in_node, node)
                    graph.ep.overlap[ue] = overlap
                    graph.ep.color[ue] = 'black'
                    simp_edge_dict[(graph.vp.id[in_node], graph.vp.id[node])] = ue                 

        for (_,v) in out_edges:
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                ve = graph.add_edge(node, simp_node_dict[v])
                graph.ep.overlap[ve] = overlap
                graph.ep.color[ve] = 'black'
                simp_edge_dict[(graph.vp.id[node], v)] = ve
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict:
                    ve = graph.add_edge(node, out_node)
                    graph.ep.overlap[ve] = overlap
                    graph.ep.color[ve] = 'black'
                    simp_edge_dict[(graph.vp.id[node], graph.vp.id[out_node])] = ve
    return contig_node_dict

def contig_node_cov_rise(graph: Graph, simp_node_dict: dict, contig_dict: dict, node_to_contig_dict: dict):
    """
    for any node that involved in one or more contigs, rise the depth if less than the sum of related contig cov
    """
    for no, cnos in node_to_contig_dict.items():
        sum_covs = numpy.sum([contig_dict[cno][2] for cno in cnos])
        node = simp_node_dict[no]
        if sum_covs > graph.vp.dp[node]:
            if DEBUG_MODE:
                print("Node: {0} dp is really low: {1} vs {2}, rise it up".format(no, graph.vp.dp[node], sum_covs))
            graph.vp.dp[node] = sum_covs
    return

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, single_iter=False):
    def expectation_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
        for (u,v),e in simp_edge_dict.items():
            u_node = simp_node_dict[u]
            v_node = simp_node_dict[v]
            flow = 0
            if u_node.out_degree() == 1 and v_node.in_degree() == 1:
                flow = max(graph.vp.dp[u_node], graph.vp.dp[v_node])
            elif u_node.out_degree() > 1 and v_node.in_degree() == 1:
                u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
                flow = max(graph.vp.dp[v_node], (graph.vp.dp[v_node]/u_out_sum)*graph.vp.dp[u_node])
            elif u_node.out_degree() == 1 and v_node.in_degree() > 1:
                v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
                flow = max(graph.vp.dp[u_node], (graph.vp.dp[u_node]/v_in_sum)*graph.vp.dp[v_node])
            else:
                u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
                v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
                flow = max((graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node])
            graph.ep.flow[e] = max(graph.ep.flow[e], flow)
            # print_edge(graph, e, "{0}".format(flow))
        return
    def maximization_node_depth(graph: Graph, simp_node_dict: dict):
        is_update = False
        for no, node in simp_node_dict.items():
            us = list(node.in_neighbors())
            u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
            in_neighbor_dp_sum = -1
            if u_out_degrees == node.in_degree():
                in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

            vs = list(node.out_neighbors())
            v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 

            out_neighbor_dp_sum = -1
            if v_in_degrees == node.out_degree():
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

    if is_DAG(graph):
        cutoff = 0.0001
    else:
        cutoff = 0.001 * len(simp_node_dict)
    sum_delta = 0
    sum_depth_before = numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()])
    
    # store previous node depth
    prev_dp_dict = {}
    for no, v in simp_node_dict.items():
        prev_dp_dict[no] = graph.vp.dp[v]
    if single_iter:
        expectation_edge_flow(graph, simp_node_dict, simp_edge_dict)
        maximization_node_depth(graph, simp_node_dict)
        return prev_dp_dict, {}, {}, 1
    print(cutoff)
    is_update = True
    while is_update:
        # E Step
        expectation_edge_flow(graph, simp_node_dict, simp_edge_dict)
        # validation
        sum_delta = 0
        for no, node in simp_node_dict.items():
            inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
            outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
            if node.in_degree() == 0 or node.out_degree() == 0:
                continue
            else:
                sum_delta += (abs(inflow - outflow))/((inflow + outflow)/2)
        if round(sum_delta, 2) < cutoff:
            break
        # M Step
        is_update = maximization_node_depth(graph, simp_node_dict)

    # final evaluation
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Ratio: ", sum_ratio, "Delta: ", sum_delta)
    for node in simp_node_dict.values():
        graph.vp.dp[node] = graph.vp.dp[node] / sum_ratio
    node_ratio_dict = {}
    for no in prev_dp_dict.keys():
        node_ratio_dict[no] = (graph.vp.dp[simp_node_dict[no]] / prev_dp_dict[no])
    curr_dp_dict = {}
    for no, v in simp_node_dict.items():
        curr_dp_dict[no] = graph.vp.dp[v]

    return prev_dp_dict, curr_dp_dict, node_ratio_dict, sum_ratio

def graph_simplification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, min_cov):
    """
    Directly remove all the vertex with coverage less than minimum coverage and related edge

    Node belongs to any contigs should not be removed
    return:
        removed_node_dict
        removed_edge_dict
    """
    print("-------------------------graph simplification----------------------")
    print("Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    removed_node_dict = {}
    removed_edge_dict = {}
    # iterate until no more node be removed from the graph
    for id, node in list(simp_node_dict.items()):
        if graph.vp.dp[node] < min_cov:
            if id in node_to_contig_dict:
                if DEBUG_MODE:
                    print("node: {0} should not be removed although with ccov: {1}".format(id, graph.vp.dp[node]))
                continue
            if DEBUG_MODE:
                print_vertex(graph, node, "Node removed by graph simplification -")

            # delete the node
            simp_node_dict.pop(id)
            graph.vp.color[node] = 'gray'
            removed_node_dict[id] = node
            total_reduce_depth = graph.vp.dp[node]
            # delete related edges
            for out_node in node.out_neighbors():
                out_id = graph.vp.id[out_node]
                if (id, out_id) in edge_to_contig_dict:
                    if DEBUG_MODE:
                        print("edge: {0} should not be removed".format((id, out_id)))
                    continue
                if (id, out_id) in simp_edge_dict:
                    e = simp_edge_dict[(id, out_id)]
                    graph.ep.color[e] = 'gray'
                    simp_edge_dict.pop((id, out_id))
                    removed_edge_dict[(id, out_id)] = e

            for in_node in node.in_neighbors():
                in_id = graph.vp.id[in_node]
                if (in_id, id) in edge_to_contig_dict:
                    if DEBUG_MODE:
                        print("edge: {0} should not be removed".format((in_id, id)))
                    continue
                if (in_id, id) in simp_edge_dict:
                    e = simp_edge_dict[(in_id, id)]
                    graph.ep.color[e] = 'gray'
                    simp_edge_dict.pop((in_id, id))
                    removed_edge_dict[(in_id, id)] = e
    print("Remain: Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    print("-------------------------graph simplification end----------------------")
    return removed_node_dict, removed_edge_dict

def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    for (u,v),e in simp_edge_dict.items():
        u_node = simp_node_dict[u]
        v_node = simp_node_dict[v]
        flow = 0
        if u_node.out_degree() == 1 and v_node.in_degree() == 1:
            flow = max(graph.vp.dp[u_node], graph.vp.dp[v_node])
        elif u_node.out_degree() > 1 and v_node.in_degree() == 1:
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
            flow = max(graph.vp.dp[v_node], (graph.vp.dp[v_node]/u_out_sum)*graph.vp.dp[u_node])
        elif u_node.out_degree() == 1 and v_node.in_degree() > 1:
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
            flow = max(graph.vp.dp[u_node], (graph.vp.dp[u_node]/v_in_sum)*graph.vp.dp[v_node])
        else:
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
            flow = max((graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node])
        graph.ep.flow[e] = round(max(graph.ep.flow[e], flow), 2)
    return

def contig_reduction(graph: Graph, contig, cno, clen, ccov, simp_node_dict: dict, simp_edge_dict: dict, min_cov):
    """
    reduce and update the graph by given contig
    """
    print("*---Contig reduction: ", cno, clen, ccov)
    next_node_index = 1
    for node in contig:
        if next_node_index >= len(contig):
            break
        adj_node = contig[next_node_index]

        u = simp_node_dict[node] if node in simp_node_dict else None
        v = simp_node_dict[adj_node] if adj_node in simp_node_dict else None
        e = simp_edge_dict[(node,adj_node)] if (node,adj_node) in simp_edge_dict else None 

        # edge may be eliminated from previous execution already
        if u == None or v == None or e == None:
            if e != None:
                graph.ep.flow[e] = 0
                graph.ep.color[e] = 'gray'
                simp_edge_dict.pop((node,adj_node))
                print_edge(graph, e, "edge removed 1st")
            continue
        if DEBUG_MODE:
            print_edge(graph, e, "current edge eval")
        # reduce the depth for involving edge
        if graph.ep.flow[e] - ccov <= min_cov:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            print_edge(graph, e, "edge removed 2nd")
        else:
            graph.ep.flow[e] = round(graph.ep.flow[e] - ccov, 2)

        # reduce the depth for involving node, gray color for forbiddened node
        if graph.vp.dp[u] - ccov <= min_cov:
            graph.vp.dp[u] = 0
            graph.vp.color[u] = 'gray'
            simp_node_dict.pop(node)
        else:
            graph.vp.dp[u] = round(graph.vp.dp[u] - ccov, 2)

        # last node in the contig, reduce its depth
        if next_node_index == len(contig) - 1:
            if graph.vp.dp[v] - ccov <= min_cov: 
                graph.vp.dp[v] = 0
                graph.vp.color[v] = 'gray'
                simp_node_dict.pop(adj_node)
            else:
                graph.vp.dp[v] = round(graph.vp.dp[v] - ccov, 2)
        
        # update edges
        if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (node,adj_node) in simp_edge_dict:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            print_edge(graph, e, "edge removed 3rd")

        next_node_index = next_node_index + 1
    return


### Legacy ...



def distance_search(graph: Graph, simp_node_dict: dict, node_usage_dict: dict, src_contig, source, sink_contig, sink, overlap: int):
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
                        if ss_path == None:
                            ss_path = cmp_path
                        else:
                            #FIXME fixed, less len, more cov, less usage
                            ss_path = sorted([cmp_path, ss_path], key=lambda p: (path_usage(graph, node_usage_dict, p), path_len(graph, p, overlap)))[0]
                    visited[v] = False
        return ss_path

    print("source: ", source, "sink: ", sink)
    s_path = None
    s_len = 0
    if source not in simp_node_dict or sink not in simp_node_dict:
        print("source/sink not found, error -1")
    elif simp_node_dict[source].out_degree() == 0:
        print("source has 0 out degree")
    elif simp_node_dict[sink].in_degree() == 0:
        print("sink has 0 in degree")
    else:
        print("start ssp")
        visited = {}
        for u in graph.vertices():
            visited[u] = False
        
        # avoid cycle
        visited[simp_node_dict[source]] = True
        # avoid double contig path
        for s in src_contig:
            visited[simp_node_dict[s]] = True
        for s in sink_contig:
            visited[simp_node_dict[s]] = True

        u = simp_node_dict[source]

        s_path = dfs_helper(graph, u, sink, visited)
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

def contig_merge(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, node_usage_dict: dict, output_file, min_cov, min_len, max_len, overlap):

    print("-----------------------contig pair-wise concatenaton--------------------")

    is_linear_strain = is_DAG(graph)
    if not is_linear_strain:
        print("graph is cyclic, cyclic strian may exist")
    else:
        print("graph is not cyclic, lienar strain may exist")

    # find pair-wise shortest path among contigs and construct the pairwise contig path matrix
    # TODO if sp found, do we also reduce its coverage?
    dist_matrix = {}
    concat_dicision = {}
    for tail_cno, [tail_contig, tail_clen, tail_ccov] in temp_contigs_dict.items():
        for head_cno, [head_contig, head_clen, head_ccov] in temp_contigs_dict.items():
            print("------------------------------------------------------")
            print("Tail Contig: ", tail_cno, " -> Head Contig: ", head_cno)
            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len:
                    print("Contig ", tail_cno, "-", "Contig ", head_cno, " with len over maxlen")
                    continue
                intersect = list(set(tail_contig) & set(head_contig))
                if intersect != []:
                    print("Contig ", tail_cno, "-", "Contig ", head_cno, " Intersection with ", len(intersect), " nodes")
                    if DEBUG_MODE:
                        print(intersect)
                    continue
            else:
                # definitely no self-to-self path.
                if is_linear_strain:
                    print("No cyclic strain possible, skip")
                    continue
            print("start path construction")
            s_path, s_len = distance_search(graph, simp_node_dict, node_usage_dict, tail_contig[:-1], tail_contig[-1], head_contig[1:], head_contig[0], overlap)
            if s_path != None:
                s_path_ids = [graph.vp.id[v] for v in s_path]
                print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])
                s_path_edge_flow = contig_flow(graph, simp_edge_dict,  s_path_ids)
                s_path_ccov = numpy.mean(s_path_edge_flow) if len(s_path_edge_flow) != 0 else 0
                dist_matrix[(tail_cno, head_cno)] = (s_path_ids, s_len, s_path_ccov)

                concat_ccov = min(head_ccov, tail_ccov, s_path_ccov) if s_path_ccov != 0.0 else min(head_ccov, tail_ccov)
                concat_len = head_clen + tail_clen - overlap if s_len == 0 else head_clen + s_len + tail_clen - overlap * 2
                print("coverage: head contig eflow: ", head_ccov, " s path eflow: ", s_path_ccov, " tail contig eflow: ", tail_ccov)
                print("potential concatenated length: ", concat_len)

                if concat_len > max_len:
                    print("exceed maximum strain length")
                else:
                    print("length satisfied, concat_ccov: ", concat_ccov)
                    concat_dicision[(tail_cno, head_cno)] = concat_ccov
            else:
                print("shortest path not found")
    print("------------------------------------------------------")
    
    # TODO re-evaluate the concat dicision, sort the concat dicision via concat_ccov, with reverse order or not, TBD
    sorted_pair = sorted(concat_dicision.items(), key=lambda x:x[1])
    concat_dicision_s = [t[0] for t in sorted_pair]

    #start concatenation
    # FIXME fix the concat dicision, test for used contig only concat n times.
    concat_strain_dict = {}
    concat_contig_dict = {}
    skip_key = set()
    used_contig = set()
    #FIXME not only used_contig, but also the nodes been used throughout the path.
    for (tail_cno, head_cno) in concat_dicision_s:
        print("------------------------------------------------------")
        print("Concatentate contigs: ", tail_cno, " <-> ", head_cno)
        [tail_contig, tail_clen, tail_ccov] = temp_contigs_dict[tail_cno]
        [head_contig, head_clen, head_ccov] = temp_contigs_dict[head_cno]

        if head_cno in used_contig:
            print("contig {0} be used from previous concatenation, skip".format(head_cno))
            continue
        if tail_cno in used_contig:
            print("contig {0} be used from previous concatenation, skip".format(tail_cno))
            continue

        if (head_cno, tail_cno) in skip_key:
            continue
        if (tail_cno, head_cno) in skip_key:
            continue
        skip_key.add((head_cno, tail_cno))
        skip_key.add((tail_cno, head_cno))

        is_linear = False
        if (head_cno, tail_cno) not in concat_dicision_s:
            print("no reverse concatenation exists, potential lienar strain")
            is_linear = True
        else: 
            print("access contig in both direction, potential cyclic strain")
            is_linear = False
        
        (s_path_ids_l, s_len_l, s_path_ccov_l) = dist_matrix[(tail_cno, head_cno)]

        concat_cno = tail_cno + "_" + head_cno
        overlap_count = 1 if s_len_l == 0 else 2
        if is_linear:
            #linear strain concatenation
            concat_clen = tail_clen + s_len_l + head_clen - overlap_count * overlap
            concat_ccov = min(v for v in [head_ccov, tail_ccov, s_path_ccov_l] if v > 0.0)
            concat_c = tail_contig + s_path_ids_l + head_contig
        else:
            #cyclic strain concatentation
            (s_path_ids_r, s_len_r, s_path_ccov_r) = dist_matrix[(head_cno, tail_cno)]
            overlap_count = overlap_count + 1 if s_len_r == 0 else overlap_count + 2

            concat_clen = tail_clen + s_len_l + head_clen + s_len_r - overlap_count * overlap
            concat_ccov = min(v for v in [head_ccov, tail_ccov, s_path_ccov_l, s_path_ccov_r] if v > 0.0)
            concat_c = tail_contig + s_path_ids_l + head_contig + s_path_ids_r
        
        # reduce the pre-defined contig's ccov by concat_ccov
        if head_ccov - concat_ccov < min_cov:
            used_contig.add(head_cno)
        else:
            temp_contigs_dict[head_cno][2] = temp_contigs_dict[head_cno][2] - concat_ccov

        if tail_ccov - concat_ccov < min_cov:
            used_contig.add(tail_cno)
        else:
            temp_contigs_dict[tail_cno][2] = temp_contigs_dict[tail_cno][2] - concat_ccov

        if concat_clen >= min_len:
            concat_strain_dict[concat_cno] = (concat_c, concat_clen, concat_ccov)
        else:
            concat_contig_dict[concat_cno] = [concat_c, concat_clen, concat_ccov]
    
    # re-append rest of contig
    for cno, [contig, clen, ccov] in temp_contigs_dict.items():
        if cno not in used_contig:
            if clen >= min_len:
                print("Full len contig directly store into strain dict", cno, clen, ccov)
                used_contig.add(cno)
                concat_strain_dict[cno] = (contig, clen, ccov)
            else:
                print("Re-appending contig to temp_contig_dict: ", cno, clen, ccov)
                concat_contig_dict[cno] = [contig, clen, ccov]
        else:
            print("contig {0} is used".format(cno))
    
    udpate_node_to_contig_dict(graph, node_to_contig_dict, simp_node_dict)
    update_edge_to_contig_dict(graph, edge_to_contig_dict, simp_edge_dict)
    # add concat contig cno to dict
    for cno, [c, _, _] in concat_contig_dict.items():
        for n in c:
            if n not in node_to_contig_dict:
                node_to_contig_dict[n] = [{cno},graph.vp.dp[simp_node_dict[n]], simp_node_dict[n]]
            else:
                node_to_contig_dict[n][0].add(cno)     
           
        for i in range(len(c)):
            c_i = c[i]
            c_i_1 = c[i+1] if (i < len(c) - 1) else None
            if c_i_1 != None:
                if (c_i, c_i_1) not in edge_to_contig_dict:
                    edge_to_contig_dict[(c_i, c_i_1)] = [{cno}, graph.ep.flow[simp_edge_dict[(c_i, c_i_1)]], simp_edge_dict[(c_i, c_i_1)]]
                else:
                    edge_to_contig_dict[(c_i, c_i_1)][0].add(cno)

    # update edge_to_contig_dict
    for (u,v), [cnos, _, _] in list(edge_to_contig_dict.items()):
        for cno in used_contig:
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove used contig cno {0} from edge: {1} cno list".format(cno, (u,v)))
        for cno in concat_strain_dict.keys():
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remov cand strain cno {0} from edge: {1} cno list".format(cno, (u,v)))

        if len(cnos) == 0:
            edge_to_contig_dict.pop((u,v))
        else:
            edge_to_contig_dict[(u,v)][0] = cnos

    # update node_to_contig_dict
    for no, [cnos, _, _] in list(node_to_contig_dict.items()):
        for cno in used_contig:
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove used contig cno {0} from node: {1} cno list".format(cno, no))
        for cno in concat_strain_dict.keys():
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove cand strain cno {0} from node: {1} cno list".format(cno, no))

        if len(cnos) == 0:
            node_to_contig_dict.pop(no)
        else:
            node_to_contig_dict[no][0] = cnos

    # reduce all the full length concated contigs as cand strain
    for cno, (c, clen, ccov) in concat_strain_dict.items():
        contig_reduction(graph, c, cno, clen, ccov, simp_node_dict, simp_edge_dict, min_cov)
        increment_node_usage_dict(node_usage_dict, c, ccov)

    for no, [cnos, dp, node] in node_to_contig_dict.items():
        if no not in simp_node_dict:
            simp_node_dict[no] = node
            sumcov = numpy.sum([concat_contig_dict[c][2] for c in cnos])
            graph.vp.dp[node] = sumcov
            graph.vp.color[node] = "black"
    # recover edge
    for (u,v), [cnos, flow, edge] in edge_to_contig_dict.items():
        if (u,v) not in simp_edge_dict:
            simp_edge_dict[(u,v)] = edge
            sumcov = numpy.sum([concat_contig_dict[c][2] for c in cnos])
            graph.ep.flow[edge] = sumcov
            graph.ep.color[edge] = "black"

    for cno, [contig, clen, ccov] in concat_contig_dict.items():
        print("------------------------------------------------------")
        print_contig(cno, clen, ccov, contig, "partial length residue contig found after concatenation")
        u = simp_node_dict[contig[0]]
        v = simp_node_dict[contig[-1]]
        if u.in_degree() == 0 and v.out_degree() == 0:
            print("no further concatentation exist")

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, output_file)

    print("--------------------contig pair-wise concatenaton end--------------------")
    return concat_strain_dict, concat_contig_dict

def simp_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    find simple edges, simple edge is the only edge between its source and sink
    """
    simple_edges = []
    out_edge = {}
    in_edge = {}
    for e in simp_edge_dict.values():
        src = e.source()
        src_out_d = len([u for u in src.out_neighbors()])
        target = e.target()
        target_in_d = len([u for u in target.in_neighbors()])
        if graph.vp.id[src] not in simp_node_dict or graph.vp.id[target] not in simp_node_dict:
            continue
        if src_out_d == 1 and target_in_d == 1:
            assert src != target
            simple_edges.append([src, target])
            in_edge[int(src)] = e
            out_edge[int(target)] = e

    # build simple paths from simple edges
    def extend_path(p):
        v = int(p[-1])
        if v in in_edge:
            p.append(in_edge[v].target())
            return extend_path(p)
        else:
            return p
    simple_paths = []
    for v, e in in_edge.items():
        if v not in out_edge:
            p = extend_path([e.source(), e.target()])
            simple_paths.append(p) 
    return simple_paths

def simple_paths_to_dict(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_node_dict: dict, overlap):
    simple_paths = simp_path(graph, simp_node_dict, simp_edge_dict)
    simp_path_dict = {}
    contig_node_ids = set(contig_node_dict.values())

    for id, p in enumerate(simple_paths):
        print("path: ", [int(graph.vp.id[u]) for u in p])
        pids = [graph.vp.id[n] for n in p]
        if contig_node_ids.intersection(set(pids)):
            print("simple path forbidden, contig is involved")
            continue
        name = "0" + str(id) + "0"
        clen = path_len(graph, p, overlap)
        cov = numpy.min([graph.vp.dp[n] for n in p])
        simp_path_dict[name] = [pids, clen, cov]
    return simp_path_dict

def path_extraction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, node_usage_dict: dict, overlap, min_cov, min_len):
    """
    extract the last mile path from the graph, with support from residue contigs
    """
    print("--------------------path extraction--------------------")

    paths_per_group_dict = {}

    graph, groups = graph_grouping(graph, simp_node_dict)
    print("number of groups:", len(groups))
    final_strain_dict = {}
    
    for gno, group in groups.items():
        print("group number: ", gno, " member: ", [graph.vp.id[v] for v in group])
        if len(group) < 2:
            # single node group
            if len(group) == 1:
                node = group[0]
                pcov = graph.vp.dp[node]
                id = graph.vp.id[node]
                plen = len(graph.vp.seq[node])
                print("Single node Path: ", int(id), "path len: ", plen, "cov: ", pcov)
                # ratio%
                usage_ratio = round((node_usage_dict[id][0]/node_usage_dict[id][1])*100, 2)
                if pcov >= min_cov and usage_ratio > 0:
                    print("accept")
                    final_strain_dict[id] = ([id], plen, pcov)
            else:
                continue
        srcs = []
        sinks = []
        isolations = []
        middles = []
        paths = []
        # classify the nodes
        for u in group:
            no = graph.vp.id[u]
            if u.in_degree() == 0 and u.out_degree() == 0:
                isolations.append(no)
            elif u.in_degree() == 0:
                srcs.append(no)
            elif u.out_degree() == 0:
                sinks.append(no)
            else:
                middles.append(no)
        #FIXME potential issue: no src/sink node
        for src in srcs:
            for sink in sinks:
                print("Path extraction")
                p, plen = distance_search(graph, simp_node_dict, node_usage_dict, [], src, [], sink, overlap)
                if p != None:
                    s_path_ids = [graph.vp.id[v] for v in p]
                    s_path_ids.append(sink)
                    s_path_ids.insert(0, src)             
                    plen = plen + len(graph.vp.seq[simp_node_dict[src]]) + len(graph.vp.seq[simp_node_dict[sink]]) - 2 * overlap
                    paths.append((s_path_ids, plen))
        paths_per_group_dict[gno] = paths

    for gno, paths in paths_per_group_dict.items():
        print("Current group: ", gno)
        for p_ids, plen in paths:
            print("------------------------------------------------------")
            print("Path: ", [int(u) for u in p_ids])
            pcov = numpy.mean([graph.vp.dp[simp_node_dict[u]] for u in p_ids if u in simp_node_dict])
            pno = str(p_ids[0]) + "_" + str(p_ids[-1])
            print("path no: {0} path len: {1} path cov: {2}".format(pno, plen, pcov))

            if pcov >= min_cov:
                print("path accepted")
                final_strain_dict[pno] = (p_ids, plen, pcov)
                contig_reduction(graph, p_ids, pno, plen, pcov, simp_node_dict, simp_edge_dict, min_cov)
                increment_node_usage_dict(node_usage_dict, p_ids, pcov)

    print("------------------------------------------------------")
    print("--------------------path extraction end--------------------")
    return final_strain_dict

def strain_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strain_dict: dict, min_cov, max_len, overlap):
    """
    Extend the strain length
    1. map all the strain back into the graph
    """

    return

if __name__ == "__main__":
    main()