#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
# from graph_tool import GraphView, _in_degree
# import graph_tool
# from graph_tool.search import dfs_iterator
# from graph_tool.topology import topological_sort, all_circuits
# from graph_tool.clustering import local_clustering

import subprocess
from graph_tool.all import Graph
from graph_tool.topology import is_DAG

import argparse

import numpy

from graph_converter import *

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = False
TEMP_DIR = "acc/"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, default=100, help=("minimum coverage for strains"))
    parser.add_argument('-minlen', '--minimum_strain_length', dest='min_len', default=4000, type=int, help=("minimum strain length"))
    parser.add_argument('-maxlen', '--maximum_strain_length', dest='max_len', default=12000, type=int, help=("maximum strain length"))
    parser.add_argument('-overlap', '--vertex_overlap', dest='overlap', default=127, type=int, help=("adjacent vertex overlap in the graph"))
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, DEBUG_MODE only')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, required=True, help='output directory')
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

    #FIXME keep track of the node usage, add an attr or var for that.
    # Read in as Level -1 graph
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict)

    # store the usage info for each node.
    node_usage = {}
    for no in simp_node_dict.keys():
        node_usage[no] = 0

    contig_dict, node_to_contig_dict, edge_to_contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_cov, args.min_len, args.overlap)
    graph_simplification(graph, simp_node_dict, simp_edge_dict, node_to_contig_dict, edge_to_contig_dict, args.min_cov)
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}pre_graph.gfa".format(TEMP_DIR))
    
    # Read in as pre graph
    pre_graph, simp_node_dict_pre, simp_edge_dict_pre = flipped_gfa_to_graph("{0}pre_graph.gfa".format(TEMP_DIR))
    assign_edge_flow(pre_graph, simp_node_dict_pre, simp_edge_dict_pre)

    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict_pre, "{0}pre_graph.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))

    # selected contig from SPAdes
    contig_dict_to_fasta(graph, contig_dict, simp_node_dict, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))

    # reduce SPAdes full-length contig as init cand strain
    cand_strains_dict, temp_contigs_dict = graph_reduction(graph, contig_dict, simp_node_dict, simp_edge_dict, node_to_contig_dict, edge_to_contig_dict, "{0}graph_L0.gfa".format(TEMP_DIR), args.min_cov, args.min_len, args.overlap)


    iter_no = 0
    strain_dict = cand_strains_dict.copy()
    iter_condition = set(temp_contigs_dict.keys())
    concat_contig_dict_iter = temp_contigs_dict
    # contig pair-wise concatenation iteration
    # with iteration stopped when no further concatenation occurred
    while True:
        print("Current iteration: ", iter_no)
        graph_iter, simp_node_dict_iter, simp_edge_dict_iter = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(iter_no)))
        assign_edge_flow(graph_iter, simp_node_dict_iter, simp_edge_dict_iter)
        graph_simplification(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, node_to_contig_dict, edge_to_contig_dict, args.min_cov)

        concat_strain_dict_iter, concat_contig_dict_iter = contig_merge(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, concat_contig_dict_iter, node_to_contig_dict, edge_to_contig_dict, "{0}graph_L{1}.gfa".format(TEMP_DIR, str(iter_no+1)), args.min_cov, args.min_len, args.max_len, args.overlap)

        contig_dict_to_fasta(graph_iter, concat_contig_dict_iter, simp_node_dict_iter, args.overlap, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(iter_no)))
        minimap_api(args.ref_file, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(iter_no)), "{0}L{1}_contigs_to_strain.paf".format(TEMP_DIR, str(iter_no)))

        strain_dict.update(concat_strain_dict_iter.copy())
        iter_no = iter_no + 1

        iter_compare = set(concat_contig_dict_iter.keys())
        if iter_condition == iter_compare:
            print("end of contig pair-wise concatenation iteration: ", iter_no)
            break
        else:
            iter_condition = iter_compare


    # print strain and minimap overlaps
    contig_dict_to_fasta(pre_graph, strain_dict, simp_node_dict_pre, args.overlap, "{0}cand_strains.fasta".format(TEMP_DIR))
    
    # further step to recover the rest of strain from the reduced graph
    #
    graph_red, simp_node_dict_red, simp_edge_dict_red = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(iter_no)))
    assign_edge_flow(graph_red, simp_node_dict_red, simp_edge_dict_red)
    graph_simplification(graph_red, simp_node_dict_red, simp_edge_dict_red, node_to_contig_dict, edge_to_contig_dict, args.min_cov)

    # contig_map_node_dict = graph_compactification(graph_red, simp_node_dict_red, simp_edge_dict_red, concat_contig_dict_iter, "{0}graph_compacted.gfa".format(TEMP_DIR), args.overlap, TEMP_DIR)
    
    # extract the last mile paths from the graph
    # graph_comp, simp_node_dict_comp, simp_edge_dict_comp = flipped_gfa_to_graph("{0}graph_compacted.gfa".format(TEMP_DIR))
    # assign_edge_flow(graph_comp, simp_node_dict_comp, simp_edge_dict_comp)

    # no reason to simplify again
    # final_strain_dict = path_extraction(graph_comp, simp_node_dict_comp, simp_edge_dict_comp, contig_map_node_dict, args.overlap, args.min_cov, args.min_len)   
    final_strain_dict = path_extraction(graph_red, simp_node_dict_red, simp_edge_dict_red, args.overlap, args.min_cov, args.min_len)   

    # with open("{0}cand_strains.fasta".format(TEMP_DIR), 'a') as fasta:
    #     for cno, (c, clen, ccov) in final_strain_dict.items():
    #         seq = path_ids_to_seq(graph_comp, c, cno, simp_node_dict_comp, args.overlap)
    #         seq += "\n"
    #         name = ">" + str(cno) + "_" + str(clen) + "_" + str(ccov) + "\n"
    #         fasta.write(name)
    #         fasta.write(seq)
    #     fasta.close()

    minimap_api(args.ref_file, "{0}cand_strains.fasta".format(TEMP_DIR), "{0}cand_strain_to_strain.paf".format(TEMP_DIR))


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
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for no, node in simp_node_dict.items():
            w = graph.vp.dp[node]
            if node.in_degree() == 1:
                for e in node.in_edges():
                    if graph.ep.flow[e] == 0.0 and (graph.vp.id[e.target()], no) in simp_edge_dict:
                        graph.ep.flow[e] = w
                        un_assigned_edge = un_assigned_edge - 1
            if node.out_degree() == 1:
                for e in node.out_edges():
                    if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[e.target()]) in simp_edge_dict:
                        graph.ep.flow[e] = w
                        un_assigned_edge = un_assigned_edge - 1
    print("un-assigned edges after init iteration : ", un_assigned_edge)

    # coverage iteration
    converage_flag = 0
    while True:          
        for no, node in simp_node_dict.items():
            in_d = node.in_degree()
            in_w = graph.vp.dp[node]
            in_e = []
            for e in node.in_edges():
                f = graph.ep.flow[e]
                if f != 0.0:
                    in_d = in_d - 1
                    in_w = in_w - f
                else:
                    in_e.append(e)

            out_d = node.out_degree()
            out_w = graph.vp.dp[node]
            out_e = []
            for e in node.out_edges():
                f = graph.ep.flow[e]
                if f != 0.0:
                    out_d = out_d - 1
                    out_w = out_w - f
                else:
                    out_e.append(e)

            if in_d == 1:
                for e in in_e:
                    if graph.ep.flow[e] == 0.0 and (graph.vp.id[e.source()], no) in simp_edge_dict:
                        graph.ep.flow[e] = in_w
                        un_assigned_edge = un_assigned_edge - 1
            if out_d == 1:
                for e in out_e:
                    if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[e.target()]) in simp_edge_dict:
                        graph.ep.flow[e] = out_w
                        un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  

    print("un-assigned edges after node-weight coverage iteration : ", un_assigned_edge)
    for (_,_), e in simp_edge_dict.items():
        if graph.ep.flow[e] == 0.0:
            u = e.source()
            v = e.target()
            u_flow_remain = graph.vp.dp[u]
            u_degree = u.out_degree()

            v_flow_remain = graph.vp.dp[v]
            v_degree = v.in_degree()
            for ue in u.out_edges():
                if graph.ep.flow[ue] != 0.0:
                    # assigned edge
                    u_degree = u_degree - 1
                    u_flow_remain = u_flow_remain - graph.ep.flow[ue]

            for ve in v.in_edges():
                if graph.ep.flow[ve] != 0.0:
                    # assigned edge
                    v_degree = v_degree - 1
                    v_flow_remain = v_flow_remain - graph.ep.flow[ve]
            u_flow_remain = u_flow_remain if u_flow_remain > 0 else 0
            v_flow_remain = v_flow_remain if v_flow_remain > 0 else 0
            if u_flow_remain == 0 or v_flow_remain == 0:
                assign_flow = (u_flow_remain / u_degree) + (v_flow_remain / v_degree)
            else:
                assign_flow = ((u_flow_remain / u_degree) + (v_flow_remain / v_degree)) / 2
            if assign_flow <= 0:
                print("low edge flow error: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                print("u_flow_remain: ", u_flow_remain, " u_degree: ", u_degree)
                print("v_flow_remain: ", v_flow_remain, " v_degree: ", v_degree)
            else:
                if DEBUG_MODE:
                    print("manual assigned edge: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                if graph.ep.flow[e] == 0.0:
                    un_assigned_edge = un_assigned_edge - 1
                    graph.ep.flow[e] = assign_flow

    print("un-assigned edges after manual assign iteration : ", un_assigned_edge)
    print("-----------------------assign edge flow end--------------------")

def contig_reduction(graph: Graph, contig, cno, clen, ccov, simp_node_dict: dict, simp_edge_dict: dict, min_cov):
    """
    reduce and update the graph by given contig
    """
    if DEBUG_MODE:
        print("*---Contig: ", cno, clen, ccov)
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
                if DEBUG_MODE:
                    print_edge(graph, e, "edge been removed due to either u or v is removed")
            else:
                if DEBUG_MODE:
                    print("edge: ", node, " -> ", adj_node, "already been removed")
            continue
        if DEBUG_MODE:
            print_edge(graph, e, "current edge eval")
        # reduce the depth for involving edge
        if graph.ep.flow[e] - ccov <= min_cov:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            if DEBUG_MODE:
                print_edge(graph, e, "edge been removed")
        else:
            graph.ep.flow[e] = graph.ep.flow[e] - ccov

        # reduce the depth for involving node, gray color for forbiddened node
        if graph.vp.dp[u] - ccov <= min_cov:
            graph.vp.dp[u] = 0
            graph.vp.color[u] = 'gray'
            simp_node_dict.pop(node)
            if DEBUG_MODE:
                print("node ", node, "been removed")
        else:
            graph.vp.dp[u] = graph.vp.dp[u] - ccov

        # last node in the contig, reduce its depth
        if next_node_index + 1 == len(contig):
            if graph.vp.dp[v] - ccov <= min_cov: 
                graph.vp.dp[v] = 0
                graph.vp.color[v] = 'gray'
                simp_node_dict.pop(adj_node)
                if DEBUG_MODE:
                    print("node ", adj_node, "been removed")
            else:
                graph.vp.dp[v] = graph.vp.dp[v] - ccov
        
        # update edges
        if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (node,adj_node) in simp_edge_dict:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            if DEBUG_MODE:
                print_edge(graph, e, "edge been removed in the final step")
        next_node_index = next_node_index + 1
    return

def graph_reduction(graph: Graph, contig_dict: dict, simp_node_dict: dict, simp_edge_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, output_file, min_cov, min_len, overlap):
    """
    reduce the node/edge weight based on existing contig found by SPAdes.
    only contig with minimum strain length satisfied be removed
    return non-removed contig dict
    """

    print("-------------------------graph reduction----------------------")
    udpate_node_to_contig_dict(node_to_contig_dict, simp_node_dict)
    update_edge_to_contig_dict(edge_to_contig_dict, simp_edge_dict)
    cand_strains_dict = {}
    temp_contigs_dict = {}
    for cno, (contig, clen, ccov) in contig_dict.items():
        if clen >= min_len:
            cand_strains_dict[cno] = (contig, clen, ccov)
            print("full-length contig found: ", cno, clen, ccov)
            contig_reduction(graph, contig, cno, clen, ccov, simp_node_dict, simp_edge_dict, min_cov)  
        else:
            temp_contigs_dict[cno] = [contig, clen, ccov]
            print("imcomplete contig found: ", cno, clen, ccov) 

    # reappend nodes
    for no, [cnos, dp, node] in list(node_to_contig_dict.items()):

        for strain_cno in cand_strains_dict.keys():
            if strain_cno in cnos:
                cnos.remove(strain_cno)

        # update cnos
        node_to_contig_dict[no][0] = cnos
        if len(cnos) == 0:
            node_to_contig_dict.pop(no)
            continue
        
        if no not in simp_node_dict and no in node_to_contig_dict:
            print("contig node {0} be removed, append it back".format(no))
            graph.vp.dp[node] = dp
            graph.vp.color[node] = "black"
            simp_node_dict[no] = node
            print_vertex(graph, node, "from graph reduction: ")

    # reappends edges
    for (u,v), [cnos, flow, edge] in list(edge_to_contig_dict.items()):

        for strain_cno in cand_strains_dict.keys():
            if strain_cno in cnos:
                cnos.remove(strain_cno)

        # update cnos
        edge_to_contig_dict[(u,v)][0] = cnos
        if len(cnos) == 0:
            edge_to_contig_dict.pop((u,v))
            continue
        
        if (u,v) not in simp_edge_dict and (u,v) in edge_to_contig_dict:
            print("contig edge {0} be removed, append it back".format((u,v)))
            graph.ep.flow[edge] = flow
            graph.ep.color[edge] = "black"
            simp_edge_dict[(u,v)] = edge

    # store output graph 
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, output_file)

    print("-----------------------graph reduction end--------------------")
    return cand_strains_dict, temp_contigs_dict

def distance_search(graph: Graph, simp_node_dict: dict, source, sink, overlap: int):
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
                if not visited[v]:
                    visited[v] = True
                    cmp_path = dfs_helper(graph, v, sink, visited)
                    if cmp_path != None:
                        # path lead to sink
                        cmp_path.insert(0, u)
                        if ss_path == None:
                            ss_path = cmp_path
                        else:
                            #FIXME fixed
                            if path_len(graph, cmp_path, overlap) > path_len(graph, ss_path, overlap) and path_cov(graph, cmp_path) >= path_cov(graph, ss_path):
                                ss_path = ss_path
                            else:
                                ss_path = cmp_path
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
        visited = {}
        for u in graph.vertices():
            visited[u] = False
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

def pairwise_contig_dist(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, overlap, max_len):
    """
    find pair-wise shortest path among contigs and construct the pairwise contig path matrix
    TODO if sp found, do we also reduce its coverage?
    """
    dist_matrix = {}
    for cno_i, [contig_i, clen_i, ccov_i] in temp_contigs_dict.items():
        for cno_j, [contig_j, clen_j, ccov_j] in temp_contigs_dict.items():
            print("------------------------------------------------------")
            print("Tail Contig: ", cno_i, " -> Head Contig: ", cno_j)
            if cno_i != cno_j:
                if clen_i + clen_j > max_len:
                    print("Contig ", cno_i, "-", "Contig ", cno_j, " with len over maxlen")
                    continue
                intersect = list(set(contig_i) & set(contig_j))
                if intersect != []:
                    print("Contig ", cno_i, "-", "Contig ", cno_j, " Intersection with ", len(intersect), " nodes")
                    if DEBUG_MODE:
                        print(intersect)
                    continue
            s_path, s_len = distance_search(graph, simp_node_dict, contig_i[-1], contig_j[0], overlap)
            if s_path != None:
                s_path_ids = [graph.vp.id[v] for v in s_path]
                print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])
                s_path_edge_flow = contig_flow(graph, simp_edge_dict,  s_path_ids)
                s_path_ccov = numpy.mean(s_path_edge_flow) if len(s_path_edge_flow) != 0 else 0
                dist_matrix[(cno_i, cno_j)] = (s_path_ids, s_len, s_path_ccov)

    print("------------------------------------------------------")
    return dist_matrix

def contig_merge(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, output_file, min_cov, min_len, max_len, overlap):

    print("-----------------------contig pair-wise concatenaton--------------------")

    if not is_DAG(graph):
        print("graph is cyclic, cyclic strian may exist")
    else:
        print("graph is not cyclic, lienar strain may exist")
    
    dist_matrix = pairwise_contig_dist(graph, simp_node_dict, simp_edge_dict, temp_contigs_dict, overlap, max_len)

    concat_dicision = {}
    for (tail_cno, head_cno), (s_path_ids, s_len, s_path_ccov) in dist_matrix.items():
        [tail_contig, tail_clen, tail_ccov] = temp_contigs_dict[tail_cno]
        [head_contig, head_clen, head_ccov] = temp_contigs_dict[head_cno]
        print("------------------------------------------------------")
        print("Tail Contig: ", tail_cno, " -> Head Contig: ", head_cno)

        if (head_cno, tail_cno) in dist_matrix:
            print("reverse concatentation exist")
        else:
            print("linear concatentation exist")
    
        print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])

        concat_ccov = min(head_ccov, tail_ccov, s_path_ccov) if s_path_ccov != 0.0 else min(head_ccov, tail_ccov)
        concat_len = head_clen + tail_clen - overlap if s_len == 0 else head_clen + s_len + tail_clen - overlap * 2
        print("coverage: head contig eflow: ", head_ccov, " s path eflow: ", s_path_ccov, " tail contig eflow: ", tail_ccov)
        print("potential concatenated length: ", concat_len)

        # decide on cancatenation
        if concat_len > max_len:
            print("exceed maximum strain length")
        else:
            print("length satisfied, concat_ccov: ", concat_ccov)
            concat_dicision[(tail_cno, head_cno)] = concat_ccov
            
            # if concat_ccov >= min_cov:
            #     print("coverage satisfied, concat dicision candidate")
            #     concat_dicision[(tail_cno, head_cno)] = concat_ccov
            # else:
            #     print("lower than minimum coverage")
        print("------------------------------------------------------")
    
    # TODO re-evaluate the concat dicision, sort the concat dicision via concat_ccov, with reverse order or not, TBD
    sorted_pair = sorted(concat_dicision.items(), key=lambda x:x[1], reverse=False)
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
            print("Re-appending contig to temp_contig_dict: ", cno)
            concat_contig_dict[cno] = [contig, clen, ccov]
        else:
            print("contig {0} is used".format(cno))
    
    udpate_node_to_contig_dict(node_to_contig_dict, simp_node_dict)
    update_edge_to_contig_dict(edge_to_contig_dict, simp_edge_dict)
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
    
    # TODO recovered node coverage should equal to the rest of contig sum cov, not full
    # recover nodes
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

def graph_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, output_file, overlap, TEMP_DIR):
    """
    Compactify the reduced De Bruijin graph
    """
    print("--------------------graph compactification--------------------")
    # find simple edges, simple edge is the only edge between its source and sink
    simple_edges = []
    out_edge = {}
    in_edge = {}
    for e in simp_edge_dict.values():
        if e.source().out_degree() == 1 and e.target().in_degree() == 1:
            assert e.source() != e.target()
            simple_edges.append([e.source(), e.target()])
            in_edge[int(e.source())] = e
            out_edge[int(e.target())] = e

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
    
    path_dict = {}
    path_src = {}
    path_sink = {}
    for p in simple_paths:
        print("path: ", [int(graph.vp.id[u]) for u in p])
        pid = str(graph.vp.id[p[0]]) + "00" + str(len(p)) + "00" + str(graph.vp.id[p[-1]])
        plen = path_len(graph, p, overlap)
        pseq = path_to_seq(graph, p, pid, overlap)
        pdp = numpy.mean([graph.vp.dp[u] for u in p])
        pkc = numpy.mean([graph.vp.kc[u] for u in p])

        # recolor the used node to gray
        for i in range(len(p)):
            u = p[i]
            v = p[i+1] if (i + 1) < len(p) else None
            graph.vp.color[u] = 'gray'
            if v != None:
                graph.vp.color[v] = 'gray'
            if u != None and v != None:
                graph.ep.color[graph.edge(u,v)] = 'gray'
        
        # add the path to graph and dict
        pv = graph.add_vertex()
        graph.vp.seq[pv] = pseq
        graph.vp.dp[pv] = pdp
        graph.vp.kc[pv] = pkc
        graph.vp.id[pv] = pid
        graph.vp.color[pv] = 'black'
        
        uid = graph.vp.id[p[0]]
        vid = graph.vp.id[p[-1]]
        path_dict[(pid, uid, vid)] = pv
        path_src[uid] = pv
        path_sink[vid] = pv
        simp_node_dict[pid] = pv
    
    # re-append path edges
    for (pid, uid, vid), pv in path_dict.items():
        u_in_neighbors = simp_node_dict[uid].in_neighbors()
        for u_neighbor in u_in_neighbors:
            un_id = graph.vp.id[u_neighbor]
            if graph.vp.color[u_neighbor] == 'black':
                e = graph.add_edge(u_neighbor, pv)
                graph.ep.overlap[e] = overlap
                graph.ep.color[e] = 'black'
                simp_edge_dict[(un_id, pid)] = e
            else:
                if un_id in path_sink:
                    u_alternative = path_sink[un_id]
                    e = graph.add_edge(u_alternative, pv)
                    graph.ep.overlap[e] = overlap
                    graph.ep.color[e] = 'black'
                    simp_edge_dict[(graph.vp.id[u_alternative], pid)] = e
                
        v_out_neighbors = simp_node_dict[vid].out_neighbors()
        for v_neighbor in v_out_neighbors:
            vn_id = graph.vp.id[v_neighbor]
            if graph.vp.color[v_neighbor] == 'black':
                e = graph.add_edge(pv, v_neighbor)
                graph.ep.overlap[e] = overlap
                graph.ep.color[e] = 'black'
                simp_edge_dict[(pid, vn_id)] = e
            else:
                if vn_id in path_src:
                    v_alternative = path_src[vn_id]
                    e = graph.add_edge(pv, v_alternative)
                    graph.ep.overlap[e] = overlap
                    graph.ep.color[e] = 'black'
                    simp_edge_dict[(pid, graph.vp.id[v_alternative])] = e
    

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, output_file)

    # map contig to graph
    gfa_to_fasta(output_file, "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    contig_dict_to_fasta(graph, contig_dict, simp_node_dict, overlap, "{0}temp_contigs.fasta".format(TEMP_DIR))
    minimap_api("{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR), "{0}temp_contigs.fasta".format(TEMP_DIR), "{0}temp_contigs_to_graph.paf".format(TEMP_DIR))
    
    contig_map_node_dict = {}
    with open("{0}temp_contigs_to_graph.paf".format(TEMP_DIR), 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            cno = splited[0]
            ccov = cno.split("_")[-1]
            # c_l = int(splited[1])
            # c_s = int(splited[2])
            # c_f = int(splited[3])
            node_no = splited[5]
            # if ((c_f - c_s) / c_l) >= 0.6:
            if node_no not in contig_map_node_dict:
                contig_map_node_dict[node_no] = []
            contig_map_node_dict[node_no].append((cno, ccov))
        paf.close()
    subprocess.check_call("rm {0}temp*".format(TEMP_DIR), shell=True)

    for node_no, cnos in contig_map_node_dict.items():
        print(node_no, cnos)
    print([int(no) for no in contig_map_node_dict.keys()])
    print("------------------graph compactification end------------------")
    return contig_map_node_dict

def graph_grouping(graph: Graph, simp_node_dict: dict, forward="", reverse="", partition_length_cut_off=0):
    """
    Maximimize graph connectivity by minimizing node with 0 in-degree or out-degree, detect and remove all the cycles.
    Out-of-date, TBD
    """
    # determine the isolated subgraphs, and assign each node with its group No, which indicates they are belong to same group
    def bfs_grouping(graph: Graph, start_node, group_no, groups):
        """
        Perform a breadth-first search and assign the group no to all the connected nodes.
        """
        groups[group_no] = []

        queue = []
        queue.append(start_node)
        while queue:
            v = queue.pop()
            graph.vp.group[v] = group_no
            groups[group_no].append(v)

            for neighbor in v.all_neighbors():
                if graph.vp.group[neighbor] == -1:
                    queue.append(neighbor)
        return graph, group_no, groups

    graph.vp.group = graph.new_vertex_property("int16_t", val=-1)
    group_no = 1
    groups = {}
    for v in simp_node_dict.values():
        # grouping
        if graph.vp.group[v] == -1:
            graph, group_no, groups = bfs_grouping(graph, v, group_no, groups)
            group_no = group_no + 1

    print("number of groups:", len(groups))
    for key, item in groups.items():
        print("group number: ", key, " member: ", [graph.vp.id[v] for v in item])

    # connect sub-graphs based on pair-end reads information
    # TODO
    return graph, groups


def path_extraction(graph_comp: Graph, simp_node_dict_comp: dict, simp_edge_dict_comp: dict, overlap, min_cov, min_len):
    """
    extract the last mile path from the graph, with support from residue contigs
    """
    print("--------------------path extraction--------------------")

    paths_per_group_dict = {}

    graph_comp, groups = graph_grouping(graph_comp, simp_node_dict_comp)
    
    final_strain_dict = {}
    
    for gno, group in groups.items():
        if len(group) < 2:
            # single node group
            if len(group) == 1:
                node = group[0]
                pcov = graph_comp.vp.dp[node]
                id = graph_comp.vp.id[node]
                plen = len(graph_comp.vp.seq[node])
                print("Single node Path: ", int(id), "path len: ", plen, "cov: ", pcov)
                if plen < min_len / 5 or pcov < min_cov:
                    print("reject")
                else:
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
            no = graph_comp.vp.id[u]
            if u.in_degree() == 0 and u.out_degree() == 0:
                isolations.append(no)
            elif u.in_degree() == 0:
                srcs.append(no)
            elif u.out_degree() == 0:
                sinks.append(no)
            else:
                middles.append(no)

        for src in srcs:
            for sink in sinks:
                p, plen = distance_search(graph_comp, simp_node_dict_comp, src, sink, overlap)
                if p != None:
                    s_path_ids = [graph_comp.vp.id[v] for v in p]
                    s_path_ids.append(sink)
                    s_path_ids.insert(0, src)             
                    plen = plen + len(graph_comp.vp.seq[simp_node_dict_comp[src]]) + len(graph_comp.vp.seq[simp_node_dict_comp[sink]]) - 2 * overlap
                    paths.append((s_path_ids, plen))
        paths_per_group_dict[gno] = paths

    for gno, paths in paths_per_group_dict.items():
        print("Current group: ", gno)
        for p_ids, plen in paths:
            print("------------------------------------------------------")
            print("Path: ", [int(u) for u in p_ids])
            print("path len: \n", plen)
            # involved_contig = [(id, contig_map_node_dict[id]) for id in p_ids if id in contig_map_node_dict]
            # for id, cs in involved_contig:
            #     print(id, cs)
            # TODO if the path involved some used contig and the contig has been exhausted many time
            # within its limited ccov, then skip this path.
            # also reduce the cov via the path cov on the nodes

            pcov = numpy.min([graph_comp.vp.dp[simp_node_dict_comp[u]] for u in p_ids])
            pno = str(p_ids[0]) + "_" + str(p_ids[-1])
            if plen >= min_len / 5 and pcov >= min_cov:
                final_strain_dict[pno] = (p_ids, plen, pcov)
    print("------------------------------------------------------")
    print("--------------------path extraction end--------------------")
    return final_strain_dict

if __name__ == "__main__":
    main()
