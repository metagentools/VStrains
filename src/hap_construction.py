#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
import subprocess
# from graph_tool import GraphView, _in_degree
# import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw
# from graph_tool.clustering import local_clustering

import argparse

import gfapy
import numpy

from graph_converter import *

usage = "Construct haplotypes using divide-and-conquer method"
author = "Runpeng Luo"

debug = False

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, default=100, help=("minimum coverage for strains"))
    parser.add_argument('-minlen', '--minimum_strain_length', dest='min_len', default=8000, type=int, help=("minimum strain length"))
    parser.add_argument('-maxlen', '--maximum_strain_length', dest='max_len', default=10000, type=int, help=("maximum strain length"))
    parser.add_argument('-overlap', '--vertex_overlap', dest='overlap', default=127, type=int, help=("adjacent vertex overlap in the graph"))
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, debug only')
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')

    ## TODO may add gfa validation
    args = parser.parse_args()
    if not args.gfa_file:
        print("gfa file is not imported")
        return 1
    
    subprocess.check_call("rm -rf acc/ && mkdir acc/", shell=True)

    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict)
    graph_simplification(graph, simp_node_dict, simp_edge_dict, args.min_cov)

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, args.min_cov, "acc/graph_L1.gfa")

    # graph_draw(graph, vprops={'text': graph.vp.id}, eprops={'text': graph.ep.flow}, output="graph.pdf", output_size=(2000,2000))
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict, "acc/graph_L1.gfa")
        
    contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_cov, args.min_len, args.overlap)

    for cno, (contig, clen, ccov) in contig_dict.items():
        print_contig(cno, clen, ccov, contig)

    cand_strains_dict, temp_contigs_dict = graph_reduction(graph, contig_dict, simp_node_dict, simp_edge_dict, "acc/graph_L2.gfa", args.min_cov, args.min_len)

    graph_L2, simp_node_dict_L2, simp_edge_dict_L2 = gfa_to_graph("acc/graph_L2.gfa", init_ori=1)
    assign_edge_flow(graph_L2, simp_node_dict_L2, simp_edge_dict_L2)
    graph_simplification(graph_L2, simp_node_dict_L2, simp_edge_dict_L2, args.min_cov)
    concat_strain_dict, concat_contig_dict = contig_classification(graph_L2, simp_node_dict_L2, simp_edge_dict_L2, temp_contigs_dict, "acc/graph_L3.gfa", args.min_cov, args.min_len, args.max_len, args.overlap)
    
    # graph_L3, simp_node_dict_L3, simp_edge_dict_L3 = gfa_to_graph("acc/graph_L3.gfa", init_ori=1)
    # assign_edge_flow(graph_L3, simp_node_dict_L3, simp_edge_dict_L3)
    # graph_simplification(graph_L3, simp_node_dict_L3, simp_edge_dict_L3, args.min_cov)
    # concat_strain_dict_2, concat_contig_dict_2 = contig_classification(graph_L3, simp_node_dict_L3, simp_edge_dict_L3, concat_contig_dict, "acc/graph_L4.gfa", args.min_cov, args.min_len, args.max_len, args.overlap)

    graph_L1, simp_node_dict_L1, simp_edge_dict_L1 = gfa_to_graph("acc/graph_L1.gfa", init_ori=1)
    assign_edge_flow(graph_L1, simp_node_dict_L1, simp_edge_dict_L1)
    graph_simplification(graph_L2, simp_node_dict_L2, simp_edge_dict_L2, args.min_cov)

    strain_dict = concat_strain_dict.copy()
    # strain_dict.update(concat_strain_dict_2)
    strain_dict.update(cand_strains_dict)

    contig_dict_to_fq(graph_L1, strain_dict, simp_node_dict_L1, args.overlap, "acc/cand_strains.fq")
    minimap_api(args.ref_file, "acc/cand_strains.fq", "acc/ref_map_cand_strain.paf")

def graph_stat(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    print("-------------------------graph stat----------------------")
    for seg_no, v in simp_node_dict.items():
        print_vertex(graph, v, "stat")
    for (_,_), e in simp_edge_dict.items():
        print_edge(graph, e, "stat")
    
    print("-----------------------graph stat end--------------------")

def graph_simplification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, min_cov):
    """
    Directly remove all the vertex with coverage less than minimum coverage and related edge
    """
    print("-------------------------graph simplification----------------------")
    print("Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    for id, node in list(simp_node_dict.items()):
        if graph.vp.dp[node] < min_cov:
            if debug:
                print_vertex(graph, node, "Node removed by graph simplification -")
            # delete the node
            simp_node_dict.pop(id)

            # delete related edges
            for out_node in node.out_neighbors():
                out_id = graph.vp.id[out_node]
                if (id, out_id) in simp_edge_dict:
                    simp_edge_dict.pop((id, out_id))
            for in_node in node.in_neighbors():
                in_id = graph.vp.id[in_node]
                if (in_id, id) in simp_edge_dict:
                    simp_edge_dict.pop((in_id, id))
    print("Remain: Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    print("-------------------------graph simplification end----------------------")
    return

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

    # converage iteration
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

    print("un-assigned edges after node-weight converage iteration : ", un_assigned_edge)
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
                if debug:
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
    if debug:
        print("*---Contig: ", cno, clen, ccov)
    adj_index = 1
    for node in contig:
        if adj_index >= len(contig):
            break
        adj_node = contig[adj_index]

        u = simp_node_dict[node] if node in simp_node_dict else None
        v = simp_node_dict[adj_node] if adj_node in simp_node_dict else None
        e = simp_edge_dict[(node,adj_node)] if (node,adj_node) in simp_edge_dict else None 

        # edge may be eliminated from previous execution already
        if u == None or v == None or e == None:
            if e != None:
                graph.ep.flow[e] = 0
                graph.ep.color[e] = 'gray'
                simp_edge_dict.pop((node,adj_node))
                if debug:
                    print_edge(graph, e, "edge been removed due to either u or v is removed")
            else:
                if debug:
                    print("edge: ", node, " -> ", adj_node, "already been removed")
            continue
        if debug:
            print_edge(graph, e, "current edge eval")
        # reduce the depth for involving edge
        if graph.ep.flow[e] - ccov <= min_cov:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            if debug:
                print_edge(graph, e, "edge been removed")
        else:
            graph.ep.flow[e] = graph.ep.flow[e] - ccov

        # reduce the depth for involving node, gray color for forbiddened node
        if graph.vp.dp[u] - ccov <= min_cov:
            graph.vp.dp[u] = 0
            graph.vp.color[u] = 'gray'
            simp_node_dict.pop(node)
            if debug:
                print("node ", node, "been removed")
        else:
            graph.vp.dp[u] = graph.vp.dp[u] - ccov

        # last node in the contig, reduce its depth
        if adj_index + 1 == len(contig):
            if graph.vp.dp[v] - ccov <= min_cov: 
                graph.vp.dp[v] = 0
                graph.vp.color[v] = 'gray'
                simp_node_dict.pop(adj_node)
                if debug:
                    print("node ", adj_node, "been removed")
            else:
                graph.vp.dp[v] = graph.vp.dp[v] - ccov
        
        # update edges
        if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (node,adj_node) in simp_edge_dict:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            if debug:
                print_edge(graph, e, "edge been removed in the final step")
        adj_index = adj_index + 1
    return

def graph_reduction(graph: Graph, contig_dict: dict, simp_node_dict: dict, simp_edge_dict: dict, output_file, min_cov, min_len):
    """
    reduce the node/edge weight based on existing contig found by SPAdes.
    only contig with minimum strain length satisfied be removed
    return non-removed contig dict
    """

    print("-------------------------graph reduction----------------------")
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
    # store level 2 graph
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, min_cov, output_file)

    print("-----------------------graph reduction end--------------------")
    return cand_strains_dict, temp_contigs_dict

def distance_search(graph: Graph, simp_node_dict: dict, source, source_contig, sink, sink_contig, overlap: int):
    """
    Compute minimal distance and its path between source node and sink node
    optimise the function with contig overlap check TODO
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
            rtn = None
            for v in u.out_neighbors():
                if not visited[v]:
                    visited[v] = True
                    cmp = dfs_helper(graph, v, sink, visited)
                    if cmp != None:
                        # path lead to sink
                        cmp.insert(0, u)
                        if rtn == None:
                            rtn = cmp
                        else:
                            rtn = rtn if path_len(graph, cmp, overlap) > path_len(graph, rtn, overlap) else cmp
                    visited[v] = False
        return rtn
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

def pairwise_contig_dist(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, overlap):
    """
    find pair-wise shortest path among contigs and construct the pairwise contig path matrix
    """
    dist_matrix = {}
    for cno_i, [contig_i, clen_i, ccov_i] in temp_contigs_dict.items():
        for cno_j, [contig_j, clen_i, ccov_i] in temp_contigs_dict.items():
            print("------------------------------------------------------")
            if cno_i == cno_j:
                continue
            print("Tail Contig: ", cno_i, " -> Head Contig: ", cno_j)
            intersect = list(set(contig_i) & set(contig_j))
            if intersect != []:
                print("Contig ", cno_i, "-", "Contig ", cno_j, " Intersection with ", len(intersect), " nodes")
                if debug:
                    print(intersect)
                continue
            s_path, s_len = distance_search(graph, simp_node_dict, contig_i[-1], contig_i, contig_j[0], contig_j, overlap)
            if s_path != None:
                s_path_ids = [graph.vp.id[v] for v in s_path]
                print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])
                s_path_edge_flow = contig_flow(graph, simp_edge_dict,  s_path_ids)
                s_path_ccov = numpy.mean(s_path_edge_flow) if len(s_path_edge_flow) != 0 else 0
                dist_matrix[(cno_i, cno_j)] = (s_path_ids, s_len, s_path_ccov)

    print("------------------------------------------------------")
    return dist_matrix

def contig_classification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, output_file, min_cov, min_len, max_len, overlap):

    print("-----------------------contig classifcation--------------------")

    if not is_DAG(graph):
        print("graph is cyclic, cyclic strian may exist")
    else:
        print("graph is not cyclic, lienar strain may exist")
    
    dist_matrix = pairwise_contig_dist(graph, simp_node_dict, simp_edge_dict, temp_contigs_dict, overlap)

    concat_dicision = []
    for (tail_cno, head_cno), (s_path_ids, s_len, s_path_ccov) in dist_matrix.items():
        [tail_contig, tail_clen, tail_ccov] = temp_contigs_dict[tail_cno]
        [head_contig, head_clen, head_ccov] = temp_contigs_dict[head_cno]
        print("------------------------------------------------------")
        print("Tail Contig: ", tail_cno, " -> Head Contig: ", head_cno)

        if ((head_cno, head_clen, head_ccov), (tail_cno, tail_clen, tail_ccov)) in dist_matrix:
            print("reverse concatentation exist")
        else:
            print("linear concatentation exist")
    
        print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])

        min_mu = min(head_ccov, tail_ccov, s_path_ccov) if s_path_ccov != 0.0 else min(head_ccov, tail_ccov)
        concat_len = head_clen + tail_clen - overlap if s_len == 0 else head_clen + s_len + tail_clen - overlap * 2
        print("coverage: head contig eflow: ", head_ccov, " s path eflow: ", s_path_ccov, " tail contig eflow: ", tail_ccov)
        print("potential concatenated length: ", concat_len)

        # decide on cancatenation
        if concat_len > max_len:
            print("exceed maximum strain length")
        else:
            print("length satisfied")
            if min_mu >= min_cov:
                print("coverage satisfied")
                concat_dicision.append((tail_cno, head_cno))
        print("------------------------------------------------------")
    
    #start concatenation
    concat_strain_dict = {}
    concat_contig_dict = {}
    skip_key = set()
    used_contig = set()
    for (tail_cno, head_cno) in concat_dicision:
        [tail_contig, tail_clen, tail_ccov] = temp_contigs_dict[tail_cno]
        [head_contig, head_clen, head_ccov] = temp_contigs_dict[head_cno]

        if (head_cno, tail_cno) in skip_key:
            continue
        if (tail_cno, head_cno) in skip_key:
            continue
        skip_key.add((head_cno, tail_cno))
        skip_key.add((tail_cno, head_cno))

        is_linear = False
        if (head_cno, tail_cno) not in concat_dicision:
            print("no reverse concatenation exists, potential lienar strain")
            is_linear = True
        else: 
            print("access contig in both direction, potential cyclic strain")
            is_linear = False
        
        (s_path_ids_l, s_len_l, s_path_ccov_l) = dist_matrix[(tail_cno, head_cno)]

        print("Concatentate contigs: ", tail_cno, " <-> ", head_cno)
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
        
    # reduce all the full length concated contigs
    for cno, (contig, clen, ccov) in concat_strain_dict.items():
        print("------------------------------------------------------")
        print_contig(cno, clen, ccov, contig, "Cand concat strain")
        contig_reduction(graph, contig, cno, clen, ccov, simp_node_dict, simp_edge_dict, min_cov)
        print("------------------------------------------------------")
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, min_cov, output_file)

    for cno, [contig, clen, ccov] in temp_contigs_dict.items():
        if cno not in used_contig and ccov >= min_cov:
            concat_contig_dict[cno] = [contig, clen, ccov]
        else:
            print("contig {0} is used or coverage is lower than threshold".format(cno))
    
    for cno, [contig, clen, ccov] in list(concat_contig_dict.items()):
        concat_contig_dict.pop(cno)
        contig_list = contig_split(graph, cno, contig, simp_node_dict, simp_edge_dict, overlap) #FIXME
        if contig_list == []:
            print("No sub contig be found for original contig: ", cno)
        else:
            print("Update sub contigs to the concat contig dict for cno: ", cno)
            for (sub_cno, sub_contig, sub_clen, sub_ccov) in contig_list:
                if sub_cno in concat_contig_dict:
                    print("sub cno: ", sub_cno, " already exist, error")
                else:
                    concat_contig_dict[sub_cno] = [sub_contig, sub_clen, sub_ccov]

    for cno, [contig, clen, ccov] in concat_contig_dict.items():
        print("------------------------------------------------------")
        print_contig(cno, clen, ccov, contig, "partial length residue contig found after concatenation")
        print("------------------------------------------------------")

    print("--------------------contig classification end--------------------")
    return concat_strain_dict, concat_contig_dict

if __name__ == "__main__":
    main()
