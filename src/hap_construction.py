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
from matplotlib.cbook import maxdict
import numpy

from graph_converter import *

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = True
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

    # Read in as Level -1 graph
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, init_ori=1)
    contig_dict, node_to_contig_dict, edge_to_contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_cov, args.min_len, args.overlap)
    
    contig_node_cov_rise(graph, contig_dict, node_to_contig_dict)

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}init_graph.gfa".format(TEMP_DIR))
    graph_init, simp_node_dict_init, simp_edge_dict_init = flipped_gfa_to_graph("{0}init_graph.gfa".format(TEMP_DIR))
    
    prev_dp_dict = coverage_rebalance_formal(graph_init, simp_node_dict_init, simp_edge_dict_init)

    graph_to_gfa(graph_init, simp_node_dict_init, simp_edge_dict_init, "{0}pre_graph.gfa".format(TEMP_DIR))
    
    # Read in as pre graph, only used for path seq extraction.
    pre_graph, simp_node_dict_pre, simp_edge_dict_pre = flipped_gfa_to_graph("{0}pre_graph.gfa".format(TEMP_DIR))
    assign_edge_flow(pre_graph, simp_node_dict_pre, simp_edge_dict_pre)

    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict_pre, "{0}pre_graph.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    # selected contig from SPAdes
    contig_dict_to_fasta(graph_init, contig_dict, simp_node_dict_init, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))

    for cno, (contig, clen, ccov) in list(contig_dict.items()):
        print("---------------------------------------------------------------")
        print_contig(cno, clen, ccov, contig)
        if len(contig) < 2:
            contig_edge_flow = [pre_graph.vp.dp[simp_node_dict_pre[contig[0]]]]
        else:
            contig_edge_flow = contig_flow(pre_graph, simp_edge_dict_pre, contig)
        print("edge flow - mean: {0}, min: {1}, max: {2}, median: {3}".format(numpy.mean(contig_edge_flow), numpy.min(contig_edge_flow), numpy.max(contig_edge_flow), numpy.median(contig_edge_flow)))
        
        contig_depths = [pre_graph.vp.dp[simp_node_dict_pre[u]] for u in contig]
        print("depth - mean: {0}, min: {1}, max: {2}, median: {3}".format(numpy.mean(contig_depths), numpy.min(contig_depths), numpy.max(contig_depths), numpy.median(contig_depths)))
        
        contig_kc = [pre_graph.vp.kc[simp_node_dict_pre[u]] for u in contig]
        print("kc - mean: {0}, min: {1}, max: {2}, median: {3}".format(numpy.mean(contig_kc), numpy.min(contig_kc), numpy.max(contig_kc), numpy.median(contig_kc)))
        
        prev_dp_sum = numpy.sum([prev_dp_dict[n] for n in contig])
        curr_dp_sum = numpy.sum([pre_graph.vp.dp[simp_node_dict_pre[n]] for n in contig])
        ratio = curr_dp_sum / prev_dp_sum
        print("Prev dp sum: ", prev_dp_sum)
        print("Curr dp sum: ", curr_dp_sum)
        print("Ratio: ", curr_dp_sum / prev_dp_sum)
        print("Prev cov: {0} vs median flow: {1}".format(ccov, numpy.median(contig_edge_flow)))
        print("normalised cov: {0} vs normalised median flow (new ccov): {1}".format(ccov/ratio, numpy.median([d/ratio for d in contig_edge_flow])))
        contig_dict[cno][2] = numpy.median([d/ratio for d in contig_edge_flow])

    # the max depth node on the connected graph must be visited by all the s-t path
    # proof:
    # if there exist a s-t path that is not involved by max dp node n
    # then there exist a node n' other than n such that dp(n') = dp(n) + path_flow(p)
    # in which n is not the max dp node, contradiction.
    #
    # the min depth node on the connected graph must be visited by exactly one s-t path 
    # store all the overlap/intersection contig case, key can never concat with any contig from the values
    
    contig_concat_plans = get_concat_plan(contig_dict, args.max_len)
    
    contig_info = []
    contig_node_dict = {}
    # reduce all the contig to a single node from the graph
    for cno, (contig, clen, ccov) in list(contig_dict.items()):
        src = contig[0]
        tgt = contig[-1]
        id = src + "00" + cno + "00" + tgt
        cseq = path_to_seq(graph_init, [simp_node_dict_init[n] for n in contig], cno, args.overlap)
        kc = numpy.median([graph_init.vp.kc[simp_node_dict_init[u]] for u in contig])
        in_edges = list((graph_init.vp.id[e.source()], src) for e in simp_node_dict_init[src].in_edges())
        out_edges = list((tgt, graph_init.vp.id[e.target()],) for e in simp_node_dict_init[tgt].out_edges())
        
        contig_reduction(pre_graph, contig, cno, clen, ccov, simp_node_dict_pre, simp_edge_dict_pre, 0)

        cv = pre_graph.add_vertex()
        pre_graph.vp.seq[cv] = cseq
        pre_graph.vp.dp[cv] = ccov
        pre_graph.vp.kc[cv] = kc
        pre_graph.vp.id[cv] = id
        pre_graph.vp.color[cv] = 'black'
        simp_node_dict_pre[id] = cv

        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])
        contig_node_dict[cno] = cv
    
    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u,_) in in_edges:
            if u in simp_node_dict_pre and (u, pre_graph.vp.id[node]) not in simp_edge_dict_pre:
                ue = pre_graph.add_edge(simp_node_dict_pre[u], node)
                pre_graph.ep.overlap[ue] = args.overlap
                pre_graph.ep.color[ue] = 'black'
                simp_edge_dict_pre[(u, pre_graph.vp.id[node])] = ue
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (pre_graph.vp.id[in_node], pre_graph.vp.id[node]) not in simp_edge_dict_pre:
                    ue = pre_graph.add_edge(in_node, node)
                    pre_graph.ep.overlap[ue] = args.overlap
                    pre_graph.ep.color[ue] = 'black'
                    simp_edge_dict_pre[(pre_graph.vp.id[in_node], pre_graph.vp.id[node])] = ue                 

        for (_,v) in out_edges:
            if v in simp_node_dict_pre and (pre_graph.vp.id[node], v) not in simp_edge_dict_pre:
                ve = pre_graph.add_edge(node, simp_node_dict_pre[v])
                pre_graph.ep.overlap[ve] = args.overlap
                pre_graph.ep.color[ve] = 'black'
                simp_edge_dict_pre[(pre_graph.vp.id[node], v)] = ve
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (pre_graph.vp.id[node], pre_graph.vp.id[out_node]) not in simp_edge_dict_pre:
                    ve = pre_graph.add_edge(node, out_node)
                    pre_graph.ep.overlap[ve] = args.overlap
                    pre_graph.ep.color[ve] = 'black'
                    simp_edge_dict_pre[(pre_graph.vp.id[node], pre_graph.vp.id[out_node])] = ve 
    
    graph_simplification(pre_graph, simp_node_dict_pre, simp_edge_dict_pre, {}, {}, args.min_cov/2)
    graph_to_gfa(pre_graph, simp_node_dict_pre, simp_edge_dict_pre, "{0}red_graph.gfa".format(TEMP_DIR))
    
    pre_graph_v2, simp_node_dict_pre_v2, simp_edge_dict_pre_v2 = flipped_gfa_to_graph("{0}red_graph.gfa".format(TEMP_DIR))
    coverage_rebalance_formal(pre_graph_v2, simp_node_dict_pre_v2, simp_edge_dict_pre_v2)
    
    graph_to_gfa(pre_graph_v2, simp_node_dict_pre_v2, simp_edge_dict_pre_v2, "{0}graph_L0.gfa".format(TEMP_DIR))

    level_no = 0
    graph_iter, simp_node_dict_iter, simp_edge_dict_iter = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no)))
    # pairwise_contig_concatenation(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, contig_node_dict, contig_concat_plans)
    # strain_dict = {}
    # iter_condition = set()
    # concat_contig_dict_iter = contig_dict.copy()
    # concat_strain_dict_iter = {}
    # contig pair-wise concatenation iteration
    # with iteration stopped when no further concatenation occurred
    # while True:
    #     print("Current iteration: ", level_no)
    #     graph_iter, simp_node_dict_iter, simp_edge_dict_iter = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no)))
    #     graph_simplification(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, {}, {}, args.min_cov/2)
    #     assign_edge_flow(graph_iter, simp_node_dict_iter, simp_edge_dict_iter)
    #     output_file = "{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no+1))
    #     # contig pair-wise concatenation
    #     concat_strain_dict_iter, concat_contig_dict_iter = contig_merge(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, concat_contig_dict_iter, node_to_contig_dict, edge_to_contig_dict, output_file, args.min_cov, args.min_len, args.max_len, args.overlap)

    #     contig_dict_to_fasta(graph_iter, concat_contig_dict_iter, simp_node_dict_iter, args.overlap, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(level_no)))
    #     minimap_api(args.ref_file, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(level_no)), "{0}L{1}_contigs_to_strain.paf".format(TEMP_DIR, str(level_no)))

    #     strain_dict.update(concat_strain_dict_iter.copy())
    #     level_no += 1

    #     iter_compare = set(concat_contig_dict_iter.keys())
    #     if iter_condition == iter_compare:
    #         print("end of contig pair-wise concatenation iteration: ", level_no)
    #         break
    #     else:
    #         iter_condition = iter_compare
    
    # # further step to recover the rest of strain from the reduced graph

    # graph_red, simp_node_dict_red, simp_edge_dict_red = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no)))
    # # graph_simplification(graph_red, simp_node_dict_red, simp_edge_dict_red, node_to_contig_dict, edge_to_contig_dict, args.min_cov)
    # assign_edge_flow(graph_red, simp_node_dict_red, simp_edge_dict_red)
    
    # # extract the last mile paths from the graph
    # final_strain_dict = path_extraction(graph_red, simp_node_dict_red, simp_edge_dict_red, node_usage_dict, args.overlap, args.min_cov, args.min_len)   
    
    # graph_to_gfa(graph_red, simp_node_dict_red, simp_edge_dict_red, "{0}final_graph.gfa".format(TEMP_DIR))

    # strain_dict.update(final_strain_dict)

    # # print strain and minimap overlaps
    # contig_dict_to_fasta(pre_graph_v2, strain_dict, simp_node_dict_pre_v2, args.overlap, "{0}cand_strains.fasta".format(TEMP_DIR))
    # contig_dict_to_path(strain_dict, "{0}cand_strains.paths".format(TEMP_DIR))
    # minimap_api(args.ref_file, "{0}cand_strains.fasta".format(TEMP_DIR), "{0}cand_strain_to_strain.paf".format(TEMP_DIR))
    
    # # print out node usage stat
    # node_usage_pair = sorted(node_usage_dict.items(), key=lambda x: x[1][1] - x[1][0], reverse=True)
    # for id, used in node_usage_pair:
    #     print("id: {0} has been used {1} times".format(id, used))

    # TODO extend the strain end length. really high chance
    # map all the final strains back into pre graphs

def get_concat_plan(contig_dict: dict, max_len):
    contig_impossible_dict = {}
    all_contig_ids = contig_dict.keys()
    for no in contig_dict.keys():
        contig_impossible_dict[no] = set()
    for tail_cno, [tail_contig, tail_clen, _] in contig_dict.items():
        for head_cno, [head_contig, head_clen, _] in contig_dict.items():
            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len:
                    contig_impossible_dict[tail_cno].add(head_cno)
                if list(set(tail_contig) & set(head_contig)) != []:
                    contig_impossible_dict[tail_cno].add(head_cno)
    contig_concat_plans = {}
    for key, item in contig_impossible_dict.items():
        ps = all_contig_ids - item
        print("cno: ", key, " can concat with following: ", ps)
        contig_concat_plans[key] = ps
    return contig_concat_plans
# def pairwise_contig_concatenation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_node_dict: dict, concat_plan: dict):
#     for 
#     return 

def contig_node_cov_rise(graph: Graph, contig_dict: dict, node_to_contig_dict: dict):
    """
    for any node that involved in one or more contigs, rise the depth if less than the sum of related contig cov
    """
    for no, [cnos, dp, node] in list(node_to_contig_dict.items()):
        sum_covs = numpy.sum([contig_dict[cno][2] for cno in cnos])
        if sum_covs > graph.vp.dp[node]:
            print("Node: {0} dp is really low: {1} vs {2}, rise it up".format(no, graph.vp.dp[node], sum_covs))
            graph.vp.dp[node] = sum_covs
            node_to_contig_dict[no][2] = sum_covs
    return

def coverage_rebalance_formal(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
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
    def maximization_node_depth(graph: Graph, simp_node_dict: dict, usage_dict: dict):
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
                usage_dict[no] += 1
                is_update = True
        return is_update
    
    cutoff = 0.001 * len(simp_node_dict)
    sum_delta = 0
    sum_depth_before = numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()])
    
    # store previous node depth
    prev_dp_dict = {}
    for no, v in simp_node_dict.items():
        prev_dp_dict[no] = graph.vp.dp[v]

    print(cutoff)
    is_update = True
    usage_dict = {}
    for no in simp_node_dict.keys():
        usage_dict[no] = 0
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
        is_update = maximization_node_depth(graph, simp_node_dict, usage_dict)

    # final evaluation
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Ratio: ", sum_ratio, "Delta: ", sum_delta)
    for node in simp_node_dict.values():
        graph.vp.dp[node] = graph.vp.dp[node] / sum_ratio
    return prev_dp_dict

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
        graph.ep.flow[e] = max(graph.ep.flow[e], flow)
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

def simp_path(graph: Graph, simp_edge_dict: dict):
    """
    find simple edges, simple edge is the only edge between its source and sink
    """
    simple_edges = []
    out_edge = {}
    in_edge = {}
    for e in simp_edge_dict.values():
        src = e.source()
        src_out_d = len([u for u in src.out_neighbors() if graph.vp.color[u] == 'black'])
        target = e.target()
        target_in_d = len([u for u in target.in_neighbors() if graph.vp.color[u] == 'black'])
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

def graph_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, output_file, overlap, TEMP_DIR):
    """
    Compactify the reduced De Bruijin graph
    """
    print("--------------------graph compactification--------------------")
    simple_paths = simp_path(graph, simp_edge_dict)
    
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
            node_no = splited[5]
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