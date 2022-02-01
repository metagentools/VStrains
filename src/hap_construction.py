#!/usr/bin/env python3

# import re
# import sys, os
# import json
import re
from graph_tool import GraphView, _in_degree
import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw
# from graph_tool.clustering import local_clustering
from graph_tool import flow
# import subprocess
import argparse

import gfapy
import numpy

from graph_converter import gfa_to_graph, graph_to_gfa, map_ref_to_graph, get_contig, contig_to_seq

usage = "Construct haplotypes using divide-and-conquer method"
author = "Runpeng Luo"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, help=("minimum coverage for strains"))
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, debug only')

    ## TODO may add gfa validation
    args = parser.parse_args()
    if not args.gfa_file:
        print("gfa file is not imported")
        return 1

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa2').from_file(filename=args.gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))
    
    graph, node_dict, edge_dict, dp_dict = gfa_to_graph(gfa)
    graph, simp_node_dict, edge_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict.copy(), 1)
    # graph, simp_node_dict = simplify_graph(graph, node_dict, pick_dict)

    assign_edge_flow(graph, simp_node_dict, edge_dict)

    contig_dict = get_contig(graph, args.contig_file, simp_node_dict, edge_dict, args.min_cov)

    graph_stat(graph, simp_node_dict, edge_dict)

    # graph_draw(graph, vprops={'text': graph.vp.id}, eprops={'text': graph.ep.flow}, output="graph.pdf", output_size=(2000,2000))
    
    # strain_dict = map_ref_to_graph(args.ref_file, graph, simp_node_dict)
    graph, simp_node_dict, edge_dict = graph_reduction(graph, contig_dict, simp_node_dict, edge_dict)
    # graph_to_gfa(graph, simp_node_dict, edge_dict)

def flip_graph_bfs(graph: Graph, node_dict: dict, edge_dict: dict, dp_dict: dict, init_ori):
    """
    Flip all the node orientation.

    return an node_dict, which only contains one orientation per node for simplicity.
    rename all the used node to positive, and forbidden the opponent node.
    """

    pick_dict = {}
    while set(dp_dict):
        seg_no = source_node_via_dp(dp_dict)
        source_pos, source_neg = node_dict[seg_no]
        graph.vp.visited[source_pos] = 0
        graph.vp.visited[source_neg] = 0
        queue = []
        queue.append([node_dict[seg_no], init_ori]) 

        while queue:
            (v_pos, v_neg), ori = queue.pop()
            dp_dict.pop(graph.vp.id[v_pos])
            
            u = None
            if ori == 1:
                u = v_pos
                pick_dict[graph.vp.id[u]] = '+'
                # print_vertex(graph, v_neg, "node to reverse")
                for e in list(v_neg.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            else:
                u = v_neg
                pick_dict[graph.vp.id[u]] = '-'
                # print_vertex(graph, v_pos, "node to reverse")
                for e in list(v_pos.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            
            graph.vp.visited[v_pos] = 1
            graph.vp.visited[v_neg] = 1
            # add further nodes into the queue TODO, all or out only
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1:
                    graph.vp.visited[adj_node] = 0
                    queue.append([node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]])
    # verify sorted graph
    print("-------------------------verify graph----------------------")
    check = True

    check = check and (len(pick_dict) == len(node_dict))
    for key, item in pick_dict.items():
        # print("pick: ", key, item)
        v_pos, v_neg = node_dict[key]
        if item == '+':
            if v_neg.in_degree() + v_neg.out_degree() > 0:
                print_vertex(graph, v_neg, "pick error found")
                check = False
        else:
            if v_pos.in_degree() + v_pos.out_degree() > 0:
                print_vertex(graph, v_pos, "pick error found")
                check = False

    for key, (v_pos, v_neg) in node_dict.items():
        if not v_pos.in_degree() + v_pos.out_degree() == 0 and not v_neg.in_degree() + v_neg.out_degree() == 0:
            check = False
            print_vertex(graph, v_pos, "erroroness node pos found")
            print_vertex(graph, v_neg, "erroroness node neg found")
            print("re-pict nodes, pick both")
            pick_dict[key] = 0
            graph.vp.id[v_pos] = graph.vp.id[v_pos] + "a"
            graph.vp.id[v_neg] = graph.vp.id[v_neg] + "b"

    check = check and len(node_dict) == len(pick_dict)
    if check: print("Graph is verified")
    print("-------------------------end verify------------------------")

    simp_node_dict = {}
    for seg_no, pick in pick_dict.items():
        if pick == '+':
            picked = node_dict[seg_no][0]
        else:
            picked = node_dict[seg_no][1]
        graph.vp.ori[picked] = 1
        simp_node_dict[seg_no] = picked
    if is_DAG(graph):
        print("The graph is acylic")
    else:
        print("The graph is cyclic")
    
    return graph, simp_node_dict, edge_dict

# def simplify_graph(graph: Graph, node_dict: dict, pick_dict: dict):
#     """
#     return an node_dict, which only contains one orientation per node for simplicity.
#     rename all the used node to positive, and forbidden the opponent.
#     """
#     simp_node_dict = {}
#     for seg_no, pick in pick_dict.items():
#         if pick == '+':
#             picked = node_dict[seg_no][0]
#         else:
#             picked = node_dict[seg_no][1]
#         graph.vp.ori[picked] = 1
#         simp_node_dict[seg_no] = picked
#     if is_DAG(graph):
#         print("The graph is acylic")
#     else:
#         print("The graph is cyclic")
#     return graph, simp_node_dict

def graph_grouping(graph: Graph, simp_node_dict: dict, forward, reverse, partition_length_cut_off):
    """
    Maximimize graph connectivity by minimizing node with 0 in-degree or out-degree, detect and remove all the cycles.
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

    group_no = 1
    groups = {}
    for v in simp_node_dict.values():
        # grouping
        if graph.vp.group[v] == -1:
            graph, group_no, groups = bfs_grouping(graph, v, group_no, groups)
            group_no = group_no + 1

    print("number of groups:", len(groups))
    # for key, item in groups.items():
    #     print("group number: ", key, " member: ", [(graph.vp.id[v],graph.vp.partition[v]) for v in item])

    # connect sub-graphs based on pair-end reads information
    # TODO
    return graph, groups

def graph_partition(graph: Graph, simp_node_dict: dict, edge_dict: dict, dp_dict: dict):
    """
    start from source node, run BFS to determine the partition for each nodes until reaching the sink node.
    partition is determined based on seg-length
    """
    None
    # queue = []
    # queue.append(src)
    # while queue:
    #     u = queue.pop()
    #     if graph.vp.id[u] == 'src':

    #         for v in u.out_neighbors():
    #             queue.append(v)
    #     elif graph.vp.id[u] == 'sink':
    #         None
    #     else:
    #         None

    # # pick a source node based on depth.
    # seg_no = source_node_via_dp(dp_dict)
    # src = simp_node_dict[seg_no]
    # graph.vp.visited[src] = 0
    # max_len_per_partition = 1000
    # max_node_per_partition = len(simp_node_dict) // 10
    # # estimate #partitions in the graph.
    # print("max len per partition: ", max_len_per_partition, ", max node per partition: ", max_node_per_partition)

    # queue = []
    # queue.append(src)
    # partition_dict = {}
    # partition_len_dict = {}
    # partition = 1
    # length = 0
    # partition_dict[partition] = []
    # while queue:
    #     u = queue.pop()
    #     # print_vertex(graph, u, "current vertex")
    #     # label with partition tag
    #     # consider all the 
    #     curr_seq = graph.vp.seq[u]
    #     curr_len = len(curr_seq)

    #     if len(partition_dict[partition]) >= max_node_per_partition:
    #         partition = partition + 1
    #         length = 0
        
    #     if partition not in partition_dict:
    #         partition_dict[partition] = []

    #     if length + curr_len > max_len_per_partition and length + curr_len > (max_len_per_partition * 1.5):
    #         # overfit to put into the current partition 
    #         next_partition = partition + 1
    #         while next_partition in partition_dict:
    #             if len(partition_dict[next_partition]) >= max_node_per_partition:
    #                 next_partition = next_partition + 1
    #             else:
    #                 break
    #         if next_partition not in partition_dict:
    #             partition_dict[next_partition] = []
    #         # print("overfit in current partition ", partition, ", jump to next partition ", next_partition)
    #         partition_dict[next_partition].append(u)
    #         graph.vp.partition[u] = next_partition
    #     else:
    #         # best fit
    #         # print("best fit in partition ", partition)
    #         partition_dict[partition].append(u)
    #         graph.vp.partition[u] = partition
        
    #     for v in u.out_neighbors():
    #         if graph.vp.partition[v] == -1:
    #             queue.append(v)
    # # for s, v in simp_node_dict.items():
    # #     # print_vertex(graph, v, "partition check")
    # #     if graph.vp.partition[v] == -1:
    # #         print("current node with seg no: ", s, " is in partition ", graph.vp.partition[v])
    # return graph

def longest_distance(graph: Graph, s):
    return

def graph_stat(graph: Graph, simp_node_dict: dict, edge_dict: dict):
    for seg_no, v in simp_node_dict.items():
        if v.in_degree() == 0 and v.out_degree() != 0:
            print_vertex(graph, v, "vertex with only out degree")
        if v.out_degree() == 0 and v.in_degree() != 0:
            print_vertex(graph, v, "vertex with only in degree")

    # for (u, v), e in edge_dict.items():
    #         print_edge(graph, e, "flow {0}".format(graph.ep.flow[e]))

def graph_dfs(graph: Graph, source):
    """
    Count the maximum depth path the source node can reach via directed edge in the graph
    """
    visited = {}
    for u in graph.vertices():
        visited[u] = False
    
    def dfs_helper(graph: Graph, u, visited):
        if u.out_degree() == 0:
            return [u]
        else:
            rtn = []
            for v in u.out_neighbors():
                # print_vertex(graph, v)
                if not visited[v]:
                    visited[v] = True
                    cmp = dfs_helper(graph, v, visited)
                    cmp.insert(0, (v, len(graph.vp.seq[v])))
                    rtn = rtn if len(cmp) < len(rtn) else cmp
            return rtn

    return dfs_helper(graph, source, visited)

def source_node_via_dp(dp_dict: dict):
    """
    return the pos-neg node with maximum depth
    """
    seg_no = max(dp_dict, key=dp_dict.get)
    print("source node id: ", seg_no, ", depth: ", dp_dict[seg_no])
    return seg_no

def swap_node_ori_name(graph: Graph, node_dict: dict, seg_no):
    """
    swap the node's orientation, and update related edges
    """
    v_pos, v_neg = node_dict[seg_no]
    graph.vp.ori[v_pos] = -1
    graph.vp.ori[v_neg] = 1
    node_dict[seg_no] = (v_neg, v_pos)

    return graph, node_dict

def reverse_edge(graph: Graph, edge, node_dict: dict, edge_dict: dict):
    """
    reverse an edge with altered orientation and direction.
    """
    tmp_s = edge.source()
    tmp_t = edge.target()
    
    edge_dict.pop((tmp_s, tmp_t))

    tmp_s_pos, tmp_s_neg = node_dict[graph.vp.id[tmp_s]]
    tmp_t_pos, tmp_t_neg = node_dict[graph.vp.id[tmp_t]]
    t = tmp_s_pos if graph.vp.ori[tmp_s] == -1 else tmp_s_neg
    s = tmp_t_pos if graph.vp.ori[tmp_t] == -1 else tmp_t_neg

    o = graph.ep.overlap[edge]
    graph.remove_edge(edge)
    e = graph.add_edge(s, t)
    graph.ep.overlap[e] = o
    edge_dict[(s, t)] = e

    return graph, e, edge_dict

def print_edge(graph, e, s=""):
    print(s, " edge: ", graph.vp.id[e.source()],graph.vp.ori[e.source()], "->", graph.vp.id[e.target()], graph.vp.ori[e.target()])

def print_vertex(graph, v, s=""):
    print(s, " vertex: ", graph.vp.id[v], ", dp: ", graph.vp.dp[v], ", ori: ", graph.vp.ori[v], ", in_degree: ", v.in_degree(), ", out_degree: ", v.out_degree())

def assign_edge_flow(graph: Graph, simp_node_dict: dict, edge_dict: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    un_assigned_edge = len(edge_dict)
    print("Total edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for node in simp_node_dict.values():
            w = graph.vp.dp[node]
            if node.in_degree() == 1:
                for e in node.in_edges():
                    graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1
            if node.out_degree() == 1:
                for e in node.out_edges():
                    graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1

    # converage iteration
    converage_flag = 0
    while True:          
        for node in simp_node_dict.values():
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
                    graph.ep.flow[e] = in_w
                    un_assigned_edge = un_assigned_edge - 1
            if out_d == 1:
                for e in out_e:
                    graph.ep.flow[e] = out_w
                    un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  

    print("un-assigned edges after node-weight converage iteration : ", un_assigned_edge)
    # def subfinder(mylist, pattern):
    #     rtn = -1
    #     for i in range(len(mylist)):
    #         if mylist[i] == pattern[0] and mylist[i:i+len(pattern)] == pattern:
    #             rtn = i
    #             break
    #     return rtn
    # deal with rest of un assigned edges
    for (u,v), e in edge_dict.items():
        uv = [graph.vp.id[u],graph.vp.id[v]]
        if graph.ep.flow[e] == 0.0:
            
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
                print("low flow error: ", uv, assign_flow)
                print("u_flow_remain: ", u_flow_remain, " u_degree: ", u_degree)
                print("v_flow_remain: ", v_flow_remain, " v_degree: ", v_degree)
            else:
                print("manual assigned edge: ", uv, assign_flow)
                un_assigned_edge = un_assigned_edge - 1
                graph.ep.flow[e] = assign_flow

    print("un-assigned edges after manual assign iteration : ", un_assigned_edge)
            # assigned = False
            # for (cno, clen, ccov), contig in contig_dict.items():
            #     ind = subfinder(contig, uv)
            #     if ind != -1:
            #         assigned = True
            #         if ind + 1 == len(contig):
            #             print("unavailable edge, error", uv, ind)
            #         elif ind + 1 == len(contig) - 1:
            #             print("ending edge in contig", uv, ind)
            #             prior_edge = edge_dict[(simp_node_dict[contig[ind - 1]], simp_node_dict[contig[ind]])]
            #             graph.ep.flow[e] = graph.ep.flow[prior_edge]
            #             # only prior edge support
            #         elif ind == 0:
            #             print("starting edge in contig", uv, ind)
            #             forward_edge = edge_dict[(simp_node_dict[contig[ind + 1]], simp_node_dict[contig[ind + 2]])]
            #             graph.ep.flow[e] = graph.ep.flow[forward_edge]
            #             # only forward edge support
            #         else:
            #             print("interior edge in contig", uv, ind)
            #             prior_edge = edge_dict[(simp_node_dict[contig[ind - 1]], simp_node_dict[contig[ind]])]
            #             forward_edge = edge_dict[(simp_node_dict[contig[ind + 1]], simp_node_dict[contig[ind + 2]])]
            #             graph.ep.flow[e] = (graph.ep.flow[prior_edge] + graph.ep.flow[forward_edge])/2
            #             # both support, pick the ave for the adjacent edges

def depth_decomposition(graph: Graph, simp_node_dict: dict):
    None

def src_sink(graph: Graph, simp_node_dict: dict, edge_dict: dict, contig_dict: dict):
    """
    add manufactured source/sink node for finding maximum flow over graph
    """
    src = graph.add_vertex()
    sink = graph.add_vertex()
    graph.vp.id[src] = 'src'
    graph.vp.id[sink] = 'sink'
    graph.vp.ori[src] = 1
    graph.vp.ori[sink] = 1
    graph.vp.color[src] = 'white'
    graph.vp.color[sink] = 'white'

    src_is_connect = False
    sink_is_connect = False
    inf_value = 2**20 # may replace with numpy inf TODO
    for u in simp_node_dict.values():
        if u.in_degree() == 0 and u.out_degree() == 0:
            print_vertex(graph, u, "isolated vertex")
            continue

        if u.in_degree() == 0:
            print_vertex(graph, u, "vertex with 0 in degree, connect to src")
            e = graph.add_edge(src, u)
            graph.ep.flow[e] = inf_value
            src_is_connect = True
        
        if u.out_degree() == 0:
            print_vertex(graph, u, "vertex with 0 out degree, connect to sink")
            e = graph.add_edge(u, sink)
            graph.ep.flow[e] = inf_value
            sink_is_connect = True
    # TODO create a src and sink node if is cyclic
    # attempt to break the cycle by removing one edge contains the maximum depth among all the edges
    if not src_is_connect or not sink_is_connect:
        # obtain all the involving node from contigs
        contig_node_set = set()
        for (contig, contig_rev) in contig_dict.values():
            for n in contig:
                contig_node_set.add(n)

        largest_circuit = max(all_circuits(graph), key=lambda c: len(c))
        circuit_nodes = [simp_node_dict[graph.vp.id[v]] for v in largest_circuit]
        for v in sorted(circuit_nodes, key=lambda n: graph.vp.dp[n], reverse=True):
            if src_is_connect and sink_is_connect:
                break
            if graph.vp.id[v] in contig_node_set:
                print("already in the contig set, never consider to break the relevant edge")
                continue
            if v.out_degree() == 1 and v.in_degree() != 0 and not sink_is_connect:
                # break the edge, connect to sink
                u = list(v.out_neighbors())[0]
                print_edge(graph, edge_dict[(v, u)], "edge be removed")
                graph.remove_edge(edge_dict[(v, u)])
                edge_dict.pop((v, u))
                e = graph.add_edge(v, sink)
                graph.ep.flow[e] = inf_value
                sink_is_connect = True

            if v.in_degree() == 1 and v.out_degree() != 0 and not src_is_connect:
                # break the edge, connect to src        
                u = list(v.in_neighbors())[0]
                print_edge(graph, edge_dict[(u, v)], "edge be removed")
                graph.remove_edge(edge_dict[(u, v)])
                edge_dict.pop((u, v))
                e = graph.add_edge(src, v)
                graph.ep.flow[e] = inf_value
                src_is_connect = True

    if not src_is_connect:
        print("source is not connected still")
    if not sink_is_connect:
        print("sink is not connected still")
    if src_is_connect and sink_is_connect:
        print("src and sink is initialised")
    if not is_DAG(graph):
        print("The graph still have cycles")
    return graph, simp_node_dict, edge_dict, src, sink

def max_flow(graph: Graph, simp_node_dict: dict, edge_dict: dict, contig_dict: dict):
    print("find max flow")
    c1 = contig_dict[("1", 9557, 1897.189396)]
    src = simp_node_dict[c1[0]]
    
    sink = simp_node_dict[c1[-1]]
    cap = graph.ep.flow
    res = flow.boykov_kolmogorov_max_flow(graph, src, sink, cap)
    res.a = cap.a - res.a  # the actual flow
    max_flow = sum(res[e] for e in sink.in_edges())
    print(max_flow)

    for i in range(len(c1)-1):
        u = simp_node_dict[c1[i]]
        v = simp_node_dict[c1[i+1]]
        e = edge_dict[(u,v)]
        f = graph.ep.flow[e]
        print_edge(graph, e, "edge flow: {0}, res flow: {1}".format(f, res[e]))

def graph_reduction(graph: Graph, contig_dict: dict, simp_node_dict: dict, edge_dict: dict):
    for (cno, clen, ccov), (contigs,edge_flow) in contig_dict.items():
        print("Contig: ", cno, clen, ccov)
    # print("-------------------------------------------------")
    #     covs = []
    #     kcs = []
    #     for v in contigs:
    #         if v not in simp_node_dict:
    #             continue
    #         covs.append(graph.vp.dp[simp_node_dict[v]])
    #         kcs.append(graph.vp.kc[simp_node_dict[v]])
    #     covs = sorted(covs)

    #     print("node coverage - min: ", min(covs), " ave: ", (sum(covs))/len(covs), " median: ", numpy.median(covs))
    #     print("kmer count total: ", sum(kcs))
    #     print("all covs: ", [(graph.vp.dp[simp_node_dict[v]], v) for v in contigs])
    #     edge_flow = []
    #     for i in range(len(contigs)-1):
    #         u = simp_node_dict[contigs[i]]
    #         v = simp_node_dict[contigs[i+1]]
    #         e = edge_dict[(u,v)]
    #         f = graph.ep.flow[e]
    #         edge_flow.append(f)
    #         print_edge(graph, e, "edge flow: {0}".format(f))
    #     print("edge flow coverage - min: ", min(edge_flow), " ave: ", (sum(edge_flow))/len(edge_flow), " median: ", numpy.median(edge_flow))
    #     for v in contigs:
    #         n = simp_node_dict[v]
    #         if n.in_degree() <= 1 and n.out_degree() <= 1:
    #             print_vertex(graph, n, "single-connect nodes with in degree: {0}, out degree: {1}".format(n.in_degree(), n.out_degree()))
    #     print("-------------------------------------------------")
        # ccov = min(edge_flow)
        j = 1
        for i in contigs:
            if j >= len(contigs):
                break
            ii = contigs[j]
            j = j + 1
            u = simp_node_dict[i] if i in simp_node_dict else None
            v = simp_node_dict[ii] if ii in simp_node_dict else None
            e = edge_dict[(u,v)] if (u,v) in edge_dict else None 


            # edge may be eliminated from previous execution
            if e == None or u == None or v == None:
                continue

            # reduce the depth for involving edge
            if graph.ep.flow[e] - ccov <= 0:
                graph.ep.flow[e] = 0
                edge_dict.pop((u,v))
                print_edge(graph, e, "edge been removed")
            else:
                graph.ep.flow[e] = graph.ep.flow[e] - ccov

            # reduce the depth for involving node, gray color for forbiddened node
            if graph.vp.dp[u] - ccov <= 0:
                graph.vp.dp[u] = 0
                graph.vp.color[u] = 'gray'
                if i in simp_node_dict:
                    simp_node_dict.pop(i)
                    print("node ", i, "been removed")
                else:
                    print("node ", i, "is removed already")
            else:
                graph.vp.dp[u] = graph.vp.dp[u] - ccov


            if graph.vp.dp[v] - ccov <= 0: 
                graph.vp.dp[v] = 0
                graph.vp.color[v] = 'gray'
                if ii in simp_node_dict:
                    simp_node_dict.pop(ii)
                    print("node ", ii, "been removed")
                else:
                    print("node ", ii, "is removed already")
            else:
                graph.vp.dp[v] = graph.vp.dp[v] - ccov
            
            # # update edges
            if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (u,v) in edge_dict:
                edge_dict.pop((u,v))

    return graph, simp_node_dict, edge_dict

if __name__ == "__main__":
    main()
