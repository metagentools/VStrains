#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
from graph_tool import GraphView, _in_degree
# import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw
# from graph_tool.clustering import local_clustering

"""
This file is used to store out of date scripts to keep a record for further usage needed.
"""

def graph_grouping(graph: Graph, simp_node_dict: dict, forward, reverse, partition_length_cut_off):
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
    #     # #print_vertex(graph, u, "current vertex")
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
    # #     # #print_vertex(graph, v, "partition check")
    # #     if graph.vp.partition[v] == -1:
    # #         print("current node with seg no: ", s, " is in partition ", graph.vp.partition[v])
    return graph

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
            #print_vertex(graph, u, "isolated vertex")
            continue

        if u.in_degree() == 0:
            #print_vertex(graph, u, "vertex with 0 in degree, connect to src")
            e = graph.add_edge(src, u)
            graph.ep.flow[e] = inf_value
            src_is_connect = True
        
        if u.out_degree() == 0:
            #print_vertex(graph, u, "vertex with 0 out degree, connect to sink")
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
                #print_edge(graph, edge_dict[(v, u)], "edge be removed")
                graph.remove_edge(edge_dict[(v, u)])
                edge_dict.pop((v, u))
                e = graph.add_edge(v, sink)
                graph.ep.flow[e] = inf_value
                sink_is_connect = True

            if v.in_degree() == 1 and v.out_degree() != 0 and not src_is_connect:
                # break the edge, connect to src        
                u = list(v.in_neighbors())[0]
                #print_edge(graph, edge_dict[(u, v)], "edge be removed")
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
    # res = flow.boykov_kolmogorov_max_flow(graph, src, sink, cap)
    # res.a = cap.a - res.a  # the actual flow
    # max_flow = sum(res[e] for e in sink.in_edges())
    # print(max_flow)

    # for i in range(len(c1)-1):
    #     u = simp_node_dict[c1[i]]
    #     v = simp_node_dict[c1[i+1]]
    #     e = edge_dict[(u,v)]
    #     f = graph.ep.flow[e]
        #print_edge(graph, e, "edge flow: {0}, res flow: {1}".format(f, res[e]))
    return


"""
assign edge flow
"""
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



"""
graph reduction
"""
    # print("-------------------------------------------------")
    #     covs = []
    #     kcs = []
    #     for v in contig:
    #         if v not in simp_node_dict:
    #             continue
    #         covs.append(graph.vp.dp[simp_node_dict[v]])
    #         kcs.append(graph.vp.kc[simp_node_dict[v]])
    #     covs = sorted(covs)

    #     print("node coverage - min: ", min(covs), " ave: ", (sum(covs))/len(covs), " median: ", numpy.median(covs))
    #     print("kmer count total: ", sum(kcs))
    #     print("all covs: ", [(graph.vp.dp[simp_node_dict[v]], v) for v in contig])
    #     edge_flow = []
    #     for i in range(len(contig)-1):
    #         u = simp_node_dict[contig[i]]
    #         v = simp_node_dict[contig[i+1]]
    #         e = edge_dict[(u,v)]
    #         f = graph.ep.flow[e]
    #         edge_flow.append(f)
    #         print_edge(graph, e, "edge flow: {0}".format(f))
    #     print("edge flow coverage - min: ", min(edge_flow), " ave: ", (sum(edge_flow))/len(edge_flow), " median: ", numpy.median(edge_flow))
    #     for v in contig:
    #         n = simp_node_dict[v]
    #         if n.in_degree() <= 1 and n.out_degree() <= 1:
    #             print_vertex(graph, n, "single-connect nodes with in degree: {0}, out degree: {1}".format(n.in_degree(), n.out_degree()))
    #     print("-------------------------------------------------")