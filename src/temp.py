#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
from this import d
from graph_tool import GraphView, _in_degree
# import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw
# from graph_tool.clustering import local_clustering

from graph_converter import *

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

def contig_reduction():
    return


def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, consistent_node: dict, inconsistent_node: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    #TODO assign non-branch path node dp first, keep the start and end node dp.
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for no, node in simp_node_dict.items():
        if no not in consistent_node:
            continue
        w = graph.vp.dp[node]
        in_d = len([n for n in node.in_neighbors() if graph.vp.color[n] == 'black'])
        node_in_edges = [e for e in node.in_edges() if graph.ep.color[e] == 'black']
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if graph.ep.flow[e] == 0.0 and (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = round((w + graph.vp.dp[src])/2, 2)
                    else:
                        graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1

        out_d = len([n for n in node.out_neighbors() if graph.vp.color[n] == 'black'])
        node_out_edges = [e for e in node.out_edges() if graph.ep.color[e] == 'black']
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = round((w + graph.vp.dp[src])/2, 2)
                    else:
                        graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1
    print("un-assigned edges after init iteration : ", un_assigned_edge)

    # coverage iteration
    converage_flag = 0
    while True:          
        for no, node in simp_node_dict.items():
            if no not in consistent_node:
                continue
            in_d = len([n for n in node.in_neighbors() if graph.vp.color[n] == 'black'])
            node_in_edges = [e for e in node.in_edges() if graph.ep.color[e] == 'black']
            in_w = graph.vp.dp[node]
            in_e = []
            for e in node_in_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    in_d = in_d - 1
                    in_w = in_w - f
                else:
                    in_e.append(e)
            if in_d == 1:
                for e in in_e:
                    if graph.ep.flow[e] == 0.0 and (graph.vp.id[e.source()], no) in simp_edge_dict:
                        if in_w <= 0:
                            print("in low edge flow error: ", graph.vp.id[e.source()], " -> ", graph.vp.id[e.target()], in_w)
                        graph.ep.flow[e] = in_w
                        un_assigned_edge = un_assigned_edge - 1

            out_d = len([n for n in node.out_neighbors() if graph.vp.color[n] == 'black'])
            node_out_edges = [e for e in node.out_edges() if graph.ep.color[e] == 'black']
            out_w = graph.vp.dp[node]
            out_e = []
            for e in node_out_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    out_d = out_d - 1
                    out_w = out_w - f
                else:
                    out_e.append(e)
            if out_d == 1:
                for e in out_e:
                    if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[e.target()]) in simp_edge_dict:
                        if out_w <= 0:
                            print_vertex(graph, node, "Current node")
                            print("out low edge flow error: ", graph.vp.id[e.source()], " -> ", graph.vp.id[e.target()], out_w)
                        graph.ep.flow[e] = out_w
                        un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  
    print("un-assigned edges after node-weight coverage iteration : ", un_assigned_edge)

    # FIXME for inconsistent node part

    for (u,v), e in simp_edge_dict.items():
        if graph.ep.flow[e] == 0.0:
            if u in inconsistent_node and v in inconsistent_node:
                continue

            u = e.source()
            v = e.target()
            u_flow_remain = graph.vp.dp[u]
            u_degree = len([n for n in u.out_neighbors() if graph.vp.color[n] == 'black'])
            u_out_edges = [e for e in u.out_edges() if graph.ep.color[e] == 'black']

            for ue in u_out_edges:
                if graph.ep.flow[ue] != 0.0:
                    # assigned edge
                    u_degree = u_degree - 1
                    u_flow_remain = u_flow_remain - graph.ep.flow[ue]

            v_flow_remain = graph.vp.dp[v]
            v_degree = len([n for n in v.in_neighbors() if graph.vp.color[n] == 'black'])
            v_in_edges = [e for e in v.in_edges() if graph.ep.color[e] == 'black']

            for ve in v_in_edges:
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
                print("manual assigned edge: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                if graph.ep.flow[e] == 0.0:
                    un_assigned_edge = un_assigned_edge - 1
                    graph.ep.flow[e] = assign_flow
    print("un-assigned edges after manual assign iteration : ", un_assigned_edge)
    u_nodes = []
    for (u,v), e in simp_edge_dict.items():
        if graph.ep.flow[e] == 0.0:
            u_nodes.append(int(u))
            u_nodes.append(int(v))
        if u == '49355':
            print_edge(graph, e, "")
    print("unassigned related node: ", u_nodes)
    print("-----------------------assign edge flow end--------------------")


def node_rebalance_via_simple_path(graph: Graph, simp_edge_dict: dict):
    # get all the simple path (no branch path)
    simple_paths = simp_path(graph, simp_edge_dict)
    for sp in simple_paths:
        if DEBUG_MODE:
            print("path: ", [int(graph.vp.id[u]) for u in sp])
        dp_mean_p = numpy.mean([graph.vp.dp[u] for u in sp])
        # rebalance the node depth
        for u in sp:
            graph.vp.dp[u] = dp_mean_p

def node_dp_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Re-balance all the node depth, such that for any given node u, sum(dp(in_neighbors)) == dp(u) == sum(dp(out_neigbors))
    where dp(u) = max(sum(dp(in_neighbors)), dp(u), sum(dp(out_neigbors)))
    plan to run BFS from the source node (how to pick?) until all the node is rebalanced.
    source node may pick from the node with 0 indegree.
    """
    def source_node_via_dp(dp_dict: dict):
        """
        return the pos-neg node with minimum depth
        """
        seg_no = max(dp_dict, key=dp_dict.get)
        print("source node id: ", seg_no, ", depth: ", dp_dict[seg_no])
        return seg_no
    print("-------------------------node depth rebalance----------------------")
    graph.vp.visited = graph.new_vertex_property("int16_t", val=-1)
    dp_dict = {}
    for no, u in simp_node_dict.items():
        dp_dict[no] = graph.vp.dp[u]
        graph.vp.visited[u] = -1

    while set(dp_dict):
        for no in dp_dict.keys():
            dp_dict[no] = graph.vp.dp[simp_node_dict[no]]
       
        seg_no = source_node_via_dp(dp_dict)
        src = simp_node_dict[seg_no]
        graph.vp.visited[src] = 0
        fifo_queue = [src]
        # BFS
        while fifo_queue:
            u = fifo_queue.pop()
            dp_dict.pop(graph.vp.id[u])

            u_dp = graph.vp.dp[u]

            in_es = [e for e in u.in_edges() if graph.ep.color[e] == 'black']
            in_flow_sum = numpy.sum([graph.ep.flow[e] for e in in_es])

            out_es = [e for e in u.out_edges() if graph.ep.color[e] == 'black']
            out_flow_sum = numpy.sum([graph.ep.flow[e] for e in out_es])

            max_dp = max(in_flow_sum, u_dp, out_flow_sum)
            graph.vp.dp[u] = max_dp
            if DEBUG_MODE:
                print("Node: {0} dp from {1} to {2}".format(graph.vp.id[u], u_dp, round(max_dp)))
                print("in {0} dp {1} out {2}".format(in_flow_sum, u_dp, out_flow_sum))

            if len(in_es) == 1:
                graph.ep.flow[in_es[0]] = max_dp
            elif len(in_es) > 1:
                in_flow_load = max_dp
                for i in range(len(in_es) - 1):
                    e = in_es[i]
                    flow = round((graph.ep.flow[e]/in_flow_sum) * max_dp)
                    graph.ep.flow[e] = flow
                    in_flow_load -= flow
                graph.ep.flow[in_es[-1]] = round(in_flow_load)
            else:
                print("No in edges, skip")

            if len(out_es) == 1:
                graph.ep.flow[out_es[0]] = max_dp
            elif len(out_es) > 1:
                out_flow_load = max_dp
                for i in range(len(out_es) - 1):
                    e = out_es[i]
                    flow = round((graph.ep.flow[e]/out_flow_sum) * max_dp)
                    graph.ep.flow[e] = flow
                    out_flow_load -= flow
                graph.ep.flow[out_es[-1]] = round(out_flow_load)
            else:
                print("No out edges, skip")

            graph.vp.visited[u] = 1

            # append neighbors into fifo queue
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1 and graph.vp.color[adj_node] == 'black':
                    graph.vp.visited[adj_node] = 0
                    fifo_queue.append(adj_node)
    print("Done...")
    print("Check node dp balance now...")
    # sum check
    false_count = 0
    for no, u in simp_node_dict.items():
        dp = graph.vp.dp[u]
        in_es = [n for n in u.in_edges() if graph.ep.color[n] == 'black']
        in_flow_sum = numpy.sum([graph.ep.flow[v] for v in in_es])

        out_es = [n for n in u.out_edges() if graph.ep.color[n] == 'black']
        out_flow_sum = numpy.sum([graph.ep.flow[v] for v in out_es])

        if len(in_es) != 0:
            if abs(in_flow_sum - dp) > 100:
                if len(out_es) == 0:
                    # edge case, manual fix
                    print("mannual fix, sink node", no)
                    graph.vp.dp[no] = max(in_flow_sum, dp)
                else:
                    print("Current node: ", no, dp)
                    false_count += 1
                    for e in u.in_edges():
                        if graph.ep.color[e] == 'black':
                            print_edge(graph, e, "in edge")  
                    print("Dp diff: ", in_flow_sum - dp)

        if len(out_es) != 0:
            if abs(out_flow_sum - dp) > 100:
                if len(in_es) == 0:
                    # edge case, manual fix
                    print("mannual fix, source node", no)
                    graph.vp.dp[no] = max(out_flow_sum, dp)
                else:
                    print("Current node: ", no, dp)
                    false_count += 1
                    for e in u.out_edges():
                        if graph.ep.color[e] == 'black':
                            print_edge(graph, e, "out edge")  
                    print("Dp diff: ", out_flow_sum - dp)             
    print("False case: ", false_count)
    print("-------------------------node depth rebalance end----------------------")
    return


# refine all the node coverage by reassign all the edge flow
    for (u,v), e in simp_edge_dict.items():
        graph.ep.flow[e] = 0.0
    
    print("Final refine")
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict, consistent_node, inconsistent_node)

    # re assign the node depth = max(curr_dp, out_neighbor_dp_sum if no branch, in_neighbor_dp_sum if no branch, out edge flow sum, in edge flow sum)
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        node = queue.pop()
        visited.append(graph.vp.id[node])
        curr_dp = graph.vp.dp[node]

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
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
        graph.vp.dp[node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[node]))
        assign_edge_flow2(graph, simp_node_dict, simp_edge_dict)
        
        for n in node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)


def assign_edge_flow2(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total unassigned edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for no, node in simp_node_dict.items():
        w = graph.vp.dp[node]
        in_d = len([n for n in node.in_neighbors()])
        node_in_edges = [e for e in node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if src.out_degree() == 1:
                    graph.ep.flow[e] = numpy.max([w, graph.ep.flow[e], graph.vp.dp[src]])
                    un_assigned_edge = un_assigned_edge - 1
        out_d = len([n for n in node.out_neighbors()])
        node_out_edges = [e for e in node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if tgt.in_degree() == 1:
                    graph.ep.flow[e] = numpy.max([w, graph.ep.flow[e], graph.vp.dp[tgt]])
                    un_assigned_edge = un_assigned_edge - 1

    # coverage iteration
    converage_flag = 0
    while True:          
        for no, node in simp_node_dict.items():
            in_d = len([n for n in node.in_neighbors()])
            node_in_edges = [e for e in node.in_edges()]
            in_w = graph.vp.dp[node]
            in_e = []
            for e in node_in_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    in_d = in_d - 1
                    in_w = in_w - f
                else:
                    in_e.append(e)
            if in_d == 1:
                for e in in_e:
                    graph.ep.flow[e] = max(in_w,graph.ep.flow[e])
                    un_assigned_edge = un_assigned_edge - 1

            out_d = len([n for n in node.out_neighbors()])
            node_out_edges = [e for e in node.out_edges()]
            out_w = graph.vp.dp[node]
            out_e = []
            for e in node_out_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    out_d = out_d - 1
                    out_w = out_w - f
                else:
                    out_e.append(e)
            if out_d == 1:
                for e in out_e:
                        graph.ep.flow[e] = max(out_w,graph.ep.flow[e])
                        un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  

    # double cross consistent node assign
    for (u,v), e in simp_edge_dict.items():

        src = e.source()
        tgt = e.target()

        u_flow_remain = graph.vp.dp[src]
        u_degree = len([n for n in src.out_neighbors()])
        u_out_edges = [e for e in src.out_edges()]
        for ue in u_out_edges:
            if graph.ep.flow[ue] != 0.0:
                # assigned edge
                u_degree = u_degree - 1
                u_flow_remain = u_flow_remain - graph.ep.flow[ue]

        v_flow_remain = graph.vp.dp[tgt]
        v_degree = len([n for n in tgt.in_neighbors()])
        v_in_edges = [e for e in tgt.in_edges()]
        for ve in v_in_edges:
            if graph.ep.flow[ve] != 0.0:
                # assigned edge
                v_degree = v_degree - 1
                v_flow_remain = v_flow_remain - graph.ep.flow[ve]

        u_flow_remain = u_flow_remain if u_flow_remain > 0 else 0
        v_flow_remain = v_flow_remain if v_flow_remain > 0 else 0
        udiv = (u_flow_remain / u_degree) if u_degree != 0 else 0
        vdiv = (v_flow_remain / v_degree) if v_degree != 0 else 0
        if udiv != 0 and vdiv != 0:
            assign_flow = (udiv + vdiv) / 2
        elif udiv != 0:
            assign_flow = udiv
        elif vdiv != 0:
            assign_flow = vdiv
        else:
            assign_flow = 0
        un_assigned_edge = un_assigned_edge - 1
        graph.ep.flow[e] = max(assign_flow,graph.ep.flow[e])
    print("-----------------------assign edge flow v2--------------------")
    return


def graph_split_by_breadth(graph: Graph, simp_node_dict: dict):
    srcs = []
    for no, node in simp_node_dict.items():
        if node.in_degree() == 0:
            srcs.append(node)
    
    queue = []
    visited = []
    srcs = sorted(srcs, key=lambda x: graph.vp.dp[x], reverse=True)
    queue.append(srcs[0])
    split_by_level = {}
    split_by_level[0] = [graph.vp.id[srcs[0]]]
    curr_level = 1
    while queue:
        node = queue.pop()
        visited.append(graph.vp.id[node])
        split_by_level[curr_level] = []
        for out_n in node.out_neighbors():
            if not graph.vp.id[out_n] in visited:
                split_by_level[curr_level].append(graph.vp.id[out_n])
                queue.append(out_n)
        curr_level += 1
    sorted_level_pair = sorted(split_by_level.items(), key = lambda x: x[0])
    allNodes = []
    for level, nodes in sorted_level_pair:
        print("Current level: ", level)
        print("Nodes: ", [int(v) for v in nodes])
        [allNodes.append((int(v))) for v in nodes]
    print("All nodes: ", allNodes)









# print("Final refine")
    for no, node in simp_node_dict.items():
        curr_dp = graph.vp.dp[node]
        in_edges = list(node.in_edges())
        in_sum = numpy.sum([graph.ep.flow[e] for e in in_edges])
        out_edges = list(node.out_edges())
        out_sum = numpy.sum([graph.ep.flow[e] for e in out_edges])
        if in_edges != [] and curr_dp > in_sum:
            for e in in_edges:
                graph.ep.flow[e] = round(curr_dp * (graph.ep.flow[e] / in_sum), 2)
        
        if out_edges != [] and curr_dp > out_sum:
            for e in out_edges:
                graph.ep.flow[e] = round(curr_dp * (graph.ep.flow[e] / out_sum), 2)

    for edge in simp_edge_dict.values():
        print_edge(graph, edge, "After refine")

    # re assign the node depth = max(curr_dp, out_neighbor_dp_sum if no branch, in_neighbor_dp_sum if no branch, out edge flow sum, in edge flow sum)
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        curr_node = queue.pop()
        visited.append(graph.vp.id[curr_node])
        curr_dp = graph.vp.dp[curr_node]

        us = list(curr_node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        in_neighbor_dp_sum = -1
        if u_out_degrees == curr_node.in_degree():
            in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

        vs = list(curr_node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        out_neighbor_dp_sum = -1
        if v_in_degrees == curr_node.out_degree():
            out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
        graph.vp.dp[curr_node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[curr_node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))

        in_d = len([n for n in curr_node.in_neighbors()])
        node_in_edges = [e for e in curr_node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[src])
                    else:
                        graph.ep.flow[e] = w

        out_d = len([n for n in curr_node.out_neighbors()])
        node_out_edges = [e for e in curr_node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[tgt] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[tgt])
                    else:
                        graph.ep.flow[e] = w

        for n in curr_node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)
    print("Next")
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        curr_node = queue.pop()
        visited.append(graph.vp.id[curr_node])
        curr_dp = graph.vp.dp[curr_node]

        us = list(curr_node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        in_neighbor_dp_sum = -1
        if u_out_degrees == curr_node.in_degree():
            in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

        vs = list(curr_node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        out_neighbor_dp_sum = -1
        if v_in_degrees == curr_node.out_degree():
            out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
        graph.vp.dp[curr_node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[curr_node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))

        in_d = len([n for n in curr_node.in_neighbors()])
        node_in_edges = [e for e in curr_node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[src])
                    else:
                        graph.ep.flow[e] = w

        out_d = len([n for n in curr_node.out_neighbors()])
        node_out_edges = [e for e in curr_node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[tgt] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[tgt])
                    else:
                        graph.ep.flow[e] = w

        for n in curr_node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)


def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, min_cov, cutoff=125):
    """
    rebalance all the node coverage to ensure flow consistency
    """
    # helper functions
    def node_expansion_both(graph: Graph, con_nodes, incon_nodes, curr_node):
        """
        from given node, BFS search all the in/out branches, if the curr node is covered by c node layer
        then return true, if any inconsistent node be found then return false
        """
        visited = []
        queue = [curr_node]
        while queue:
            node = queue.pop()
            id = graph.vp.id[node]
            visited.append(id)
            if id in con_nodes:
                #skip
                None
            elif id in incon_nodes:
                # print("From src node {0} reach in consistent node {1} return false".format(graph.vp.id[curr_node], id))
                return False
            else:
                for u in node.all_neighbors():
                    if graph.vp.id[u] not in visited:
                        queue.append(u)
        return True
    
    def subfix(graph: Graph, curr_node, no_changes_in, consistent_node: dict, inconsistent_node: dict, simp_node_dict: dict):
        no_changes = no_changes_in
        for edge in curr_node.all_edges():
            u = graph.vp.id[edge.source()]
            v = graph.vp.id[edge.target()]

            if u in inconsistent_node and v in inconsistent_node:
                None
            elif v in inconsistent_node:
                # u is consistent
                u_node = simp_node_dict[u]
                all_incon_dp = 0
                incon_edges = []
                out_flow_remain = graph.vp.dp[u_node]
                for out_edge in u_node.out_edges():
                    tgt = out_edge.target()
                    if graph.ep.flow[out_edge] != 0.0:
                        out_flow_remain -= graph.ep.flow[out_edge]
                    else:
                        all_incon_dp += graph.vp.dp[tgt]
                        incon_edges.append(out_edge)
                for incon_edge in incon_edges:
                    tgt = incon_edge.target()
                    maxflow =  max(graph.ep.flow[incon_edge], round((graph.vp.dp[tgt]/all_incon_dp) * out_flow_remain, 2))
                    if maxflow != graph.ep.flow[incon_edge]:
                        graph.ep.flow[incon_edge] = maxflow
                        no_changes = False
            elif u in inconsistent_node:
                # v is consistent
                v_node = simp_node_dict[v]
                all_incon_dp = 0
                incon_edges = []
                in_flow_remain = graph.vp.dp[v_node]
                for in_edge in v_node.in_edges():
                    src = in_edge.source()
                    if graph.ep.flow[in_edge] != 0.0:
                        in_flow_remain -= graph.ep.flow[in_edge]
                    else:
                        all_incon_dp += graph.vp.dp[src]
                        incon_edges.append(in_edge)
                for incon_edge in incon_edges:
                    src = incon_edge.source()
                    maxflow = numpy.max([graph.ep.flow[incon_edge], round((graph.vp.dp[src]/all_incon_dp) * in_flow_remain, 2)])
                    if maxflow != graph.ep.flow[incon_edge]:
                        graph.ep.flow[incon_edge] = maxflow
                        no_changes = False     
            else:
                None

        if curr_id in consistent_node:
            w = graph.vp.dp[curr_node]
            node_in_edges = list(curr_node.in_edges())
            if len(node_in_edges) < 1:
                None
            elif len(node_in_edges) == 1:
                e = node_in_edges[0]
                src = e.source()
                if src.out_degree() == 1:
                    if graph.vp.id[src] in consistent_node:
                        maxflow = numpy.max([graph.ep.flow[e], w, graph.vp.dp[src]])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
                    else:
                        maxflow = numpy.max([graph.ep.flow[e], w])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
            else:
                # mutliple in edges
                total_dp = numpy.sum([graph.vp.dp[u] for u in curr_node.in_neighbors()])
                for e in node_in_edges:
                    src = e.source()
                    maxflow = numpy.max([graph.ep.flow[e], round((graph.vp.dp[src]/total_dp) * w)])
                    if graph.ep.flow[e] != maxflow:
                        graph.ep.flow[e] = maxflow
                        no_changes = False                            

            node_out_edges = list(curr_node.out_edges())
            if len(node_out_edges) < 1:
                None
            elif len(node_out_edges) == 1:
                e = node_out_edges[0]
                tgt = e.target()
                if tgt.in_degree() == 1:
                    if graph.vp.id[tgt] in consistent_node:
                        maxflow = numpy.max([graph.ep.flow[e], w, graph.vp.dp[tgt]])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
                    else:
                        maxflow = numpy.max([graph.ep.flow[e], w])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
            else:
                # mutliple in edges
                total_dp = numpy.sum([graph.vp.dp[u] for u in curr_node.out_neighbors()])
                for e in node_out_edges:
                    tgt = e.target()
                    maxflow = numpy.max([graph.ep.flow[e], round((graph.vp.dp[tgt]/total_dp) * w)])
                    if graph.ep.flow[e] != maxflow:
                        graph.ep.flow[e] = maxflow
                        no_changes = False
        return no_changes
    #### Helper functions
    print("-------------------------coverage rebalance----------------------")
    print("Total node: ", len(simp_node_dict))
    consistent_count = 0

    consistent_node = {}
    inconsistent_node = {}
    untracked_node = {}

    # select the most confident nodes
    for no, node in simp_node_dict.items():
        dp = graph.vp.dp[node]
        untracked = False

        us = list(node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        out_consistent = False
        if u_out_degrees == 0:
            out_consistent = False
        else:
            if u_out_degrees == node.in_degree():
                u_sum = numpy.sum([graph.vp.dp[u] for u in us])
                out_consistent = abs(u_sum - dp) < cutoff
            else:
                out_consistent = True
                untracked = True

        vs = list(node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        in_consistent = False
        if v_in_degrees == 0:
            in_consistent = False
        else:
            if v_in_degrees == node.out_degree():
                v_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                in_consistent = abs(v_sum - dp) < cutoff
            else:
                in_consistent = True
                untracked = True

        if in_consistent and out_consistent:
            if untracked:
                untracked_node[graph.vp.id[node]] = node
            else:
                consistent_node[graph.vp.id[node]] = node
                consistent_count += 1
        else:
            inconsistent_node[graph.vp.id[node]] = node
    
    print("Total consistent node after most confident extraction: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node.keys()])
    print([int(n) for n in inconsistent_node.keys()]) 
    print([int(n) for n in untracked_node.keys()])

    # deal with untracked node
    # from current untracked node, BFS to both end, for every branch, stop until non-untracked
    # node be found, if found a in-consistent node during the search, then keep it untracked/mark it
    # inconsistent, otherwise mark it consistent.
    for no, node in list(untracked_node.items()):
        untracked_node.pop(no)
        if node_expansion_both(graph, consistent_node, inconsistent_node, node):
            dp = graph.vp.dp[node]

            us = list(node.in_neighbors())
            u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
            out_consistent = False
            out_amb = False
            if u_out_degrees == 0:
                out_consistent = False
            else:
                if u_out_degrees == node.in_degree():
                    u_sum = numpy.sum([graph.vp.dp[u] for u in us])
                    out_consistent = abs(u_sum - dp) < cutoff
                else:
                    out_amb = True

            vs = list(node.out_neighbors())
            v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
            in_consistent = False
            in_amb = False
            if v_in_degrees == 0:
                in_consistent = False
            else:
                if v_in_degrees == node.out_degree():
                    v_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                    in_consistent = abs(v_sum - dp) < cutoff
                else:
                    in_amb = True
            
            if out_amb and in_amb:
                # no way to determine its consistency
                inconsistent_node[no] = node
            else:
                if in_consistent and out_consistent:
                    consistent_node[no] = node
                    consistent_count += 1
                elif in_consistent and out_amb:
                    consistent_node[no] = node
                    consistent_count += 1
                elif out_consistent and in_amb:
                    consistent_node[no] = node
                    consistent_count += 1
                else:
                    inconsistent_node[no] = node
        else:
            inconsistent_node[no] = node
    print("Total consistent node after untracked classification: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node])
    print([int(n) for n in inconsistent_node]) 
    print([int(n) for n in untracked_node])

    # assign edge flow only to consistent node related edge
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict, consistent_node, inconsistent_node)
    
    sorted_consistent_node = sorted(consistent_node.items(), key=lambda x: graph.vp.dp[x[1]], reverse=True)
    # fix the inconsistent node depth, and also try to 
    # increment all the imbalance consistent node
    i = 0
    break_point = consistent_count
    no_changes = False
    cyclic = not is_DAG(graph)

    while len(inconsistent_node) != 0 or not no_changes:
        no_changes = True
        no, node = sorted_consistent_node[i]
        print("Current iteration: ", i, "Total consistent node: ", consistent_count)
        print([int(n) for n in inconsistent_node])
        print("Current solid node: ", no)

        queue = [node]
        visited = []
        while queue:
            curr_node = queue.pop()
            curr_id = graph.vp.id[curr_node]
            curr_dp = graph.vp.dp[curr_node]
            visited.append(curr_id)
            # fix the curr node if is in consistent
            if curr_id in inconsistent_node:
                fix = True
                vtotal = 0
                for e in curr_node.in_edges():
                    if graph.ep.flow[e] != 0.0:
                        vtotal += graph.ep.flow[e]
                    else:
                        # there exist an unassigned edge around the curr node
                        fix = False
                        break
                # perform fix
                if fix:
                    graph.vp.dp[curr_node] = vtotal
                    if DEBUG_MODE:
                        print(curr_id, "inconsistent node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))
                    inconsistent_node.pop(curr_id)
                    consistent_node[curr_id] = curr_node
                    consistent_count += 1
                    no_changes = False
            else:
                us = list(curr_node.in_neighbors())
                u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
                in_neighbor_dp_sum = -1
                if u_out_degrees == curr_node.in_degree():
                    in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

                vs = list(curr_node.out_neighbors())
                v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
                out_neighbor_dp_sum = -1
                if v_in_degrees == curr_node.out_degree():
                    out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                
                in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
                out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
                maxdp = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
                if curr_dp != maxdp:
                    graph.vp.dp[curr_node] = maxdp
                    no_changes = False
            
            no_changes = subfix(graph, curr_node, no_changes, consistent_node, inconsistent_node, simp_node_dict)
            # append all the out neighbors
            for n in curr_node.out_neighbors():
                if graph.vp.id[n] not in visited:
                    queue.append(n)

        queue = [node]
        visited = []
        while queue:
            curr_node = queue.pop()
            curr_id = graph.vp.id[curr_node]
            curr_dp = graph.vp.dp[curr_node]
            visited.append(curr_id)
            # fix the curr node if is in consistent
            if curr_id in inconsistent_node:
                fix = True
                vtotal = 0
                for e in curr_node.out_edges():
                    if graph.ep.flow[e] != 0.0:
                        vtotal += graph.ep.flow[e]
                    else:
                        # there exist an unassigned edge around the curr node
                        fix = False
                        break
                # perform fix
                if fix:
                    graph.vp.dp[curr_node] = vtotal
                    if DEBUG_MODE:
                        print(curr_id, "inconsistent node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))
                    inconsistent_node.pop(curr_id)
                    consistent_node[curr_id] = curr_node
                    consistent_count += 1
                    no_changes = False
            else:
                us = list(curr_node.in_neighbors())
                u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
                in_neighbor_dp_sum = -1
                if u_out_degrees == curr_node.in_degree():
                    in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

                vs = list(curr_node.out_neighbors())
                v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
                out_neighbor_dp_sum = -1
                if v_in_degrees == curr_node.out_degree():
                    out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                
                in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
                out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
                maxdp = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
                if curr_dp != maxdp:
                    graph.vp.dp[curr_node] = maxdp
                    no_changes = False
            
            no_changes = subfix(graph, curr_node, no_changes, consistent_node, inconsistent_node, simp_node_dict)
            # append all the in neighbors
            for n in curr_node.in_neighbors():
                if graph.vp.id[n] not in visited:
                    queue.append(n)        
        i += 1
        if i >= len(sorted_consistent_node):
            if break_point == consistent_count and (no_changes or cyclic):
                # no more changes
                break
            else:
                i = 0
                break_point = consistent_count

    print("Total consistent node after solid fix: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node])
    print([int(n) for n in inconsistent_node]) 
    print([int(n) for n in untracked_node])
    if DEBUG_MODE:
        for edge in simp_edge_dict.values():
            print_edge(graph, edge, "")
    return

######


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




# def get_concat_plan(contig_dict: dict, max_len):
#     contig_impossible_dict = {}
#     all_contig_ids = contig_dict.keys()
#     for no in contig_dict.keys():
#         contig_impossible_dict[no] = set()
#     for tail_cno, [tail_contig, tail_clen, _] in contig_dict.items():
#         for head_cno, [head_contig, head_clen, _] in contig_dict.items():
#             if tail_cno != head_cno:
#                 if tail_clen + head_clen > max_len:
#                     contig_impossible_dict[tail_cno].add(head_cno)
#                 if list(set(tail_contig) & set(head_contig)) != []:
#                     contig_impossible_dict[tail_cno].add(head_cno)
#     contig_concat_plans = {}
#     for key, item in contig_impossible_dict.items():
#         ps = all_contig_ids - item
#         print("cno: ", key, " can concat with following: ", ps)
#         contig_concat_plans[key] = ps
#     return contig_concat_plans



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