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