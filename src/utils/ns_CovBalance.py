#!/usr/bin/env python3

from graph_tool.all import Graph
import numpy

from utils.ns_Utilities import *
from utils.ns_IO import graph_to_gfa, flipped_gfa_to_graph

def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    for (u,v), e in simp_edge_dict.items():

        if graph.ep.color[e] != 'black':
            #forbidden edge
            continue

        u_node = simp_node_dict[u]
        u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])

        v_node = simp_node_dict[v]
        v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])

        graph.ep.flow[e] = numpy.mean([(graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], 
            (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node]])
    return

def coverage_rebalance_s(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    
    # break the cycle first
    removed_edges = []
    if not graph_is_DAG(graph, simp_node_dict):
        removed_edges = cyclic_to_dag(graph, simp_node_dict, simp_edge_dict)
    # incident node insertion
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict)
    print("add incident nodes")
    incident_vs = []
    for edge in set(graph.edges()):
        u = edge.source()
        uid = graph.vp.id[u]
        v = edge.target()
        vid = graph.vp.id[v]
        node = graph_add_vertex(graph, simp_node_dict, uid+"*"+vid, graph.ep.flow[edge], "*")
        incident_vs.append((u, node, v, graph.ep.overlap[edge]))
        graph_remove_edge(graph, simp_edge_dict, uid, vid)
        graph.remove_edge(edge)
        graph_add_edge(graph, simp_edge_dict, u, uid, node, graph.vp.id[node], 0, 0)
        graph_add_edge(graph, simp_edge_dict, node, graph.vp.id[node], v, vid, 0, 0)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth has been stored in the dict
    # ratio: normalised balanced node depth / previous node depth
    coverage_rebalance_ave(graph, simp_node_dict, simp_edge_dict)

    print("remove incident nodes")
    for (src, node, tgt, overlap) in sorted(incident_vs, key=lambda t: t[1], reverse=True):
        for edge in set(node.all_edges()):
            u = edge.source()
            uid = graph.vp.id[u]
            v = edge.target()
            vid = graph.vp.id[v]
            graph_remove_edge(graph, simp_edge_dict, uid, vid)
            graph.remove_edge(edge)
        graph_add_edge(graph, simp_edge_dict, src, graph.vp.id[src], tgt, graph.vp.id[tgt], overlap, graph.vp.dp[node])
        graph_remove_vertex(graph, simp_node_dict, graph.vp.id[node])
        graph.remove_vertex(node)

    for (uid,vid,o) in removed_edges:
        graph_add_edge(graph, simp_edge_dict, simp_node_dict[uid], uid, simp_node_dict[vid], vid, o)
    print("done")
    return

def coverage_rebalance_ave(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    print("coverage rebalance..")
    # set cutoff delta
    cutoff = 0.00001 * len(simp_node_dict)
    print("cutoff coverage-unbalance rate: ", cutoff)
    # store previous node depth
    prev_dp_dict = {}
    for no, v in simp_node_dict.items():
        prev_dp_dict[no] = graph.vp.dp[v]

    # EM optimization
    sum_delta = 0
    sum_depth_before = numpy.sum([dp for dp in prev_dp_dict.values()])
    while True:
        # ****************************-E Step
        for (u,v),e in simp_edge_dict.items():
            if graph.ep.color[e] != 'black':
                #forbidden edge
                continue

            u_node = simp_node_dict[u]
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])

            v_node = simp_node_dict[v]
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])

            graph.ep.flow[e] = numpy.mean([(graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], 
                (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node]])
        # ****************************-Validation
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
                deltas.append((no, (graph.vp.dp[node] * abs(inflow - outflow))/dominator))

        sum_delta = sum([k[1] for k in deltas]) / sum([graph.vp.dp[node] for node in graph.vertices()])
        # print("sum delta: ", sum_delta, "worst delta: ", sorted(deltas, key=lambda p: p[1], reverse=True)[:10])
        if sum_delta < cutoff:
            print("final delta: ", sum_delta)
            break
        # rescal 
        sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
        if sum_ratio > 2 or sum_ratio < 0.5:
            print("ratio doubled: ", sum_ratio, "rescalling..")
            for node in graph.vertices():
                if node.in_degree() > 0 or node.out_degree() > 0:
                    graph.vp.dp[node] = graph.vp.dp[node]/sum_ratio
            for edge in graph.edges():
                graph.ep.flow[edge] = graph.ep.flow[edge]/sum_ratio
        # ****************************-M Step
        for no, node in simp_node_dict.items():
            inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
            outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
            if node.in_degree() == 0 and node.out_degree() == 0:
                graph.vp.dp[node] = graph.vp.dp[node]
            elif node.in_degree() == 0 and node.out_degree() != 0:
                graph.vp.dp[node] = outflow
            elif node.in_degree() != 0 and node.out_degree() == 0:
                graph.vp.dp[node] = inflow
            else:
                graph.vp.dp[node] = numpy.mean([inflow, outflow])
    # final evaluation
    sum_ratio = (sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)

    print("scaled ratio, normalizing: ", sum_ratio)
    for node in simp_node_dict.values():
        # only scale the connected nodes
        if node.in_degree() > 0 or node.out_degree() > 0:
            graph.vp.dp[node] = graph.vp.dp[node] / sum_ratio
    for edge in graph.edges():
        graph.ep.flow[edge] = graph.ep.flow[edge] / sum_ratio
    print("done")
    return