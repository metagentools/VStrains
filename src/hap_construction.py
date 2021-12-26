#!/usr/bin/env python3

import re
import sys, os
import json
from gfapy import gfa
from gfapy.gfa import Gfa
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort
import subprocess
import argparse

import gfapy
from numpy import sin

usage = "Construct haplotypes using divide-and-conquer method"
author = "Runpeng Luo"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')

    ## TODO may add gfa validation
    args = parser.parse_args()
    if not args.gfa_file:
        print("gfa file is not imported")
        return 1

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa1').from_file(filename=args.gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))
    
    graph, node_dict, edge_dict = gfa_to_graph(gfa)
    refine_graph(graph, node_dict, edge_dict, args.forward, args.reverse)


def gfa_to_graph(gfa: gfapy.Gfa):
    """
    Convert gfa to graph
    Nodes: segment with corresponding 
    """
    graph = Graph(directed=True)

    vprop_seq = graph.new_vertex_property("string", val="")
    vprop_abd = graph.new_vertex_property("double")
    vprop_kc = graph.new_vertex_property("int32_t")
    vprop_id = graph.new_vertex_property("string", val="UD")
    vprop_visited = graph.new_vertex_property("int16_t", val=0)
    vprop_ori = graph.new_vertex_property("int16_t") # 1 = +, -1 = -
    vprop_group = graph.new_vertex_property("int16_t", val=-1)

    graph.vp.seq = vprop_seq
    graph.vp.abd = vprop_abd
    graph.vp.kc = vprop_kc
    graph.vp.id = vprop_id
    graph.vp.visited = vprop_visited
    graph.vp.ori = vprop_ori
    graph.vp.group = vprop_group

    eprop_overlap = graph.new_edge_property("int", val=0)
    eprop_visited = graph.new_edge_property("int", val=0)

    graph.ep.overlap = eprop_overlap
    graph.ep.visited = eprop_visited

    # source = graph.add_vertex()
    # sink   = graph.add_vertex()
    # graph.vp.id[source] = "source"
    # graph.vp.id[sink] = "sink"

    # graph.list_properties()
    # S
    node_dict = {}
    edge_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [line_type, seg_no, seg, dp, kc] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        kc_float = float(kc.split(":")[2])
        v_pos = graph.add_vertex()
        v_neg = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.abd[v_pos] = dp_float
        graph.vp.kc[v_pos] = kc_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.group[v_pos] = -1
        
        graph.vp.seq[v_neg] = seg
        graph.vp.abd[v_neg] = dp_float
        graph.vp.kc[v_neg] = kc_float
        graph.vp.id[v_neg] = seg_no
        graph.vp.ori[v_neg] = -1
        graph.vp.group[v_neg] = -1

        node_dict[seg_no] = (v_pos, v_neg)
    # L
    for edge in gfa.edges:
        [line_type, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
        u_pos, u_neg = node_dict[seg_no_l]
        v_pos, v_neg = node_dict[seg_no_r]
        u = u_pos if ori_l == '+' else u_neg
        v = v_pos if ori_r == '+' else v_neg
        e = graph.add_edge(source=u, target=v)
        # gfa format check
        assert overlap_len[-1] == 'M'
        graph.ep.overlap[e] = int(overlap_len[:-1])
        edge_dict[(u,v)] = e
    
    # P
    for path in gfa.paths:
        [line_type, path_no, seg_names, seg_overlap] = str(path).split("\t")
        #TODO


    return graph, node_dict, edge_dict


def refine_graph(graph: Graph, node_dict: dict, edge_dict: dict, forward, reverse):
    """
    Maximimize graph connectivity by minimizing node with 0 in-degree or out-degree, detect and remove all the cycles.
    """
    # determine the isolated subgraphs, and assign each node with its group No, which indicates they are belong to same group
    # assign non-connected node with group No 0.

    def bfs_grouping(graph: Graph, start_node, group_no, groups):
        """
        Perform a breadth-first search and assign the group no to all the connected nodes.
        """
        group_no = group_no + 1
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

    group_no = 0
    groups = {}
    for (v_pos, v_neg) in node_dict.values():
        if graph.vp.group[v_pos] == -1:
            graph, group_no, groups = bfs_grouping(graph, v_pos, group_no, groups)
        if graph.vp.group[v_neg] == -1:
            graph, group_no, groups = bfs_grouping(graph, v_neg, group_no, groups)

    print("groups:")
    for key, item in groups.items():
        print("group number: ", key, " member: ", [(graph.vp.id[v], '+' if graph.vp.ori[v] == 1 else '-') for v in item])

    # connect sub-graphs based on pair-end reads information
    
    

    # check if the graph is cyclic, break cycle if true, then assign edge from s&t node to nodes.
    # DFS search, walk through all paths, for all cycle paths, break the cycle by removing the lowest aboundance nodes.

    # start from source node, run BFS to determine the level for each nodes until reaching the sink node.
    # level is determined based on seg-length, 
    return {}


if __name__ == "__main__":
    main()
