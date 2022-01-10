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
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')

    ## TODO may add gfa validation
    args = parser.parse_args()
    if not args.gfa_file:
        print("gfa file is not imported")
        return 1

    # sys.setrecursionlimit(1000000000)

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa1').from_file(filename=args.gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))
    
    graph, node_dict, edge_dict, dp_dict = gfa_to_graph(gfa)
    # for node in graph.vertices():
    #     print_vertex(graph, node, "node summary: ")
    graph, node_dict, edge_dict, pick_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict)
    graph_to_gfa(graph, node_dict, edge_dict, pick_dict)
    # graph, node_dict, edge_dict, groups = refine_graph(graph, node_dict, edge_dict, args.forward, args.reverse, 10)


def gfa_to_graph(gfa: gfapy.Gfa):
    """
    Convert gfa to graph
    Nodes: segment with corresponding 
    """
    graph = Graph(directed=True)

    vprop_seq = graph.new_vertex_property("string", val="")
    vprop_dp = graph.new_vertex_property("double")
    vprop_kc = graph.new_vertex_property("int32_t")
    vprop_id = graph.new_vertex_property("string", val="UD")
    vprop_visited = graph.new_vertex_property("int16_t", val=0)
    vprop_ori = graph.new_vertex_property("int16_t") # 1 = +, -1 = -
    vprop_group = graph.new_vertex_property("int16_t", val=-1)

    graph.vp.seq = vprop_seq
    graph.vp.dp = vprop_dp
    graph.vp.kc = vprop_kc
    graph.vp.id = vprop_id
    graph.vp.visited = vprop_visited
    graph.vp.ori = vprop_ori
    graph.vp.group = vprop_group

    eprop_overlap = graph.new_edge_property("int", val=0)
    eprop_visited = graph.new_edge_property("int", val=0)
    eprop_flow = graph.new_edge_property("float", val=0.0)

    graph.ep.overlap = eprop_overlap
    graph.ep.visited = eprop_visited
    graph.ep.flow = eprop_flow

    # graph.list_properties()
    # S
    node_dict = {}
    dp_dict = {}
    edge_dict = {}
    # matrix_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [line_type, seg_no, seg, dp, kc] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        kc_float = float(kc.split(":")[2])
        v_pos = graph.add_vertex()
        v_neg = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.dp[v_pos] = dp_float
        graph.vp.kc[v_pos] = kc_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.group[v_pos] = -1
        graph.vp.visited[v_pos] = -1
        
        graph.vp.seq[v_neg] = reverse_seq(seg)
        graph.vp.dp[v_neg] = dp_float
        graph.vp.kc[v_neg] = kc_float
        graph.vp.id[v_neg] = seg_no
        graph.vp.ori[v_neg] = -1
        graph.vp.group[v_neg] = -1
        graph.vp.visited[v_neg] = -1

        node_dict[seg_no] = (v_pos, v_neg)
        dp_dict[seg_no] = dp_float
        # matrix_dict[seg_no] = [[[],[]],
        #                         [[],[]]]
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
        graph.ep.flow[e] = min(graph.vp.dp[u], graph.vp.dp[v])
        edge_dict[(u,v)] = e
        
    # P
    for path in gfa.paths:
        [line_type, path_no, seg_names, seg_overlap] = str(path).split("\t")
        #TODO


    return graph, node_dict, edge_dict, dp_dict


def refine_graph(graph: Graph, node_dict: dict, edge_dict: dict, forward, reverse, level_length_cut_off):
    """
    Maximimize graph connectivity by minimizing node with 0 in-degree or out-degree, detect and remove all the cycles.
    """
    # determine the isolated subgraphs, and assign each node with its group No, which indicates they are belong to same group
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
    groups[0] = []
    for (v_pos, v_neg) in node_dict.values():
        # grouping
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
    # level is determined based on seg-length
    return graph, node_dict, edge_dict, groups

# def divide_and_conquer(graph, base_case_length)


def flip_graph_bfs(graph: Graph, node_dict: dict, edge_dict: dict, dp_dict: dict):
    pick_dict = {}
    while set(dp_dict):
        seg_no = source_node(dp_dict)
        source_pos, source_neg = node_dict[seg_no]
        graph.vp.visited[source_pos] = 0
        graph.vp.visited[source_neg] = 0
        queue = []
        queue.append([node_dict[seg_no], -1])

        while queue:
            (v_pos, v_neg), ori = queue.pop()
            dp_dict.pop(graph.vp.id[v_pos])
            
            u = None
            if ori == 1:
                u = v_pos
                pick_dict[graph.vp.id[u]] = 1
                # print_vertex(graph, v_neg, "node to reverse")
                for e in list(v_neg.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
                    None
            else:
                u = v_neg
                pick_dict[graph.vp.id[u]] = -1
                # print_vertex(graph, v_pos, "node to reverse")
                for e in list(v_pos.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
                    None
            
            graph.vp.visited[v_pos] = 1
            graph.vp.visited[v_neg] = 1
            # add further nodes into the queue TODO, all or out only
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1:
                    graph.vp.visited[adj_node] = 0
                    queue.append([node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]])

    for key, item in pick_dict.items():
        print("pick: ", key, item)
    # verify sorted graph
    print("-------------------------verify graph----------------------")
    check = True
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
    # print all the path
    # for key, item in pick_dict.items():
    #     print("-----------------------------------------------------------")
    #     print("longest path {0} {1} can reach:".format(key, item))
    #     path = graph_dfs(graph, node_dict[key][0 if item == 1 else 1])
    #     [print_vertex(graph, v) for v in path]
    #     print("---------------------------{0}-----------------------------".format(len(path)))
    return graph, node_dict, edge_dict, pick_dict

def is_ambiguous(graph: Graph, node):
    return True

def graph_dfs(graph: Graph, source):
    """
    Count the maximum depth the source can reach via directed edge in the graph
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
                    cmp.insert(0, v)
                    rtn = rtn if len(cmp) < len(rtn) else cmp
            return rtn

    return dfs_helper(graph, source, visited)


def source_node(dp_dict: dict):
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
    edge_dict[(tmp_s, tmp_t)] = None
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

def reverse_seq(seq: str):
    return ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[x] for x in seq[::-1])

def graph_to_gfa(graph: Graph, node_dict: dict, edge_dict: dict, pick_dict: dict):
    subprocess.check_call("rm -rf {0} && mkdir {0} && echo "" > {1}".format(
    "acc", "acc/acc_graph.gfa"), shell=True)
    
    def to_ori(i):
        return '+' if i == 1 else '-'
    with open("acc/acc_graph.gfa", 'w') as gfa:
        for seg_no, (v_pos, v_neg) in node_dict.items():
            if pick_dict[seg_no] == 1:
                pick_nodes = [v_pos]
            elif pick_dict[seg_no] == 0:
                pick_nodes = [v_pos, v_neg]
            else:
                pick_nodes = [v_neg]
            # print_vertex(graph, pick_node, "picked node")
            for n in pick_nodes:
                gfa.write("S\t{0}\t{1}\tDP:f:{2}\tKC:i:{3}\n".format
                (graph.vp.id[n], graph.vp.seq[n], graph.vp.dp[n], graph.vp.kc[n]))
        for (u, v), e in edge_dict.items():
            if e != None:
                # print_edge(graph, e, "picked edge")
                gfa.write("L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format
                (graph.vp.id[u], to_ori(graph.vp.ori[u]), 
                graph.vp.id[v], to_ori(graph.vp.ori[v]), 
                graph.ep.overlap[e]))
        gfa.close()
    return 0

if __name__ == "__main__":
    main()
