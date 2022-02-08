#!/usr/bin/env python3

# import re
# import sys, os
# import json
# import re
import subprocess
from tkinter.messagebox import NO
from turtle import ht
from graph_tool import GraphView, _in_degree
# import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw
# from graph_tool.clustering import local_clustering

import argparse

import gfapy
import numpy

from graph_converter import gfa_to_graph, graph_to_gfa, map_ref_to_graph, get_contig, contig_to_seq, simplify_edge_dict

usage = "Construct haplotypes using divide-and-conquer method"
author = "Runpeng Luo"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, default=500, help=("minimum coverage for strains"))
    parser.add_argument('-minlen', '--minimum_strain_length', dest='min_len', default=8000, type=int, help=("minimum strain length"))
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

    graph, node_dict, edge_dict, dp_dict = gfa_to_graph(args.gfa_file)

    graph, simp_node_dict, simp_edge_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict.copy(), 1)

    assign_edge_flow(graph, simp_node_dict, simp_edge_dict)

    contig_dict = get_contig(graph, args.contig_file, simp_node_dict, simp_edge_dict, args.min_cov)

    graph_stat(graph, simp_node_dict, simp_edge_dict)

    # graph_draw(graph, vprops={'text': graph.vp.id}, eprops={'text': graph.ep.flow}, output="graph.pdf", output_size=(2000,2000))
    if args.ref_file:
        strain_dict = map_ref_to_graph(args.ref_file, graph, simp_node_dict)

    graph_reduction(graph, contig_dict, simp_node_dict, simp_edge_dict, args.min_cov, args.min_len)

    contig_classification(contig_dict, args.min_cov, args.min_len)

def flip_graph_bfs(graph: Graph, node_dict: dict, edge_dict: dict, dp_dict: dict, init_ori):
    """
    Flip all the node orientation.

    return an node_dict, which only contains one orientation per node for simplicity.
    rename all the used node to positive, and forbidden the opponent node.
    """
    print("-------------------------flip graph orientation----------------------")
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

    simp_edge_dict = simplify_edge_dict(edge_dict)
    print("-------------------------flip graph orientation end------------------")
    return graph, simp_node_dict, simp_edge_dict

def longest_distance(graph: Graph, s):
    return

def graph_stat(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    print("-------------------------graph stat----------------------")
    for seg_no, v in simp_node_dict.items():
        print_vertex(graph, v, "stat")
    for (_,_), e in simp_edge_dict.items():
        print_edge(graph, e, "stat")
    
    print("-----------------------graph stat end--------------------")


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
    
    edge_dict.pop((graph.vp.id[tmp_s], graph.vp.ori[tmp_s], graph.vp.id[tmp_t], graph.vp.ori[tmp_t]))

    tmp_s_pos, tmp_s_neg = node_dict[graph.vp.id[tmp_s]]
    tmp_t_pos, tmp_t_neg = node_dict[graph.vp.id[tmp_t]]
    s = tmp_t_pos if graph.vp.ori[tmp_t] == -1 else tmp_t_neg
    t = tmp_s_pos if graph.vp.ori[tmp_s] == -1 else tmp_s_neg

    o = graph.ep.overlap[edge]
    graph.remove_edge(edge)
    e = graph.add_edge(s, t)
    graph.ep.overlap[e] = o
    edge_dict[(graph.vp.id[s], graph.vp.ori[s], graph.vp.id[t], graph.vp.ori[t])] = e

    return graph, e, edge_dict

def print_edge(graph, e, s=""):
    print(s, " edge: ", graph.vp.id[e.source()],graph.vp.ori[e.source()], "->", graph.vp.id[e.target()], graph.vp.ori[e.target()], graph.ep.flow[e], graph.ep.color[e])

def print_vertex(graph, v, s=""):
    print(s, " vertex: ", graph.vp.id[v], ", dp: ", graph.vp.dp[v], ", ori: ", graph.vp.ori[v], ", in_degree: ", v.in_degree(), ", out_degree: ", v.out_degree(), graph.vp.color[v])

def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total edges: ", un_assigned_edge)
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
                print("manual assigned edge: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                un_assigned_edge = un_assigned_edge - 1
                graph.ep.flow[e] = assign_flow

    print("un-assigned edges after manual assign iteration : ", un_assigned_edge)
    print("-----------------------assign edge flow end--------------------")

def graph_reduction(graph: Graph, contig_dict: dict, simp_node_dict: dict, simp_edge_dict: dict, min_cov, min_len):
    """
    reduce the node/edge weight based on existing contig found by SPAdes.
    """

    print("-------------------------graph reduction----------------------")
    def contig_reduction(graph: Graph, contig, ccov, clen, simp_node_dict: dict, simp_edge_dict: dict, min_cov, min_len):
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
                    print_edge(graph, e, "edge been removed due to either u or v is removed")
                else:
                    print("edge: ", node, " -> ", adj_node, "already been removed")
                continue
            

            print_edge(graph, e, "current edge eval")

            # reduce the depth for involving edge
            if graph.ep.flow[e] - ccov <= min_cov:
                graph.ep.flow[e] = 0
                graph.ep.color[e] = 'gray'
                simp_edge_dict.pop((node,adj_node))
                print_edge(graph, e, "edge been removed")
            else:
                graph.ep.flow[e] = graph.ep.flow[e] - ccov

            # reduce the depth for involving node, gray color for forbiddened node
            if graph.vp.dp[u] - ccov <= min_cov:
                graph.vp.dp[u] = 0
                graph.vp.color[u] = 'gray'
                simp_node_dict.pop(node)
                print("node ", node, "been removed")
            else:
                graph.vp.dp[u] = graph.vp.dp[u] - ccov

            # last node in the contig, reduce its depth
            if adj_index + 1 == len(contig):
                if graph.vp.dp[v] - ccov <= min_cov: 
                    graph.vp.dp[v] = 0
                    graph.vp.color[v] = 'gray'
                    simp_node_dict.pop(adj_node)
                    print("node ", adj_node, "been removed")
                else:
                    graph.vp.dp[v] = graph.vp.dp[v] - ccov
            
            # update edges
            if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (node,adj_node) in simp_edge_dict:
                graph.ep.flow[e] = 0
                graph.ep.color[e] = 'gray'
                simp_edge_dict.pop((node,adj_node))
                print_edge(graph, e, "edge been removed in the final step")
            adj_index = adj_index + 1
        return
    
    for (cno, clen, ccov), (contig,_) in contig_dict.items():
        if clen >= min_len:
            contig_reduction(graph, contig, ccov, clen, simp_node_dict, simp_edge_dict, min_cov, min_len)    
    # store level 2 graph
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, min_cov, "acc/graph_L2.gfa")

    for (cno, clen, ccov), (contig,_) in contig_dict.items():
        if clen < min_len:
            # keep the head&tail of the contig
            # contig = contig[1:-1]
            contig_reduction(graph, contig, ccov, clen, simp_node_dict, simp_edge_dict, min_cov, min_len)
    # store level 3 graph
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, min_cov, "acc/graph_L3.gfa")

    # for (cno, clen, ccov), (contig,_) in contig_dict.items():
    #     # if clen < min_len and cno == "5":
    #         # keep the head&tail of the contig
    #         contig_head = contig[0:2]
    #         contig_reduction(graph, contig_head, ccov, clen, simp_node_dict, simp_edge_dict, min_cov, min_len)
    #         contig_tail = contig[-2:]
    #         contig_reduction(graph, contig_tail, ccov, clen, simp_node_dict, simp_edge_dict, min_cov, min_len)
    # # store level 4 graph
    # graph_to_gfa(graph, simp_node_dict, simp_edge_dict, min_cov, "acc/graph_L4.gfa")
    print("-----------------------graph reduction end--------------------")
    return 0

def contig_classification(contig_dict: dict, min_cov, min_len):

    print("-----------------------contig classifcation--------------------")

    graph_L2, node_dict_L2, edge_dict_L2, dp_dict_L2  = gfa_to_graph("acc/graph_L2.gfa")
    graph_L2, simp_node_dict_L2, simp_edge_dict_L2 = flip_graph_bfs(graph_L2, node_dict_L2, edge_dict_L2, dp_dict_L2, 1)
    
    graph_L3, node_dict_L3, edge_dict_L3, dp_dict_L3 = gfa_to_graph("acc/graph_L3.gfa")
    graph_L3, simp_node_dict_L3, simp_edge_dict_L3 = flip_graph_bfs(graph_L3, node_dict_L3, edge_dict_L3, dp_dict_L3, 1)

    cand_strains = {}
    temp_contigs = {}
    for (cno, clen, ccov), (contig, flow) in contig_dict.items():
        if clen > min_len:
            cand_strains[(cno, clen, ccov)] = (contig, flow)
            print("full-length contig found: ", cno, clen, ccov)
        else:
            temp_contigs[(cno, clen, ccov)] = (contig, flow)
    h_list = []
    t_list = []
    for (cno, clen, ccov), (contig, flow) in temp_contigs.items():
        head = contig[0]
        tail = contig[-1]
        h_list.append(int(head))
        t_list.append(int(tail))
        print("contig: ", cno, clen, ccov, " head: ", head, " tail: ", tail)
    print("head list: ", h_list)
    print("tail list: ", t_list)


    # find pair-wise distance among contigs
    l = len(h_list)
    dist_list = [[None for i in range(l)] for j in range(l)]
    for i in range(l):
        for j in range(l):
            if i == j:
                continue
            dist_list[i][j] = (h_list[i], t_list[j])
            # compte distance between ith head to tth tail


    # dist_matrix = numpy.array(dist_list, dtype=object)
    # print(dist_matrix)
    print(dist_list)
    print("--------------------contig classification end--------------------")
    return

# def build_residue_graph(graph: Graph, simp_node_dict: dict, edge_dict: dict):
#     for u in graph.vertices():

#     for e in graph.edges():
        

if __name__ == "__main__":
    main()
