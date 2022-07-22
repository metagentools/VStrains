#!/usr/bin/env python3
from utils.ns_IO import flipped_gfa_to_graph, gfa_to_graph, graph_to_gfa
import numpy
# from utils.ns_Utilities import graph_add_edge, graph_add_vertex, graph_remove_edge, graph_remove_vertex, print_edge, print_vertex
# from utils.ns_Utilities import graph_is_DAG
import sys
from graph_tool import Graph
from graph_tool.topology import all_circuits

if __name__ == '__main__':
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph("graph_c.gfa")
    # print(graph_is_DAG(graph, simp_node_dict))
    # circ = all_circuits(graph, True)
    # print(circ)
    graph = Graph(directed=True)
    u = graph.add_vertex()
    v = graph.add_vertex()
    print(graph.edge(u, v))