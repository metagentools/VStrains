#!/usr/bin/env python3
import logging
from utils.vsAware_Path import label_filter_dfs, bellman_ford
from utils.vsAware_Path import set_filter
from utils.vsAware_Path import set_edge_weight
from utils.vsAware_IO import flipped_gfa_to_graph, gfa_to_graph, graph_to_gfa
import numpy

import sys
from utils.vsAware_Utilities import *
from utils.vsAware_Split import *
from graph_tool.topology import label_components, shortest_path

if __name__ == "__main__":
    logger = logging.getLogger("log")
    logger.setLevel(logging.DEBUG)
    consoleHeader = logging.StreamHandler()
    consoleHeader.setLevel(logging.DEBUG)
    consoleHeader.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(consoleHeader)

    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(
        "acc_5hiv_labmix/gfa/final_graph.gfa", logger
    )
    logger.debug("graph is dag? " + str(graph_is_DAG(graph, simp_node_dict)))
    graph.ep.keep = graph.new_edge_property("boolean", True)
    contig = [
        "102*A&1395&289&288&287&286&285*A&316*A&317*A&319*A&320*A&325&270*A*A*B&327*A*B&306*A*B&305*A*B&304*A*B&303*A*B&309*B&310*B&311*B&312*B&313*B&314*B&400&241&240*B&239*B&238*B&237*0&236*B&235*B&404&345&346&1426&71&70&69&1429&353&352*A&354*A&355*A&356*A&414*A&224*A&223*A&222*A&416*A&378*A&377*A&376*A&422*A&213*A*A&212*A*A&423*A*A"
    ]
    s = None
    t = None
    if len(contig) == 1:
        # split the single node contig, zip it again later
        ns = graph.add_vertex()
        nt = simp_node_dict[contig[0]]
        graph.vp.id[ns] = str(graph.vp.id[nt]) + "temp"
        graph.vp.color[ns] = 'black'
        for e in list(nt.in_edges()):
            ie = graph.add_edge(e.source(), ns)
            graph.ep.color[ie] = 'black'
            simp_edge_dict[(graph.vp.id[e.source()], graph.vp.id[ns])] = ie
            simp_edge_dict.pop((graph.vp.id[e.source()], graph.vp.id[e.target()]))
            graph.remove_edge(e)

        s = nt
        t = ns
        contig.insert(0, graph.vp.id[ns])
        simp_node_dict[graph.vp.id[ns]] = ns
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, logger, "graph_f0.gfa")

    for e in graph.edges():
        graph.ep.keep[e] = True
    # interior edges
    for i in range(1, len(contig) - 1):
        for e in simp_node_dict[contig[i]].all_edges():
            graph.ep.keep[e] = False

    if len(contig) > 1:
        # first contig node edges
        for oe in simp_node_dict[contig[0]].out_edges():
            print_edge(graph,oe, logger, "oe")
            graph.ep.keep[oe] = False
        # last contig node edges
        for ie in simp_node_dict[contig[-1]].in_edges():
            print_edge(graph, ie, logger, "ie")
            graph.ep.keep[ie] = False
    graph.ep.eval = graph.new_edge_property("int", 0)
    # forward cycle
    label_filter_dfs(graph, contig, False)
    # reverse cycle
    label_filter_dfs(graph, contig, True)
    # graph.set_edge_filter(graph.ep.keep)

    # print(reachable(graph, s, t))
    true_s = [v for v in graph.vertices() if graph.vp.id[v] == contig[-1]][0]
    true_t = [v for v in graph.vertices() if graph.vp.id[v] == contig[0]][0]
    logger.debug("reachable s->t {0}".format(reachable(graph, true_s, true_t)))
    sp_vlist, sp_score = bellman_ford(graph, true_s, true_t, logger)
    
    for e in graph.edges():
        if not graph.ep.keep[e]:
            # print_edge(graph, e, logger)
            graph_remove_edge(
                graph, simp_edge_dict, graph.vp.id[e.source()], graph.vp.id[e.target()]
            )
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, logger, "graph_f1.gfa")


    # # test
    # comp, hist = label_components(graph, None, True, False)
    # comp_dict = {}
    # for node in graph.vertices():
    #     if comp[node] not in comp_dict:
    #         comp_dict[comp[node]] = []
    #     comp_dict[comp[node]].append(graph.vp.id[node])

    # for label, sccs in comp_dict.items():
    #     logger.debug("SCC-" + str(label) + " " + list_to_string(sccs))
    # # end test
