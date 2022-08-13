#!/usr/bin/env python3
import logging
from utils.vsAware_Path import label_filter_dfs
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
    # contig = ['77&76*B', '78', '100&13*1&101&83*1&82&81*1&91&92*B&112*C&1*B*C&0*6&115*A&72', '25', '24']
    # contig = ['97&96*A&17&98&83*0&104&9*1&8&7*0&6*A', '5', '113&0*0&115*B&28*A&29&30&70&32&33&34*1&35&36*0&67&38*0&66&42*0&43&44*0&45&46&47&48&62']
    # contig = ['102*C&101*B&100&99&98&97*A&1400&293&295&296&298*A&328*A&268*0&329&330&336&261*B&260*B&259*B&258*B&257*B*B&256*0&255*A&341&401&237&402']
    contig = [
        "693*A&692&691&1213&706*0&705*B&710&711&1207&741&740&739&738&737&736*A&1200*B&759*B*B&1195*4&769*A*A*B*A*A&770*A*A*B*A*A&771*A*A*B*A*A&1191*A*A*B*A*A&787*A*A*B*A*A&786*A*B*A*A&785*A*B*A*A&784*A*B*A*A&783*A*B*A*A&782*A*B*A*A&781*A*B*A*A&780*A*B*A*A&779*A*B*A*A&778*A*B*A*A&1188*A*B*A*A&798*A*B*A*A&797*A*B*A*A&1180*A*B*A*A&836*A*B*A*A&835*A*B*A*A&834*A*B*A*A&1163*A*B*A*A&858*B*A*A&857*B*A*A&856*B*A*A&855*B*A*A&854*B*A*A&1152*B*A*A&904*B*A*A&903*B*A*A&901*B*A*A&900*B*A*A&899*B*A*A&898*B*A*A&897*A*A&896*A*A&895*A*A&1138*A*A&967*A&966*A&968*A&969*A&970*A&971*A&972*A&1380*A&116*A&1381*A&997*A&993*A&992*A&991*A&1106*A&1092&1020*A&1019*A&1018*A&1017*A&1096*A&1047*A&1048*A&1087*A&1387*A&109*A&108*A&1389*A&1082*A&1086*A*A&1393*B*A&103*3&102*B*A&101*A*A&1396*A&283*A&282*A&281*A&279*A&278*A&277&321&322&1403&93*A&92&91&90&1407&273&272&324&300&299&298*B&328*B&268&267&266&265&264&263&262&261*A&260*A&259*A&258*A&257*B*A*A&256&255*B&254&253*B&342*B&343*B&357*B&405*C&233*B*C&232*B*C&231&230*B&229*B&228*B&227*B&226*B&411*B&412*B&1433&64&63&62&61&1437&366*B"
    ]

    for e in graph.edges():
        graph.ep.keep[e] = True
    # interior edges
    for i in range(1, len(contig) - 1):
        u = simp_node_dict[contig[i - 1]]
        k = simp_node_dict[contig[i]]
        v = simp_node_dict[contig[i + 1]]
        for e in k.all_edges():
            if (e.source() == u and e.target() == k) or (
                e.source() == k and e.target() == v
            ):
                continue
            graph.ep.keep[e] = False
    # # first contig node edges
    # for oe in simp_node_dict[contig[0]].out_edges():
    #     graph.ep.keep[oe] = False
    # # last contig node edges
    # for ie in simp_node_dict[contig[-1]].in_edges():
    #     graph.ep.keep[ie] = False

    cs = simp_node_dict[contig[0]]
    ct = simp_node_dict[contig[-1]]
    can_reach = reachable(graph, simp_node_dict, ct, ct)
    can_reach_rev = reachable(graph, simp_node_dict, cs, cs)

    for e in graph.edges():
        if not graph.ep.keep[e]:
            print_edge(graph, e, logger)
            graph_remove_edge(
                graph, simp_edge_dict, graph.vp.id[e.source()], graph.vp.id[e.target()]
            )
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, logger, "graph_f1.gfa")

    visited = label_filter_dfs(graph, contig, False)
    # reverse cycle
    visited_rev = label_filter_dfs(graph, contig, True)

    for ino in contig[1:-1]:
        for e in simp_node_dict[ino].all_edges():
            graph.ep.keep[e] = False

    for e in graph.edges():
        if not graph.ep.keep[e] and graph.ep.color[e] == "black":
            print_edge(graph, e, logger)
            graph_remove_edge(
                graph, simp_edge_dict, graph.vp.id[e.source()], graph.vp.id[e.target()]
            )
    print(
        "forward: ",
        can_reach,
    )
    print(
        "reverse: ",
        can_reach_rev,
    )
    # print(path_to_id_string(graph, retrieve_cycle(graph)[0]))
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, logger, "graph_f2.gfa")

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
