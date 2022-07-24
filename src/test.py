#!/usr/bin/env python3
import logging
from utils.vsAware_IO import flipped_gfa_to_graph, gfa_to_graph, graph_to_gfa
import numpy

import sys
from utils.vsAware_Utilities import *
from utils.vsAware_Split import *

if __name__ == "__main__":
    logger = logging.getLogger("log")
    logger.setLevel(logging.DEBUG)
    consoleHeader = logging.StreamHandler()
    consoleHeader.setLevel(logging.DEBUG)
    consoleHeader.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(consoleHeader)

    graph, simp_node_dict, simp_edge_dict = gfa_to_graph("graph_c.gfa", logger)
    logger.debug("graph is dag? " + str(graph_is_DAG(graph, simp_node_dict)))
    if not graph_is_DAG(graph, simp_node_dict):
        circs = list(all_circuits(graph, unique=False))
        for cir in circs:
            is_par_trivial = False
            for v in cir:
                node = graph.vertex(v)
                if node.in_degree() > 1 and node.out_degree() > 1:
                    is_par_trivial = False
                    break
                elif node.in_degree() > 1 or node.out_degree() > 1:
                    is_par_trivial = True
            # if is_par_trivial:
            #     node =

    logger.debug(circs)

    graph_split_trivial(graph, simp_node_dict, simp_edge_dict, logger)

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, logger, "split_graphc.gfa")
