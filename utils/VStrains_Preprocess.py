#!/usr/bin/env python3

from logging import Logger
import subprocess
from graph_tool.all import Graph

import numpy
import matplotlib.pyplot as plt

from utils.VStrains_Utilities import *


def reindexing(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Reindex the nodes, with idx-node_id mappings
    """
    idx_mapping = {}
    idx_node_dict = {}
    idx_edge_dict = {}
    idx = 0
    for no, node in simp_node_dict.items():
        if graph.vp.color[node] == "black":
            idx_mapping[no] = str(idx)
            graph.vp.id[node] = str(idx)
            idx_node_dict[str(idx)] = node
            idx += 1
    for (u, v), e in simp_edge_dict.items():
        if (
            graph.ep.color[e] == "black"
            and graph.vp.color[e.source()] == "black"
            and graph.vp.color[e.target()] == "black"
        ):
            idx_edge_dict[(idx_mapping[u], idx_mapping[v])] = e
    return graph, idx_node_dict, idx_edge_dict, idx_mapping


def threshold_estimation(graph: Graph, logger: Logger, temp_dir):
    dps = [graph.vp.dp[node] for node in graph.vertices()]
    # handle edge case, when the graph contains uniform coverage
    if max(dps) == min(dps):
        return 0.00
    regions, bins = numpy.histogram(
        dps, bins=int((max(dps) - min(dps)) // (0.05 * numpy.median(dps)))
    )
    pidx, _ = max(list(enumerate(regions)), key=lambda p: p[1])
    ratio = 0.00
    if pidx == 0:
        ratio = 0.05
        # global peak belongs to first filter region, find maximum peak range, bound by 25% Median
        for i in range(0, 4):
            if i >= len(regions):
                logger.warning(
                    "histogram is not properly set, reset cutoff to default (0.05*M)"
                )
                ratio = 0.05
                break
            if regions[i] > regions[i + 1]:
                ratio += 0.05
            else:
                break
    threshold = ratio * numpy.median(dps)
    plt.figure(figsize=(128, 64))
    for b in bins:
        plt.axvline(b, color="blue")
    plt.hist(x=dps, bins=len(dps))
    plt.axvline(threshold, color="r")
    plt.title("node coverage bar plot")
    plt.xticks(numpy.arange(min(dps), max(dps) + 1, 50.0))
    plt.savefig("{0}{1}".format(temp_dir, "/tmp/bar_plot.png"))
    return threshold


def graph_simplification(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
    min_cov,
):
    """
    Directly remove all the vertex with coverage less than minimum coverage and related edge

    Node belongs to any contigs should not be removed
    return:
        removed_node_dict
        removed_edge_dict
    """
    logger.info("graph simplification")
    logger.debug(
        "Total nodes: "
        + str(len(simp_node_dict))
        + " Total edges: "
        + str(len(simp_edge_dict))
    )
    node_to_contig_dict = {}
    edge_to_contig_dict = {}
    if contig_dict != None:
        node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    # iterate until no more node be removed from the graph
    for id, node in list(simp_node_dict.items()):
        if graph.vp.dp[node] <= min_cov:
            if id in node_to_contig_dict:
                continue

            graph_remove_vertex(graph, simp_node_dict, id, printout=False)

            for e in set(node.all_edges()):
                uid = graph.vp.id[e.source()]
                vid = graph.vp.id[e.target()]
                if (uid, vid) in edge_to_contig_dict:
                    continue
                if (uid, vid) in simp_edge_dict:
                    graph_remove_edge(graph, simp_edge_dict, uid, vid, printout=False)

    logger.debug(
        "Remain nodes: "
        + str(len(simp_node_dict))
        + " Total edges: "
        + str(len(simp_edge_dict))
    )
    logger.info("done")
    return


# ------------------------------------LEGACY------------------------------------#
def paths_from_src(graph: Graph, simp_node_dict: dict, self_node, src, maxlen):
    """
    retrieve all the path from src node to any node
    within maxlen restriction, in straight direction
    """

    def dfs_rev(graph: Graph, u, curr_path: list, maxlen, visited, all_path):
        visited[u] = True
        curr_path.append(u)
        curr_len = path_len(graph, curr_path)
        if curr_len >= maxlen:
            all_path.append(list(curr_path))
        else:
            for v in u.out_neighbors():
                if not visited[v]:
                    dfs_rev(graph, v, curr_path, maxlen, visited, all_path)
        curr_path.pop(-1)
        visited[u] = False
        return

    visited = {}
    for u in graph.vertices():
        if graph.vp.id[u] not in simp_node_dict:
            visited[u] = True
        else:
            visited[u] = False
    visited[self_node] = True
    all_path = []
    dfs_rev(graph, src, [], maxlen, visited, all_path)
    return all_path


def paths_to_tgt(graph: Graph, simp_node_dict: dict, self_node, tgt, maxlen):
    """
    retrieve all the path from any node to tgt node
    within maxlen restriction, in reverse direction
    """

    def dfs_rev(graph: Graph, v, curr_path: list, maxlen, visited, all_path):
        visited[v] = True
        curr_path.insert(0, v)
        curr_len = path_len(graph, curr_path)
        if curr_len >= maxlen:
            all_path.append(list(curr_path))
        else:
            for u in v.in_neighbors():
                if not visited[u]:
                    dfs_rev(graph, u, curr_path, maxlen, visited, all_path)
        curr_path.pop(0)
        visited[v] = False
        return

    visited = {}
    for u in graph.vertices():
        if graph.vp.id[u] not in simp_node_dict:
            visited[u] = True
        else:
            visited[u] = False
    visited[self_node] = True
    all_path = []
    dfs_rev(graph, tgt, [], maxlen, visited, all_path)
    return all_path


def tip_removal_s(
    graph: Graph,
    simp_node_dict: dict,
    contig_dict: dict,
    logger: Logger,
    tempdir,
    accept_rate=0.99,
):
    if not graph_is_DAG(graph, simp_node_dict):
        logger.info("Graph is Cyclic, tip removal start..")
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(
                graph, simp_node_dict, logger, tempdir, accept_rate
            )
        for cno, [contig, _, ccov] in list(contig_dict.items()):
            if not all([no in simp_node_dict for no in contig]):
                subcontigs = []
                curr_contig = []
                addLast = False
                for no in contig:
                    if no in simp_node_dict:
                        addLast = True
                        curr_contig.append(no)
                    else:
                        addLast = False
                        if curr_contig != []:
                            subcontigs.append(curr_contig[:])
                        curr_contig = []
                if addLast:
                    subcontigs.append(curr_contig[:])

                contig_dict.pop(cno)
                for i, subc in enumerate(subcontigs):
                    sublen = path_len(graph, [simp_node_dict[c] for c in subc])
                    contig_dict[cno + "^" + str(i)] = [subc, sublen, ccov]
    else:
        logger.info("Graph is DAG, tip removal skipped.")
    logger.info("done")
    return


def tip_removal(
    graph: Graph, simp_node_dict: dict, logger: Logger, tempdir, accept_rate
):
    """
    retrieve all the source/tail simple path, and merge them into adjacent neighbor path if possible

    the collapse step can be done before node depeth rebalance, since it only regards to
    matching score within node seq len

    if is the case, then spades contig may also be modified.
    """

    def remove_tip(graph: Graph, simp_node_dict: dict, from_node, to_path):
        """
        collapse the node with the given path, increment given path depth, remove related information
        about the node.
        """
        graph.vp.color[from_node] = "gray"
        pending_dp = graph.vp.dp[from_node]
        for node in to_path:
            graph.vp.dp[node] += pending_dp
        simp_node_dict.pop(graph.vp.id[from_node])
        for e in from_node.all_edges():
            graph.ep.color[e] = "gray"
        logger.debug(
            path_to_id_string(
                graph,
                to_path,
                "Tip Node {0} collapsed to path".format(graph.vp.id[from_node]),
            )
        )
        return

    def cand_collapse_path(graph: Graph, from_node, to_paths, temp_dir):
        """
        use minimap2 -c to evaluation the node-path similarity, sort based on matching score in DESC order

        return: the most similar path if there exist a path with score >= accept rate, else return None
        """
        ref_loc = "{0}/ref.fa".format(temp_dir)
        query_loc = "{0}/query.fa".format(temp_dir)
        overlap_loc = "{0}/overlap.paf".format(temp_dir)
        subprocess.check_call(
            "touch {0}; echo > {0}; touch {1}; echo > {1}".format(ref_loc, query_loc),
            shell=True,
        )

        id_path_dict = {}
        for id, path in list(enumerate(to_paths)):
            id_path_dict[id] = path

        # retrieve all the path information and save into ref.fa
        with open(ref_loc, "w") as ref_file:
            for id, path in id_path_dict.items():
                name = ">" + str(id) + "\n"
                seq = path_to_seq(graph, path, id) + "\n"
                ref_file.write(name)
                ref_file.write(seq)
            ref_file.close()

        # save from node info to query.fa
        with open(query_loc, "w") as query_file:
            name = ">" + graph.vp.id[from_node] + "\n"
            seq = path_to_seq(graph, [from_node], name) + "\n"
            query_file.write(name)
            query_file.write(seq)
            query_file.close()

        # minimap to obtain matching score for all node-path
        id_evalscore = {}
        minimap_api(ref_loc, query_loc, overlap_loc)
        with open(overlap_loc, "r") as overlap_file:
            for Line in overlap_file:
                splited = (Line[:-1]).split("\t")
                path_no = int(splited[5])
                nmatch = int(splited[9])
                nblock = int(splited[10])
                if path_no not in id_evalscore:
                    id_evalscore[path_no] = [nmatch / nblock]
                else:
                    id_evalscore[path_no].append(nmatch / nblock)
            overlap_file.close()

        # remove temp file
        subprocess.check_call(
            "rm {0}; rm {1}; rm {2}".format(ref_loc, query_loc, overlap_loc), shell=True
        )

        id_evalscore_sum = []
        for id, scores in id_evalscore.items():
            mean_score = numpy.mean(scores) if len(scores) != 0 else 0
            id_evalscore_sum.append((id, mean_score))

        best_match = sorted(id_evalscore_sum, key=lambda t: t[1], reverse=True)
        logger.debug("Tip Node: " + str(graph.vp.id[from_node]) + str(best_match))
        if len(best_match) == 0:
            return None
        elif best_match[0][1] >= accept_rate:
            return id_path_dict[best_match[0][0]]
        else:
            return None

    is_removed = True
    # get all the source simple path
    src_nodes = []
    tgt_nodes = []
    isolated_node = []
    for node in simp_node_dict.values():
        if node.in_degree() + node.out_degree() == 0:
            isolated_node.append(node)
        elif node.in_degree() == 0:
            src_nodes.append(node)
        elif node.out_degree() == 0:
            tgt_nodes.append(node)
        else:
            None

    # src node collapse
    src_nodes = sorted(src_nodes, key=lambda x: graph.vp.dp[x])
    for src in src_nodes:
        src_len = path_len(graph, [src])
        potential_paths = []
        # path retrieve
        for out_branch in src.out_neighbors():
            if graph.vp.id[out_branch] not in simp_node_dict:
                continue
            # print("current out branch: ", graph.vp.id[out_branch])
            for in_tgt in out_branch.in_neighbors():
                if graph.vp.id[in_tgt] == graph.vp.id[src]:
                    # coincidence path
                    continue
                if graph.vp.id[in_tgt] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                # print("current in tgt: ", graph.vp.id[in_tgt])
                potential_paths.extend(
                    paths_to_tgt(graph, simp_node_dict, src, in_tgt, src_len)
                )
        cand_path = cand_collapse_path(graph, src, potential_paths, tempdir)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, src, cand_path)
            is_removed = False

    # target node collapse
    tgt_nodes = sorted(tgt_nodes, key=lambda x: graph.vp.dp[x])
    for tgt in tgt_nodes:
        tgt_len = path_len(graph, [tgt])
        potential_paths = []
        # path retrieve
        for in_branch in tgt.in_neighbors():
            if graph.vp.id[in_branch] not in simp_node_dict:
                continue
            # print("current in branch: ", graph.vp.id[in_branch])
            for out_src in in_branch.out_neighbors():
                if graph.vp.id[out_src] == graph.vp.id[tgt]:
                    # coincidence path
                    continue
                if graph.vp.id[out_src] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                # print("current out src: ", graph.vp.id[out_src])
                potential_paths.extend(
                    paths_from_src(graph, simp_node_dict, tgt, out_src, tgt_len)
                )
        cand_path = cand_collapse_path(graph, tgt, potential_paths, tempdir)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, tgt, cand_path)
            is_removed = False
    return is_removed
