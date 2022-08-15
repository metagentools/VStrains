#!/usr/bin/env python3

from logging import Logger
from graph_tool.all import Graph
from graph_tool.draw import graph_draw
from graph_tool.topology import all_circuits
import subprocess
import sys

import numpy
from functools import reduce


def map_ref_to_graph(
    ref_file,
    simp_node_dict: dict,
    graph_file,
    logger: Logger,
    store_mapping=False,
    output_file="overlap.paf",
    fasta_file="temp_gfa_to_fasta.fasta",
):
    """
    map reference strain to the graph, debug only
    """
    if not ref_file:
        logger.debug("No ref file imported")
        return -1
    with open(graph_file, "r") as gfa:
        with open(fasta_file, "w") as fasta:
            for Line in gfa:
                splited = (Line[:-1]).split("\t")
                if splited[0] == "S":
                    fasta.write(">{0}\n{1}\n".format(splited[1], splited[2]))
            fasta.close()
        gfa.close()

    subprocess.check_call(
        "minimap2 {0} {1} -c > {2}".format(ref_file, fasta_file, output_file),
        shell=True,
    )

    strain_dict = {}
    with open(output_file, "r") as paf:
        for Line in paf:
            splited = (Line[:-1]).split("\t")
            seg_no = str(splited[0])
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = str(splited[5])
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if seg_no not in simp_node_dict:
                continue
            if (nmatch / nblock) == 1:
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = []
                strain_dict[ref_no].append(seg_no)
        paf.close()

    if not store_mapping:
        subprocess.check_call(
            "rm {0}; rm {1}".format(output_file, fasta_file), shell=True
        )

    logger.debug("strain dict mapping")
    for seg_no, strains in strain_dict.items():
        logger.debug("strains: " + str(seg_no) + " Path: " + list_to_string(strains))
        logger.debug("-------------------")
    return strain_dict


def map_ref_to_contig(contig_dict: dict, logger: Logger, paf_file):
    logger.debug("map ref to contig")
    strain_dict = {}
    with open(paf_file, "r") as paf:
        for Line in paf:
            splited = (Line[:-1]).split("\t")
            seg_no = str(splited[0].split("_")[0])
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = str(splited[5])
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if seg_no not in contig_dict:
                continue
            if (nmatch / nblock) >= 0.99:
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = set()
                strain_dict[ref_no].add(seg_no)
        paf.close()

    for sno, cnos in strain_dict.items():
        logger.debug("--------------------------------->")
        logger.debug(
            "contig-strains: "
            + str(sno)
            + "Count: "
            + str(len(cnos))
            + " - Contigs: "
            + str(
                [
                    (cno1, clen, ccov)
                    for cno1, [_, clen, ccov] in contig_dict.items()
                    if cno1 in cnos
                ]
            )
        )
        rel_nodes = []
        for cno in cnos:
            rel_nodes.extend(contig_dict[cno][0])
        logger.debug("related nodes in contig: " + list_to_string(rel_nodes))


def minimap_api(ref_file, fasta_file, output_file):
    subprocess.check_call(
        "minimap2 {0} {1} -c > {2}".format(ref_file, fasta_file, output_file),
        shell=True,
    )
    return


def trim_contig_dict(
    graph: Graph, simp_node_dict: dict, contig_dict: dict, logger: Logger
):
    logger.info("trim contig..")
    for cno, [contig, _, ccov] in list(contig_dict.items()):
        new_contig = list(dict.fromkeys(contig))
        contig_dict[cno] = [
            new_contig,
            path_len(graph, [simp_node_dict[no] for no in new_contig]),
            ccov,
        ]
    logger.info("done")
    return contig_dict


def contig_map_node(contig_dict: dict):
    node_to_contig_dict = {}
    edge_to_contig_dict = {}
    for cno, (c, _, _) in contig_dict.items():
        for n in c:
            if n not in node_to_contig_dict:
                node_to_contig_dict[n] = {cno}
            else:
                node_to_contig_dict[n].add(cno)
        for i in range(len(c)):
            c_i = c[i]
            c_i_1 = c[i + 1] if (i < len(c) - 1) else None
            if c_i_1 != None:
                if (c_i, c_i_1) not in edge_to_contig_dict:
                    edge_to_contig_dict[(c_i, c_i_1)] = {cno}
                else:
                    edge_to_contig_dict[(c_i, c_i_1)].add(cno)
    return node_to_contig_dict, edge_to_contig_dict


def contig_cov_fix(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
):
    """
    if found a single node contig that also not appearing in the simp_node_dict, check if is mapped to split contig
    """
    for cno, [contig, clen, _] in list(contig_dict.items()):
        newccov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
        contig_dict[cno][2] = newccov
        if cno in contig_dict:
            if logger != None:
                print_contig(cno, clen, contig_dict[cno][2], contig, logger)
    return


def graph_reduction_c(graph: Graph, cand_path, usage_dict: dict, cand_cov):
    """
    reduce the graph coverage based on given path and cov,
    only applied after udp be deployed in the graph
    """
    for i in range(len(cand_path)):
        graph.vp.dp[cand_path[i]] -= cand_cov
        usage_dict[graph.vp.id[cand_path[i]]] += 1

    for i in range(len(cand_path) - 1):
        e = graph.edge(cand_path[i], cand_path[i + 1])
        # print(e, graph.vp.id[cand_path[i]], graph.vp.id[cand_path[i+1]])
        graph.ep.flow[e] -= cand_cov


def contig_dict_remapping(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    id_mapping: dict,
    prev_ids: list,
    logger: Logger,
):
    """
    Update the contig nodes to mapped nodes.
    """

    def map_contig_tree(cno, contig, id_mappingP: dict):
        paths = []
        if len(id_mappingP[contig[0]]) == 0:
            paths = [[contig[0]]]
        else:
            paths = [[s] for s in id_mappingP[contig[0]]]
        for i in range(1, len(contig)):
            acc_paths = []
            next = contig[i]
            for p in paths:
                last = p[-1]
                if len(id_mappingP[next]) == 0:
                    if (last, next) in simp_edge_dict:
                        acc_paths.append((p + [next]))
                else:
                    for nextm in id_mappingP[next]:
                        if (last, nextm) in simp_edge_dict:
                            acc_paths.append((p + [nextm]))
            paths = acc_paths
        # print("----------------> curr contig tree mapping: ", cno, " mapped count: ", len(paths))
        # for p in paths:
        #     print(list_to_string(p))
        return paths

    def merge_id(id_mapping_r: dict, curr_set: set, myid):
        if len(curr_set) == 0:
            return set([myid])
        else:
            rtn_set = set()
            for id in curr_set:
                rtn_set = rtn_set.union(merge_id(id_mapping_r, id_mapping_r[id], id))
            return rtn_set

    logger.info("contig resolution..")
    # id_mapping merging, recursive merge down.
    red_id_mapping = {}

    for id in prev_ids:
        all_set = merge_id(id_mapping, id_mapping[id], id)
        red_id_mapping[id] = all_set
        logger.debug("Node {0} maps to {1}".format(id, all_set))

    for cno, (contig, _, _) in list(contig_dict.items()):
        logger.debug("---------------------------------------------")
        logger.debug(
            "Current mapping contig: " + str(cno) + ", " + list_to_string(contig)
        )
        paths = map_contig_tree(cno, contig, red_id_mapping)
        # split the contig tree to avoid the ambiguity variation
        if len(paths) < 1:
            logger.debug("error, contig missed: " + str(cno) + str(contig))
        elif len(paths) == 1:
            if paths[0] == contig:
                logger.debug("single mapping, keep original")
            else:
                logger.debug("single mapping, replace" + list_to_string(paths[0]))
                contig_dict.pop(cno)
                subcov = path_cov(graph, simp_node_dict, simp_edge_dict, paths[0])
                contig_dict[cno] = [
                    paths[0],
                    path_len(graph, [simp_node_dict[no] for no in paths[0]]),
                    subcov,
                ]
        else:
            contig_dict.pop(cno)
            logger.debug(
                "multi mapping for the current contig: whole contig is ambiguous mapping, keep the intersection reduced one only"
                + cno
            )
            final_path = reduce(lambda a, b: [i for i in a if i in b], paths)
            if len(final_path) > 0:
                # at least one node
                logger.debug("selected mapped contig: " + str(final_path))
                sublen = path_len(graph, [simp_node_dict[no] for no in final_path])
                subcov = path_cov(graph, simp_node_dict, simp_edge_dict, final_path)
                contig_dict[cno] = [final_path, sublen, subcov]
    logger.info("done")
    return


def simp_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    find simple edges, simple edge is the only edge between its source and sink
    """
    simple_edges = []
    out_edge = {}
    in_edge = {}
    for e in simp_edge_dict.values():
        src = e.source()
        target = e.target()
        if (
            graph.vp.id[src] not in simp_node_dict
            or graph.vp.id[target] not in simp_node_dict
        ):
            continue
        if src.out_degree() == 1 and target.in_degree() == 1:
            if src != target:
                simple_edges.append([src, target])
                in_edge[src] = e
                out_edge[target] = e

    # build simple paths from simple edges
    def extend_path(p):
        v = p[-1]
        if v in in_edge:
            p.append(in_edge[v].target())
            return extend_path(p)
        else:
            return p

    simple_paths = []
    for v, e in in_edge.items():
        if v not in out_edge:
            p = extend_path([e.source(), e.target()])
            simple_paths.append(p)
    return simple_paths


def simple_paths_to_dict(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    simple_paths = simp_path(graph, simp_node_dict, simp_edge_dict)
    simp_path_dict = {}
    for id, p in enumerate(simple_paths):
        pids = [graph.vp.id[n] for n in p]
        name = str(id)
        clen = path_len(graph, p)
        cov = numpy.max([graph.vp.dp[n] for n in p])
        simp_path_dict[name] = [pids, clen, cov]
        # print("Simple PATH: ", list_to_string(pids))
    return simp_path_dict


def simp_path_compactification(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
):
    """
    reduce all the contig to a single node, and keep all the potential src/tgt edge.

    1. reduce the coverage for each involving node by the amount of contig cov
    2. reconnect end-to-end nodes to the contig node
    """
    logger.info("simple path compactification..")
    simp_path_dict = simple_paths_to_dict(graph, simp_node_dict, simp_edge_dict)

    graph_backup = graph.copy()
    simp_node_dict_backup = simp_node_dict.copy()

    node_to_simp_node = {}
    for id in simp_node_dict.keys():
        node_to_simp_node[id] = id

    contig_info = []
    # reduce all the simple path to a single node from the graph
    for cno, (contig, _, ccov) in list(simp_path_dict.items()):
        src = contig[0]
        tgt = contig[-1]
        id = contig[0]
        for nid in contig[1:]:
            id += "&" + nid
        cseq = path_to_seq(
            graph_backup, [simp_node_dict_backup[n] for n in contig], cno
        )
        in_edges = list(
            (graph_backup.vp.id[e.source()], src, graph_backup.ep.overlap[e])
            for e in simp_node_dict_backup[src].in_edges()
        )
        out_edges = list(
            (tgt, graph_backup.vp.id[e.target()], graph_backup.ep.overlap[e])
            for e in simp_node_dict_backup[tgt].out_edges()
        )

        for i in range(len(contig)):
            no = contig[i]
            node_to_simp_node[no] = id
            graph_remove_vertex(graph, simp_node_dict, no, printout=False)
            if i != len(contig) - 1:
                graph_remove_edge(
                    graph, simp_edge_dict, contig[i], contig[i + 1], printout=False
                )
        cv = graph_add_vertex(graph, simp_node_dict, id, ccov, cseq, printout=False)
        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])

    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u, v, o) in in_edges:
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                graph_add_edge(
                    graph,
                    simp_edge_dict,
                    simp_node_dict[u],
                    u,
                    node,
                    graph.vp.id[node],
                    o,
                    printout=False,
                )

            for [_, tgt, _, in_node, _, _] in contig_info:
                if (
                    tgt == u
                    and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict
                ):
                    graph_add_edge(
                        graph,
                        simp_edge_dict,
                        in_node,
                        graph.vp.id[in_node],
                        node,
                        graph.vp.id[node],
                        o,
                        printout=False,
                    )

        for (u, v, o) in out_edges:
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                graph_add_edge(
                    graph,
                    simp_edge_dict,
                    node,
                    graph.vp.id[node],
                    simp_node_dict[v],
                    v,
                    o,
                    printout=False,
                )

            for [src, _, _, out_node, _, _] in contig_info:
                if (
                    src == v
                    and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict
                ):
                    graph_add_edge(
                        graph,
                        simp_edge_dict,
                        node,
                        graph.vp.id[node],
                        out_node,
                        graph.vp.id[out_node],
                        o,
                        printout=False,
                    )
    # fix the contig, with simple path be concated
    for cno, (contig, _, ccov) in list(contig_dict.items()):
        new_contig = []
        for no in contig:
            if node_to_simp_node[no] == no:
                new_contig.append(no)
            else:
                if len(new_contig) == 0:
                    new_contig.append(node_to_simp_node[no])
                else:
                    if node_to_simp_node[no] != new_contig[-1]:
                        new_contig.append(node_to_simp_node[no])
        logger.debug(
            "cno: {0} from {1} to {2}".format(
                cno, list_to_string(contig), list_to_string(new_contig)
            )
        )
        contig_dict[cno] = [
            new_contig,
            path_len(graph, [simp_node_dict[no] for no in new_contig]),
            ccov,
        ]
    logger.info("done")
    return


def contig_low_cov_removal(contig_dict: dict, logger: Logger, threshold):
    for cno in list(contig_dict.keys()):
        if contig_dict[cno][2] <= threshold:
            logger.debug(
                "remove low coverage contig: "
                + str(cno)
                + " with cov: "
                + str(contig_dict[cno][2])
            )
            contig_dict.pop(cno)


def contig_dup_removed_s(contig_dict: dict, logger: Logger):
    logger.info("drop duplicated contigs..")
    dup_contig_ids = set()
    for cno1 in contig_dict.keys():
        contig1, _, _ = contig_dict[cno1]
        for cno2 in contig_dict.keys():
            if (
                cno1 not in dup_contig_ids
                and cno2 not in dup_contig_ids
                and cno1 != cno2
            ):
                contig2, _, _ = contig_dict[cno2]
                # use set equality to avoid cyclic contig
                intersect = set(contig1).intersection(set(contig2))
                if len(intersect) == len(contig1) and len(intersect) == len(contig2):
                    # duplicated
                    dup_contig_ids.add(cno2)
                elif len(intersect) == len(contig1):
                    # contig1 is subcontig of contig2
                    dup_contig_ids.add(cno1)
                elif len(intersect) == len(contig2):
                    # contig2 is subcontig of contig1
                    dup_contig_ids.add(cno2)
    for cno in dup_contig_ids:
        contig_dict.pop(cno)
    logger.debug("duplicated contigs: " + str(dup_contig_ids))
    logger.info("done")
    return contig_dict


def concat_overlap_contig(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
):
    def self_loop(contig):
        return (contig[-1], contig[0]) in simp_edge_dict

    logger.info("concat overlapped contig..")
    contig_overlap_dict = {}
    for key in contig_dict.keys():
        contig_overlap_dict[key] = []

    for cno, [contig, _, _] in contig_dict.items():
        for cno2, [contig2, _, _] in contig_dict.items():
            if cno == cno2:
                continue
            if self_loop(contig) or self_loop(contig2):
                continue
            isParallel, intersects, status = check_contig_intersection(contig, contig2)
            if not isParallel:
                if status in ["f", "d"]:
                    # forward/double e2e overlap
                    contig_overlap_dict[cno].append((cno2, intersects))
                elif status == "n":
                    if (
                        simp_node_dict[contig2[0]]
                        in simp_node_dict[contig[-1]].out_neighbors()
                        and simp_node_dict[contig[0]]
                        in simp_node_dict[contig2[-1]].out_neighbors()
                    ):
                        # contig <-> contig2, circular touch
                        contig_overlap_dict[cno].append((cno2, []))
    logger.debug("--contig overlap info: " + str(contig_overlap_dict))
    overlap_graph = Graph(directed=True)
    overlap_graph.vp.id = overlap_graph.new_vertex_property("string")
    node_dict = dict()
    concat_dict = dict()
    for cno in contig_overlap_dict.keys():
        node = overlap_graph.add_vertex()
        overlap_graph.vp.id[node] = cno
        node_dict[cno] = node
    for cno, cno2s in contig_overlap_dict.items():
        u = node_dict[cno]
        for (cno2, intersects) in cno2s:
            v = node_dict[cno2]
            edge = overlap_graph.add_edge(u, v)
            concat_dict[(cno, cno2)] = intersects
    if not graph_is_DAG_simp(overlap_graph, node_dict):
        # contig overlap graph contains circuits
        # all_circuits can only be called if graph is not dag
        circuits = list(all_circuits(overlap_graph, unique=True))
        for k in range(len(circuits)):
            cyc = circuits[k]
            logger.debug(
                "current cyc: "
                + str([overlap_graph.vp.id[overlap_graph.vertex(v)] for v in cyc])
            )
            unique_cyc = True
            for j in range(len(circuits)):
                cyc2 = circuits[j]
                if k != j and len(set(cyc).intersection(set(cyc2))) > 0:
                    unique_cyc = False
            for i in range(len(cyc)):
                u = overlap_graph.vertex(cyc[i])
                v = overlap_graph.vertex(cyc[(i + 1) % len(cyc)])
                cyc_e = overlap_graph.edge(u, v)

                for e in set(u.out_edges()):
                    if e != cyc_e or not unique_cyc:
                        overlap_graph.remove_edge(e)
            if unique_cyc:
                s = overlap_graph.vertex(cyc[0])
                t = overlap_graph.vertex(cyc[1])
                overlap_graph.remove_edge(overlap_graph.edge(s, t))
    has_del = True
    while has_del:
        has_del = False
        sorted_nodes = sorted(overlap_graph.vertices(), reverse=True)
        for node in sorted_nodes:
            if (node.in_degree() == 0 and node.out_degree() == 0) or (
                node.in_degree() > 1 or node.out_degree() > 1
            ):
                for edge in set(node.all_edges()):
                    overlap_graph.remove_edge(edge)
                overlap_graph.remove_vertex(node)
                has_del = True
    srcs = filter(lambda n: n.in_degree() == 0, overlap_graph.vertices())
    for src in srcs:
        contig_path = []
        curr = src
        while curr != None:
            contig_path.append(curr)
            curr = list(curr.out_neighbors())[0] if curr.out_degree() == 1 else None
        concat_contig = []
        cnos = ""
        logger.debug(
            "contig path: " + str([overlap_graph.vp.id[k] for k in contig_path])
        )
        for ind, contig_node in enumerate(contig_path):
            ccno = overlap_graph.vp.id[contig_node]
            contig, _, _ = contig_dict.pop(ccno)
            if ind < len(contig_path) - 1:
                cnos += ccno + "&"
                vid = overlap_graph.vp.id[contig_path[ind + 1]]
                intersect = concat_dict[(ccno, vid)]
                if intersect != []:
                    if intersect.count(None) > 0:
                        cut = list(reversed(intersect)).index(None)
                        contig = contig[:-cut]
                    else:
                        logger.debug("error: " + str(contig) + str(intersect))
                        sys.exit(1)
            else:
                cnos += ccno
            concat_contig.extend(contig)
        logger.debug("concat end-to-end overlap contig: " + str(cnos))
        logger.debug("after concat: " + str(concat_contig))
        concat_len = path_len(graph, [simp_node_dict[id] for id in concat_contig])
        concat_cov = path_cov(graph, simp_node_dict, simp_edge_dict, concat_contig)
        contig_dict[cnos] = [concat_contig, concat_len, concat_cov]
    logger.info("done")
    return


def check_contig_intersection(contig, contig2):
    """
    check if contig1 and 2 are overlapped end-to-end or intersect as parallel
    return true if two contigs are parallel, false if overlap end-to-end
    direction:
    'o': overlap
    'f': forward
    'b': backward
    'd': double, both forward and backward
    """
    # intersection region is proper subset for both contig1 and contig2
    # check intersection
    intersect = set(contig).intersection(set(contig2))
    if len(intersect) <= 0:
        return False, None, "n"

    if len(intersect) == len(contig) or len(intersect) == len(contig2):
        return True, None, "o"

    # determine intermediate intersection
    intersect_maps = [intersect.__contains__(c) for c in contig]
    prev_false_index = intersect_maps.index(False)
    for j in range(prev_false_index + 1, len(intersect_maps)):
        if not intersect_maps[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                return True, None, "o"

    # determine intermediate intersection
    intersect_maps2 = [intersect.__contains__(c) for c in contig2]
    prev_false_index = intersect_maps2.index(False)
    for j in range(prev_false_index + 1, len(intersect_maps2)):
        if not intersect_maps2[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                return True, None, "o"

    # determine direction
    if contig[0] == contig2[0]:
        return True, None, "o"
    if contig[-1] == contig2[-1]:
        return True, None, "o"

    intersect_path = [n if intersect_maps[i] else None for i, n in enumerate(contig)]
    direction = None
    if intersect_maps[0]:
        direction = "b"
    if intersect_maps[-1]:
        direction = "f" if direction == None else "d"
    return False, intersect_path, direction


def strain_repeat_resol(
    graph: Graph,
    simp_node_dict: dict,
    strain_dict: dict,
    contig_info: dict,
    copy_contig_dict: dict,
    logger,
):
    logger.info("resolving repeat nodes in strain..")
    # if contig is fully aligned to the strain, then map the repeat to the strain
    for sno, [strain, _, scov] in list(strain_dict.items()):
        cnos = set()
        subids = []
        for id in strain:
            for iid in str(id).split("&"):
                if iid.find("*") != -1:
                    iid = iid[: iid.find("*")]
                subids.append(iid)
        for cno, [contig, _, _] in copy_contig_dict.items():
            if set(contig).issubset(set(subids)):
                cnos.add(cno)

        repeat_dec = dict.fromkeys(subids, 1)
        # cnos: all related contig id for current strain
        for cno in cnos:
            (_, repeat_dict) = contig_info[cno]
            for no, rpc in repeat_dict.items():
                repeat_dec[no] = max(repeat_dec[no], rpc)
        strain_r = []
        [strain_r.extend([id] * repeat_dec[id]) for id in subids]
        strain_dict[sno] = [
            strain_r,
            path_len(graph, [simp_node_dict[no] for no in strain_r]),
            scov,
        ]
    logger.info("done")
    return


def path_len(graph: Graph, path):
    """
    Find length of the linear path.
    """
    lens = sum([len(graph.vp.seq[u]) for u in path])
    for i in range(len(path) - 1):
        u = path[i]
        v = path[i + 1]
        e = graph.edge(u, v)
        if e != None:
            lens -= graph.ep.overlap[graph.edge(u, v)]
    return lens


def path_cov(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, path):
    """
    Compute the coverage for the path
    """
    eflow = contig_flow(graph, simp_edge_dict, path)
    if len(eflow) < 1:
        # only one vertex
        return graph.vp.dp[simp_node_dict[path[0]]]
    else:
        return min(eflow)


def contig_edges(contig):
    """
    contig edges
    """
    edges = []
    if len(contig) < 2:
        return edges
    for i in range(len(contig) - 1):
        edges.append((contig[i], contig[i + 1]))

    return edges


def contig_flow(graph: Graph, edge_dict: dict, contig):
    """
    edge flow for the contig
    """
    edge_flow = []
    if len(contig) < 2:
        return edge_flow
    for i in range(len(contig) - 1):
        e = edge_dict[(contig[i], contig[i + 1])]
        f = graph.ep.flow[e]
        edge_flow.append(f)

    return edge_flow


def path_ids_to_seq(graph: Graph, path_ids: list, path_name, simp_node_dict: dict):
    seq = ""
    for i in range(len(path_ids)):
        u = simp_node_dict[path_ids[i]]
        if i == len(path_ids) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            e = graph.edge(u, simp_node_dict[path_ids[i + 1]])
            overlap_len = graph.ep.overlap[e] if e != None else 0
            if overlap_len == 0:
                seq = seq + graph.vp.seq[u]
            else:
                seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq


def path_to_seq(graph: Graph, path: list, path_name):
    seq = ""
    for i in range(len(path)):
        u = path[i]
        if i == len(path) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            overlap_len = graph.ep.overlap[graph.edge(u, path[i + 1])]
            if overlap_len == 0:
                seq = seq + graph.vp.seq[u]
            else:
                seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq


def graph_stat(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    print("-------------------------graph stat----------------------")
    for v in simp_node_dict.values():
        print_vertex(graph, v, "stat")
    for e in simp_edge_dict.values():
        print_edge(graph, e, "stat")

    print("-----------------------graph stat end--------------------")


def graph_add_vertex(
    graph: Graph,
    simp_node_dict: dict,
    id,
    dp,
    seq,
    s="add vertex",
    color="black",
    printout=False,
):
    node = graph.add_vertex()
    graph.vp.id[node] = id
    graph.vp.dp[node] = dp
    graph.vp.seq[node] = seq
    graph.vp.color[node] = color
    simp_node_dict[id] = node
    if printout:
        print_vertex(graph, node, s)
    return node


def graph_remove_vertex(
    graph, simp_node_dict: dict, id, s="remove vertex", color="gray", printout=False
):
    node = simp_node_dict[id]
    graph.vp.color[node] = color
    simp_node_dict.pop(id)
    if printout:
        print_vertex(graph, node, s)
    return node


def graph_add_edge(
    graph: Graph,
    simp_edge_dict: dict,
    src,
    src_id,
    tgt,
    tgt_id,
    overlap,
    flow=0,
    s="add edge",
    color="black",
    printout=False,
):
    edge = graph.add_edge(src, tgt)
    graph.ep.overlap[edge] = overlap
    graph.ep.color[edge] = color
    graph.ep.flow[edge] = flow
    simp_edge_dict[(src_id, tgt_id)] = edge
    if printout:
        print_edge(graph, edge, s)
    return edge


def graph_remove_edge(
    graph: Graph,
    simp_edge_dict: dict,
    src_id,
    tgt_id,
    s="remove edge",
    color="gray",
    printout=False,
):
    edge = simp_edge_dict.pop((src_id, tgt_id))
    graph.ep.color[edge] = color
    if printout:
        print_edge(graph, edge, s)
    return edge


def draw_graph_api(graph: Graph, output_file):
    output_size = 120 * (graph.num_edges() + graph.num_vertices())
    vsize = 30
    graph_draw(
        g=graph,
        output=output_file,
        bg_color="white",
        vertex_size=vsize,
        output_size=(output_size, output_size),
    )


def reverse_seq(seq: str):
    return "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[x] for x in reversed(seq))


def print_edge(graph, e, logger: Logger, s=""):
    logger.debug(
        s
        + " edge: "
        + str(graph.vp.id[e.source()])
        + " -> "
        + str(graph.vp.id[e.target()])
        + " "
        + str(graph.ep.flow[e])
        + " "
        + str(graph.ep.color[e])
    )


def print_vertex(graph, v, logger: Logger, s=""):
    logger.debug(
        s
        + " vertex: "
        + graph.vp.id[v]
        + ", dp: "
        + str(graph.vp.dp[v])
        + ", in_degree: "
        + str(v.in_degree())
        + ", out_degree: "
        + str(v.out_degree())
        + graph.vp.color[v]
    )


def print_contig(cno, clen, ccov, contig, logger: Logger, s=""):
    logger.debug(
        s
        + " Contig: "
        + cno
        + ", length: "
        + str(clen)
        + ", cov: "
        + str(ccov)
        + "Path: "
        + list_to_string(contig)
    )


def list_to_string(ids: list, s=""):
    string = s + " - "
    for id in ids:
        string += str(id) + ", "
    return string[:-2] if len(string) >= 2 else ""


def path_to_id_string(graph: Graph, path, s=""):
    return list_to_string([graph.vp.id[node] for node in path], s)


def add_global_source_sink(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    create a global source node & sink node, and concat all the curr src/sink nodes to the
    global source/sink node
    """
    # find all the srcs-targets
    src_nodes = [node for node in graph.vertices() if node.in_degree() == 0]
    tgt_nodes = [node for node in graph.vertices() if node.out_degree() == 0]

    global_src = graph.add_vertex()
    graph.vp.id[global_src] = "global_src"
    graph.vp.dp[global_src] = 0
    graph.vp.color[global_src] = "black"
    simp_node_dict[graph.vp.id[global_src]] = global_src

    for src in src_nodes:
        e = graph.add_edge(global_src, src)
        graph.ep.flow[e] = graph.vp.dp[src]
        graph.ep.color[e] = "black"
        graph.ep.overlap[e] = 0
        graph.vp.dp[global_src] += graph.ep.flow[e]
        simp_edge_dict[(graph.vp.id[global_src], graph.vp.id[src])] = e

    global_sink = graph.add_vertex()
    graph.vp.id[global_sink] = "global_sink"
    graph.vp.dp[global_sink] = 0
    graph.vp.color[global_sink] = "black"
    simp_node_dict[graph.vp.id[global_sink]] = global_sink

    for tgt in tgt_nodes:
        e = graph.add_edge(tgt, global_sink)
        graph.ep.flow[e] = graph.vp.dp[tgt]
        graph.ep.color[e] = "black"
        graph.ep.overlap[e] = 0
        graph.vp.dp[global_sink] += graph.ep.flow[e]
        simp_edge_dict[(graph.vp.id[tgt], graph.vp.id[global_sink])] = e
    return global_src, global_sink


#######################################################################################################
#######################################PATH FINDING ALGORITHM##########################################
#######################################################################################################


def graph_is_DAG_simp(graph: Graph, simp_node_dict: dict):
    """
    check if the graph is a DAG, advanced to check all the isolated subgraph by all mean
    graphtool is_DAG() may not work if the graph is not connected as several parts
    """

    def isCyclicUtil(graph: Graph, v, visited, recStack):
        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        recStack[v] = True

        # Recur for all neighbours
        # if any neighbour is visited and in
        # recStack then graph is cyclic
        for e in v.out_edges():
            neighbour = e.target()
            if not visited[neighbour]:
                if isCyclicUtil(graph, neighbour, visited, recStack):
                    return True
            elif recStack[neighbour]:
                return True

        # The node needs to be poped from
        # recursion stack before function ends
        recStack[v] = False
        return False

    # init
    visited = dict.fromkeys(simp_node_dict.values(), False)
    recStack = dict.fromkeys(simp_node_dict.values(), False)
    for node in simp_node_dict.values():
        visited[node] = False
        recStack[node] = False
    for node in simp_node_dict.values():
        if not visited[node]:
            if isCyclicUtil(graph, node, visited, recStack):
                return False
    return True


def graph_is_DAG(graph: Graph, simp_node_dict: dict):
    """
    check if the graph is a DAG, advanced to check all the isolated subgraph by all mean
    graphtool is_DAG() may not work if the graph is not connected as several parts
    """

    def isCyclicUtil(graph: Graph, v, visited, recStack):
        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        recStack[v] = True

        # Recur for all neighbours
        # if any neighbour is visited and in
        # recStack then graph is cyclic
        for e in v.out_edges():
            if graph.ep.color[e] != "black":
                continue
            neighbour = e.target()
            if not visited[neighbour]:
                if isCyclicUtil(graph, neighbour, visited, recStack):
                    return True
            elif recStack[neighbour]:
                return True

        # The node needs to be poped from
        # recursion stack before function ends
        recStack[v] = False
        return False

    # init
    visited = {}
    recStack = {}
    for node in simp_node_dict.values():
        if graph.vp.color[node] == "black":
            visited[node] = False
            recStack[node] = False
        else:
            visited[node] = True
            recStack[node] = True
    for node in simp_node_dict.values():
        if not visited[node]:
            if isCyclicUtil(graph, node, visited, recStack):
                return False
    return True


def retrieve_cycle(graph: Graph, n=1):
    """
    retrieve a cycle, if any, else return None, sequential
    """
    cycles = []

    def processDFSTree(graph: Graph, stack: list, visited: list, n):
        for out_e in stack[-1].out_edges():
            if graph.ep.color[out_e] != "black":
                continue
            if n == 0:
                return n
            next = out_e.target()
            if visited[next] == "instack":
                # revisit a visited node, cycle
                n -= 1
                store_cycle(stack, next)
            elif visited[next] == "unvisited":
                visited[next] = "instack"
                stack.append(next)
                n = processDFSTree(graph, stack, visited, n)
        visited[stack[-1]] = "done"
        stack.pop()
        return n

    def store_cycle(stack: list, next):
        stack = stack[stack.index(next) :]
        # print("Cycle: ", list_to_string([graph.vp.id[node] for node in stack]))
        cycles.append(stack[:])

    visited = dict.fromkeys(list(graph.vertices()), "unvisited")

    for v in graph.vertices():
        if visited[v] == "unvisited":
            stack = [v]
            visited[v] = "instack"
            n = processDFSTree(graph, stack, visited, n)
            if n == 0:
                break

    return cycles if len(cycles) > 0 else None


def cyclic_to_dag(
    graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, logger: Logger
):
    """
    convert graph to dag by delete minimum coverage not-in-contig edge in cycle until reduced to dag, legacy
    """

    def remove_edge(fst, snd):
        logger.debug(
            "removing edge: {0} -> {1} to reduce a cycle".format(
                graph.vp.id[fst], graph.vp.id[snd]
            )
        )
        e = graph.edge(fst, snd)
        graph.ep.color[e] = "gray"
        removed_edges.append((graph.vp.id[fst], graph.vp.id[snd], graph.ep.overlap[e]))
        return fst, snd

    removed_edges = []
    logger.debug("Turn cyclic graph to dag..")
    if graph_is_DAG(graph, simp_node_dict):
        logger.debug("graph is dag already, skip")
    else:
        while not graph_is_DAG(graph, simp_node_dict):
            cycle = retrieve_cycle(graph)[0]
            max_node = max(cycle, key=lambda v: graph.vp.dp[v])
            prev_node = cycle[(cycle.index(max_node) - 1) % len(cycle)]
            next_node = cycle[(cycle.index(max_node) + 1) % len(cycle)]
            if graph.vp.dp[prev_node] < graph.vp.dp[next_node]:
                remove_edge(prev_node, max_node)
            else:
                remove_edge(max_node, next_node)
    for (uid, vid, _) in removed_edges:
        e = simp_edge_dict.pop((uid, vid))
        graph.remove_edge(e)
    logger.debug("done")
    return removed_edges


def reachable(graph: Graph, src, tgt):
    """
    determine whether src can possibly reach the tgt
    """
    visited = dict.fromkeys(graph.vertices(), False)

    count_down = 1 if src != tgt else 2
    queue = [src]
    while queue:
        curr = queue.pop()
        visited[curr] = True
        if curr == tgt:
            count_down -= 1
            if count_down == 0:
                return True
            else:
                visited[curr] = False

        for oute in curr.out_edges():
            out = oute.target()
            if not visited[out]:
                queue.append(out)
    return False
