#!/usr/bin/env python3

from graph_tool.all import Graph
import numpy

from utils.vsAware_Utilities import *
from utils.vsAware_CovBalance import assign_edge_flow, coverage_rebalance_s
from utils.vsAware_IO import graph_to_gfa, flipped_gfa_to_graph
from utils.vsAware_IO import contig_dict_to_fasta, contig_dict_to_path

__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, vsAware Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"


def iterated_graph_split(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
    tempdir,
    b0,
    b1,
    threshold,
    ref_file,
):
    BOUND_ITER = len(simp_node_dict) ** 2
    it = 0

    grapha = graph
    simp_node_dicta = simp_node_dict
    simp_edge_dicta = simp_edge_dict
    total_removed_branch_nt = 0
    total_removed_branch_t = 0
    iterCount = "A"
    num_split = 0
    trivial_split_count = 0
    while it < BOUND_ITER:
        it += 1
        # trivial branch split
        prev_ids = list(simp_node_dicta.keys())
        trivial_split_count, id_mapping = graph_split_trivial(
            grapha, simp_node_dicta, simp_edge_dicta, logger
        )

        graph_to_gfa(
            grapha,
            simp_node_dicta,
            simp_edge_dicta,
            logger,
            "{0}/gfa/split_graph_L{1}1.gfa".format(tempdir, iterCount),
        )
        graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph(
            "{0}/gfa/split_graph_L{1}1.gfa".format(tempdir, iterCount), logger
        )
        assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)

        # fix trivial splitted contig
        contig_dict_remapping(
            graphb,
            simp_node_dictb,
            simp_edge_dictb,
            contig_dict,
            id_mapping,
            prev_ids,
            logger,
        )

        contig_dict_to_path(contig_dict, "{0}/tmp/inte_contigs.paths".format(tempdir))
        contig_dict_to_fasta(
            graphb,
            simp_node_dictb,
            contig_dict,
            "{0}/tmp/inte_contigs.fasta".format(tempdir),
        )
        minimap_api(
            ref_file,
            "{0}/tmp/inte_contigs.fasta".format(tempdir),
            "{0}/paf/inte_contigs_to_strain.paf".format(tempdir),
        )
        map_ref_to_contig(
            contig_dict, logger, "{0}/paf/inte_contigs_to_strain.paf".format(tempdir)
        )

        # non-trivial branch split
        num_split = graph_splitting(
            graphb,
            simp_node_dictb,
            simp_edge_dictb,
            contig_dict,
            logger,
            b0,
            b1,
            threshold,
            ref_file
        )
        graph_to_gfa(
            graphb,
            simp_node_dictb,
            simp_edge_dictb,
            logger,
            "{0}/gfa/split_graph_L{1}2.gfa".format(tempdir, iterCount),
        )
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph(
            "{0}/gfa/split_graph_L{1}2.gfa".format(tempdir, iterCount), logger
        )

        simp_path_compactification(
            graphc, simp_node_dictc, simp_edge_dictc, contig_dict, None, None, logger
        )

        graph_to_gfa(
            graphc,
            simp_node_dictc,
            simp_edge_dictc,
            logger,
            "{0}/gfa/split_graph_L{1}3.gfa".format(tempdir, iterCount),
        )
        grapha, simp_node_dicta, simp_edge_dicta = flipped_gfa_to_graph(
            "{0}/gfa/split_graph_L{1}3.gfa".format(tempdir, iterCount), logger
        )
        assign_edge_flow(grapha, simp_node_dicta, simp_edge_dicta)

        if trivial_split_count == None:
            logger.warning("trigger iteration split process edge case, skip")
            break

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            iterCount = chr(ord(iterCount) + 1)
            contig_cov_fix(grapha, simp_node_dicta, simp_edge_dicta, contig_dict, None)
        else:
            break
    logger.debug(
        "Total non-trivial branches removed: "
        + str(total_removed_branch_nt)
        + " total trivial branches removed: "
        + str(total_removed_branch_t)
    )
    coverage_rebalance_s(grapha, simp_node_dicta, simp_edge_dicta, logger)
    graph_to_gfa(
        grapha,
        simp_node_dicta,
        simp_edge_dicta,
        logger,
        "{0}/gfa/final_graph.gfa".format(tempdir),
    )
    grapho, simp_node_dicto, simp_edge_dicto = flipped_gfa_to_graph(
        "{0}/gfa/final_graph.gfa".format(tempdir), logger
    )
    assign_edge_flow(grapho, simp_node_dicto, simp_edge_dicto)

    # remove duplicated contigs
    # contig_dup_removed_s(contig_dict, logger)
    # concat overlapped contigs
    concat_overlap_contig(grapho, simp_node_dicto, simp_edge_dicto, contig_dict, logger)
    trim_contig_dict(grapho, simp_node_dicto, contig_dict, logger)
    contig_cov_fix(grapho, simp_node_dicto, simp_edge_dicto, contig_dict, logger)

    return grapho, simp_node_dicto, simp_edge_dicto


def graph_split_trivial(
    graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, logger: Logger
):
    """
    Split the graph, for any (0|1)->N, N->(0|1) branch, split by forking the 1 edge to N edge.
    """
    logger.info("graph trivial split..")

    BOUND_ITER = len(simp_node_dict) ** 2
    has_split = True
    trivial_split_count = 0
    id_mapping = {}
    for id in simp_node_dict.keys():
        id_mapping[id] = set()
    while has_split and trivial_split_count < BOUND_ITER:
        has_split = False
        for id in list(simp_node_dict.keys()):
            node = simp_node_dict[id]
            if graph.vp.color[node] != "black":
                continue
            if id not in id_mapping:
                id_mapping[id] = set()
            ines = [ue for ue in node.in_edges() if graph.ep.color[ue] == "black"]
            outes = [ve for ve in node.out_edges() if graph.ep.color[ve] == "black"]
            # if len(ines) == 0 and len(outes) > 1:
            #     logger.debug(id + " current node is source node, split left")
            #     graph.vp.color[node] = "gray"
            #     # create len(outes) subnodes
            #     s = "A"
            #     for i in range(len(outes)):
            #         oute = outes[i]
            #         tgt = oute.target()
            #         snode = graph_add_vertex(
            #             graph,
            #             simp_node_dict,
            #             id + "*" + chr(ord(s) + i),
            #             graph.ep.flow[oute],
            #             graph.vp.seq[node],
            #         )
            #         graph.ep.color[oute] = "gray"
            #         sedge = graph_add_edge(
            #             graph,
            #             simp_edge_dict,
            #             snode,
            #             tgt,
            #             graph.ep.overlap[oute],
            #             graph.ep.flow[oute],
            #         )
            #         simp_node_dict[graph.vp.id[snode]] = snode
            #         simp_edge_dict[
            #             (graph.vp.id[sedge.source()], graph.vp.id[sedge.target()])
            #         ] = sedge
            #         id_mapping[id].add(graph.vp.id[snode])
            #     has_split = True
            #     trivial_split_count += 1
            # elif len(ines) > 1 and len(outes) == 0:
            #     logger.debug(id + " current node is sink node, split right")
            #     graph.vp.color[node] = "gray"
            #     # create len(ines) subnodes
            #     s = "A"
            #     for i in range(len(ines)):
            #         ine = ines[i]
            #         src = ine.source()
            #         snode = graph_add_vertex(
            #             graph,
            #             simp_node_dict,
            #             id + "*" + chr(ord(s) + i),
            #             graph.ep.flow[ine],
            #             graph.vp.seq[node],
            #         )
            #         graph.ep.color[ine] = "gray"
            #         sedge = graph_add_edge(
            #             graph,
            #             simp_edge_dict,
            #             src,
            #             snode,
            #             graph.ep.overlap[ine],
            #             graph.ep.flow[ine],
            #         )
            #         simp_node_dict[graph.vp.id[snode]] = snode
            #         simp_edge_dict[
            #             (graph.vp.id[sedge.source()], graph.vp.id[sedge.target()])
            #         ] = sedge
            #         id_mapping[id].add(graph.vp.id[snode])
            #     has_split = True
            #     trivial_split_count += 1
            # el
            if len(ines) == 1 and len(outes) > 1:
                logger.debug(id + " split left")
                graph.vp.color[node] = "gray"
                ine = ines[0]
                src = ine.source()
                graph.ep.color[ine] = "gray"
                s = "A"
                for i in range(len(outes)):
                    oute = outes[i]
                    tgt = oute.target()
                    snode = graph_add_vertex(
                        graph,
                        simp_node_dict,
                        id + "*" + chr(ord(s) + i),
                        graph.ep.flow[oute],
                        graph.vp.seq[node],
                    )
                    graph.ep.color[oute] = "gray"
                    sedge_out = graph_add_edge(
                        graph,
                        simp_edge_dict,
                        snode,
                        tgt,
                        graph.ep.overlap[oute],
                        graph.ep.flow[oute],
                    )
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[
                        (
                            graph.vp.id[sedge_out.source()],
                            graph.vp.id[sedge_out.target()],
                        )
                    ] = sedge_out

                    sedge_in = graph_add_edge(
                        graph,
                        simp_edge_dict,
                        src,
                        snode,
                        graph.ep.overlap[ine],
                        graph.ep.flow[oute],
                    )
                    simp_edge_dict[
                        (graph.vp.id[sedge_in.source()], graph.vp.id[sedge_in.target()])
                    ] = sedge_in
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            elif len(ines) > 1 and len(outes) == 1:
                logger.debug(id + " split right")
                graph.vp.color[node] = "gray"
                oute = outes[0]
                tgt = oute.target()
                graph.ep.color[oute] = "gray"
                s = "A"
                for i in range(len(ines)):
                    ine = ines[i]
                    src = ine.source()
                    snode = graph_add_vertex(
                        graph,
                        simp_node_dict,
                        id + "*" + chr(ord(s) + i),
                        graph.ep.flow[ine],
                        graph.vp.seq[node],
                    )
                    graph.ep.color[ine] = "gray"
                    sedge_in = graph_add_edge(
                        graph,
                        simp_edge_dict,
                        src,
                        snode,
                        graph.ep.overlap[ine],
                        graph.ep.flow[ine],
                    )
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[
                        (graph.vp.id[sedge_in.source()], graph.vp.id[sedge_in.target()])
                    ] = sedge_in

                    sedge_out = graph_add_edge(
                        graph,
                        simp_edge_dict,
                        snode,
                        tgt,
                        graph.ep.overlap[oute],
                        graph.ep.flow[ine],
                    )
                    simp_edge_dict[
                        (
                            graph.vp.id[sedge_out.source()],
                            graph.vp.id[sedge_out.target()],
                        )
                    ] = sedge_out
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            else:
                None
    if trivial_split_count >= BOUND_ITER:
        logger.warning("Strange topology detected, exit trivial split immediately")
        return None, id_mapping
    else:
        logger.debug("No of trivial branch be removed: " + str(trivial_split_count))
        logger.info("done")
        return trivial_split_count, id_mapping

def best_aln_score(graph: Graph, ori, strain, ref_file, temp_dir):
    fname = "{0}/temp_{1}.fa".format(temp_dir, ori)
    pafname = "{0}/temp_{1}_aln.paf".format(temp_dir, ori)
    subprocess.check_call("echo \"\" > {0}".format(fname), shell=True)
    with open(fname, "w") as f:
        f.write(">{0}\n".format(ori))
        f.write("{0}\n".format(path_to_seq(graph, strain, "")))
        f.close()
    minimap_api(ref_file, fname, pafname)
    subprocess.check_call("rm {0}".format(fname), shell=True)
    best_aln = []
    with open(pafname, "r") as paf:
        for line in paf.readlines():
            splited = line[:-1].split("\t")
            if len(splited) < 12:
                continue
            best_aln.append([splited[0], int(splited[10]), splited[5], int(splited[10]) - int(splited[9])])
        paf.close()
    subprocess.check_call("rm {0}".format(pafname), shell=True)
    return best_aln

def graph_splitting(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
    b0,
    b1,
    threshold,
    ref_file
):
    """
    n-n branch splitting, # FIXME add more restrict rule to avoid false positive split
    """
    logger.info("graph non-trivial split..")
    split_branches = []
    node_to_contig_dict, _ = contig_map_node(contig_dict)
    no_mapping = {}
    for no, node in list(simp_node_dict.items()):
        ine = [e for e in node.in_edges() if graph.ep.color[e] == "black"]
        oute = [e for e in node.out_edges() if graph.ep.color[e] == "black"]
        support_contigs = []
        # non-trivial branch
        if len(ine) > 1 and len(oute) > 1:
            threshold2 = 2 * abs(b0 + b1 * graph.vp.dp[node])
            support_contigs = (
                [cno for cno in node_to_contig_dict[no]]
                if no in node_to_contig_dict
                else []
            )
            print_vertex(
                graph,
                node,
                logger,
                "---------- branch node, support by contig {0}".format(support_contigs),
            )
            logger.debug("local threshold: " + str(threshold2))
            ine_usage = {}
            for ie in ine:
                ine_usage[ie] = 0
            oute_usage = {}
            for oe in oute:
                oute_usage[oe] = 0
            cproduct = []
            for ie in ine:
                for oe in oute:
                    delta = abs(graph.ep.flow[ie] - graph.ep.flow[oe])
                    if delta <= threshold2:
                        cproduct.append(
                            (
                                ie,
                                oe,
                                numpy.mean([graph.ep.flow[ie], graph.ep.flow[oe]]),
                                delta,
                            )
                        )
                        ine_usage[ie] += 1
                        oute_usage[oe] += 1

            if cproduct != []:
                for i, (ie, oe, subcov, delta) in enumerate(cproduct):
                    prev_node = ie.source()
                    prev_no = graph.vp.id[prev_node]
                    next_node = oe.target()
                    next_no = graph.vp.id[next_node]
                    if (prev_no, no) not in simp_edge_dict or (
                        no,
                        next_no,
                    ) not in simp_edge_dict:
                        # logger.debug("edge has been removed already")
                        continue
                    logger.debug("Delta: " + str(delta))
                    involved_contigs = []
                    is_conf = False
                    cross_talk = False
                    for cno in support_contigs:
                        contig, clen, ccov = contig_dict[cno]
                        if len(contig) < 3:
                            # logger.debug("unnecessary contig (len < 3) check: ", contig)
                            if prev_no in contig or next_no in contig:
                                involved_contigs.append(cno)
                        elif prev_no in contig and next_no in contig:
                            is_conf = True
                            involved_contigs.append(cno)
                            print_contig(cno, clen, ccov, contig, logger, "support contig ")
                        elif prev_no in contig and next_no not in contig:
                            logger.debug("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        elif prev_no not in contig and next_no in contig:
                            logger.debug("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        else:
                            None
                    if cross_talk:
                        None
                        logger.debug("- current branch split forbidden, cross talk {0} -> {1} -> {2}".format(prev_no, no, next_no))
                    elif not is_conf and (ine_usage[ie] > 1 or oute_usage[oe] > 1):
                        None
                        logger.debug("- current branch split forbidden, no supporting contig path, ambigous split {0} -> {1} -> {2}".format(prev_no, no, next_no))
                    else:
                        logger.debug("- branch split performed - {0}:{1} - {2} - {3}:{4} - (is_conf: {5})".format(prev_no, graph.ep.flow[ie], no, next_no, graph.ep.flow[oe], is_conf))
                        # map to reference to check correctness
                        ref_l = best_aln_score(graph, "L", [prev_node], ref_file, ".")
                        ref_r = best_aln_score(graph, "R", [next_node], ref_file, ".")
                        best_aln_l = [ref for [_, l, ref, nm] in ref_l if nm == 0 and l == len(graph.vp.seq[prev_node])]
                        best_aln_r = [ref for [_, l, ref, nm] in ref_r if nm == 0 and l == len(graph.vp.seq[next_node])]
                        its_refs = set(best_aln_l).intersection(set(best_aln_r))
                        logger.debug("aln left: " + str(best_aln_l))
                        logger.debug("aln right: " + str(best_aln_r))
                        if len(its_refs) != 0:
                            logger.debug("perfect split")
                        else:
                            logger.debug("potential incorrect split")

                        split_branches.append(no)

                        subid = no + "*" + str(i)
                        sub_node = graph_add_vertex(
                            graph, simp_node_dict, subid, subcov, graph.vp.seq[node]
                        )

                        if no not in no_mapping:
                            no_mapping[no] = []
                        no_mapping[no].append(str(subid))

                        graph.vp.dp[node] -= subcov

                        graph_remove_edge(graph, simp_edge_dict, prev_no, no)
                        graph_remove_edge(graph, simp_edge_dict, no, next_no)
                        graph_add_edge(
                            graph,
                            simp_edge_dict,
                            prev_node,
                            sub_node,
                            graph.ep.overlap[ie],
                            graph.ep.flow[ie],
                        )
                        graph_add_edge(
                            graph,
                            simp_edge_dict,
                            sub_node,
                            next_node,
                            graph.ep.overlap[oe],
                            graph.ep.flow[oe],
                        )

                        for icno in involved_contigs:
                            if contig_dict[icno][0].count(no) == 1:
                                contig_dict[icno][0][
                                    contig_dict[icno][0].index(no)
                                ] = subid
                                support_contigs.remove(icno)
                            else:
                                logger.warning(
                                    "contig error, previous node {0} is not in contig or multiple occurance in contig {1}, potential bug".format(
                                        no, icno
                                    )
                                )

    # fix single node contig
    for cno, [contig, _, _] in list(contig_dict.items()):
        if len(contig) <= 1 and contig[0] not in simp_node_dict:
            contig_dict.pop(cno)
            # logger.debug("isolated contig node {0}, prepared to pop, check any replacement".format(cno))
            if contig[0] in no_mapping:
                # logger.debug("mapping: {0} -> {1}".format(contig[0], no_mapping[contig[0]]))
                s = "A"
                for i, mapped in enumerate(no_mapping[contig[0]]):
                    if mapped in simp_node_dict:
                        mapped_node = simp_node_dict[mapped]
                        contig_dict[cno + chr(ord(s) + i)] = [
                            [mapped],
                            path_len(graph, [mapped_node]),
                            graph.vp.dp[mapped_node],
                        ]
                        # logger.debug("mapped cno: ", cno+chr(ord(s) + i), mapped)

    # remove all the isolated low cov node&edge not in contig
    node_to_contig_dict, _ = contig_map_node(contig_dict)
    for node in list(graph.vertices()):
        if (
            graph.vp.id[node] not in node_to_contig_dict
            and graph.vp.dp[node] <= threshold
        ):
            alle = [e for e in set(node.all_edges()) if graph.ep.color[e] == "black"]
            graph_remove_vertex(
                graph, simp_node_dict, graph.vp.id[node], "remove isolated low cov node"
            )
            for e in alle:
                graph_remove_edge(
                    graph,
                    simp_edge_dict,
                    graph.vp.id[e.source()],
                    graph.vp.id[e.target()],
                )

    logger.debug("No of branch be removed: " + str(len(set(split_branches))))
    logger.debug("Split branches: " + list_to_string(set(split_branches)))
    logger.info("done")
    return len(set(split_branches))
