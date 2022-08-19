#!/usr/bin/env python3

from graph_tool.all import Graph
from utils.vsAware_Utilities import *

__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, vsAware Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"


def bellman_ford(graph: Graph, src, tgt, logger: Logger):
    logger.debug(
        "bellman ford, s: {0}, t: {1}".format(graph.vp.id[src], graph.vp.id[tgt])
    )
    if src == tgt:
        logger.error("co-incident source-target: " + str(graph.vp.id[tgt]))
        sys.exit(1)

    dist = dict.fromkeys(graph.vertices(), sys.maxsize)
    parent = dict.fromkeys(graph.vertices(), None)
    dist[src] = 0
    for _ in range(1, graph.num_vertices()):
        for e in graph.edges():
            if not graph.ep.keep[e]:
                continue
            u = e.source()
            v = e.target()
            if dist[u] != sys.maxsize and dist[u] + graph.ep.eval[e] < dist[v]:
                dist[v] = dist[u] + graph.ep.eval[e]
                parent[v] = u
    if dist[tgt] == sys.maxsize or parent[tgt] == None:
        logger.error("target is not reachable: " + str(graph.vp.id[tgt]))
        sys.exit(1)
    else:
        for e in graph.edges():
            if not graph.ep.keep[e]:
                continue
            u = e.source()
            v = e.target()
            if (
                dist[u] != sys.maxsize
                and dist[v] != sys.maxsize
                and dist[u] + graph.ep.eval[e] < dist[v]
            ):
                logger.error("negative cycle exist in connected path")
                sys.exit(1)
        path = []
        v = tgt
        itercount = 0
        score = 0
        while v != src:
            itercount += 1
            assert itercount < graph.num_vertices()
            path.append(v)
            score += graph.ep.eval[graph.edge(parent[v], v)]
            v = parent[v]
        path.append(src)
        path.reverse()
        return path, score


def cond_contig_check(
    usage_dict: dict,
    strain_dict: dict,
    cno,
    clen,
    ccov,
    contig,
    min_clen,
    delta,
    threshold,
    logger: Logger,
):
    logger.debug(
        "----------------------------------------------------------------------------"
    )
    logger.debug(
        "Contig: "
        + str(cno)
        + "len: "
        + str(clen)
        + "->bound: ["
        + str(ccov - delta)
        + ","
        + str(ccov)
        + ","
        + str(ccov + delta)
        + "], delta: "
        + str(delta)
    )
    logger.debug("***contig path: " + list_to_string(contig))
    if ccov <= threshold or (all([usage_dict[k] != 0 for k in contig])):
        logger.debug(
            "current contig {0} is used previously {1}".format(ccov, threshold)
        )
        return False
    if clen < min_clen:
        logger.debug(
            "current contig {0} is too short: {1}bp (unreliable) for s-t extension".format(
                cno, clen
            )
        )
        strain_dict["A" + cno] = [contig, clen, ccov]
        return False
    return True


def final_strain_check(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    usage_dict: dict,
    strain_dict: dict,
    strain: list,
    threshold,
    b0,
    b1,
    cno,
    ccov,
    logger: Logger,
):
    if strain == []:
        return
    logger.debug(path_to_id_string(graph, strain, "strain"))
    score = []
    for flow in contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in strain]):
        delta = 3 * abs(b0 + b1 * flow)
        diff = flow - ccov
        if diff < -delta:
            score.append("P4")
        elif diff >= -delta and diff <= delta:
            score.append("P1")
        elif diff > delta and diff <= 2 * delta:
            score.append("P3")
        elif diff > 2 * delta:
            score.append("P2")
    logger.debug("related edges score: " + str(score))
    ccov = path_cov(
        graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in strain]
    )
    plen = path_len(graph, strain)
    logger.debug("len: {0}, ccov: {1}".format(plen, ccov))

    # filter low coverage strain
    if ccov > threshold:
        logger.debug("cand strain found")
        strain_dict["A" + cno] = [[graph.vp.id[n] for n in strain], plen, ccov]
        graph_reduction_c(graph, strain, usage_dict, ccov)
        contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict, None)
    else:
        logger.debug("low cov strain, removed")
    return


def get_similarity(graph: Graph, contig_dict: dict, logger: Logger, b0, b1):
    contig_similarity = {}
    for (cno, [_, _, ccov]) in contig_dict.items():
        delta = 3 * abs(b0 + b1 * ccov)
        lb = ccov - delta
        ub = ccov + delta
        similar_count = 0
        for e in graph.edges():
            if graph.ep.flow[e] >= lb and graph.ep.flow[e] <= ub:
                similar_count += 1
        contig_similarity[cno] = similar_count
    logger.debug("contig histogram generated")
    return contig_similarity


def label_filter_dfs(graph: Graph, contig: list, is_rev=False):
    def dfs(graph: Graph, stack: list, visited: dict):
        if not is_rev:
            for oe in list(stack[-1].out_edges()):
                if not graph.ep.keep[oe]:
                    continue
                next = oe.target()
                if visited[next] == "instack":
                    graph.ep.keep[oe] = False
                    # print(path_to_id_string(graph, stack[stack.index(next):], "cyc"))
                elif visited[next] == "unvisited":
                    visited[next] = "instack"
                    stack.append(next)
                    dfs(graph, stack, visited)
        else:
            for ie in list(stack[-1].in_edges()):
                if not graph.ep.keep[ie]:
                    continue
                prev = ie.source()
                if visited[prev] == "instack":
                    graph.ep.keep[ie] = False
                elif visited[prev] == "unvisited":
                    visited[prev] = "instack"
                    stack.append(prev)
                    dfs(graph, stack, visited)
        visited[stack[-1]] = "done"
        stack.pop()
        return

    s = [v for v in graph.vertices() if graph.vp.id[v] == contig[0]][0]
    t = [v for v in graph.vertices() if graph.vp.id[v] == contig[-1]][0]
    visited = dict.fromkeys(list(graph.vertices()), "unvisited")
    stack = [s] if is_rev else [t]
    visited[stack[0]] = "instack"
    dfs(graph, stack, visited)
    return


def set_filter(graph: Graph, simp_node_dict: dict, contig: list):
    for e in graph.edges():
        graph.ep.keep[e] = True
    # interior edges
    for i in range(1, len(contig) - 1):
        for e in simp_node_dict[contig[i]].all_edges():
            graph.ep.keep[e] = False
    if len(contig) > 1:
        # first contig node edges
        for oe in simp_node_dict[contig[0]].out_edges():
            graph.ep.keep[oe] = False
        # last contig node edges
        for ie in simp_node_dict[contig[-1]].in_edges():
            graph.ep.keep[ie] = False

    # forward cycle
    label_filter_dfs(graph, contig, False)
    # reverse cycle
    label_filter_dfs(graph, contig, True)

    return


def set_edge_weight(graph: Graph, ccov, b0, b1):
    """
    set edge weight for path finding
    """
    for e in graph.edges():
        diff = graph.ep.flow[e] - ccov
        delta = 3 * abs(b0 + b1 * graph.ep.flow[e])
        if diff < -delta:
            # P4 worst
            graph.ep.eval[e] = (-diff) / delta
        elif diff >= -delta and diff <= delta:
            # P1 best
            alen = len(str(graph.vp.id[e.source()]).split("_")) + len(
                str(graph.vp.id[e.target()]).split("_")
            )
            # negative weight, guarantee selection
            graph.ep.eval[e] = -(alen / (abs(diff) + 1)) * delta
        elif diff > delta and diff <= 2 * delta:
            # P3
            graph.ep.eval[e] = (diff - delta) / delta
        else:
            # P2
            graph.ep.eval[e] = 0
    return


def extract_cand_path(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
    b0,
    b1,
    threshold,
    min_clen=600,
):
    logger.info("contig path extending..")
    global_src, global_sink = add_global_source_sink(
        graph, simp_node_dict, simp_edge_dict
    )
    is_dag = graph_is_DAG(graph, simp_node_dict)
    graph.ep.eval = graph.new_edge_property("double")
    graph.ep.keep = graph.new_edge_property("boolean")

    strain_dict = {}
    contig_similarity = {}
    usage_dict = dict.fromkeys(simp_node_dict.keys(), 0)

    self_contig_dict = {}
    linear_contig_dict = {}
    cyclic_contig_dict = {}
    in_tip_contig_dict = {}
    out_tip_contig_dict = {}
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        if ccov <= threshold:
            continue
        s = simp_node_dict[contig[-1]]
        t = simp_node_dict[contig[0]]
        if (
            s in global_sink.in_neighbors() and t in global_src.out_neighbors()
        ) or t in s.out_neighbors():
            self_contig_dict[cno] = contig_dict.pop(cno)
        else:
            if not is_dag:
                if reachable(graph, s, t):
                    cyclic_contig_dict[cno] = contig_dict.pop(cno)
                    logger.debug("cyclic contig - " + str(cno))
                else:
                    cs = simp_node_dict[contig[0]]
                    ct = simp_node_dict[contig[-1]]
                    in_tip = cs != ct and cs in global_src.out_neighbors()
                    out_tip = cs != ct and ct in global_sink.in_neighbors()
                    can_reach = reachable(graph, ct, ct)
                    can_reach_rev = reachable(graph, cs, cs)

                    assert not (
                        can_reach and can_reach_rev
                    )  # impossible case, unless contig is not connected

                    if can_reach:
                        if in_tip:
                            in_tip_contig_dict[cno] = contig_dict.pop(cno)
                            logger.debug("in tip contig - " + str(cno))
                        else:
                            logger.error(
                                list_to_string(
                                    contig,
                                    "forward cycle, but not in tip, error contig",
                                )
                            )
                            sys.exit(1)
                    elif can_reach_rev:
                        if out_tip:
                            out_tip_contig_dict[cno] = contig_dict.pop(cno)
                            logger.debug("out tip contig - " + str(cno))
                        else:
                            logger.error(
                                list_to_string(
                                    contig,
                                    "backward cycle, but not out tip, error contig",
                                )
                            )
                            sys.exit(1)
                    else:
                        linear_contig_dict[cno] = contig_dict.pop(cno)
                        logger.debug("linear contig - " + str(cno))
            else:
                linear_contig_dict[cno] = contig_dict.pop(cno)
                logger.debug("linear contig - " + str(cno))

    # primary combination, self cycle and s-t path contig first
    logger.debug("=====Primary Step: direct path extraction=====")
    for cno, [contig, clen, ccov] in self_contig_dict.items():
        # s-t path contig
        pcov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
        logger.debug(
            "cand strain found, st path/self cycle contig, retrieve first, cno: "
            + cno
            + " ccov: "
            + str(pcov)
            + " clen: "
            + str(clen)
        )
        strain_dict["A" + cno] = [contig, clen, pcov]
        graph_reduction_c(graph, [simp_node_dict[n] for n in contig], usage_dict, pcov)

    logger.debug("=====Secondary Step: interior cyclic contig connection=====")
    forced = dict.fromkeys(cyclic_contig_dict.keys(), False)
    while len(cyclic_contig_dict.keys()) > 0:
        contig_similarity = get_similarity(graph, cyclic_contig_dict, logger, b0, b1)
        sorted_cno = sorted(
            cyclic_contig_dict.keys(),
            key=lambda k: contig_similarity[k],
            reverse=True,
        )

        cno = sorted_cno[0]
        for scno in sorted_cno:
            if not forced[scno]:
                cno = scno
                break

        contig, clen, ccov = cyclic_contig_dict.pop(cno)
        delta = 3 * abs(b0 + b1 * ccov)
        keep = cond_contig_check(
            usage_dict,
            strain_dict,
            cno,
            clen,
            ccov,
            contig,
            min_clen,
            delta,
            threshold,
            logger,
        )
        if not keep:
            continue
        set_edge_weight(graph, ccov, b0, b1)
        # current contig is a internal cyclic contig
        s = simp_node_dict[contig[-1]]
        t = simp_node_dict[contig[0]]
        strain = []
        ns = None
        prev_edges = []
        if len(contig) == 1:
            # split the single node contig, zip it again later
            ns = graph.add_vertex()
            nt = simp_node_dict[contig[0]]
            graph.vp.id[ns] = str(graph.vp.id[nt]) + "temp"
            for e in list(nt.in_edges()):
                prev_edges.append(
                    (
                        e.source(),
                        graph.ep.overlap[e],
                        graph.ep.color[e],
                        graph.ep.flow[e],
                    )
                )
                ie = graph.add_edge(e.source(), ns)
                graph.ep.eval[ie] = graph.ep.eval[e]
                graph.remove_edge(e)
            s = nt
            t = ns
            contig.insert(0, graph.vp.id[ns])
            simp_node_dict[graph.vp.id[ns]] = ns

        set_filter(graph, simp_node_dict, contig)

        sp_vlist, sp_score = bellman_ford(graph, s, t, logger)

        if ns != None:
            # re-zip split node
            for e in list(ns.in_edges()):
                graph.remove_edge(e)
            graph.remove_vertex(ns)
            for (src, o, c, f) in prev_edges:
                ie = graph.add_edge(src, s)
                graph.ep.overlap[ie] = o
                graph.ep.color[ie] = c
                graph.ep.flow[ie] = f
                simp_edge_dict[(graph.vp.id[src], graph.vp.id[s])] = ie
            contig.pop(0)
            simp_node_dict.pop(graph.vp.id[ns])
        logger.debug(
            path_to_id_string(graph, sp_vlist, "found path score: {0}".format(sp_score))
        )

        if sp_score > 0:
            if not forced[cno]:
                logger.debug("current path is positive weighted, not optimal, re-queue")
                forced = True
                continue
            else:
                logger.debug(
                    "current path is positive weighted, not optimal, forced output"
                )
        else:
            # one <= 0 path found, may benefit others
            forced = {k: False for k in forced}

        strain.extend([simp_node_dict[n] for n in contig])
        strain.extend(sp_vlist[1:-1])

        final_strain_check(
            graph,
            simp_node_dict,
            simp_edge_dict,
            cyclic_contig_dict,
            usage_dict,
            strain_dict,
            strain,
            threshold,
            b0,
            b1,
            cno,
            ccov,
            logger,
        )

    logger.debug("=====Secondary Step: in tip contig path extraction=====")
    for cno, [contig, clen, ccov] in in_tip_contig_dict.items():
        # in tip path contig
        pcov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
        logger.debug(
            "cand strain found, in tip contig, cno: "
            + cno
            + " ccov: "
            + str(pcov)
            + " clen: "
            + str(clen)
        )
        strain_dict["A" + cno] = [contig, clen, pcov]
        graph_reduction_c(graph, [simp_node_dict[n] for n in contig], usage_dict, pcov)

    logger.debug("=====Secondary Step: out tip contig path extraction=====")
    for cno, [contig, clen, ccov] in out_tip_contig_dict.items():
        # in tip path contig
        pcov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
        logger.debug(
            "cand strain found, out tip contig, cno: "
            + cno
            + " ccov: "
            + str(pcov)
            + " clen: "
            + str(clen)
        )
        strain_dict["A" + cno] = [contig, clen, pcov]
        graph_reduction_c(graph, [simp_node_dict[n] for n in contig], usage_dict, pcov)

    logger.debug("=====Secondary Step: linear contig path extraction=====")
    forced = dict.fromkeys(linear_contig_dict.keys(), False)
    while len(linear_contig_dict.keys()) > 0:
        contig_similarity = get_similarity(graph, linear_contig_dict, logger, b0, b1)
        sorted_cno = sorted(
            linear_contig_dict.keys(),
            key=lambda k: contig_similarity[k],
            reverse=True,
        )

        cno = sorted_cno[0]
        for scno in sorted_cno:
            if not forced[scno]:
                cno = scno
                break

        contig, clen, ccov = linear_contig_dict.pop(cno)
        delta = 3 * abs(b0 + b1 * ccov)
        keep = cond_contig_check(
            usage_dict,
            strain_dict,
            cno,
            clen,
            ccov,
            contig,
            min_clen,
            delta,
            threshold,
            logger,
        )
        if not keep:
            continue
        set_filter(graph, simp_node_dict, contig)
        set_edge_weight(graph, ccov, b0, b1)
        cs = simp_node_dict[contig[0]]
        ct = simp_node_dict[contig[-1]]

        strain = []
        sphead_vlist, sphead_score = bellman_ford(graph, global_src, cs, logger)
        sptail_vlist, sptail_score = bellman_ford(graph, ct, global_sink, logger)

        logger.debug(
            path_to_id_string(
                graph, sphead_vlist, "found gs-s path score: {0}".format(sphead_score)
            )
        )

        logger.debug(
            path_to_id_string(
                graph, sptail_vlist, "found t-gt path score: {0}".format(sptail_score)
            )
        )

        if sphead_score + sptail_score > 0:
            if not forced[cno]:
                logger.debug("current path is positive weighted, not optimal, re-queue")
                forced = True
                continue
            else:
                logger.debug(
                    "current path is positive weighted, not optimal, forced output"
                )
        else:
            # one <= 0 path found, may benefit others
            forced = {k: False for k in forced}

        strain.extend(sphead_vlist[1:-1])
        strain.extend([simp_node_dict[n] for n in contig])
        strain.extend(sptail_vlist[1:-1])

        final_strain_check(
            graph,
            simp_node_dict,
            simp_edge_dict,
            linear_contig_dict,
            usage_dict,
            strain_dict,
            strain,
            threshold,
            b0,
            b1,
            cno,
            ccov,
            logger,
        )

    for v in [global_sink, global_src]:
        graph.remove_vertex(v)
    logger.info("done")
    return strain_dict
