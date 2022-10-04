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
        logger.warning("co-incident source-target: " + str(graph.vp.id[tgt]))
        return [], 0

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
        logger.warning("target is not reachable: " + str(graph.vp.id[tgt]))
        return [], 0
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
                logger.warning("negative cycle exist in connected path")
                return [], 0
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


def label_filter_dfs(graph: Graph, contig: list, is_rev=False):
    def dfs(graph: Graph, stack: list, visited: dict):
        if not is_rev:
            for oe in list(stack[-1].out_edges()):
                if not graph.ep.keep[oe]:
                    continue
                next = oe.target()
                if visited[next] == "instack":
                    graph.ep.keep[oe] = False
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

def set_edge_weight(graph: Graph, edge_precord: dict, ccov, min_ccov, b0, b1):
    """
    set edge weight for path finding
    """
    status = {}
    delta = 2 * abs(b0 + b1 * ccov)
    worst_score = graph.num_edges()
    for e in graph.edges():
        precord = edge_precord[(graph.vp.id[e.source()], graph.vp.id[e.target()])]
        s = 0
        if precord == 5: # previously set to P1 or P4 and be used, worst case
            graph.ep.eval[e] = worst_score
            s = 5
        else:
            diff = graph.ep.flow[e] - ccov
            if diff < -delta:
                # P4 worst
                graph.ep.eval[e] = 2
                s = 4
            elif diff >= -delta and diff <= delta:
                # P1 best
                graph.ep.eval[e] = -1
                s = 1
            elif diff > delta and diff <= max(min_ccov - abs(b0 + b1 * min_ccov), delta):
                # P3
                graph.ep.eval[e] = 1
                s = 3
            else:
                # P2
                graph.ep.eval[e] = 0
                s = 2
        status[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = s
    return status


def cyc_contig_split(graph: Graph, simp_node_dict: dict, contig: list):
    s = simp_node_dict[contig[-1]]
    t = simp_node_dict[contig[0]]
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
    return s, t, ns, prev_edges

def cyc_contig_rezip(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, prev_edges: list, contig: list, s, ns):
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

def print_curr_contig(cno, clen, contig, ccov, b0, b1, logger: Logger):
    logger.debug(
        "Contig: "
        + str(cno)
        + " len: "
        + str(clen)
        + "->bound: ["
        + str(ccov - 2 * abs(b0 + b1 * ccov))
        + ","
        + str(ccov)
        + ","
        + str(ccov + 2 * abs(b0 + b1 * ccov))
        + "]"
        + ", 1x delta: "
        + str(abs(b0 + b1 * ccov))
        + list_to_string(contig)
    )

def best_aln_score(graph: Graph, cno, clen, ccov, strain, ref_file):
    fname = "temp.fa"
    pafname = "temp_aln.paf"
    subprocess.check_call("echo "" > {0}".format(fname), shell=True)
    with open(fname, "w") as f:
        f.write(">{0}_{1}_{2}\n".format(cno, clen, ccov))
        f.write("{0}\n".format(path_to_seq(graph, strain, cno)))
        f.close()
    minimap_api(ref_file, fname, pafname)
    subprocess.check_call("rm {0}".format(fname), shell=True)
    best_aln = None
    with open(pafname, "r") as paf:
        best_aln = paf.readline()[:-1].split("\t")
        paf.close()
    subprocess.check_call("rm {0}".format(pafname), shell=True)
    return best_aln

def count_amb(graph: Graph, status_dict: dict, s_strain: list, t_strain: list):
    amb_count = 0
    for i in range(len(s_strain)-1):
        u = s_strain[i]
        uid = graph.vp.id[u]
        v = s_strain[i+1]
        vid = graph.vp.id[v]
        if status_dict[(uid, vid)] == 1:
            # P1, check if any ambiguous P1
            for in_e in v.in_edges():
                wid = graph.vp.id[in_e.source()]
                if wid != uid:
                    if status_dict[(wid, vid)] == 1:
                        amb_count += 1
            for out_e in u.out_edges():
                wid = graph.vp.id[out_e.target()]
                if wid != vid:
                    if status_dict[(uid, wid)] == 1:
                        amb_count += 1
    for i in range(len(t_strain)-1):
        u = t_strain[i]
        uid = graph.vp.id[u]
        v = t_strain[i+1]
        vid = graph.vp.id[v]
        if status_dict[(uid, vid)] == 1:
            # P1, check if any ambiguous P1
            for in_e in v.in_edges():
                wid = graph.vp.id[in_e.source()]
                if wid != uid:
                    if status_dict[(wid, vid)] == 1:
                        amb_count += 1
            for out_e in u.out_edges():
                wid = graph.vp.id[out_e.target()]
                if wid != vid:
                    if status_dict[(uid, wid)] == 1:
                        amb_count += 1
    return amb_count

def extract_cand_path(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    logger: Logger,
    b0,
    b1,
    threshold,
    ref_file,
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
    edge_precord = dict.fromkeys(simp_edge_dict.keys(), 0)
    usage_dict = dict.fromkeys(simp_node_dict.keys(), 0)

    for e in graph.edges():
        if graph.ep.flow[e] <= threshold:
            edge_precord[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = 5

    self_contig_dict = {}
    linear_contig_dict = {}
    cyclic_contig_dict = {}
    tip_store = {}
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        s = simp_node_dict[contig[-1]]
        t = simp_node_dict[contig[0]]
        if (
            s in global_sink.in_neighbors() and t in global_src.out_neighbors()
        ) or t in s.out_neighbors():
            self_contig_dict[cno] = contig_dict.pop(cno)
            logger.debug("s-t contig - " + str(cno))
        else:
            if not is_dag:
                if reachable(graph, s, t):
                    cyclic_contig_dict[cno] = contig_dict.pop(cno)
                    logger.debug("cyclic contig - " + str(cno))
                else:
                    ss_reachable = reachable(graph, t, t)
                    tt_reachable = reachable(graph, s, s)
                    if ss_reachable or tt_reachable:
                        tip = []
                        ori = "T" if ss_reachable else "H"
                        for n in contig[::-1] if ss_reachable else contig:
                            if not reachable(graph, simp_node_dict[n], simp_node_dict[n]):
                                if ss_reachable:
                                    tip.insert(0, simp_node_dict[n])
                                else:
                                    tip.append(simp_node_dict[n])
                            else:
                                tip_store[cno] = (tip, ori, simp_node_dict[n])
                                break
                        new_contig = contig[:-len(tip)] if ss_reachable else contig[len(tip):]
                        new_clen = path_len(graph, [simp_node_dict[v] for v in new_contig])                       
                        cyclic_contig_dict[cno] = [new_contig, new_clen, ccov]
                        print_contig(cno, clen, ccov, new_contig, logger, "tip-removed: {0}, convert to cyclic contig".format(tip))
                    else:
                        linear_contig_dict[cno] = contig_dict.pop(cno)
                        logger.debug("linear contig - " + str(cno))
            else:
                linear_contig_dict[cno] = contig_dict.pop(cno)
                logger.debug("linear contig - " + str(cno))

    # primary combination, self cycle and s-t path contig first
    logger.debug("=====Primary Step: direct path extraction=====")
    for cno, [contig, clen, ccov] in sorted(self_contig_dict.items(), key=lambda item: item[1][2]):
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
    while len(cyclic_contig_dict.keys()) > 0:
        for cno, [contig, clen, ccov] in list(cyclic_contig_dict.items()):
            if ccov <= threshold or (all([usage_dict[k] != 0 for k in contig])):
                logger.debug("Remove contig be used previously " + str(cno))
                cyclic_contig_dict.pop(cno)
        cases = {}
        min_ccov = min([cov for [_, _, cov] in cyclic_contig_dict.values()])
        logger.debug("*Find the most optimal contig among all possible cases, current minimum contig cov: " + str(min_ccov))
        for cno, [contig, clen, ccov] in list(cyclic_contig_dict.items()):
            logger.debug("-----------------------------------------------------")
            # find a path
            print_curr_contig(cno, clen, contig, ccov, b0, b1, logger)

            status_dict = set_edge_weight(graph, edge_precord, ccov, min_ccov, b0, b1)
            allPi = {1: set(), 2: set(), 3: set(), 4: set(), 5: set()}
            for (uid, vid), s in status_dict.items():
                allPi[s] = allPi[s].union({uid, vid})             
            logger.debug(list_to_string(list(allPi[1]), "All P1"))
            logger.debug(list_to_string(list(allPi[2]), "All P2"))
            logger.debug(list_to_string(list(allPi[3]), "All P3"))
            logger.debug(list_to_string(list(allPi[4]), "All P4"))
            logger.debug(list_to_string(list(allPi[5]), "All P5"))

            s, t, ns, prev_edges = cyc_contig_split(graph, simp_node_dict, contig)
            set_filter(graph, simp_node_dict, contig)
            sp_vlist, sp_score = bellman_ford(graph, s, t, logger)
            cyc_contig_rezip(graph, simp_node_dict, simp_edge_dict, prev_edges, contig, s, ns)
            
            strain = [simp_node_dict[n] for n in contig] + sp_vlist[1:-1]

            # evaluation
            amb_count = count_amb(graph, status_dict, sp_vlist[1:-1], [])
            raw_score = sum([graph.ep.eval[graph.edge(strain[i], strain[(i+1)%len(strain)])] for i in range(len(strain))])
            score = clen * (raw_score + amb_count)
            best_aln = best_aln_score(graph, cno, clen, ccov, strain, ref_file)
            status = ["P" + str(status_dict[(graph.vp.id[strain[i]], graph.vp.id[strain[(i+1)%len(strain)]])]) for i in range(len(strain))]
            cases[cno] = [contig, clen, ccov, amb_count, score, status, status_dict, strain, int(best_aln[9])/int(best_aln[10]), best_aln]

            logger.debug(
                path_to_id_string(graph, sp_vlist, "found path score: {0}".format(sp_score))
            )
            logger.debug("Ambiguous P1 count: " + str(amb_count))
            logger.debug("raw score: " + str(raw_score))
            logger.debug("total scalled score (+amb count): " + str(score))
            logger.debug("match rate: " + str(int(best_aln[9])/int(best_aln[10])))
            logger.debug("related edges status: " + str(status))
            logger.debug(path_to_id_string(graph, strain, "strain"))

        if len(cases) == 0:
            continue
        (optimal_cno, [_, _, _, optimal_amb_count, optimal_score, optimal_status, optimal_status_dict, optimal_strain, optimal_rate, optimal_aln]) = min(cases.items(), key=lambda t: t[1][4])
        (ref_cno, [_, _, _, ref_amb_count, ref_score, ref_status, ref_status_dict, ref_strain, match_rate, ref_aln]) = max(cases.items(), key=lambda t: t[1][7])
        cyclic_contig_dict.pop(optimal_cno)
        logger.debug("-----------------------------------------------------")
        logger.debug(str(optimal_aln))
        logger.debug("greedy score: " + str(optimal_score))
        logger.debug("greedy match rate: " + str(optimal_rate))
        logger.debug("greedy amb count: " + str(optimal_amb_count))
        logger.debug("greedy related edges status: " + str(optimal_status))
        logger.debug(path_to_id_string(graph, optimal_strain, "greedy optimal strain from contig {0}".format(optimal_cno)))
        logger.debug("-----------------------------------------------------")
        logger.debug(str(ref_aln))
        logger.debug("reference score: " + str(ref_score))
        logger.debug("reference match rate: " + str(match_rate))
        logger.debug("reference amb count: " + str(ref_amb_count))
        logger.debug("reference related edges status: " + str(ref_status))
        logger.debug(path_to_id_string(graph, ref_strain, "reference optimal strain from contig {0} in this iteration".format(ref_cno)))
        
        if optimal_strain == []:
            continue
        pcov = path_cov(
            graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in optimal_strain]
        )
        logger.debug("pcov: {0}".format(pcov))

        # filter low coverage strain
        if pcov > threshold:
            logger.debug("cand strain found")
            if optimal_cno in tip_store:
                (tip, ori, neighbor) = tip_store[optimal_cno]
                logger.debug("current cand strain is coming from a tip contig, maximize extension, tip: {0}".format(tip))
                neighbor_at = optimal_strain.index(neighbor)
                if ori == "T":
                    if neighbor_at == len(optimal_strain) - 1:
                        optimal_strain.extend(tip)
                    else:
                       optimal_strain = optimal_strain[neighbor_at:] + optimal_strain[:neighbor_at] + tip
                else:
                    if neighbor_at == 0:
                        optimal_strain = tip + optimal_strain
                    else:
                        optimal_strain = tip + optimal_strain[neighbor_at:] + optimal_strain[:neighbor_at]
                logger.debug(path_to_id_string(graph, optimal_strain, "cand strain after tip concat"))
                plen = path_len(graph, optimal_strain)             
            strain_dict["A" + optimal_cno] = [[graph.vp.id[n] for n in optimal_strain], plen, pcov]
            graph_reduction_c(graph, optimal_strain, usage_dict, pcov)
            contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict, None)
            
            for i in range(len(optimal_strain)):
                uid, vid = (graph.vp.id[optimal_strain[i]], graph.vp.id[optimal_strain[(i+1)%len(optimal_strain)]])
                if (uid, vid) in edge_precord:
                    edge_precord[(uid, vid)] = max({1: 5, 2: 0, 3: 0, 4: 5, 5: 5}[
                        optimal_status_dict[(uid, vid)]], edge_precord[(uid, vid)])
        else:
            logger.warning("low cov strain, removed")

    logger.debug("=====Secondary Step: linear contig path extraction=====")
    while len(linear_contig_dict.keys()) > 0:
        cases = {}
        min_ccov = min([cov for [_, _, cov] in linear_contig_dict.values()])
        logger.debug("*Find the most optimal contig among all possible cases, current minimum contig cov: " + str(min_ccov))
        for cno, [contig, clen, ccov] in list(linear_contig_dict.items()):
            logger.debug("-----------------------------------------------------")
            # find a path
            print_curr_contig(cno, clen, contig, ccov, b0, b1, logger)

            status_dict = set_edge_weight(graph, edge_precord, ccov, min_ccov, b0, b1)
            allPi = {1: set(), 2: set(), 3: set(), 4: set(), 5: set()}
            for (uid, vid), s in status_dict.items():
                allPi[s] = allPi[s].union({uid, vid})             
            logger.debug(list_to_string(list(allPi[1]), "All P1"))
            logger.debug(list_to_string(list(allPi[2]), "All P2"))
            logger.debug(list_to_string(list(allPi[3]), "All P3"))
            logger.debug(list_to_string(list(allPi[4]), "All P4"))
            logger.debug(list_to_string(list(allPi[5]), "All P5"))

            set_filter(graph, simp_node_dict, contig)
            sphead_vlist, sphead_score = bellman_ford(graph, global_src, simp_node_dict[contig[0]], logger)
            sptail_vlist, sptail_score = bellman_ford(graph, simp_node_dict[contig[-1]], global_sink, logger)

            strain = sphead_vlist[1:-1] + [simp_node_dict[n] for n in contig] + sptail_vlist[1:-1]

            # evaluation
            amb_count = count_amb(graph, status_dict, sphead_vlist[1:-1], sptail_vlist[1:-1])
            raw_score = sum([graph.ep.eval[graph.edge(strain[i], strain[i+1])] for i in range(len(strain) - 1)])
            score = clen * (raw_score + amb_count)
            best_aln = best_aln_score(graph, cno, clen, ccov, strain, ref_file)
            status = ["P" + str(status_dict[(graph.vp.id[strain[i]], graph.vp.id[strain[i+1]])]) for i in range(len(strain) - 1)]
            cases[cno] = [contig, clen, ccov, amb_count, score, status, status_dict, strain, int(best_aln[9])/int(best_aln[10]), best_aln]
            
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
            logger.debug("Ambiguous P1 count: " + str(amb_count))
            logger.debug("raw score: " + str(raw_score))
            logger.debug("total scalled score (+amb count): " + str(score))
            logger.debug("match rate: " + str(int(best_aln[9])/int(best_aln[10])))
            logger.debug("related edges status: " + str(status))
            logger.debug(path_to_id_string(graph, strain, "strain"))
        
        if len(cases) == 0:
            continue
        (optimal_cno, [_, _, _, optimal_amb_count, optimal_score, optimal_status, optimal_status_dict, optimal_strain, optimal_rate, optimal_aln]) = min(cases.items(), key=lambda t: t[1][4])
        (ref_cno, [_, _, _, ref_amb_count, ref_score, ref_status, ref_status_dict, ref_strain, match_rate, ref_aln]) = max(cases.items(), key=lambda t: t[1][8])
        linear_contig_dict.pop(optimal_cno)
        logger.debug("-----------------------------------------------------")
        logger.debug(str(optimal_aln))
        logger.debug("greedy score: " + str(optimal_score))
        logger.debug("greedy match rate: " + str(optimal_rate))
        logger.debug("greedy amb count: " + str(optimal_amb_count))
        logger.debug("greedy related edges status: " + str(optimal_status))
        logger.debug(path_to_id_string(graph, optimal_strain, "greedy optimal strain from contig {0}".format(optimal_cno)))
        logger.debug("-----------------------------------------------------")
        logger.debug(str(ref_aln))
        logger.debug("reference score: " + str(ref_score))
        logger.debug("reference match rate: " + str(match_rate))
        logger.debug("reference amb count: " + str(ref_amb_count))
        logger.debug("reference related edges status: " + str(ref_status))
        logger.debug(path_to_id_string(graph, ref_strain, "reference optimal strain from contig {0} in this iteration".format(ref_cno)))
        
        if optimal_strain == []:
            continue
        pcov = path_cov(
            graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in ref_strain]
        )
        plen = path_len(graph, ref_strain)
        logger.debug("plen: {0}, pcov: {1}".format(plen, pcov))

        # filter low coverage strain
        if ccov > threshold:
            logger.debug("cand strain found")
            strain_dict["A" + ref_cno] = [[graph.vp.id[n] for n in ref_strain], plen, pcov]
            graph_reduction_c(graph, ref_strain, usage_dict, pcov)
            contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict, None)
            
            for i in range(len(ref_strain) - 1):
                uid, vid = (graph.vp.id[ref_strain[i]], graph.vp.id[ref_strain[i+1]])
                edge_precord[(uid, vid)] = max({1: 5, 2: 0, 3: 0, 4: 5, 5: 5}[
                    optimal_status_dict[(uid, vid)]], edge_precord[(uid, vid)])
        else:
            logger.warning("low cov strain, removed")

        # for cno, [contig, clen, ccov] in list(linear_contig_dict.items()):
        #     if ccov <= threshold or (all([usage_dict[k] != 0 for k in contig])):
        #         logger.debug("Remove contig be used previously " + str(cno))
        #         linear_contig_dict.pop(cno)

    for v in [global_sink, global_src]:
        graph.remove_vertex(v)

    for sno, [strain, slen, scov] in strain_dict.items():
        print_contig(sno, slen, scov, strain, logger, "Strain")
    
    logger.info("done")
    return strain_dict