from utils.VStrains_Utilities import *
from utils.VStrains_IO import store_reinit_graph
import matplotlib.pyplot as plt
import numpy


__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, VStrains Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"


def link_split(
    sec_comb: list,
    kept_link: dict,
    in_usage: dict,
    in_capacity: dict,
    out_usage: dict,
    out_capacity: dict,
    logger,
):
    """update split plan using paired end & single end information"""
    logger.debug("attempt to split via paired end information")
    sorted_sec_comb = sorted(sec_comb, key=lambda x: x[2], reverse=True)
    for (uid, wid, pe) in sorted_sec_comb:
        if pe <= 0:
            break
        logger.debug("-----SEC LINK {0} -> {1} PE: {2}".format(uid, wid, pe))
        logger.debug("Capacity: {0} -> {1}".format(in_capacity[uid], out_capacity[wid]))
        logger.debug("- distinct compatiable case, added")
        in_usage[uid] += 1
        out_usage[wid] += 1
        kept_link[(uid, wid)] = ((in_capacity[uid] + out_capacity[wid]) / 2, pe)
    return


def cov_split(
    us: list,
    ws: list,
    pe_info: dict,
    sec_comb: list,
    kept_link: dict,
    in_usage: dict,
    in_capacity: dict,
    out_usage: dict,
    out_capacity: dict,
    logger,
):
    """update split plan using coverage information"""
    logger.debug("attempt to split via coverage information")
    logger.debug(
        "align paired end/single end information first (if any) to isolated nodes"
    )
    sorted_sec_comb = sorted(sec_comb, key=lambda x: x[2], reverse=True)
    for (uid, wid, pe) in sorted_sec_comb:
        if pe <= 0:
            break
        if in_usage[uid] > 0 or out_usage[wid] > 0:
            continue
        logger.debug("-----SEC LINK {0} -> {1} PE: {2}-----".format(uid, wid, pe))
        logger.debug("Capacity: {0} -> {1}".format(in_capacity[uid], out_capacity[wid]))
        logger.debug("- link [ > 0] supported case, added")
        in_usage[uid] += 1
        out_usage[wid] += 1
        kept_link[(uid, wid)] = ((in_capacity[uid] + out_capacity[wid]) / 2, pe)

    logger.debug("obtain best match via coverage similarity")
    for uid in us:
        if in_usage[uid] > 0:
            continue
        opt_ws = sorted(ws, key=lambda wwid: abs(in_capacity[uid] - out_capacity[wwid]))
        wid = opt_ws[0]
        opt_us = sorted(us, key=lambda uuid: abs(in_capacity[uuid] - out_capacity[wid]))
        if opt_us[0] == uid and out_usage[wid] == 0 and (uid, wid) not in kept_link:
            delta = 2 * abs(in_capacity[uid] - out_capacity[wid])
            logger.debug(
                "Found coverage best match: {0} -> {1} with cov: {2}, {3}, checking delta bound: {4}".format(
                    uid, wid, in_capacity[uid], out_capacity[wid], delta
                )
            )
            if (
                abs(in_capacity[opt_us[1]] - out_capacity[wid]) <= delta
                or abs(in_capacity[uid] - out_capacity[opt_ws[1]]) <= delta
            ):
                logger.debug("ambiguous matching, skip")
            else:
                logger.debug("added")
                in_usage[uid] += 1
                out_usage[wid] += 1
                kept_link[(uid, wid)] = (
                    (in_capacity[uid] + out_capacity[wid]) / 2,
                    pe_info[(min(uid, wid), max(uid, wid))],
                )
    return


def balance_split(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    pe_info: dict,
    logger: Logger,
    ref_file: str,
    temp_dir: str,
    count_id: int,
    threshold,
    is_prim: bool,
):
    logger.info(
        "balance split using contigs&paired end links&coverage information.. isPrim: {0}".format(
            is_prim
        )
    )
    correct_X = []
    correct_Y = []
    false_error_X = []
    false_error_Y = []
    error_X = []
    error_Y = []
    error_text = []
    cut = 100

    # detect all non-trivial branches right now
    non_trivial_branches = get_non_trivial_branches(graph, simp_node_dict)
    split_branches = []
    node_to_contig_dict, _ = contig_map_node(contig_dict)
    for no, node in non_trivial_branches.items():
        us = [
            graph.vp.id[e.source()]
            for e in node.in_edges()
            if graph.ep.color[e] == "black"
        ]
        ws = [
            graph.vp.id[e.target()]
            for e in node.out_edges()
            if graph.ep.color[e] == "black"
        ]
        logger.debug("---------------------------------------------")
        logger.debug(
            "current non trivial branch: {0}, in-degree: {1}, out-degree: {2}".format(
                no, len(us), len(ws)
            )
        )

        # authenticate if split-able
        if any([pe_info[(uid, uid)] == None for uid in us]) or any(
            [pe_info[(wid, wid)] == None for wid in ws]
        ):
            logger.debug(
                "current non-trivial branch: {0} is related to current iteration, split later".format(
                    no
                )
            )
            continue
        if not is_non_trivial(graph, node):
            logger.debug(
                "current non-trivial branch: {0} is not non-trivial, potential bug".format(
                    no
                )
            )
            continue
        if len(us) != len(ws):
            logger.debug("Not N-N split, skip")
            continue

        # check if link-split
        split_via_link = True

        # not perform link-split if any leaf is from a splitted node
        for id in us + ws:
            singles = id.split("&")
            if all([single.count("*") > 0 for single in singles]):
                logger.debug(
                    "leaf:{0} is total branch nodes, no link information, skip link split".format(
                        id
                    )
                )
                split_via_link = False
                break

        # not perform link-split if no combination has link information
        if all(
            [pe_info[(min(uid, wid), max(uid, wid))] == 0 for uid in us for wid in ws]
        ):
            logger.debug(
                "current branch node too long, no link information, skip link split"
            )
            split_via_link = False

        # add contig supports
        support_contigs = node_to_contig_dict.get(no, [])
        con_info = {}
        for cno in support_contigs:
            [contig, clen, ccov] = contig_dict[cno]
            loc = contig.index(no)
            if loc > 0 and loc < len(contig) - 1:
                con_info[(contig[loc - 1], contig[loc + 1])] = con_info.get(
                    (contig[loc - 1], contig[loc + 1]), []
                )
                con_info[(contig[loc - 1], contig[loc + 1])].append((cno, clen, ccov))
            print_contig(
                cno,
                clen,
                round(ccov, 2),
                contig[max(loc - 1, 0) : loc + 2],
                logger,
                "support contig",
            )

        # debug only
        # obtain perfect split via reference
        expect_link = []
        ref_pair_dict = {}
        ref_all_dict = {}
        if ref_file:
            lrefs = set()
            rrefs = set()
            error_nos = set()
            for uid in us:
                for wid in ws:
                    u = simp_node_dict[uid]
                    w = simp_node_dict[wid]
                    ref_l = best_aln_score(graph, "L", [u], ref_file, temp_dir)
                    best_ref_l = [
                        ref
                        for [_, l, ref, nm] in ref_l
                        if nm == 0 and l == len(graph.vp.seq[u])
                    ]
                    ref_r = best_aln_score(graph, "R", [w], ref_file, temp_dir)
                    best_ref_r = [
                        ref
                        for [_, l, ref, nm] in ref_r
                        if nm == 0 and l == len(graph.vp.seq[w])
                    ]
                    lrefs = lrefs.union(best_ref_l)
                    rrefs = rrefs.union(best_ref_r)
                    ref_pair_dict[(uid, wid)] = set(best_ref_l).intersection(
                        set(best_ref_r)
                    )
                    ref_all_dict[(uid, wid)] = set(
                        [ref for [_, _, ref, nm] in ref_l if nm < 5]
                    ).union(set([ref for [_, _, ref, nm] in ref_r if nm < 5]))
                    if len(ref_pair_dict[(uid, wid)]) > 0:
                        expect_link.append((uid, wid))
                    if len(best_ref_l) == 0:
                        error_nos.add(uid)
                    if len(best_ref_r) == 0:
                        error_nos.add(wid)
            sym_diff = lrefs.symmetric_difference(rrefs)
            if len(sym_diff) > 0:
                logger.debug(
                    "Current branch have force mismatch connection for following strains: {0}".format(
                        sym_diff
                    )
                )
        # debug only

        kept_link = {}
        sec_comb = []
        # init node usage for current branch
        in_usage = dict.fromkeys(us, 0)
        in_capacity = {}
        for uid in us:
            in_capacity[uid] = graph.ep.flow[simp_edge_dict[(uid, no)]]

        out_usage = dict.fromkeys(ws, 0)
        out_capacity = {}
        for wid in ws:
            out_capacity[wid] = graph.ep.flow[simp_edge_dict[(no, wid)]]

        # align contig link first, and update status
        logger.debug("align contig link first")
        for uid in us:
            for wid in ws:
                logger.debug("---------------------")
                u = simp_node_dict[uid]
                w = simp_node_dict[wid]
                curr_pe = pe_info[(min(uid, wid), max(uid, wid))]

                logger.debug("{0} -> {1} PE: {2}".format(uid, wid, curr_pe))
                logger.debug(
                    "cov info: {0}[{1}] -> {2}[{3}]".format(
                        graph.ep.flow[graph.edge(u, node)],
                        pe_info[(min(uid, no), max(uid, no))],
                        graph.ep.flow[graph.edge(node, w)],
                        pe_info[(min(no, wid), max(no, wid))],
                    )
                )
                if ref_file:
                    logger.debug(
                        "intersect reference: {0}".format(ref_pair_dict[(uid, wid)])
                    )
                    # potential incorrect matching, but supported by links
                    if len(ref_pair_dict[(uid, wid)]) == 0 and curr_pe > 0:
                        logger.debug("False Positive case, WARN")
                accept = False
                if (uid, wid) in con_info:
                    logger.debug(
                        "current link supported by contig: {0}, added".format(
                            con_info[(uid, wid)]
                        )
                    )
                    accept = True
                if uid == wid:
                    logger.debug(
                        "current link is a self link: {0}, potential cyclic strain, added".format(
                            uid
                        )
                    )
                    accept = True

                if accept:
                    in_usage[uid] += 1
                    out_usage[wid] += 1
                    kept_link[(uid, wid)] = (
                        (in_capacity[uid] + out_capacity[wid]) / 2,
                        curr_pe,
                    )
                else:
                    logger.debug("current link is secondary choice, process later")
                    sec_comb.append((uid, wid, curr_pe))
        if is_prim:
            if split_via_link:
                link_split(
                    sec_comb,
                    kept_link,
                    in_usage,
                    in_capacity,
                    out_usage,
                    out_capacity,
                    logger,
                )
        else:
            # secondary split, via link first, then coverage
            cov_split(
                us,
                ws,
                pe_info,
                sec_comb,
                kept_link,
                in_usage,
                in_capacity,
                out_usage,
                out_capacity,
                logger,
            )
        if not (
            all([u == 1 for u in in_usage.values()])
            and all([v == 1 for v in out_usage.values()])
        ):
            logger.debug("->Not satisfy N-N split, skip: {0}".format(kept_link))
            continue
        worst_pair_diff = max(
            [
                abs(in_capacity[uid] - out_capacity[wid])
                for (uid, wid) in kept_link.keys()
            ]
        )
        if worst_pair_diff > 4 * threshold:
            logger.debug(
                "worst pair coverage diff greater than 4 delta: {0} > {1}, too uneven, skip: {2}".format(
                    worst_pair_diff, 4 * threshold, kept_link
                )
            )
            continue
        logger.debug("->perform split, all kept links: {0}".format(kept_link))
        if ref_file:
            logger.debug("->expected links: {0}".format(expect_link))
            if set(kept_link) != set(expect_link):
                logger.debug("Incorrect split")
            else:
                logger.debug("Correct split")

        split_branches.append(no)
        link2subs = {}
        counter = 0
        for (uid, wid), (sub_flow, pe) in kept_link.items():
            logger.debug("--------> {0} - {1}".format(uid, wid))
            # debug only
            if ref_file:
                if len(ref_pair_dict[(uid, wid)]) != 0:
                    logger.debug("best pair")
                    if pe <= cut:
                        correct_X.append(pe)
                        correct_Y.append(sub_flow)
                        if pe < 5:
                            logger.debug(
                                "correct node with 0 pest {0}->{1}->{2}, with branch size: {3}".format(
                                    uid, no, wid, len(graph.vp.seq[node])
                                )
                            )
                else:
                    is_graph_error = False
                    if uid in error_nos:
                        logger.debug(
                            "src: {0} is incorrect graph erroroness node, no optimal ref".format(
                                uid
                            )
                        )
                        is_graph_error = True
                    if wid in error_nos:
                        logger.debug(
                            "tgt: {0} is incorrect graph erroroness node, no optimal ref".format(
                                wid
                            )
                        )
                        is_graph_error = True
                    if len(ref_all_dict[(uid, wid)].intersection(sym_diff)) > 0:
                        is_graph_error = True
                    if is_graph_error:
                        if pe <= cut:
                            false_error_X.append(pe)
                            false_error_Y.append(sub_flow)
                        logger.debug("false positive error pair")
                    else:
                        if pe <= cut:
                            error_X.append(pe)
                            error_Y.append(sub_flow)
                            error_text.append("{0}:{1}:{2}".format(uid, wid, pe))
                        logger.debug("error pair")
            # debug only
            # perform split
            sub_id = no + "*" + str(counter)
            counter += 1
            sub_node = graph_add_vertex(
                graph, simp_node_dict, sub_id, sub_flow, graph.vp.seq[node]
            )

            graph_add_edge(
                graph,
                simp_edge_dict,
                simp_node_dict[uid],
                sub_node,
                graph.ep.overlap[simp_edge_dict[(uid, no)]],
                sub_flow,
            )

            graph_add_edge(
                graph,
                simp_edge_dict,
                sub_node,
                simp_node_dict[wid],
                graph.ep.overlap[simp_edge_dict[(no, wid)]],
                sub_flow,
            )
            link2subs[(uid, wid)] = sub_id

        # keep track of related contig record
        for cno in support_contigs:
            curr_contig, clen, ccov = contig_dict.pop(cno)
            branch_ind = curr_contig.index(no)
            uid = curr_contig[branch_ind - 1] if branch_ind > 0 else None
            wid = (
                curr_contig[branch_ind + 1]
                if branch_ind < len(curr_contig) - 1
                else None
            )
            if uid != None and wid != None:
                # unique mapping
                curr_contig[branch_ind] = link2subs[(uid, wid)]
                contig_dict[cno] = [curr_contig, clen, ccov]
            elif uid == None and wid == None:
                for sub_id in link2subs.values():
                    # all possible contigs
                    contig_dict[cno + "$" + str(sub_id.split("*")[-1])] = [
                        [sub_id],
                        len(graph.vp.seq[simp_node_dict[sub_id]]),
                        graph.vp.dp[simp_node_dict[sub_id]],
                    ]
            elif uid != None and wid == None:
                for (uid2, _), sub_id in link2subs.items():
                    if uid == uid2:
                        curr_contig[branch_ind] = sub_id
                        contig_dict[cno + "$" + str(sub_id.split("*")[-1])] = [
                            list(curr_contig),
                            clen,
                            ccov,
                        ]
            else:
                for (_, wid2), sub_id in link2subs.items():
                    if wid == wid2:
                        curr_contig[branch_ind] = sub_id
                        contig_dict[cno + "$" + str(sub_id.split("*")[-1])] = [
                            list(curr_contig),
                            clen,
                            ccov,
                        ]

        # remove related edges and vertex, update contig tracker
        for uid in us:
            graph_remove_edge(graph, simp_edge_dict, uid, no)
        for wid in ws:
            graph_remove_edge(graph, simp_edge_dict, no, wid)
        graph_remove_vertex(graph, simp_node_dict, no)
        node_to_contig_dict, _ = contig_map_node(contig_dict)

        # update link info
        for (uid, wid), sub_id in link2subs.items():
            for nno in simp_node_dict.keys():
                pe_info[(min(sub_id, nno), max(sub_id, nno))] = None
        for (pu, pv) in list(pe_info.keys()):
            if pu == no or pv == no:
                # out of date
                pe_info.pop((min(pu, pv), max(pu, pv)))
    # final step, assign all the none val pe link to 0
    for k in pe_info.keys():
        if pe_info[k] == None:
            pe_info[k] = 0
    logger.debug("No of branch be removed: " + str(len(set(split_branches))))
    logger.debug("Split branches: " + list_to_string(set(split_branches)))
    logger.info("done")

    # plot the data
    if ref_file:
        _, (ax1) = plt.subplots(1, 1, figsize=(32, 32))
        ax1.scatter(correct_X, correct_Y, color="red", s=100, label="Correct")
        ax1.scatter(
            false_error_X, false_error_Y, color="blue", s=100, label="False-Positive"
        )
        ax1.scatter(error_X, error_Y, color="green", marker="^", s=100, label="Error")

        for index in range(len(error_X)):
            ax1.text(error_X[index], error_Y[index], error_text[index], size=10)

        ax1.set_xlabel("PE")
        ax1.set_ylabel("FLOW")
        ax1.set_title("Scatter Plot - flow vs pe")
        ax1.legend()
        plt.yticks(numpy.arange(0, 500, 10))
        plt.xticks(numpy.arange(0, cut + 1, 1))
        plt.savefig(
            "{0}{1}".format(temp_dir, "/tmp/scatter_plot_pest_{0}.png".format(count_id))
        )

    return len(set(split_branches))


def trivial_split(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    pe_info: dict,
    logger: Logger,
):
    """
    Split the graph, for any (0|1)->N, N->(0|1) branch, split by forking the 1 edge to N edge.
    """
    logger.info("graph trivial split on NT related vertices..")
    # detect all non-trivial branches right now
    non_trivial_branches = get_non_trivial_branches(graph, simp_node_dict)
    trivial_split_count = 0
    id_mapping = {}
    for id in simp_node_dict.keys():
        id_mapping[id] = set()

    for ntno, ntnode in non_trivial_branches.items():
        if graph.vp.color[ntnode] != "black":
            continue
        logger.debug("Current involving NT branch: {0}".format(ntno))
        for inode in set(ntnode.in_neighbors()):
            if graph.vp.color[inode] != "black":
                continue
            ino = graph.vp.id[inode]
            if ino not in id_mapping:
                id_mapping[ino] = set()
            ines = [ue for ue in inode.in_edges() if graph.ep.color[ue] == "black"]
            outes = [ve for ve in inode.out_edges() if graph.ep.color[ve] == "black"]
            if len(ines) > 1 and len(outes) == 1:
                # n to 1
                logger.debug("{0}, n->1 split right".format(ino))
                graph.vp.color[inode] = "gray"
                graph.ep.color[graph.edge(inode, ntnode)] = "gray"
                s = "A"
                for i in range(len(ines)):
                    ine = ines[i]
                    src = ine.source()
                    snode = graph_add_vertex(
                        graph,
                        simp_node_dict,
                        ino + "*" + chr(ord(s) + i),
                        graph.ep.flow[ine],
                        graph.vp.seq[inode],
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
                        ntnode,
                        graph.ep.overlap[graph.edge(inode, ntnode)],
                        graph.ep.flow[ine],
                    )
                    simp_edge_dict[
                        (
                            graph.vp.id[sedge_out.source()],
                            graph.vp.id[sedge_out.target()],
                        )
                    ] = sedge_out
                    id_mapping[ino].add(graph.vp.id[snode])
                    for nno in simp_node_dict.keys():
                        pe_info[
                            (min(graph.vp.id[snode], nno), max(graph.vp.id[snode], nno))
                        ] = None
                trivial_split_count += 1
                # update link information
                for (pu, pv) in list(pe_info.keys()):
                    if pu == ino or pv == ino:
                        # out of date
                        pe_info.pop((min(pu, pv), max(pu, pv)))

        for onode in set(ntnode.out_neighbors()):
            if graph.vp.color[onode] != "black":
                continue
            ono = graph.vp.id[onode]
            if ono not in id_mapping:
                id_mapping[ono] = set()
            ines = [ue for ue in onode.in_edges() if graph.ep.color[ue] == "black"]
            outes = [ve for ve in onode.out_edges() if graph.ep.color[ve] == "black"]
            if len(ines) == 1 and len(outes) > 1:
                # 1 to n
                logger.debug("{0}, 1->n split left".format(ono))
                graph.vp.color[onode] = "gray"
                graph.ep.color[graph.edge(ntnode, onode)] = "gray"
                s = "A"
                for i in range(len(outes)):
                    oute = outes[i]
                    tgt = oute.target()
                    snode = graph_add_vertex(
                        graph,
                        simp_node_dict,
                        ono + "*" + chr(ord(s) + i),
                        graph.ep.flow[oute],
                        graph.vp.seq[onode],
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
                        ntnode,
                        snode,
                        graph.ep.overlap[graph.edge(ntnode, onode)],
                        graph.ep.flow[oute],
                    )
                    simp_edge_dict[
                        (graph.vp.id[sedge_in.source()], graph.vp.id[sedge_in.target()])
                    ] = sedge_in
                    id_mapping[ono].add(graph.vp.id[snode])
                    for nno in simp_node_dict.keys():
                        pe_info[
                            (min(graph.vp.id[snode], nno), max(graph.vp.id[snode], nno))
                        ] = None
                trivial_split_count += 1
                # update link information
                for (pu, pv) in list(pe_info.keys()):
                    if pu == ono or pv == ono:
                        # out of date
                        pe_info.pop((min(pu, pv), max(pu, pv)))
    for k in pe_info.keys():
        if pe_info[k] == None:
            pe_info[k] = 0
    logger.debug(
        "Total split-ted trivial branch count: {0}".format(trivial_split_count)
    )
    return trivial_split_count, id_mapping


def global_trivial_split(
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


def edge_cleaning(
    graph: Graph, simp_edge_dict: dict, contig_dict: dict, pe_info: dict, logger: Logger
):
    """
    Detect the crossing edges and select the confident edges only.
    """
    un_assigned_edge = graph.num_edges()
    assigned = dict.fromkeys(
        [(graph.vp.id[e.source()], graph.vp.id[e.target()]) for e in graph.edges()],
        False,
    )
    _, edge_to_contig_dict = contig_map_node(contig_dict)
    logger.debug("Total edges: " + str(un_assigned_edge))
    # converage iteration
    converage_flag = 0
    while True:
        for node in graph.vertices():
            in_d = node.in_degree()
            in_e = []
            for e in node.in_edges():
                if assigned[(graph.vp.id[e.source()], graph.vp.id[e.target()])]:
                    in_d = in_d - 1
                else:
                    in_e.append(e)

            out_d = node.out_degree()
            out_e = []
            for e in node.out_edges():
                if assigned[(graph.vp.id[e.source()], graph.vp.id[e.target()])]:
                    out_d = out_d - 1
                else:
                    out_e.append(e)

            if in_d == 1:
                assigned[
                    (graph.vp.id[in_e[0].source()], graph.vp.id[in_e[0].target()])
                ] = True
                un_assigned_edge = un_assigned_edge - 1
            if out_d == 1:
                assigned[
                    (graph.vp.id[out_e[0].source()], graph.vp.id[out_e[0].target()])
                ] = True
                un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge

    logger.debug(
        "un-assigned edges after node-weight coverage iteration : {0}".format(
            un_assigned_edge
        )
    )
    for (u, v) in assigned.keys():
        if not assigned[(u, v)]:
            logger.debug(
                "***cross un-assigned edge: {0} -> {1}, with paired end link {2}".format(
                    u, v, pe_info[(min(u, v), max(u, v))]
                )
            )
            if (u, v) in edge_to_contig_dict:
                logger.debug(
                    "support contig: {0}, force assign".format(
                        edge_to_contig_dict[(u, v)]
                    )
                )
                assigned[(u, v)] = True
            else:
                logger.debug("support contig: None")
    for (u, v) in assigned.keys():
        if not assigned[(u, v)]:
            force_assign = True
            for (w, z) in assigned.keys():
                if (u == w or v == z) and assigned[(w, z)]:
                    force_assign = False
                    break
            if not force_assign:
                graph.remove_edge(simp_edge_dict.pop((u, v)))
                logger.debug(
                    "intersect unsupported edge: {0} -> {1}, removed".format(u, v)
                )
            else:
                logger.debug("disjoint unsupported edge: {0} -> {1}, kept".format(u, v))
    return assigned


def iter_graph_disentanglement(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    contig_dict: dict,
    pe_info: dict,
    ref_file: str,
    logger: Logger,
    threshold,
    temp_dir,
):
    BOUND_ITER = len(simp_node_dict) ** 2
    it = 0
    total_removed_branch = 0
    num_split = 0
    iterCount = "A"
    for is_prim in [True, False]:  # False
        do_trivial_split = True
        while it < BOUND_ITER:
            num_split = balance_split(
                graph,
                simp_node_dict,
                simp_edge_dict,
                contig_dict,
                pe_info,
                logger,
                ref_file,
                temp_dir,
                it,
                threshold,
                is_prim,
            )
            graph, simp_node_dict, simp_edge_dict = store_reinit_graph(
                graph,
                simp_node_dict,
                simp_edge_dict,
                logger,
                "{0}/gfa/split_graph_L{1}d.gfa".format(temp_dir, iterCount),
            )
            simp_path_compactification(
                graph, simp_node_dict, simp_edge_dict, contig_dict, pe_info, logger
            )
            graph, simp_node_dict, simp_edge_dict = store_reinit_graph(
                graph,
                simp_node_dict,
                simp_edge_dict,
                logger,
                "{0}/gfa/split_graph_L{1}dc.gfa".format(temp_dir, iterCount),
            )

            if num_split > 0:
                do_trivial_split = True
            else:
                if do_trivial_split:
                    # trivial split nt branch related cases FIXME
                    prev_ids = list(simp_node_dict.keys())
                    trivial_split_count, id_mapping = trivial_split(
                        graph, simp_node_dict, simp_edge_dict, pe_info, logger
                    )
                    logger.debug("my id mapping: {0}".format(id_mapping))
                    graph, simp_node_dict, simp_edge_dict = store_reinit_graph(
                        graph,
                        simp_node_dict,
                        simp_edge_dict,
                        logger,
                        "{0}/gfa/split_graph_L{1}dct.gfa".format(temp_dir, iterCount),
                    )

                    contig_dict_remapping(
                        graph,
                        simp_node_dict,
                        simp_edge_dict,
                        contig_dict,
                        id_mapping,
                        prev_ids,
                        logger,
                    )
                    simp_path_compactification(
                        graph,
                        simp_node_dict,
                        simp_edge_dict,
                        contig_dict,
                        pe_info,
                        logger,
                    )
                    graph, simp_node_dict, simp_edge_dict = store_reinit_graph(
                        graph,
                        simp_node_dict,
                        simp_edge_dict,
                        logger,
                        "{0}/gfa/split_graph_L{1}dctd.gfa".format(temp_dir, iterCount),
                    )

            contig_dup_removed_s(contig_dict, logger)
            trim_contig_dict(graph, simp_node_dict, contig_dict, logger)
            # analysis
            if ref_file:
                map_ref_to_graph(
                    ref_file,
                    simp_node_dict,
                    "{0}/gfa/split_graph_L{1}dc.gfa".format(temp_dir, iterCount),
                    logger,
                    True,
                    "{0}/paf/node_to_ref_{1}.paf".format(temp_dir, iterCount),
                    "{0}/tmp/temp_gfa_to_fasta_{1}.fasta".format(temp_dir, iterCount),
                )
            # analysis
            total_removed_branch += num_split
            it += 1
            iterCount = chr(ord(iterCount) + 1)
            if num_split == 0:
                if do_trivial_split:
                    do_trivial_split = False
                else:
                    break

    logger.debug("Total non-trivial branches removed: " + str(total_removed_branch))
    non_trivial_branches = get_non_trivial_branches(graph, simp_node_dict)
    logger.debug(
        list_to_string(
            non_trivial_branches.keys(),
            "non-trivial branches ({0}) left after paired-end&single-strand links".format(
                len(non_trivial_branches)
            ),
        )
    )

    graph, simp_node_dict, simp_edge_dict = store_reinit_graph(
        graph,
        simp_node_dict,
        simp_edge_dict,
        logger,
        "{0}/gfa/split_graph_final.gfa".format(temp_dir),
    )
    return graph, simp_node_dict, simp_edge_dict


def best_aln_score(graph: Graph, ori, strain, ref_file, temp_dir):
    fname = "{0}/temp_{1}.fa".format(temp_dir, ori)
    pafname = "{0}/temp_{1}_aln.paf".format(temp_dir, ori)
    subprocess.check_call('echo "" > {0}'.format(fname), shell=True)
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
            best_aln.append(
                [
                    splited[0],
                    int(splited[10]),
                    splited[5],
                    int(splited[10]) - int(splited[9]),
                ]
            )
        paf.close()
    subprocess.check_call("rm {0}".format(pafname), shell=True)
    return best_aln
