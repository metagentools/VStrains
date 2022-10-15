#!/usr/bin/env python3


from graph_tool.all import Graph
from utils.VStrains_Utilities import *
from utils.VStrains_Decomposition import get_non_trivial_branches, global_trivial_split
from utils.VStrains_IO import store_reinit_graph


__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, VStrains Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"

def best_matching(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, pe_info: dict, logger: Logger):
    full_link = {}
    non_trivial_branches = get_non_trivial_branches(graph, simp_node_dict)
    node_to_contig_dict, _ = contig_map_node(contig_dict)
    for no, node in non_trivial_branches.items():
        us = [graph.vp.id[src] for src in node.in_neighbors()]
        ws = [graph.vp.id[tgt] for tgt in node.out_neighbors()]
        logger.debug("---------------------------------------------")
        logger.debug("current non trivial branch: {0}, in-degree: {1}, out-degree: {2}".format(no, len(us), len(ws)))
        # add contig supports
        support_contigs = node_to_contig_dict.get(no, [])
        con_info = {}
        for cno in support_contigs:
            [contig, clen, ccov] = contig_dict[cno]
            loc = contig.index(no)
            if loc > 0 and loc < len(contig) - 1:
                con_info[(contig[loc-1], contig[loc+1])] = con_info.get((contig[loc-1], contig[loc+1]), [])
                con_info[(contig[loc-1], contig[loc+1])].append((cno,clen,ccov))
            print_contig(cno, clen, round(ccov, 2), contig[max(loc-1,0):loc+2], logger, "support contig")        
        kept_link = {}
        sec_comb = []
        # init node usage for current branch
        in_usage = dict.fromkeys(us, 0)
        out_usage = dict.fromkeys(ws, 0)

        # align contig link first, and update status
        logger.debug("align contig link first")
        for uid in us:
            for wid in ws:
                logger.debug("---------------------")
                u = simp_node_dict[uid]
                w = simp_node_dict[wid]
                curr_pe = pe_info[(min(uid, wid), max(uid, wid))]

                logger.debug("{0} -> {1} PE: {2}".format(uid, wid, curr_pe))
                logger.debug("cov info: {0}[{1}] -> {2}[{3}]".
                format(graph.ep.flow[graph.edge(u,node)], pe_info[(min(uid, no), max(uid, no))],
                        graph.ep.flow[graph.edge(node, w)], pe_info[(min(no, wid), max(no, wid))]))
                accept = False
                if (uid, wid) in con_info:
                    logger.debug("current link supported by contig: {0}, added".format(con_info[(uid,wid)]))
                    accept = True
                if uid == wid:
                    logger.debug("current link is a self link: {0}, potential cyclic strain, added".format(uid))
                    accept = True

                if accept:
                    in_usage[uid] += 1
                    out_usage[wid] += 1
                    kept_link[(uid, wid)] = curr_pe
                else:
                    logger.debug("current link is secondary choice, process later")
                    sec_comb.append((uid, wid, curr_pe))


        logger.debug("align paired end/single end information first (if any) to isolated nodes")
        sorted_sec_comb = sorted(sec_comb, key=lambda x: x[2], reverse=True)
        for (uid, wid, pe) in sorted_sec_comb:
            if pe > 0:
                logger.debug("-----SEC LINK {0} -> {1} PE: {2}-----".format(uid, wid, pe))
                logger.debug("- link [ > 0] supported case, added")
                in_usage[uid] += 1
                out_usage[wid] += 1
                kept_link[(uid, wid)] = pe
        full_link[no] = kept_link
    return full_link

# extend contigs on both end, until a non distinct extension
def contig_extension(graph: Graph, simp_node_dict: dict, contig: list, ccov, full_link: dict, logger: Logger, threshold):
    visited = dict.fromkeys(simp_node_dict.keys(), False)
    for no in contig[1:-1]:
        visited[no] = True
    final_path = []
    final_path.extend([simp_node_dict[no] for no in contig][1:-1])

    curr = simp_node_dict[contig[-1]]
    logger.debug("c-t extension")
    while curr != None and not visited[graph.vp.id[curr]]:
        visited[graph.vp.id[curr]] = True
        final_path.append(curr)
        out_branches = list([n for n in curr.out_neighbors()])
        if len(out_branches) == 0:
            curr = None
            logger.debug("Reach the end")
        elif len(out_branches) == 1:
            curr = out_branches[0]
            logger.debug("direct extending.. {0}".format(graph.vp.id[curr]))
        else:
            f_assigned = False
            if graph.vp.id[curr] in full_link and len(final_path) > 1:
                logger.debug("Curr is Branch")
                curr_links = [simp_node_dict[wid] for (uid, wid) in full_link[graph.vp.id[curr]].keys() if uid == graph.vp.id[final_path[-2]]]
                if len(curr_links) == 1:
                    # curr = curr_links[0]
                    # logger.debug("single link next: {0}".format(graph.vp.id[curr]))
                    if graph.vp.dp[curr_links[0]] - ccov <= -2 * threshold:
                        curr = None
                        logger.debug("{0} single link < 2delta, use coverage".format(graph.vp.id[curr_links[0]]))
                    else:
                        curr = curr_links[0]
                        logger.debug("single link next: {0}".format(graph.vp.id[curr]))
                elif len(curr_links) > 1:
                    logger.debug("Ambiguous, stop extension")
                    curr = None
                else:
                    logger.debug("No link in here, use coverage information")
                    f_assigned = True
            else:
                curr = None
                logger.debug("Not in full link or len of path <= 1")
            if f_assigned:
                in_branches = list([n for n in curr.in_neighbors()])
                if len(final_path) > 1 and len(in_branches) > 0:
                    curru = final_path[-2]
                    opt_ws = sorted(out_branches, key=lambda ww: abs(graph.vp.dp[curru] - graph.vp.dp[ww]))
                    bestw = opt_ws[0]
                    opt_us = sorted(in_branches, key=lambda uu: abs(graph.vp.dp[bestw] - graph.vp.dp[uu]))
                    if opt_us[0] == curru:
                        delta = max(2 * abs(graph.vp.dp[curru] - graph.vp.dp[bestw]),threshold)
                        if len(opt_us) > 1 and abs(graph.vp.dp[opt_us[1]] - graph.vp.dp[bestw]) <= delta:
                            logger.debug("ambiguous best matching, stop extension")
                            continue
                        if len(opt_ws) > 1 and abs(graph.vp.dp[curru] - graph.vp.dp[opt_ws[1]]) <= delta:
                            logger.debug("ambiguous best matching, stop extension")
                            continue
                        logger.debug("best matching")
                        curr = bestw
                    else:
                        logger.debug("Not best match")
                        curr = None
                else:
                    curr = None
                    logger.debug("No Link + Not trivial, stop extension")
            if curr == None:
                single_bests = sorted([(onode,graph.vp.dp[onode]) for onode in out_branches], key=lambda tp: tp[1], reverse=True)
                logger.debug("Try last bit: 1st: {0}, 2nd: {1}, delta: {2}, cov: {3}".
                    format((graph.vp.id[single_bests[0][0]], single_bests[0][1]), 
                    (graph.vp.id[single_bests[1][0]], single_bests[1][1]), threshold, ccov))
                if single_bests[0][1] - ccov > -threshold and single_bests[1][1] - ccov <= -threshold:
                    logger.debug("Last bit succ")
                    curr = single_bests[0][0]
                else:
                    logger.debug("Last bit fail")
    unode = simp_node_dict[contig[0]]
    if len(contig) == 1 and final_path[-1] not in unode.in_neighbors():
        visited[contig[0]] = False
        final_path.pop(0)
    curr = unode
    logger.debug("s-c extension")
    while curr != None and not visited[graph.vp.id[curr]]:
        visited[graph.vp.id[curr]] = True
        final_path.insert(0, curr)
        in_branches = list([n for n in curr.in_neighbors()])
        if len(in_branches) == 0:
            curr = None
            logger.debug("Reach the end")
        elif len(in_branches) == 1:
            curr = in_branches[0]
            logger.debug("direct extending.. {0}".format(graph.vp.id[curr]))
        else:
            f_assigned = False
            if graph.vp.id[curr] in full_link and len(final_path) > 1:
                logger.debug("Curr is Branch")
                curr_links = [simp_node_dict[uid] for (uid, wid) in full_link[graph.vp.id[curr]].keys() if wid == graph.vp.id[final_path[1]]]
                if len(curr_links) == 1:
                    # curr = curr_links[0]
                    # logger.debug("single link next: {0}".format(graph.vp.id[curr]))
                    if graph.vp.dp[curr_links[0]] - ccov <= -2 * threshold:
                        curr = None
                        logger.debug("{0} single link < 2delta, use coverage".format(graph.vp.id[curr_links[0]]))
                    else:
                        curr = curr_links[0]
                        logger.debug("prev: {0}".format(graph.vp.id[curr]))
                elif len(curr_links) > 1:
                    logger.debug("Ambiguous, stop extension")
                    curr = None
                else:
                    logger.debug("No link in here, use coverage information")
                    f_assigned = True
            else:
                curr = None
                logger.debug("Not in full link or len of path <= 1")
            if f_assigned:
                out_branches = list([n for n in curr.out_neighbors()])
                if len(final_path) > 1 and len(out_branches) > 0:
                    currw = final_path[1]
                    opt_us = sorted(in_branches, key=lambda uu: abs(graph.vp.dp[currw] - graph.vp.dp[uu]))
                    bestu = opt_us[0]
                    opt_ws = sorted(out_branches, key=lambda ww: abs(graph.vp.dp[bestu] - graph.vp.dp[ww]))
                    if opt_ws[0] == currw:
                        delta = max(2 * abs(graph.vp.dp[currw] - graph.vp.dp[bestu]), threshold)
                        if len(opt_us) > 1 and abs(graph.vp.dp[opt_us[1]] - graph.vp.dp[currw]) <= delta:
                            logger.debug("ambiguous best matching, stop extension")
                            continue
                        if len(opt_ws) > 1 and abs(graph.vp.dp[bestu] - graph.vp.dp[opt_ws[1]]) <= delta:
                            logger.debug("ambiguous best matching, stop extension")
                            continue
                        logger.debug("best matching")
                        curr = bestu
                    else:
                        logger.debug("Not best match")
                        curr = None
                else:
                    logger.debug("No Link + Not trivial, stop extension")
                    curr = None
            if curr == None:
                single_bests = sorted([(inode,graph.vp.dp[inode]) for inode in in_branches], key=lambda tp: tp[1], reverse=True)
                logger.debug("Try last bit: 1st: {0}, 2nd: {1}, delta: {2}, cov: {3}".
                    format((graph.vp.id[single_bests[0][0]], single_bests[0][1]), 
                    (graph.vp.id[single_bests[1][0]], single_bests[1][1]), threshold, ccov))
                if single_bests[0][1] - ccov > -threshold and single_bests[1][1] - ccov <= -threshold:
                    logger.debug("Last bit succ")
                    curr = single_bests[0][0]
                else:
                    logger.debug("Last bit fail")
    return final_path

def final_extension(graph: Graph, simp_node_dict: dict, contig: list, full_link: dict, logger: Logger):
    visited = dict.fromkeys(simp_node_dict.keys(), False)
    for no in contig[1:-1]:
        visited[no] = True
    curr = simp_node_dict[contig[-1]]
    final_path = []
    final_path.extend([simp_node_dict[no] for no in contig][1:-1])
    # from curr to the tail, or to the non-extendable end
    logger.debug("c-t extension")
    while curr != None and not visited[graph.vp.id[curr]]:
        visited[graph.vp.id[curr]] = True
        final_path.append(curr)
        out_branches = list([n for n in curr.out_neighbors()])
        if len(out_branches) == 0:
            curr = None
            logger.debug("Reach the end")
        elif len(out_branches) == 1:
            curr = out_branches[0]
            logger.debug("direct extending.. {0}".format(graph.vp.id[curr]))
        else:
            if graph.vp.id[curr] in full_link and len(final_path) > 1:
                logger.debug("Curr is Branch")
                curr_links = [simp_node_dict[wid] for (uid, wid) in full_link[graph.vp.id[curr]].keys() if uid == graph.vp.id[final_path[-2]]]
                if len(curr_links) == 1:
                    curr = curr_links[0]
                    logger.debug("single link next: {0}".format(graph.vp.id[curr]))
                else:
                    logger.debug("No/more link in here, end entension")
                    curr = None
            else:
                curr = None
                logger.debug("Not in full link or len of path <= 1")

    unode = simp_node_dict[contig[0]]
    if len(contig) == 1 and final_path[-1] not in unode.in_neighbors():
        visited[contig[0]] = False
        final_path.pop(0)
    curr = unode
    # from head to the curr, or to the non-extendable end
    logger.debug("s-c extension")
    while curr != None and not visited[graph.vp.id[curr]]:
        visited[graph.vp.id[curr]] = True
        final_path.insert(0, curr)
        in_branches = list([n for n in curr.in_neighbors()])
        if len(in_branches) == 0:
            curr = None
            logger.debug("Reach the end")
        elif len(in_branches) == 1:
            curr = in_branches[0]
            logger.debug("direct extending.. {0}".format(graph.vp.id[curr]))
        else:
            if graph.vp.id[curr] in full_link and len(final_path) > 1:
                logger.debug("Curr is Branch")
                curr_links = [simp_node_dict[uid] for (uid, wid) in full_link[graph.vp.id[curr]].keys() if wid == graph.vp.id[final_path[1]]]
                if len(curr_links) == 1:
                    curr = curr_links[0]
                    logger.debug("single link next: {0}".format(graph.vp.id[curr]))
                else:
                    logger.debug("No/more link in here, end extension")
                    curr = None
            else:
                curr = None
                logger.debug("Not in full link or len of path <= 1")
    return final_path

def get_bubble_nodes(simp_node_dict: dict, contig: list):
    bubbles = []
    for no in contig:
        if simp_node_dict[no].in_degree() == 1 and simp_node_dict[no].out_degree() == 1:
            bubbles.append(simp_node_dict[no])
    return bubbles

def reduce_graph(graph: Graph, simp_node_dict: dict, usages: dict, full_link: dict, logger: Logger, path, pcov, threshold):
    del_nodes_ids = []
    for node in path:
        usages[graph.vp.id[node]] += 1
        graph.vp.dp[node] -= pcov
        if graph.vp.dp[node] <= threshold:
            del_nodes_ids.append(graph.vp.id[node])
            graph.vp.color[node] = "gray"
            usages.pop(graph.vp.id[node])
    logger.debug(list_to_string(del_nodes_ids, "invalid nodes"))
    for links in full_link.values():
        for (uid, wid) in list(links.keys()):
            if graph.vp.color[simp_node_dict[uid]] != "black" or graph.vp.color[simp_node_dict[wid]] != "black":
                links.pop((uid, wid))
                logger.debug("[D]{0}, {1}".format(uid, wid))

def reduce_id_simple(id_l: list):
    ids = []
    for id in id_l:
        for iid in id.split("&"):
            if iid.find("*") != -1:
                ids.append(iid[:iid.find("*")])
            else:
                ids.append(iid)
    return ids

def path_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, full_link: dict, pe_info: dict, logger: Logger, threshold, temp_dir):
    logger.debug("-------------------------PATH Extension, delta: {0}".format(threshold))
    usages = dict.fromkeys(simp_node_dict.keys(), 0) # record the usage of each nodes
    strain_dict = {}
    rid = 1
    sno2ids = dict()
    while len(contig_dict) > 0:
        # perform trivial split
        prev_ids = list(simp_node_dict.keys())
        trivial_split_count, id_mapping = global_trivial_split(
            graph, simp_node_dict, simp_edge_dict, logger
        )
        graph, simp_node_dict, simp_edge_dict = store_reinit_graph(graph, simp_node_dict, simp_edge_dict, logger, 
            "{0}/gfa/graph_S{1}.gfa".format(temp_dir, rid))
        red_id_mapping = contig_dict_remapping(
            graph,
            simp_node_dict,
            simp_edge_dict,
            contig_dict,
            id_mapping,
            prev_ids,
            logger,
        )
        # update links
        for no in list(full_link.keys()):
            if no not in simp_node_dict:
                full_link.pop(no)
            else:
                kept_link = full_link.pop(no)
                node = simp_node_dict[no]
                for (uid, wid), pe in list(kept_link.items()):
                    # if len(red_id_mapping[uid]) != 1 or len(red_id_mapping[wid]) != 1:
                    #     kept_link.pop((uid, wid))
                    # else:
                    #     kept_link[(list(red_id_mapping[uid])[0], list(red_id_mapping[wid])[0])] = pe
                    kept_link.pop((uid, wid))
                    if len(red_id_mapping[uid]) == 1 or len(red_id_mapping[wid]) == 1:
                        for uuid in red_id_mapping[uid]:
                            for wwid in red_id_mapping[wid]:
                                if (uuid, wwid) not in kept_link and (simp_node_dict[uuid] in node.in_neighbors()) and (simp_node_dict[wwid] in node.out_neighbors()):
                                    kept_link[(uuid, wwid)] = pe
                full_link[no] = kept_link
        # update usages
        for no, u in list(usages.items()):
            usages.pop(no)
            for new_no in red_id_mapping[no]:
                usages[new_no] = u
        ############################
        # get longest contig
        (longest_cno, [contig, clen, ccov]) = max(contig_dict.items(), key=lambda tp: tp[1][1])
        contig_dict.pop(longest_cno)
        if all(usages[cn] > 0 for cn in contig):
            print_contig(longest_cno, clen, ccov, contig, logger, "-----> Used previously")
            continue
        if any(graph.vp.color[simp_node_dict[no]] == "gray" for no in contig):
            print_contig(longest_cno, clen, ccov, contig, logger, "-----> Some node low cov, skip")
            continue
        
        cbubbles = get_bubble_nodes(simp_node_dict, contig)
        bbl_cov = numpy.median([graph.vp.dp[node] for node in cbubbles]) if len(cbubbles) != 0 else ccov
        print_contig(longest_cno, clen, bbl_cov, contig, logger, "-----> Current extending contig: org ccov: {0}, use min {1}".format(ccov, min(ccov, bbl_cov)))

        path = contig_extension(graph, simp_node_dict, contig, min(ccov, bbl_cov), full_link, logger, threshold)
        pno = "A" + str(rid)
        plen = path_len(graph, path)
        path_ids = [graph.vp.id[n] for n in path]
        sno2ids[pno] = []
        for pid in path_ids:
            if pid in sno2ids:
                sno2ids[pno].extend(sno2ids[pid])
            else:
                sno2ids[pno].append(pid)
        pbubbles = get_bubble_nodes(simp_node_dict, path_ids)
        bbl_pcov = numpy.median([graph.vp.dp[node] for node in pbubbles]) if len(pbubbles) != 0 else ccov
        pcov = min([ccov, bbl_pcov, bbl_cov])
        logger.debug(path_to_id_string(graph, path, "---*extended from contig {0}".format(longest_cno)))
        logger.debug("name: {0}, plen: {1}, pcov: {2}, bubble cov: {3}".format(pno, plen, pcov, bbl_pcov))
        strain_dict[pno] = [sno2ids[pno], plen, pcov]
        for pid in path_ids:
            if pid in strain_dict:
                strain_dict.pop(pid)
        path_ins = [n for n in path[0].in_neighbors()]
        path_outs = [n for n in path[-1].out_neighbors()]
        if len(path_ins) == 0 and len(path_outs) == 0:
            # both end st
            logger.debug("st isolated, add to strain")
            reduce_graph(graph, simp_node_dict, usages, full_link, logger, path, pcov, threshold)
        elif len(path_ins) != 0 and len(path_outs) == 0:
            logger.debug("left connected, wait")
            reduce_graph(graph, simp_node_dict, usages, full_link, logger, path[1:], pcov, threshold) 
            pnode = graph_add_vertex(graph, simp_node_dict, pno, pcov, path_to_seq(graph, path[1:], pno))
            graph_add_edge(graph, simp_edge_dict, path[0], pnode, graph.ep.overlap[graph.edge(path[0], path[1])], pcov)
            usages[pno] = 0
        elif len(path_ins) == 0 and len(path_outs) != 0:
            logger.debug("right connected, wait")
            reduce_graph(graph, simp_node_dict, usages, full_link, logger, path[:-1], pcov, threshold) 
            pnode = graph_add_vertex(graph, simp_node_dict, pno, pcov, path_to_seq(graph, path[:-1], pno))
            graph_add_edge(graph, simp_edge_dict, pnode, path[-1], graph.ep.overlap[graph.edge(path[-2], path[-1])], pcov)
            usages[pno] = 0
        else:
            logger.debug("both connected, wait")
            reduce_graph(graph, simp_node_dict, usages, full_link, logger, path[1:-1], pcov, threshold)
            if len(path[1:-1]) > 0:
                pnode = graph_add_vertex(graph, simp_node_dict, pno, pcov, path_to_seq(graph, path[1:-1], pno))
                graph_add_edge(graph, simp_edge_dict, path[0], pnode, graph.ep.overlap[graph.edge(path[0], path[1])], pcov)
                graph_add_edge(graph, simp_edge_dict, pnode, path[-1], graph.ep.overlap[graph.edge(path[-2], path[-1])], pcov)
                usages[pno] = 0 
        
        graph, simp_node_dict, simp_edge_dict = store_reinit_graph(graph, simp_node_dict, simp_edge_dict, logger, 
            "{0}/gfa/graph_S{1}post.gfa".format(temp_dir, rid))
        for cno in list(contig_dict.keys()):
            delete = False
            for no in contig_dict[cno][0]:
                if no not in simp_node_dict:
                    delete = True
            if delete:
                contig_dict.pop(cno)
        rid += 1

    # remove trivial split multiple nodes
    seq_dict = {}
    for node in graph.vertices():
        if graph.vp.seq[node] not in seq_dict:
            seq_dict[graph.vp.seq[node]] = []
        seq_dict[graph.vp.seq[node]].append(node)

    for _, sp_nodes in seq_dict.items():
        if len(sp_nodes) > 1:
            sorted_sp_nodes = sorted(sp_nodes, key=lambda vnode: graph.vp.dp[vnode], reverse=True)
            for vnode in sorted_sp_nodes[1:]:
                graph_remove_vertex(graph, simp_node_dict, graph.vp.id[vnode])
                usages.pop(graph.vp.id[vnode])
    graph, simp_node_dict, simp_edge_dict = store_reinit_graph(graph, simp_node_dict, simp_edge_dict, logger, 
        "{0}/gfa/graph_S_final.gfa".format(temp_dir))
    # assign link information
    final_link_info = {}
    for node in graph.vertices():
        for node2 in graph.vertices():
            if node > node2:
                continue

            nid1s = reduce_id_simple([graph.vp.id[node]]) if graph.vp.id[node][0] != "A" else reduce_id_simple(sno2ids[graph.vp.id[node].split("*")[0]])
            nid2s = reduce_id_simple([graph.vp.id[node2]]) if graph.vp.id[node2][0] != "A" else reduce_id_simple(sno2ids[graph.vp.id[node2].split("*")[0]])
            kpair = (min(graph.vp.id[node], graph.vp.id[node2]), max(graph.vp.id[node], graph.vp.id[node2]))
            final_link_info[kpair] = 0
            for id1 in nid1s:
                for id2 in nid2s:
                    inner_kpair = (min(id1, id2), max(id1, id2))
                    final_link_info[kpair] += pe_info[inner_kpair]

    nt_branches = get_non_trivial_branches(graph, simp_node_dict)
    final_links = {}
    for no, node in nt_branches.items():
        final_links[no] = {}
        us = [graph.vp.id[src] for src in node.in_neighbors()]
        ws = [graph.vp.id[tgt] for tgt in node.out_neighbors()]
        logger.debug("---------------------------------------------")
        logger.debug("current non trivial branch: {0}, in-degree: {1}, out-degree: {2}".format(no, len(us), len(ws)))
        combs = []
        in_usage = dict.fromkeys(us, 0)

        out_usage = dict.fromkeys(ws, 0)
        for uid in us:
            for wid in ws:
                combs.append((uid, wid, final_link_info[(min(uid, wid), max(uid, wid))]))
        sorted_comb = sorted(combs, key=lambda x: x[2], reverse=True)
        for (uid, wid, lf) in sorted_comb:
            logger.debug("---------------------")
            if lf > 0 and in_usage[uid] == 0 and out_usage[wid] == 0:
                logger.debug("-----SEC LINK {0} -> {1} LINK: {2}-----".format(uid, wid, lf))
                logger.debug("- unique link [ > 0] supported case, added")
                final_links[no][(uid, wid)] = lf
                in_usage[uid] += 1
                out_usage[wid] += 1

    # add all the nodes that not be used in contig extension to final resulting sets
    for node in sorted(graph.vertices(), key=lambda nd: len(graph.vp.seq[nd]), reverse=True):
        if len(graph.vp.seq[node]) <= 600:
            break
        if usages[graph.vp.id[node]] == 0:
            logger.debug("Extend from free node: {0}".format(graph.vp.id[node]))
            ccov = graph.vp.dp[node]
            path = final_extension(graph, simp_node_dict, [graph.vp.id[node]], final_links, logger)
            pno = "N" + str(rid)
            plen = path_len(graph, path)
            path_ids = [graph.vp.id[n] for n in path]
            pids = []
            for pid in path_ids:
                if pid in sno2ids:
                    pids.extend(sno2ids[pid])
                else:
                    pids.append(pid)
            for pid in path_ids:
                if pid in strain_dict:
                    strain_dict.pop(pid)
            pbubbles = get_bubble_nodes(simp_node_dict, path_ids)
            pcov = numpy.median([graph.vp.dp[node] for node in pbubbles]) if len(pbubbles) != 0 else graph.vp.dp[node]
            logger.debug(path_to_id_string(graph, path, "---*extended from free node {0}".format(graph.vp.id[node])))
            logger.debug("name: {0}, plen: {1}, pcov: {2}".format(pno, plen, pcov))
            strain_dict[pno] = [pids, plen, pcov]
            for node in path:
                usages[graph.vp.id[node]] += 1
            rid += 1
    for sno, [_, _, scov] in list(strain_dict.items()):
        if scov <= 2*threshold:
            strain_dict.pop(sno)
    return strain_dict, usages