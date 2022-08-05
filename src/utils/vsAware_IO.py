#!/usr/bin/env python3

from logging import Logger
from graph_tool.all import Graph
import gfapy
import subprocess
import sys
import re

from utils.vsAware_Utilities import path_len, reverse_seq, print_vertex, path_ids_to_seq


def init_graph():
    graph = Graph(directed=True)
    graph.vp.seq = graph.new_vertex_property("string", val="")
    graph.vp.dp = graph.new_vertex_property("double")
    graph.vp.id = graph.new_vertex_property("string", val="UD")
    graph.vp.color = graph.new_vertex_property("string")

    graph.ep.overlap = graph.new_edge_property("int", val=0)
    graph.ep.flow = graph.new_edge_property("double", val=0.0)
    graph.ep.color = graph.new_edge_property("string")

    return graph


def gfa_to_graph(gfa_file, logger: Logger, init_ori=1):
    """
    Convert assembly graph gfa file to graph
    Nodes: segment with corresponding
    """

    logger.info("Parsing GFA format graph")
    gfa = gfapy.Gfa().from_file(filename=gfa_file)
    logger.info(
        "Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version)
    )

    graph = init_graph()
    graph.vp.visited = graph.new_vertex_property("int16_t", val=0)
    graph.vp.ori = graph.new_vertex_property("int16_t")  # 1 = +, -1 = -

    graph.ep.visited = graph.new_edge_property("int", val=0)

    # S
    node_dict = {}
    dp_dict = {}
    edge_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [t, seg_no, seg] = (str(line).split("\t"))[:3]
        tags = (str(line).split("\t"))[3:]
        dp = [tag for tag in tags if tag.startswith("dp") or tag.startswith("DP")][0]
        # gfa format check
        assert t == "S" and (dp[:2] == "DP" or dp[:2] == "dp")
        dp_float = float(dp.split(":")[2])
        v_pos = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.dp[v_pos] = dp_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.visited[v_pos] = -1
        graph.vp.color[v_pos] = "black"

        v_neg = graph.add_vertex()
        graph.vp.seq[v_neg] = reverse_seq(seg)
        graph.vp.dp[v_neg] = dp_float
        graph.vp.id[v_neg] = seg_no
        graph.vp.ori[v_neg] = -1
        graph.vp.visited[v_neg] = -1
        graph.vp.color[v_neg] = "black"

        node_dict[seg_no] = (v_pos, v_neg)
        dp_dict[seg_no] = dp_float
    # L
    for edge in gfa.edges:
        [t, seg_no_l, ori_l, seg_no_r, ori_r] = (str(edge).split("\t"))[:5]
        tags = (str(edge).split("\t"))[5:]
        overlap_len = [tag for tag in tags if tag.endswith("m") or tag.endswith("M")][0]
        # gfa format check
        assert t == "L" and overlap_len[-1] == "M"

        u_pos, u_neg = node_dict[seg_no_l]
        v_pos, v_neg = node_dict[seg_no_r]
        u = u_pos if ori_l == "+" else u_neg
        v = v_pos if ori_r == "+" else v_neg

        if (seg_no_l, graph.vp.ori[u], seg_no_r, graph.vp.ori[v]) in edge_dict:
            logger.error(
                "parallel edge found, invalid case in assembly graph, please double-check the assembly graph format"
            )
            logger.error("Pipeline aborted")
            sys.exit(1)

        if seg_no_l == seg_no_r:
            graph.vp.seq[u] = str.lower(graph.vp.seq[u])
            graph.vp.seq[v] = str.lower(graph.vp.seq[v])
            continue

        e = graph.add_edge(source=u, target=v)
        graph.ep.overlap[e] = int(overlap_len[:-1])
        graph.ep.color[e] = "black"

        edge_dict[(seg_no_l, graph.vp.ori[u], seg_no_r, graph.vp.ori[v])] = e

    graph, simp_node_dict, simp_edge_dict = flip_graph_bfs(
        graph, node_dict, edge_dict, dp_dict, logger, init_ori
    )
    red_graph, red_node_dict, red_edge_dict = reduce_graph(
        graph, simp_node_dict, simp_edge_dict
    )
    return red_graph, red_node_dict, red_edge_dict


def flip_graph_bfs(
    graph: Graph,
    node_dict: dict,
    edge_dict: dict,
    dp_dict: dict,
    logger: Logger,
    init_ori=1,
):
    """
    Flip all the node orientation.

    return an node_dict, which only contains one orientation per node for simplicity.
    rename all the used node to positive, and forbidden the opponent node.
    """

    def source_node_via_dp(dp_dict: dict):
        """
        return the pos-neg node with maximum depth
        """
        return max(dp_dict, key=dp_dict.get)

    def reverse_edge(graph: Graph, edge, node_dict: dict, edge_dict: dict):
        """
        reverse an edge with altered orientation and direction.
        """
        tmp_s = edge.source()
        tmp_t = edge.target()

        edge_dict.pop(
            (
                graph.vp.id[tmp_s],
                graph.vp.ori[tmp_s],
                graph.vp.id[tmp_t],
                graph.vp.ori[tmp_t],
            )
        )

        tmp_s_pos, tmp_s_neg = node_dict[graph.vp.id[tmp_s]]
        tmp_t_pos, tmp_t_neg = node_dict[graph.vp.id[tmp_t]]
        s = tmp_t_pos if graph.vp.ori[tmp_t] == -1 else tmp_t_neg
        t = tmp_s_pos if graph.vp.ori[tmp_s] == -1 else tmp_s_neg

        o = graph.ep.overlap[edge]
        graph.remove_edge(edge)
        e = graph.add_edge(s, t)
        graph.ep.overlap[e] = o
        edge_dict[
            (graph.vp.id[s], graph.vp.ori[s], graph.vp.id[t], graph.vp.ori[t])
        ] = e

        return graph, e, edge_dict

    logger.info("flip graph orientation..")
    pick_dict = {}
    while set(dp_dict):
        seg_no = source_node_via_dp(dp_dict)
        source_pos, source_neg = node_dict[seg_no]
        graph.vp.visited[source_pos] = 0
        graph.vp.visited[source_neg] = 0
        fifo_queue = [[node_dict[seg_no], init_ori]]

        while fifo_queue:
            (v_pos, v_neg), ori = fifo_queue.pop()
            dp_dict.pop(graph.vp.id[v_pos])

            u = None
            if ori == 1:
                u = v_pos
                pick_dict[graph.vp.id[u]] = "+"
                # print_vertex(graph, v_neg, "node to reverse pos")
                for e in set(v_neg.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            else:
                u = v_neg
                pick_dict[graph.vp.id[u]] = "-"
                # print_vertex(graph, v_pos, "node to reverse neg")
                for e in set(v_pos.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")

            graph.vp.visited[v_pos] = 1
            graph.vp.visited[v_neg] = 1
            # add further nodes into the fifo_queue
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1:
                    vpos, vneg = node_dict[graph.vp.id[adj_node]]
                    graph.vp.visited[vpos] = 0
                    graph.vp.visited[vneg] = 0
                    # print("appending node {0} to queue".format(graph.vp.id[adj_node]))
                    fifo_queue.append(
                        [node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]]
                    )

    # verify sorted graph
    logger.info("final verifying graph..")
    assert len(pick_dict) == len(node_dict)
    for key, item in list(pick_dict.items()):
        v_pos, v_neg = node_dict[key]
        if item == "+":
            # FIXME split v_neg to a new node
            if v_neg.in_degree() + v_neg.out_degree() > 0:
                print_vertex(
                    graph, v_neg, logger, "pick ambiguous found, pick both, split node"
                )
                pick_dict[key] = "t"
        else:
            # FIXME split v_neg to a new node
            if v_pos.in_degree() + v_pos.out_degree() > 0:
                print_vertex(
                    graph, v_pos, logger, "pick ambiguous found, pick both, split node"
                )
                pick_dict[key] = "t"
    logger.info("Graph is verified")

    simp_node_dict = {}
    for seg_no, pick in pick_dict.items():
        if pick == "+":
            simp_node_dict[seg_no] = node_dict[seg_no][0]
        elif pick == "-":
            simp_node_dict["-" + seg_no] = node_dict[seg_no][1]
            graph.vp.id[node_dict[seg_no][1]] = "-" + seg_no
        else:
            simp_node_dict[seg_no] = node_dict[seg_no][0]
            graph.vp.id[node_dict[seg_no][0]] = seg_no
            simp_node_dict["-" + seg_no] = node_dict[seg_no][1]
            graph.vp.id[node_dict[seg_no][1]] = "-" + seg_no

    simp_edge_dict = {}
    for e in edge_dict.values():
        simp_edge_dict[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
    logger.info("done")
    return graph, simp_node_dict, simp_edge_dict


def reduce_graph(unsimp_graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    graph = init_graph()
    red_node_dict = {}
    red_edge_dict = {}

    for no, node in simp_node_dict.items():
        v = graph.add_vertex()
        graph.vp.seq[v] = unsimp_graph.vp.seq[node]
        graph.vp.dp[v] = unsimp_graph.vp.dp[node]
        graph.vp.id[v] = unsimp_graph.vp.id[node]
        graph.vp.color[v] = "black"
        red_node_dict[no] = v

    for (u, v), e in simp_edge_dict.items():
        source = red_node_dict[u]
        sink = red_node_dict[v]

        re = graph.add_edge(source, sink)
        graph.ep.overlap[re] = unsimp_graph.ep.overlap[e]
        graph.ep.flow[re] = unsimp_graph.ep.flow[e]
        graph.ep.color[re] = "black"
        red_edge_dict[(u, v)] = re

    return graph, red_node_dict, red_edge_dict


def flipped_gfa_to_graph(gfa_file, logger: Logger):
    """
    read flipped gfa format graph in.
    """
    logger.debug("Parsing GFA format graph")
    gfa = gfapy.Gfa().from_file(filename=gfa_file)
    logger.debug(
        "Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version)
    )

    graph = init_graph()
    red_node_dict = {}
    red_edge_dict = {}

    # S
    for line in gfa.segments:
        [_, seg_no, seg, dp] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        v = graph.add_vertex()
        graph.vp.seq[v] = seg
        graph.vp.dp[v] = dp_float
        graph.vp.id[v] = seg_no
        graph.vp.color[v] = "black"
        red_node_dict[seg_no] = v
    # L
    for edge in gfa.edges:
        [_, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
        source = red_node_dict[seg_no_l]
        sink = red_node_dict[seg_no_r]

        assert overlap_len[-1] == "M" and ori_l == ori_r
        re = graph.add_edge(source, sink)
        graph.ep.overlap[re] = int(overlap_len[:-1])
        graph.ep.color[re] = "black"
        red_edge_dict[(seg_no_l, seg_no_r)] = re

    return graph, red_node_dict, red_edge_dict


def graph_to_gfa(
    graph: Graph, simp_node_dict: dict, edge_dict: dict, logger: Logger, filename
):
    """
    store the swapped graph in simplifed_graph.
    """
    subprocess.check_call("touch {0}; echo > {0}".format(filename), shell=True)

    with open(filename, "w") as gfa:
        for v in simp_node_dict.values():
            if graph.vp.color[v] == "black":
                name = graph.vp.id[v]
                gfa.write(
                    "S\t{0}\t{1}\tDP:f:{2}\n".format(
                        name, graph.vp.seq[v], graph.vp.dp[v]
                    )
                )

        for (u, v), e in edge_dict.items():
            node_u = simp_node_dict[u] if u in simp_node_dict else None
            node_v = simp_node_dict[v] if v in simp_node_dict else None

            if node_u == None or node_v == None:
                continue
            if graph.vp.color[node_u] != "black" or graph.vp.color[node_v] != "black":
                continue
            if graph.ep.color[e] != "black":
                continue
            gfa.write(
                "L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format(
                    u, "+", v, "+", graph.ep.overlap[e]
                )
            )
        gfa.close()
    logger.info(filename + " is stored..")
    return 0


def is_valid(p: list, idx_mapping: dict, simp_node_dict: dict, simp_edge_dict: dict):
    if len(p) == 0:
        return False
    if len(p) == 1:
        if p[0] not in idx_mapping:
            return False
        if idx_mapping[p[0]] not in simp_node_dict:
            return False
        return True
    for i in range(len(p) - 1):
        if p[i] not in idx_mapping or p[i + 1] not in idx_mapping:
            return False
        mu = idx_mapping[p[i]]
        mv = idx_mapping[p[i + 1]]
        if mu not in simp_node_dict:
            return False
        if mv not in simp_node_dict:
            return False
        if (mu, mv) not in simp_edge_dict:
            return False
    return True


def spades_paths_parser(
    graph: Graph,
    simp_node_dict: dict,
    simp_edge_dict: dict,
    idx_mapping: dict,
    logger: Logger,
    path_file,
    min_len=250,
    min_cov=0,
    at_least_one_edge=False,
):
    """
    Map SPAdes's contig to the graph, return all the suitable contigs.
    """

    def get_paths(fd, path):
        subpaths = []
        total_nodes = 0
        while path.endswith(";\n"):
            subpath = str(path[:-2]).split(",")
            subpath = list(
                map(
                    lambda v: str(v[:-1]) if v[-1] == "+" else "-" + str(v[:-1]),
                    subpath,
                )
            )
            subpathred = list(dict.fromkeys(subpath))
            # validity check
            if is_valid(subpathred, idx_mapping, simp_node_dict, simp_edge_dict):
                subpath = list(map(lambda v: idx_mapping[v], subpath))
                subpaths.append(subpath)
                total_nodes += len(subpath)
            path = fd.readline()

        subpath = path.rstrip().split(",")
        subpath = list(
            map(lambda v: str(v[:-1]) if v[-1] == "+" else "-" + str(v[:-1]), subpath)
        )
        subpathred = list(dict.fromkeys(subpath))
        # validity check
        if is_valid(subpathred, idx_mapping, simp_node_dict, simp_edge_dict):
            subpath = list(map(lambda v: idx_mapping[v], subpath))
            subpaths.append(subpath)
            total_nodes += len(subpath)

        return subpaths, total_nodes

    logger.info("parsing SPAdes .paths file..")
    contig_dict = {}
    contig_info = {}
    try:
        with open(path_file, "r") as contigs_file:
            name = contigs_file.readline()
            path = contigs_file.readline()

            while name != "" and path != "":
                (cno, clen, ccov) = re.search(
                    "%s(.*)%s(.*)%s(.*)" % ("NODE_", "_length_", "_cov_"), name.strip()
                ).group(1, 2, 3)
                subpaths, total_nodes = get_paths(contigs_file, path)

                name_r = contigs_file.readline()
                path_r = contigs_file.readline()
                (cno_r, clen_r, ccov_r) = re.search(
                    "%s(.*)%s(.*)%s(.*)%s" % ("NODE_", "_length_", "_cov_", "'"),
                    name_r.strip(),
                ).group(1, 2, 3)
                subpaths_r, total_nodes_r = get_paths(contigs_file, path_r)

                if not (cno == cno_r and clen == clen_r and ccov == ccov_r):
                    raise BaseException

                # next contig group
                name = contigs_file.readline()
                path = contigs_file.readline()

                # pick one direction only
                (segments, total_n) = max(
                    [(subpaths, total_nodes), (subpaths_r, total_nodes_r)],
                    key=lambda t: t[1],
                )

                # filter contig
                if segments == []:
                    continue
                if float(ccov) <= min_cov:
                    continue
                if int(clen) < min_len and total_n < 2:
                    continue
                for i, subpath in enumerate(segments):
                    repeat_dict = {}
                    for k in subpath:
                        if k not in repeat_dict:
                            repeat_dict[k] = 1
                        else:
                            repeat_dict[k] += 1
                    subpath = list(dict.fromkeys(subpath))

                    if at_least_one_edge and len(subpath) == 1:
                        if (
                            simp_node_dict[subpath[0]].in_degree() == 0
                            and simp_node_dict[subpath[-1]].out_degree() == 0
                        ):
                            # isolated contig
                            continue
                    if len(segments) != 1:
                        contig_dict[cno + "$" + str(i)] = [
                            subpath,
                            path_len(graph, [simp_node_dict[id] for id in subpath]),
                            ccov,
                        ]
                        contig_info[cno + "$" + str(i)] = (None, repeat_dict)
                    else:
                        contig_dict[cno] = [subpath, int(clen), ccov]
                        contig_info[cno] = (None, repeat_dict)

            contigs_file.close()
    except BaseException as err:
        logger.error(
            err,
            "\nPlease make sure the correct SPAdes contigs .paths file is provided.",
        )
        logger.error("Pipeline aborted")
        sys.exit(1)
    logger.debug(str(contig_dict))
    logger.debug(str(contig_info))
    logger.info("done")
    return contig_dict, contig_info


def contig_dict_to_fasta(
    graph: Graph, simp_node_dict: dict, contig_dict: dict, output_file
):
    """
    Store contig dict into fastq file
    """
    subprocess.check_call("touch {0}; echo > {0}".format(output_file), shell=True)

    with open(output_file, "w") as fasta:
        for cno, (contig, clen, ccov) in sorted(
            contig_dict.items(), key=lambda x: x[1][1], reverse=True
        ):
            contig_name = (
                ">" + str(cno) + "_" + str(clen) + "_" + str(round(ccov, 2)) + "\n"
            )
            seq = path_ids_to_seq(graph, contig, contig_name, simp_node_dict) + "\n"
            fasta.write(contig_name)
            fasta.write(seq)
        fasta.close()


def contig_dict_to_path(
    contig_dict: dict, output_file, id_mapping: dict = None, keep_original=False
):
    """
    Store contig dict into paths file
    """
    subprocess.check_call("touch {0}; echo > {0}".format(output_file), shell=True)
    rev_id_mapping = {}
    if id_mapping != None:
        for id, map in id_mapping.items():
            rev_id_mapping[map] = id
    with open(output_file, "w") as paths:
        for cno, (contig, clen, ccov) in sorted(
            contig_dict.items(), key=lambda x: x[1][1], reverse=True
        ):
            contig_name = "NODE_" + str(cno) + "_" + str(clen) + "_" + str(ccov) + "\n"
            path_ids = ""
            for id in contig:
                if keep_original:
                    for iid in str(id).split("&"):
                        if iid.find("*") != -1:
                            rid = rev_id_mapping[iid[: iid.find("*")]]
                        else:
                            rid = rev_id_mapping[iid]
                        if rid[0] == "-":
                            rid = rid[1:] + "-"
                        path_ids += rid + ","
                else:
                    path_ids += str(id) + ","
            path_ids = path_ids[:-1] + "\n"
            paths.write(contig_name)
            paths.write(path_ids)
        paths.close()
