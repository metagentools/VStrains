#!/usr/bin/env python3

from graph_tool.all import Graph
from graph_tool.draw import graph_draw
import gfapy
import subprocess
import sys

import matplotlib.pyplot as plt
import seaborn
import pandas

import numpy


def gfa_to_graph(gfa_file, overlap, init_ori=1):
    """
    Convert assembly graph gfa file to graph
    Nodes: segment with corresponding 
    """

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa2').from_file(filename=gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))

    graph = Graph(directed=True)
    graph.vp.seq = graph.new_vertex_property("string", val="")
    graph.vp.dp = graph.new_vertex_property("double")
    graph.vp.kc = graph.new_vertex_property("int32_t")
    graph.vp.id = graph.new_vertex_property("string", val="UD")
    graph.vp.visited = graph.new_vertex_property("int16_t", val=0)
    graph.vp.ori = graph.new_vertex_property("int16_t") # 1 = +, -1 = -
    graph.vp.group = graph.new_vertex_property("int16_t", val=-1)
    graph.vp.partition = graph.new_vertex_property("int16_t", val=0)
    graph.vp.color = graph.new_vertex_property("string")

    graph.ep.overlap = graph.new_edge_property("int", val=0)
    graph.ep.visited = graph.new_edge_property("int", val=0)
    graph.ep.flow = graph.new_edge_property("double", val=0.0)
    graph.ep.color = graph.new_edge_property("string")

    # S
    node_dict = {}
    dp_dict = {}
    edge_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [_, seg_no, seg, dp, kc] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        kc_float = float(kc.split(":")[2])
        v_pos = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.dp[v_pos] = dp_float
        graph.vp.kc[v_pos] = kc_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.group[v_pos] = -1
        graph.vp.visited[v_pos] = -1
        graph.vp.partition[v_pos] = -1
        graph.vp.color[v_pos] = 'black'

        v_neg = graph.add_vertex()
        graph.vp.seq[v_neg] = reverse_seq(seg)
        graph.vp.dp[v_neg] = dp_float
        graph.vp.kc[v_neg] = kc_float
        graph.vp.id[v_neg] = seg_no
        graph.vp.ori[v_neg] = -1
        graph.vp.group[v_neg] = -1
        graph.vp.visited[v_neg] = -1
        graph.vp.partition[v_neg] = -1
        graph.vp.color[v_neg] = 'black'
        

        node_dict[seg_no] = (v_pos, v_neg)
        dp_dict[seg_no] = dp_float
    # L
    for edge in gfa.edges:
        [_, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
        u_pos, u_neg = node_dict[seg_no_l]
        v_pos, v_neg = node_dict[seg_no_r]
        u = u_pos if ori_l == '+' else u_neg
        v = v_pos if ori_r == '+' else v_neg
        e = graph.add_edge(source=u, target=v)
        # gfa format check
        assert overlap_len[-1] == 'M'
        graph.ep.overlap[e] = int(overlap_len[:-1])
        assert graph.ep.overlap[e] == overlap
        graph.ep.color[e] = 'black'

        edge_dict[(seg_no_l, graph.vp.ori[u], seg_no_r, graph.vp.ori[v])] = e
        
    # P
    # for path in gfa.paths:
    #     [line_type, path_no, seg_names, seg_overlap] = str(path).split("\t")
    graph, simp_node_dict, simp_edge_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict, init_ori)
    red_graph, red_node_dict, red_edge_dict = reduce_graph(graph, simp_node_dict, simp_edge_dict)
    return red_graph, red_node_dict, red_edge_dict

def reduce_graph(unsimp_graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    graph = Graph(directed=True)

    graph.vp.seq = graph.new_vertex_property("string", val="")
    graph.vp.dp = graph.new_vertex_property("double")
    graph.vp.kc = graph.new_vertex_property("int32_t")
    graph.vp.id = graph.new_vertex_property("string", val="UD")
    graph.vp.color = graph.new_vertex_property("string")

    graph.ep.overlap = graph.new_edge_property("int", val=0)
    graph.ep.flow = graph.new_edge_property("float", val=0.0)
    graph.ep.color = graph.new_edge_property("string")

    red_node_dict = {}
    red_edge_dict = {}

    for no, node in simp_node_dict.items():
        v = graph.add_vertex()
        graph.vp.seq[v] = unsimp_graph.vp.seq[node]
        graph.vp.dp[v] = unsimp_graph.vp.dp[node]
        graph.vp.kc[v] = unsimp_graph.vp.kc[node]
        graph.vp.id[v] = unsimp_graph.vp.id[node]
        graph.vp.color[v] = 'black'
        red_node_dict[no] = v
    
    for (u,v), e in simp_edge_dict.items():
        source = red_node_dict[u]
        sink = red_node_dict[v]

        re = graph.add_edge(source, sink)
        graph.ep.overlap[re] = unsimp_graph.ep.overlap[e]
        graph.ep.flow[re] = unsimp_graph.ep.flow[e]
        graph.ep.color[re] = 'black'
        red_edge_dict[(u,v)] = re
    
    return graph, red_node_dict, red_edge_dict

def flipped_gfa_to_graph(gfa_file):
    """
    read flipped gfa format graph in.
    """
    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa2').from_file(filename=gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))

    graph = Graph(directed=True)
    graph.vp.seq = graph.new_vertex_property("string", val="")
    graph.vp.dp = graph.new_vertex_property("double")
    graph.vp.kc = graph.new_vertex_property("int32_t")
    graph.vp.id = graph.new_vertex_property("string", val="UD")
    graph.vp.color = graph.new_vertex_property("string")

    graph.ep.overlap = graph.new_edge_property("int", val=0)
    graph.ep.flow = graph.new_edge_property("float", val=0.0)
    graph.ep.color = graph.new_edge_property("string")

    red_node_dict = {}
    red_edge_dict = {}

    # S
    for line in gfa.segments:
        [_, seg_no, seg, dp, kc] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        kc_float = float(kc.split(":")[2])
        v = graph.add_vertex()
        graph.vp.seq[v] = seg
        graph.vp.dp[v] = dp_float
        graph.vp.kc[v] = kc_float
        graph.vp.id[v] = seg_no
        graph.vp.color[v] = 'black'
        red_node_dict[seg_no] = v
    # L
    for edge in gfa.edges:
        [_, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
        source = red_node_dict[seg_no_l]
        sink = red_node_dict[seg_no_r]

        assert overlap_len[-1] == 'M' and ori_l == ori_r
        re = graph.add_edge(source, sink)
        graph.ep.overlap[re] = int(overlap_len[:-1])
        graph.ep.color[re] = 'black'
        red_edge_dict[(seg_no_l,seg_no_r)] = re
    
    return graph, red_node_dict, red_edge_dict

def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    for (u,v), e in simp_edge_dict.items():

        if graph.ep.color[e] != 'black':
            #forbidden edge
            continue

        u_node = simp_node_dict[u]
        u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])

        v_node = simp_node_dict[v]
        v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])

        flow = max((graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node])
        graph.ep.flow[e] = max(graph.ep.flow[e], flow)
    return

def graph_simplification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, min_cov):
    """
    Directly remove all the vertex with coverage less than minimum coverage and related edge

    Node belongs to any contigs should not be removed
    return:
        removed_node_dict
        removed_edge_dict
    """
    print("-------------------------graph simplification----------------------")
    print("Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    removed_node_dict = {}
    removed_edge_dict = {}
    # iterate until no more node be removed from the graph
    for id, node in list(simp_node_dict.items()):
        if graph.vp.dp[node] < min_cov:
            if id in node_to_contig_dict:
                continue

            graph_remove_vertex(graph, simp_node_dict, id, printout=False)
            removed_node_dict[id] = node

            for e in node.all_edges():
                uid = graph.vp.id[e.source()]
                vid = graph.vp.id[e.target()]
                if (uid, vid) in edge_to_contig_dict:
                    continue 
                if (uid, vid) in simp_edge_dict:
                    graph_remove_edge(graph, simp_edge_dict, uid, vid, printout=False)
                    removed_edge_dict[(uid, vid)] = e

    print("Remain: Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    print("-------------------------graph simplification end----------------------")
    return removed_node_dict, removed_edge_dict

def graph_to_gfa(graph: Graph, simp_node_dict: dict, edge_dict: dict, filename):
    """
    store the swapped graph in simplifed_graph.
    """
    subprocess.check_call("touch {0}".format(
    filename), shell=True)

    with open(filename, 'w') as gfa:
        for v in simp_node_dict.values():
            if graph.vp.color[v] == 'black':
                name = graph.vp.id[v]
                gfa.write("S\t{0}\t{1}\tDP:f:{2}\tKC:i:{3}\n".format
                (name, graph.vp.seq[v], graph.vp.dp[v], graph.vp.kc[v]))

        for (u,v), e in edge_dict.items():
            node_u = simp_node_dict[u] if u in simp_node_dict else None
            node_v = simp_node_dict[v] if v in simp_node_dict else None

            if node_u == None or node_v == None:
                continue
            if graph.vp.color[node_u] != 'black' or graph.vp.color[node_v] != 'black':
                continue
            if graph.ep.color[e] != 'black':
                continue
            gfa.write("L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format(u, "+", v, "+", graph.ep.overlap[e]))
        gfa.close()
    print(filename, " is stored..")
    return 0

def flip_graph_bfs(graph: Graph, node_dict: dict, edge_dict: dict, dp_dict: dict, init_ori=1):
    """
    Flip all the node orientation.

    return an node_dict, which only contains one orientation per node for simplicity.
    rename all the used node to positive, and forbidden the opponent node.
    """
    def source_node_via_dp(dp_dict: dict):
        """
        return the pos-neg node with maximum depth
        """
        seg_no = max(dp_dict, key=dp_dict.get)
        return seg_no

    def reverse_edge(graph: Graph, edge, node_dict: dict, edge_dict: dict):
        """
        reverse an edge with altered orientation and direction.
        """
        tmp_s = edge.source()
        tmp_t = edge.target()
        
        edge_dict.pop((graph.vp.id[tmp_s], graph.vp.ori[tmp_s], graph.vp.id[tmp_t], graph.vp.ori[tmp_t]))

        tmp_s_pos, tmp_s_neg = node_dict[graph.vp.id[tmp_s]]
        tmp_t_pos, tmp_t_neg = node_dict[graph.vp.id[tmp_t]]
        s = tmp_t_pos if graph.vp.ori[tmp_t] == -1 else tmp_t_neg
        t = tmp_s_pos if graph.vp.ori[tmp_s] == -1 else tmp_s_neg

        o = graph.ep.overlap[edge]
        graph.remove_edge(edge)
        e = graph.add_edge(s, t)
        graph.ep.overlap[e] = o
        edge_dict[(graph.vp.id[s], graph.vp.ori[s], graph.vp.id[t], graph.vp.ori[t])] = e

        return graph, e, edge_dict
    print("flip graph orientation..")
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
                pick_dict[graph.vp.id[u]] = '+'
                # print_vertex(graph, v_neg, "node to reverse")
                for e in list(v_neg.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            else:
                u = v_neg
                pick_dict[graph.vp.id[u]] = '-'
                # print_vertex(graph, v_pos, "node to reverse")
                for e in list(v_pos.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            
            graph.vp.visited[v_pos] = 1
            graph.vp.visited[v_neg] = 1
            # add further nodes into the fifo_queue TODO, all or out only
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1:
                    graph.vp.visited[adj_node] = 0
                    fifo_queue.append([node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]])

    # verify sorted graph
    print("verify graph..")
    check = True

    assert len(pick_dict) == len(node_dict)
    for key, item in pick_dict.items():
        v_pos, v_neg = node_dict[key]
        if item == '+':
            if v_neg.in_degree() + v_neg.out_degree() > 0:
                print_vertex(graph, v_neg, "pick error found")
                check = False
        else:
            if v_pos.in_degree() + v_pos.out_degree() > 0:
                print_vertex(graph, v_pos, "pick error found")
                check = False

    for key, (v_pos, v_neg) in node_dict.items():
        if not v_pos.in_degree() + v_pos.out_degree() == 0 and not v_neg.in_degree() + v_neg.out_degree() == 0:
            check = False
            print_vertex(graph, v_pos, "erroroness node pos found")
            print_vertex(graph, v_neg, "erroroness node neg found")
            print("re-pict nodes, pick both")
            pick_dict[key] = 0
            graph.vp.id[v_pos] = graph.vp.id[v_pos] + "a"
            graph.vp.id[v_neg] = graph.vp.id[v_neg] + "b"

    assert check
    print("Graph is verified")

    simp_node_dict = {}
    for seg_no, pick in pick_dict.items():
        if pick == '+':
            picked = node_dict[seg_no][0]
        else:
            picked = node_dict[seg_no][1]
        graph.vp.ori[picked] = 1
        simp_node_dict[seg_no] = picked

    simp_edge_dict = simplify_edge_dict(edge_dict)
    print("done")
    return graph, simp_node_dict, simp_edge_dict

def map_ref_to_graph(ref_file, simp_node_dict: dict, graph_file, store_mapping=False, output_file="overlap.paf", fasta_file="temp_gfa_to_fasta.fasta"):
    """
    map reference strain to the graph, debug only
    """
    if not ref_file:
        print("No ref file imported")
        return -1
    with open(graph_file, 'r') as gfa:
        with open(fasta_file, 'w') as fasta:
            for Line in gfa:
                splited = Line.split('\t')
                if splited[0] == 'S':
                    quality = "B"*len(splited[2])
                    fasta.write("@{0}\n{1}\n+\n{2}\n".format(splited[1],splited[2],quality))
            fasta.close()
        gfa.close()
    
    subprocess.check_call("minimap2 {0} {1} -c > {2}".format(ref_file, fasta_file, output_file), shell=True)

    strain_dict = {}
    with open(output_file, 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
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
            if ((nmatch/nblock) == 1):
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = []
                strain_dict[ref_no].append(seg_no)
        paf.close()
        
    if not store_mapping:
        subprocess.check_call("rm {0}; rm {1}".format(output_file, fasta_file), shell=True)
    
    print("strain dict mapping")
    for seg_no, strains in strain_dict.items():
        print("strains: ", seg_no, " Path: ", list_to_string(strains))
        print("-------------------")
    return strain_dict

def map_ref_to_contig(contig_dict: dict, paf_file):
    print("map ref to contig")
    strain_dict = {}
    with open(paf_file, 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            seg_no = str(splited[0].split('_')[0])
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = str(splited[5])
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if seg_no not in contig_dict:
                continue
            if ((nmatch/nblock) >= 0.999):
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = set()
                strain_dict[ref_no].add(seg_no)
        paf.close()

    for sno, cnos in strain_dict.items():
        print("--------------------------------->")
        print("contig-strains: ", sno, "Count: ", len(cnos), "-Contigs: ", [(cno1, clen, ccov) for cno1, [_, clen, ccov] in contig_dict.items() if cno1 in cnos])
        rel_nodes = []
        for cno in cnos:
            rel_nodes.extend(contig_dict[cno][0])
        print("related nodes in contig: ", list_to_string(rel_nodes))

def gfa_to_fasta(gfa_file, fasta_file):
    if not gfa_file:
        print("No gfa file imported")
        return -1
    subprocess.check_call("touch {0}".format(fasta_file), shell=True)
    with open(gfa_file, 'r') as gfa:
        with open(fasta_file, 'w') as fasta:
            for Line in gfa:
                splited = Line.split('\t')
                if splited[0] == 'S':
                    fasta.write(">{0}\n{1}\n".format(splited[1],splited[2]))
            fasta.close()
        gfa.close()

def minimap_api(ref_file, fasta_file, output_file):
    subprocess.check_call("minimap2 {0} {1} -c > {2}".format(
        ref_file, fasta_file, output_file), shell=True)
    return  

def trim_contig_dict(graph: Graph, simp_node_dict: dict, contig_dict: dict, overlap):
    for cno, [contig, _, ccov] in list(contig_dict.items()):
        involved_node_set = set()
        new_contig = []
        for id in contig:
            if id not in involved_node_set:
                new_contig.append(id)
                involved_node_set.add(id)
        contig_dict[cno] = [new_contig, path_len(graph, [simp_node_dict[no] for no in new_contig], overlap), ccov]
    return contig_dict

def contig_dict_to_fasta(graph: Graph, simp_node_dict: dict, contig_dict: dict, overlap_len, output_file):
    """
    Store contig dict into fastq file
    """
    subprocess.check_call("touch {0}".format(
    output_file), shell=True)

    with open(output_file, 'w') as fasta:
        for cno, (contig, clen, ccov) in contig_dict.items():
            contig_name = ">" + str(cno) + "_" + str(clen) + "_" + str(round(ccov, 2)) + "\n"
            seq = path_ids_to_seq(graph, contig, contig_name, simp_node_dict, overlap_len) + "\n"
            fasta.write(contig_name)
            fasta.write(seq)
        fasta.close()

def contig_dict_to_path(contig_dict: dict, output_file):
    """
    Store contig dict into paths file
    """
    subprocess.check_call("touch {0}".format(output_file), shell=True)
    with open(output_file, 'w') as paths:
        for cno, (contig, clen, ccov) in sorted(contig_dict.items(), key=lambda x: x[1][2]):
            contig_name = "NODE_" + str(cno) + "_" + str(clen) + "_" + str(ccov) + "\n"
            path_ids = ""
            for id in contig:
                path_ids += str(id) + ","
            path_ids = path_ids[:-1]  + "\n"
            paths.write(contig_name)
            paths.write(path_ids)
        paths.close()

def get_contig(contig_file, simp_node_dict: dict, simp_edge_dict: dict, min_len=250):
    """
    Map SPAdes's contig to the graph, return all the contigs with length >= 250 or minimum 2 nodes in contig.
    """
    print("processing contigs..")
    if not contig_file:
        print("contig file not imported")
        return -1
    contig_dict = {}
    with open(contig_file, 'r') as contigs_file:
        while True:
            name = contigs_file.readline()
            seg_nos = contigs_file.readline()
            name_r = contigs_file.readline()
            seg_nos_r = contigs_file.readline()

            #TODO
            if seg_nos.find(';') != -1:
                print("find gaps in contig file, TODO")
                continue

            if not name or not seg_nos or not name_r or not seg_nos_r: 
                break

            split_name = name.split('_')
            cno = str(split_name[1])
            clen = int(split_name[3])
            ccov = float(split_name[5])

            # contig from both orientation
            contigs = [n[:-1] if n[-1] in ['-', '+'] else n[:-2] for n in seg_nos.split(',')]
            contigs_rev = [n[:-1] if n[-1] in ['-', '+'] else n[:-2] for n in seg_nos_r.split(',')]
            contig_len = len(contigs)

            # contig filter
            # use as less in-confident contigs as possible.
            if clen < min_len and len(contigs) < 2:
                continue

            if contig_len > 1:
                i = 0
                c = []
                pick = False
                while not pick:
                    e1 = (contigs[i],contigs[i+1])
                    e1_r = (contigs_rev[i], contigs_rev[i+1])
                    i = i + 1
                    if e1 not in simp_edge_dict and e1_r not in simp_edge_dict:
                        print("edge is not exist in both direction, skip contig: ", cno)
                        break
                    elif e1 not in simp_edge_dict:
                        c = contigs_rev[:]
                        pick = True
                        # print("pick forward side for contig: ", cno)
                    elif e1_r not in simp_edge_dict:
                        c = contigs[:]
                        pick = True
                        # print("pick reverse side for contig: ", cno)
                    else:
                        print("both direction edge, error edge case, skip contig: ", cno)
                        break
                    
                    if not pick and i == contig_len - 1:
                        # still not pick until last edge
                        print("all the edge is removed, no pick until last point, skip contig: ", cno)
                        break

                if not pick:
                    # whole contig is reduced already, no split chance, potential error
                    continue
                else:
                    contig_dict[cno] = [c, clen, ccov]
                    for i in range(len(c)):
                        c_i = c[i]
                        c_i_1 = c[i+1] if (i < len(c) - 1) else None
                        if c_i not in simp_node_dict:
                            print("node {0} not in contig {1}, error".format(c_i, cno))

                        if c_i_1 != None and c_i_1 not in simp_node_dict:
                            print("node {0} not in contig {1}, error".format(c_i_1, cno))
                        
                        if c_i_1 != None and (c_i, c_i_1) not in simp_edge_dict:
                            print("edge {0} not in contig {1}, error".format((c_i, c_i_1), cno))
            else:
                c = contigs
                if c[0] in simp_node_dict:
                    contig_dict[cno] = [c, clen, ccov]
        contigs_file.close()
    print("done")
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
            c_i_1 = c[i+1] if (i < len(c) - 1) else None
            if c_i_1 != None:
                if (c_i, c_i_1) not in edge_to_contig_dict:
                    edge_to_contig_dict[(c_i, c_i_1)] = {cno}
                else:
                    edge_to_contig_dict[(c_i, c_i_1)].add(cno)
    return node_to_contig_dict, edge_to_contig_dict

def clear_flow(graph: Graph, simp_edge_dict: dict):
    for _, e in simp_edge_dict.items():
        graph.ep.flow[e] = 0.0

def contig_dict_fix(graph: Graph, simp_node_dict: dict, contig_dict: dict, overlap):
    """
    fix the contig dict, reassign the contig coverage to minimum used **edge flow** coverage,
    we not use the minimum used **node** since the min node may still be shared by multiple contigs.
    however, min edge flow may still not be the true contig coverage

    """
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
                sublen = path_len(graph, [simp_node_dict[c] for c in subc], overlap)
                contig_dict[cno + "^" + str(i)] = [subc, sublen, ccov]
    return

def contig_cov_fix(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, printout=False):
    """
    if found a single node contig that also not appearing in the simp_node_dict, check if is mapped to split contig
    """
    for cno, [contig, clen, _] in list(contig_dict.items()):
        if printout:
            print("---------------------------------------------------------------")
        newccov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
        contig_dict[cno][2] = newccov
        if cno in contig_dict:
            if printout:
                print_contig(cno, clen, contig_dict[cno][2], contig)
                eflow = contig_flow(graph, simp_edge_dict, contig)
                if len(contig) > 1:
                    print(eflow)
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
        e = graph.edge(cand_path[i], cand_path[i+1])
        # print(e, graph.vp.id[cand_path[i]], graph.vp.id[cand_path[i+1]])
        graph.ep.flow[e] -= cand_cov

def is_splitter_branch(node):
    return node.in_degree() > 1 and node.out_degree() > 1
def is_trivial_branch(node):
    return (node.in_degree() == 1 and node.out_degree() > 1) or ((node.in_degree() > 1 and node.out_degree() == 1))

def graph_splitting(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, b0, b1, threshold, overlap):
    """
    n-n branch splitting, # FIXME add more restrict rule to avoid false positive split
    """
    print("-------------------------graph split----------------------")
    split_branches = []
    node_to_contig_dict, _ = contig_map_node(contig_dict)
    no_mapping = {}
    for no, node in list(simp_node_dict.items()):
        ine = [e for e in node.in_edges() if graph.ep.color[e] == 'black']
        oute = [e for e in node.out_edges() if graph.ep.color[e] == 'black']
        support_contigs = []
        # non-trivial branch
        if len(ine) > 1 and len(oute) > 1:
            threshold2 = 2*abs(b0 + b1*graph.vp.dp[node])
            support_contigs = [cno for cno in node_to_contig_dict[no]] if no in node_to_contig_dict else []
            # print_vertex(graph, node, "---------- branch node, support by contig {0}".format(support_contigs))
            # print("threshold2: ", threshold2)
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
                        cproduct.append((ie, oe, numpy.mean([graph.ep.flow[ie], graph.ep.flow[oe]]), delta))
                        ine_usage[ie] += 1
                        oute_usage[oe] += 1

            if cproduct != []:
                for i, (ie, oe, subcov, delta) in enumerate(cproduct):
                    # print("---------------")
                    # print_edge(graph, ie, "in")
                    # print_edge(graph, oe, "out")
                    prev_node = ie.source()
                    prev_no = graph.vp.id[prev_node]
                    next_node = oe.target()
                    next_no = graph.vp.id[next_node]
                    if ((prev_no, no) not in simp_edge_dict or (no, next_no) not in simp_edge_dict):
                        # print("edge has been removed already")
                        continue
                    # print("Delta: ", delta)
                    involved_contigs = []
                    is_conf = False
                    cross_talk = False
                    for cno in support_contigs:
                        contig, clen, ccov = contig_dict[cno]
                        if len(contig) < 3:
                            # print("unnecessary contig (len < 3) check: ", contig)
                            if prev_no in contig or next_no in contig:
                                involved_contigs.append(cno)
                        elif prev_no in contig and next_no in contig:
                            is_conf = True
                            involved_contigs.append(cno)
                            # print_contig(cno, clen, ccov, contig, "support contig ")
                        elif prev_no in contig and next_no not in contig:
                            # print("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        elif prev_no not in contig and next_no in contig:
                            # print("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        else:
                            None
                    involved_contig_sumcov = sum([contig_dict[cno][2] for cno in involved_contigs])
                    # print("involved contig cov sum: {0}, triple-node mean flow: {1}".format(involved_contig_sumcov, subcov))
                    # print("in edge usage: {0}, out edge usage: {1}".format(ine_usage[ie], oute_usage[oe]))
                    if cross_talk:
                        None
                        # print("- current branch split forbidden, cross talk", prev_no, no, next_no)
                    elif not is_conf and (ine_usage[ie] > 1 or oute_usage[oe] > 1):
                        None
                        # print("- current branch split forbidden, no supporting contig path, ambigous split", prev_no, no, next_no)
                    else:
                        # print("- branch split performed")
                        split_branches.append(no)

                        subid = "X" + no + "X" + str(i)
                        sub_node = graph_add_vertex(graph, simp_node_dict, subid, subcov, graph.vp.seq[node], graph.vp.kc[node])

                        if no not in no_mapping:
                            no_mapping[no] = []
                        no_mapping[no].append(str(subid))

                        graph.vp.dp[node] -= subcov

                        graph_remove_edge(graph, simp_edge_dict, prev_no, no)
                        graph_remove_edge(graph, simp_edge_dict, no, next_no)
                        graph_add_edge(graph, simp_edge_dict, prev_node, prev_no, sub_node, subid, graph.ep.overlap[ie],
                        graph.ep.flow[ie])
                        graph_add_edge(graph, simp_edge_dict, sub_node, subid, next_node, next_no, graph.ep.overlap[oe],
                        graph.ep.flow[oe])

                        for icno in involved_contigs:
                            if contig_dict[icno][0].count(no) == 1:
                                contig_dict[icno][0][contig_dict[icno][0].index(no)] = subid
                                support_contigs.remove(icno)
                            else:
                                print("contig error, previous node {0} is not in contig or multiple occurance in contig {1}, potential bug".format(no, icno))
    
    # remove all the isolated low cov node&edge
    for node in list(graph.vertices()):
        alle = [e for e in node.all_edges() if graph.ep.color[e] == 'black']
        if graph.vp.dp[node] < threshold:
            graph_remove_vertex(graph, simp_node_dict, graph.vp.id[node], "remove isolated low cov node")
            for e in alle:
                graph_remove_edge(graph, simp_edge_dict, graph.vp.id[e.source()], graph.vp.id[e.target()])
    
    # fix single node contig
    for cno, [contig, clen, _] in list(contig_dict.items()):
        if len(contig) <= 1 and contig[0] not in simp_node_dict:
            contig_dict.pop(cno)
            # print("isolated contig node {0}, prepared to pop, check any replacement".format(cno))
            if contig[0] in no_mapping:
                # print("mapping: {0} -> {1}".format(contig[0], no_mapping[contig[0]]))
                s = 'A'
                for i, mapped in enumerate(no_mapping[contig[0]]):
                    if mapped in simp_node_dict:
                        mapped_node = simp_node_dict[mapped]
                        contig_dict[cno+chr(ord(s) + i)] = [[mapped], path_len(graph, [mapped_node], overlap), graph.vp.dp[mapped_node]]
                        # print("mapped cno: ", cno+chr(ord(s) + i), mapped)  
    
    print("No of branch be removed: ", len(set(split_branches)))
    print("Split branches: ", list_to_string(set(split_branches)))
    return len(set(split_branches))

def graph_split_trivial(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Split the graph, for any (0|1)->N, N->(0|1) branch, split by forking the 1 edge to N edge.
    """
    print("-------------------------graph trivial split----------------------")
    has_split = True
    trivial_split_count = 0
    id_mapping = {}
    for id in simp_node_dict.keys():
        id_mapping[id] = set()
    while has_split:
        has_split = False
        for id in list(simp_node_dict.keys()):
            node = simp_node_dict[id]
            if graph.vp.color[node] != 'black':
                continue
            if id not in id_mapping:
                id_mapping[id] = set()
            ines = [ue for ue in node.in_edges() if graph.ep.color[ue] == 'black']
            outes = [ve for ve in node.out_edges() if graph.ep.color[ve] == 'black']
            if len(ines) == 0 and len(outes) > 1:
                # print(id, len(ines), len(outes), "current node is source node, split left")
                graph.vp.color[node] = 'gray'
                # create len(outes) subnodes
                s = 'A'
                for i in range(len(outes)):
                    oute = outes[i]
                    tgt = oute.target()
                    snode = graph_add_vertex(graph, simp_node_dict, id+chr(ord(s) + i), graph.ep.flow[oute], graph.vp.seq[node], graph.vp.kc[node])
                    graph.ep.color[oute] = 'gray'
                    sedge = graph_add_edge(graph, simp_edge_dict, snode, graph.vp.id[snode], 
                    tgt, graph.vp.id[tgt], graph.ep.overlap[oute], graph.ep.flow[oute])
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[(graph.vp.id[sedge.source()], graph.vp.id[sedge.target()])] = sedge
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            elif len(ines) > 1 and len(outes) == 0:
                # print(id, len(ines), len(outes), "current node is sink node, split right")
                graph.vp.color[node] = 'gray'
                # create len(ines) subnodes
                s = 'A'
                for i in range(len(ines)):
                    ine = ines[i]
                    src = ine.source()
                    snode = graph_add_vertex(graph, simp_node_dict, id+chr(ord(s) + i), graph.ep.flow[ine], graph.vp.seq[node], graph.vp.kc[node])
                    graph.ep.color[ine] = 'gray'
                    sedge = graph_add_edge(graph, simp_edge_dict, src, graph.vp.id[src], 
                    snode, graph.vp.id[snode], graph.ep.overlap[ine], graph.ep.flow[ine])
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[(graph.vp.id[sedge.source()], graph.vp.id[sedge.target()])] = sedge
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            elif len(ines) == 1 and len(outes) > 1:
                # print(id, len(ines), len(outes), "split left")
                graph.vp.color[node] = 'gray'
                ine = ines[0]
                src = ine.source()
                graph.ep.color[ine] = 'gray'
                s = 'A'
                for i in range(len(outes)):
                    oute = outes[i]
                    tgt = oute.target()
                    snode = graph_add_vertex(graph, simp_node_dict, id+chr(ord(s) + i), graph.ep.flow[oute], graph.vp.seq[node], graph.vp.kc[node])
                    graph.ep.color[oute] = 'gray'
                    sedge_out = graph_add_edge(graph, simp_edge_dict, snode, graph.vp.id[snode], 
                    tgt, graph.vp.id[tgt], graph.ep.overlap[oute], graph.ep.flow[oute])
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[(graph.vp.id[sedge_out.source()], graph.vp.id[sedge_out.target()])] = sedge_out

                    sedge_in = graph_add_edge(graph, simp_edge_dict, src, graph.vp.id[src],
                    snode, graph.vp.id[snode], graph.ep.overlap[ine], graph.ep.flow[oute])
                    simp_edge_dict[(graph.vp.id[sedge_in.source()], graph.vp.id[sedge_in.target()])] = sedge_in
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            elif len(ines) > 1 and len(outes) == 1:
                # print(id, len(ines), len(outes), "split right")   
                graph.vp.color[node] = 'gray'
                oute = outes[0]
                tgt = oute.target()
                graph.ep.color[oute] = 'gray'
                s = 'A'
                for i in range(len(ines)):
                    ine = ines[i]
                    src = ine.source()
                    snode = graph_add_vertex(graph, simp_node_dict, id+chr(ord(s) + i), graph.ep.flow[ine], graph.vp.seq[node], graph.vp.kc[node])
                    graph.ep.color[ine] = 'gray'
                    sedge_in = graph_add_edge(graph, simp_edge_dict, src, graph.vp.id[src], 
                    snode, graph.vp.id[snode], graph.ep.overlap[ine], graph.ep.flow[ine])
                    simp_node_dict[graph.vp.id[snode]] = snode
                    simp_edge_dict[(graph.vp.id[sedge_in.source()], graph.vp.id[sedge_in.target()])] = sedge_in

                    sedge_out = graph_add_edge(graph, simp_edge_dict, snode, graph.vp.id[snode],
                    tgt, graph.vp.id[tgt], graph.ep.overlap[oute], graph.ep.flow[ine])
                    simp_edge_dict[(graph.vp.id[sedge_out.source()], graph.vp.id[sedge_out.target()])] = sedge_out
                    id_mapping[id].add(graph.vp.id[snode])
                has_split = True
                trivial_split_count += 1
            else:
                None
                # print(id, len(ines), len(outes), ", unable to split, skip")
    print("No of trivial branch be removed: ", trivial_split_count)
    return trivial_split_count, id_mapping

def contig_dict_remapping(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, id_mapping: dict, prev_ids: list, overlap):
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
                        acc_paths.append((p+[next]))
                else:
                    for nextm in id_mappingP[next]:
                        if (last, nextm) in simp_edge_dict:
                            acc_paths.append((p+[nextm]))
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
    print("contig resolution..")
    # id_mapping merging, recursive merge down.
    red_id_mapping = {}

    for id in prev_ids:
        all_set = merge_id(id_mapping, id_mapping[id], id)
        red_id_mapping[id] = all_set
        # print("Node {0} maps to {1}".format(id, all_set))

    for cno, (contig, clen, ccov) in list(contig_dict.items()):
        # print("---------------------------------------------")
        # print("Current mapping contig: ", cno, list_to_string(contig))
        paths = map_contig_tree(cno, contig, red_id_mapping)
        # split the contig tree to avoid the ambiguity variation
        if len(paths) < 1:
            print("error, contig missed: ", cno, contig)
        elif len(paths) == 1:
            if paths[0] == contig:
                None
                # print("single mapping, keep original")
            else:
                # print("single mapping, replace", list_to_string(paths[0]))
                contig_dict.pop(cno)
                subcov = path_cov(graph, simp_node_dict, simp_edge_dict, paths[0])
                contig_dict[cno] = [paths[0], path_len(graph, [simp_node_dict[no] for no in paths[0]], overlap), subcov]
        else:
            contig_dict.pop(cno)
            # print("multi mapping for the current contig: whole contig is ambiguous mapping", cno)
            for i, path in enumerate(paths):
                dupcno = cno + "^" + str(i)
                if dupcno in contig_dict:
                    print("dup cno: ", dupcno, " already exist, error")
                else:
                    subcov = path_cov(graph, simp_node_dict, simp_edge_dict, path)
                    contig_dict[dupcno] = [path, path_len(graph, [simp_node_dict[no] for no in path], overlap), subcov]
                    # print("duplicated mapped contig: ", dupcno, list_to_string(path))
    print("done")
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
        if graph.vp.id[src] not in simp_node_dict or graph.vp.id[target] not in simp_node_dict:
            continue
        if src.out_degree() == 1 and target.in_degree() == 1:
            assert src != target
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

def simple_paths_to_dict(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, overlap):
    simple_paths = simp_path(graph, simp_node_dict, simp_edge_dict)
    simp_path_dict = {}
    for id, p in enumerate(simple_paths):
        pids = [graph.vp.id[n] for n in p]
        name = str(id)
        clen = path_len(graph, p, overlap)
        cov = numpy.max([graph.vp.dp[n] for n in p])
        simp_path_dict[name] = [pids, clen, cov]
        # print("Simple PATH: ", list_to_string(pids))
    return simp_path_dict

def simp_path_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
    """
    reduce all the contig to a single node, and keep all the potential src/tgt edge.

    1. reduce the coverage for each involving node by the amount of contig cov
    2. reconnect end-to-end nodes to the contig node
    """

    simp_path_dict = simple_paths_to_dict(graph, simp_node_dict, simp_edge_dict, overlap)

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
            id += "_" + nid
        cseq = path_to_seq(graph_backup, [simp_node_dict_backup[n] for n in contig], cno, overlap)
        kc = numpy.median([graph_backup.vp.kc[simp_node_dict_backup[u]] for u in contig])
        in_edges = list((graph_backup.vp.id[e.source()], src) for e in simp_node_dict_backup[src].in_edges())
        out_edges = list((tgt, graph_backup.vp.id[e.target()],) for e in simp_node_dict_backup[tgt].out_edges())
        

        for i in range(len(contig)):
            no = contig[i]
            node_to_simp_node[no] = id
            graph_remove_vertex(graph, simp_node_dict, no, printout=False)
            if i != len(contig) - 1:
                graph_remove_edge(graph, simp_edge_dict, contig[i], contig[i+1], printout=False)
        cv = graph_add_vertex(graph, simp_node_dict, id, ccov, cseq, kc, printout=False)
        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])
    
    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u,v) in in_edges:
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                graph_add_edge(graph, simp_edge_dict, simp_node_dict[u], u, node, graph.vp.id[node], overlap, printout=False)
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict:
                    graph_add_edge(graph, simp_edge_dict, in_node, graph.vp.id[in_node], node, graph.vp.id[node], overlap, printout=False)            

        for (u,v) in out_edges:
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                graph_add_edge(graph, simp_edge_dict, node, graph.vp.id[node], simp_node_dict[v], v, overlap, printout=False)
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict:
                    graph_add_edge(graph, simp_edge_dict, node, graph.vp.id[node], out_node, graph.vp.id[out_node], overlap, printout=False)
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
        # print("cno: {0} from {1} to {2}".format(cno, list_to_string(contig), list_to_string(new_contig)))
        contig_dict[cno] = [new_contig, path_len(graph, [simp_node_dict[no] for no in new_contig], overlap), ccov]
    
    return

def coincident_min_node(graph: Graph, simp_edge_dict: dict, cno1, contig1: list, cno2, contig2: list):
    """
    return true if contig1 and contig2 's min cov node aligned to same relative index, known both contigs are not equal
    """
    def min_cov_edge(graph: Graph, simp_edge_dict: dict, contig: list):
        """
        return index (i), where contig[i,i+1] is the min edge along the contig
        """
        flows = contig_flow(graph, simp_edge_dict, contig)
        minflow = sys.maxsize
        minIndex = -1
        for i in range(len(contig) - 1):
            if flows[i] < minflow:
                minflow = flows[i]
                minIndex = i
        return minIndex
    intersect = set(contig1).intersection(set(contig2))
    min_ind_1 = min_cov_edge(graph, simp_edge_dict, contig1)
    min_ind_2 = min_cov_edge(graph, simp_edge_dict, contig2)
    # edge_eq = contig1[min_ind_1:min_ind_1+2] == contig2[min_ind_2:min_ind_2+2] 
    if len(intersect) <= 0:
        return (None,None)
    elif len(intersect) == len(contig1) and len(intersect) == len(contig2):
        # total duplicated contig
        # print("duplicated")
        return (cno2, cno1)
    elif len(intersect) == len(contig1):
        # print(min_ind_1, min_ind_2)
        # contig1 is covered by contig2
        if len(contig1) == 1:
            if contig1[0] in contig2[min_ind_2:min_ind_2+2]:
                return (cno1, cno2)
        else:
            if min_ind_2 == contig2.index(contig1[0]) + min_ind_1:
                return (cno1, cno2)
    elif len(intersect) == len(contig2):
        # contig2 is covered by contig1
        # print(min_ind_1, min_ind_2)
        if len(contig2) == 1:
            if contig2[0] in contig1[min_ind_1:min_ind_1+2]:
                return (cno2, cno1)
        else:
            if min_ind_1 == contig1.index(contig2[0]) + min_ind_2:
                return (cno2, cno1)
    return (None,None)

def contig_dup_removed(graph: Graph, simp_edge_dict: dict, contig_dict: dict):
    print("drop duplicated contigs..")
    dup_contig_ids = set()
    for cno in contig_dict.keys():
        contig, _, _ = contig_dict[cno]       
        for cno2 in contig_dict.keys():
            if cno not in dup_contig_ids and cno2 not in dup_contig_ids and cno != cno2:
                contig2, _, _ = contig_dict[cno2]
                # use set equality to avoid cyclic contig
                (del_cno,keep_cno) = coincident_min_node(graph, simp_edge_dict, cno, contig, cno2, contig2)
                if del_cno != None:
                    # print("coincident min edge, del covered contig: ", del_cno, contig_dict[del_cno][0], " from long contig: ", keep_cno, contig_dict[keep_cno][0])
                    dup_contig_ids.add(del_cno)
    for cno in dup_contig_ids:
        contig_dict.pop(cno)
    print("done")
    return contig_dict

def contig_dup_removed_s(contig_dict: dict):
    print("drop duplicated contigs..")
    dup_contig_ids = set()
    for cno1 in contig_dict.keys():
        contig1, _, _ = contig_dict[cno1]       
        for cno2 in contig_dict.keys():
            if cno1 not in dup_contig_ids and cno2 not in dup_contig_ids and cno1 != cno2:
                contig2, _, _ = contig_dict[cno2]
                # use set equality to avoid cyclic contig
                intersect = set(contig1).intersection(set(contig2))
                if len(intersect) == len(contig1) and len(intersect) == len(contig2):
                    # duplicated
                    # print("duplicated contig: ", cno1, contig1, cno2, contig2)
                    dup_contig_ids.add(cno2)
                elif len(intersect) == len(contig1):
                    # contig1 is subcontig of contig2
                    # print("covered contig: ", cno1, contig1, cno2, contig2)
                    dup_contig_ids.add(cno1)
                elif len(intersect) == len(contig2):
                    # contig2 is subcontig of contig1
                    # print("covered contig: ", cno1, contig1, cno2, contig2)
                    dup_contig_ids.add(cno2)
    for cno in dup_contig_ids:
        contig_dict.pop(cno)
    print("done")
    return contig_dict

def concat_overlap_contig(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
    print("concat overlapped contig..")
    overlaps = {}
    for cno in contig_dict.keys():
        overlaps[cno] = ([], None)
    for cno in list(contig_dict.keys()):
        for cno2 in list(contig_dict.keys()):
            if cno == cno2:
                continue
            isParallel, intersects, status = check_contig_intersection(cno, contig_dict[cno][0], cno2, contig_dict[cno2][0])
            if not isParallel:
                if status == 'f':
                    # forward e2e overlap
                    if path_len(graph, overlaps[cno][0], overlap) < path_len(graph, [simp_node_dict[k] for k in intersects], overlap):
                        overlaps[cno] = (intersects, cno2)
                elif status == 'd':
                    if cno not in [c for _, c in overlaps[cno2]]:
                        if path_len(graph, overlaps[cno][0], overlap) < path_len(graph, [simp_node_dict[k] for k in intersects], overlap):
                            overlaps[cno] = (intersects, cno2)
    
    # print(overlaps)
    for cno in list(contig_dict.keys()):
        if cno in overlaps and overlaps[cno][1] != None:
            i, c = overlaps.pop(cno)
            concat_plan = [(None, cno)]
            while c != None:
                concat_plan.append((i, c))
                i, c = overlaps.pop(c)
            concat_contig = []
            cnos = ""
            print("plan: ", concat_plan)
            for ind, (intersect, ccno) in enumerate(concat_plan):
                contig, _, _ = contig_dict.pop(ccno)
                s = ccno + "c"
                if ind != 0:
                    contig = contig[len(intersect):]
                    s = ccno
                concat_contig.extend(contig)
                cnos += s
            print("concat e2e contig: ", cnos, "after concat: ", concat_contig)
            contig_dict[cnos] = [concat_contig, path_len(graph, [simp_node_dict[k] for k in concat_contig], overlap), path_cov(graph, simp_node_dict, simp_edge_dict, concat_contig)]

    return

def check_contig_intersection(cno, contig, cno2, contig2):
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
    # print("intersect check: {0} vs {1}, count: {2}".format(cno, cno2, len(intersect)))
    if len(intersect) <= 0:
        return False, intersect, 'n'

    if len(intersect) == len(contig):
        # print("{0} is covered by {1}".format(cno, cno2))
        return True, intersect, 'o'

    if len(intersect) == len(contig2):
        # print("{0} is covered by {1}".format(cno2, cno))
        return True, intersect, 'o'

    intermediate_nodes_index = [False for _ in contig]
    for i in [contig.index(e) for e in intersect]:
        intermediate_nodes_index[i] = True
    # print("cno: {0} intersect at: {1}".format(cno, list(enumerate(intermediate_nodes_index))))
    if not intermediate_nodes_index[0] and not intermediate_nodes_index[-1]:
        # print("intersection in the middle")
        return True, intersect, 'o'
    prev_false_index = intermediate_nodes_index.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index)):
        if not intermediate_nodes_index[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                # print("intersection in the middle")
                return True, intersect, 'o'

    intermediate_nodes_index2 = [False for _ in contig2]
    for i in [contig2.index(e) for e in intersect]:
        intermediate_nodes_index2[i] = True
    # print("cno2: {0} intersect at: {1}".format(cno2, list(enumerate(intermediate_nodes_index2))))
    if not intermediate_nodes_index2[0] and not intermediate_nodes_index2[-1]:
        # print("intersection in the middle")
        return True, intersect, 'o'
    prev_false_index = intermediate_nodes_index2.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index2)):
        if not intermediate_nodes_index2[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                # print("intersection in the middle")
                return True, intersect, 'o'
    direction = None
    if intermediate_nodes_index[0]:
        direction = 'b'
    if intermediate_nodes_index[-1]:
        direction = 'f' if direction == None else 'd'
    # print("overlap end-to-end")
    return False, intersect, direction

def path_len(graph: Graph, path, overlap):
    """
    Find length of the linear path.
    """
    lens = [len(graph.vp.seq[u]) for u in path]
    return sum(lens) - overlap * (len(lens) - 1) if len(lens) > 0 else 0

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
    for i in range(len(contig)-1):
        edges.append((contig[i], contig[i+1]))
        
    return edges

def contig_flow(graph: Graph, edge_dict: dict, contig):
    """
    edge flow for the contig
    """
    edge_flow = []
    if len(contig) < 2:
        return edge_flow
    for i in range(len(contig)-1):
        e = edge_dict[(contig[i],contig[i+1])]
        f = graph.ep.flow[e]
        edge_flow.append(f)

    return edge_flow

def path_ids_to_seq(graph: Graph, path_ids: list, path_name, simp_node_dict: dict, overlap_len):
    seq = ""
    for i in range(len(path_ids)):
        u = simp_node_dict[path_ids[i]]
        if i == len(path_ids) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq

def path_to_seq(graph: Graph, path: list, path_name, overlap_len):
    seq = ""
    for i in range(len(path)):
        u = path[i]
        if i == len(path) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq

def reverse_seq(seq: str):
    return ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[x] for x in seq[::-1])

def simplify_edge_dict(edge_dict: dict):
    simp_edge_dict = {}
    for (u, _, v, _), e in edge_dict.items():
        simp_edge_dict[(u,v)] = e
    return simp_edge_dict

def swap_node_ori_name(graph: Graph, node_dict: dict, seg_no):
    """
    swap the node's orientation, and update related edges
    """
    v_pos, v_neg = node_dict[seg_no]
    graph.vp.ori[v_pos] = -1
    graph.vp.ori[v_neg] = 1
    node_dict[seg_no] = (v_neg, v_pos)

    return graph, node_dict

def graph_stat(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    print("-------------------------graph stat----------------------")
    for v in simp_node_dict.values():
        print_vertex(graph, v, "stat")
    for e in simp_edge_dict.values():
        print_edge(graph, e, "stat")
    
    print("-----------------------graph stat end--------------------")

def graph_to_dict(graph: Graph):
    simp_node_dict = {}
    simp_edge_dict = {}
    for node in graph.vertices():
        simp_node_dict[graph.vp.id[node]] = node
    for edge in graph.edges():
        simp_edge_dict[(graph.vp.id[edge.source()], graph.vp.id[edge.target()])] = edge
    return simp_node_dict, simp_edge_dict

def graph_color_other_to_gray(graph: Graph, simp_node_dict: dict, ids: set):
    gray_ids = set(simp_node_dict.keys()).difference(ids)
    for id in gray_ids:
        graph.vp.color[simp_node_dict[id]] = 'gray'
    return 

def graph_add_vertex(graph: Graph, simp_node_dict: dict, id, dp, seq, kc, s="add vertex", color='black', printout=False):
    node = graph.add_vertex()
    graph.vp.id[node] = id
    graph.vp.dp[node] = dp
    graph.vp.seq[node] = seq
    graph.vp.kc[node] = kc
    graph.vp.color[node] = color
    simp_node_dict[id] = node
    if printout:
        print_vertex(graph, node, s)
    return node

def graph_remove_vertex(graph, simp_node_dict: dict, id, s="remove vertex", color='gray', printout=False):
    node = simp_node_dict[id]
    graph.vp.color[node] = color
    simp_node_dict.pop(id)
    if printout:
        print_vertex(graph, node, s)
    return node

def graph_add_edge(graph: Graph, simp_edge_dict: dict, src, src_id, tgt, tgt_id, overlap, flow=0, s="add edge", color='black', printout=False):
    edge = graph.add_edge(src, tgt)
    graph.ep.overlap[edge] = overlap
    graph.ep.color[edge] = color
    graph.ep.flow[edge] = flow
    simp_edge_dict[(src_id, tgt_id)] = edge
    if printout:
        print_edge(graph, edge, s)
    return edge

def graph_remove_edge(graph: Graph, simp_edge_dict: dict, src_id, tgt_id, s="remove edge", color='gray', printout=False):
    edge = simp_edge_dict.pop((src_id, tgt_id))
    graph.ep.color[edge] = color
    if printout:
        print_edge(graph, edge, s)
    return edge       

def draw_graph_api(graph: Graph, output_file):
    output_size = 120 * (graph.num_edges() + graph.num_vertices())
    vsize= 30
    esize = 30
    graph_draw(g=graph, output=output_file, bg_color="white", 
    vertex_size=vsize,
    output_size=(output_size, output_size))  

def print_edge(graph, e, s=""):
    print(s, " edge: ", graph.vp.id[e.source()], "->", graph.vp.id[e.target()], graph.ep.flow[e], graph.ep.color[e])

def print_vertex(graph, v, s=""):
    print(s, " vertex: ", graph.vp.id[v], ", dp: ", graph.vp.dp[v], ", kc: ", graph.vp.kc[v], ", in_degree: ", v.in_degree(), ", out_degree: ", v.out_degree(), graph.vp.color[v])

def print_contig(cno, clen, ccov, contig, s=""):
    print(s, " Contig: ", cno, ", length: ", clen, ", cov: ", ccov, "Path: ", list_to_string(contig))

def list_to_string(ids: list, s=""):
    string = s + " - "
    for id in ids:
        string += str(id) + ", "
    return string[:-2] if len(string) >= 2 else ""

def path_to_id_string(graph: Graph, path, s=""):
    return list_to_string([graph.vp.id[node] for node in path], s)

def reindexing(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, extended_contig_dict: dict):
    idx_mapping = {}
    idx_node_dict = {}
    idx_edge_dict = {}
    idx_contig_dict = {}
    idx = 0
    for no, node in simp_node_dict.items():
        if graph.vp.color[node] == 'black':
            idx_mapping[no] = str(idx)
            graph.vp.id[node] = str(idx)
            idx_node_dict[str(idx)] = node
            print("Node: {0} maps to {1}, seqlen: {2}".format(list_to_string(str(no).split('_')), idx, len(graph.vp.seq[node])))
            idx += 1
    for (u, v), e in simp_edge_dict.items():
        if graph.ep.color[e] == 'black' and graph.vp.color[e.source()] == 'black' and graph.vp.color[e.target()] == 'black':
            idx_edge_dict[(idx_mapping[u], idx_mapping[v])] = e

    for cno, [contig, clen, ccov] in extended_contig_dict.items():
        idx_contig_dict[cno] = [[idx_mapping[no] for no in contig], clen, ccov]
        print("indexed contig: ", cno, list_to_string(idx_contig_dict[cno][0]))

    return graph, idx_node_dict, idx_edge_dict, idx_contig_dict

def add_global_source_sink(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, overlap, store_dict=False):
    # find all the srcs-targets
    src_nodes = [node for node in graph.vertices() if node.in_degree() == 0]
    tgt_nodes = [node for node in graph.vertices() if node.out_degree() == 0]
    
    # create a global source node & sink node, and concat all the curr src to the global source node
    global_src = graph.add_vertex()
    graph.vp.id[global_src] = 'global_src'
    graph.vp.dp[global_src] = 0
    graph.vp.color[global_src] = 'black'
    if store_dict:
        simp_node_dict[graph.vp.id[global_src]] = global_src
    for src in src_nodes:
        e = graph.add_edge(global_src, src)
        graph.ep.flow[e] = graph.vp.dp[src]
        graph.ep.color[e] = 'black'
        graph.ep.overlap[e] = overlap
        graph.vp.dp[global_src] += graph.ep.flow[e]
        if store_dict:
            simp_edge_dict[(graph.vp.id[global_src], graph.vp.id[src])] = e
    # print("srcs:  ", list_to_string([graph.vp.id[src] for src in src_nodes]))

    global_sink = graph.add_vertex()
    graph.vp.id[global_sink] = 'global_sink'
    graph.vp.dp[global_sink] = 0
    graph.vp.color[global_sink] = 'black'
    if store_dict:
        simp_node_dict[graph.vp.id[global_sink]] = global_sink
    for tgt in tgt_nodes:
        e = graph.add_edge(tgt, global_sink)
        graph.ep.flow[e] = graph.vp.dp[tgt]
        graph.ep.color[e] = 'black'
        graph.ep.overlap[e] = overlap
        graph.vp.dp[global_sink] += graph.ep.flow[e]
        if store_dict:
            simp_edge_dict[(graph.vp.id[tgt], graph.vp.id[global_sink])] = e
    # print("sinks: ", list_to_string([graph.vp.id[tgt] for tgt in tgt_nodes]))
    return global_src, global_sink

def remove_global_source_sink(graph: Graph, gs, gt):
    del_v = [gs, gt]
    for v in reversed(sorted(del_v)):
        graph.remove_vertex(v)
    return

#######################################################################################################
#######################################PATH FINDING ALGORITHM##########################################
#######################################################################################################

def graph_is_DAG(graph: Graph, simp_node_dict: dict):
    """
    check if the graph is a DAG, advanced to check all the isolated subgraph by all mean
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
            if graph.ep.color[e] != 'black':
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
    
    #init
    visited = {}
    recStack = {}
    for node in simp_node_dict.values():
        if graph.vp.color[node] == 'black':
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

def reachable(graph: Graph, simp_node_dict: dict, src, tgt):
    """
    determine whether src can possibly reach the tgt
    """
    visited = {}
    for no in simp_node_dict.keys():
        visited[no] = False

    count_down = 1 if src != tgt else 2
    queue = [src]
    while queue:
        curr = queue.pop()
        visited[graph.vp.id[curr]] = True
        if curr == tgt:
            count_down -= 1
            if count_down == 0:
                return True
            else:
                visited[graph.vp.id[curr]] = False

        for oute in curr.out_edges():
            out = oute.target()
            if not visited[graph.vp.id[out]]:
                queue.append(out)
    return False

def paths_from_src(graph: Graph, simp_node_dict: dict, self_node, src, overlap, maxlen):
    """
    retrieve all the path from src node to any node 
    within maxlen restriction, in straight direction
    """
    def dfs_rev(graph: Graph, u, curr_path: list, maxlen, visited, all_path):
        visited[u] = True
        curr_path.append(u)
        curr_len = path_len(graph, curr_path, overlap)
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

def paths_to_tgt(graph: Graph, simp_node_dict: dict, self_node, tgt, overlap, maxlen):
    """
    retrieve all the path from any node to tgt node
    within maxlen restriction, in reverse direction
    """
    def dfs_rev(graph: Graph, v, curr_path: list, maxlen, visited, all_path):
        visited[v] = True
        curr_path.insert(0, v)
        curr_len = path_len(graph, curr_path, overlap)
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