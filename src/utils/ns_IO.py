#!/usr/bin/env python3

from graph_tool.all import Graph
from graph_tool.draw import graph_draw
import gfapy
import subprocess
import sys

from utils.ns_Utilities import reverse_seq, print_vertex, path_ids_to_seq, print_edge

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
def gfa_to_graph(gfa_file, init_ori=1):
    """
    Convert assembly graph gfa file to graph
    Nodes: segment with corresponding 
    """

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa().from_file(filename=gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))

    graph = init_graph()
    graph.vp.visited = graph.new_vertex_property("int16_t", val=0)
    graph.vp.ori = graph.new_vertex_property("int16_t") # 1 = +, -1 = -

    graph.ep.visited = graph.new_edge_property("int", val=0)

    # S
    node_dict = {}
    dp_dict = {}
    edge_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [t, seg_no, seg] = (str(line).split("\t"))[:3]
        tags = (str(line).split("\t"))[3:]
        dp = [tag for tag in tags if tag.startswith('dp') or tag.startswith("DP")][0]
        assert dp[:2] == 'DP' or dp[:2] == 'dp'
        assert t == 'S'
        dp_float = float(dp.split(":")[2])
        v_pos = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.dp[v_pos] = dp_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.visited[v_pos] = -1
        graph.vp.color[v_pos] = 'black'

        v_neg = graph.add_vertex()
        graph.vp.seq[v_neg] = reverse_seq(seg)
        graph.vp.dp[v_neg] = dp_float
        graph.vp.id[v_neg] = seg_no
        graph.vp.ori[v_neg] = -1
        graph.vp.visited[v_neg] = -1
        graph.vp.color[v_neg] = 'black'
        

        node_dict[seg_no] = (v_pos, v_neg)
        dp_dict[seg_no] = dp_float
    # L
    for edge in gfa.edges:
        [t, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = (str(edge).split("\t"))[:6]
        assert t == 'L' or t == 'C' or t == 'E'
        if seg_no_l == seg_no_r:
            print("self-loop edge, skip: ", seg_no_l, ori_l, "->", seg_no_r, ori_r)
            continue
        u_pos, u_neg = node_dict[seg_no_l]
        v_pos, v_neg = node_dict[seg_no_r]
        u = u_pos if ori_l == '+' else u_neg
        v = v_pos if ori_r == '+' else v_neg
        e = graph.add_edge(source=u, target=v)
        # gfa format check
        assert overlap_len[-1] == 'M'
        graph.ep.overlap[e] = int(overlap_len[:-1])
        graph.ep.color[e] = 'black'

        edge_dict[(seg_no_l, graph.vp.ori[u], seg_no_r, graph.vp.ori[v])] = e
        
    # P
    # for path in gfa.paths:
    #     [line_type, path_no, seg_names, seg_overlap] = str(path).split("\t")
    graph, simp_node_dict, simp_edge_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict, init_ori)
    red_graph, red_node_dict, red_edge_dict = reduce_graph(graph, simp_node_dict, simp_edge_dict)
    return red_graph, red_node_dict, red_edge_dict

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
        return max(dp_dict, key=dp_dict.get)

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
                # print_vertex(graph, v_neg, "node to reverse pos")
                for e in set(v_neg.all_edges()):
                    graph, r_e, edge_dict = reverse_edge(graph, e, node_dict, edge_dict)
                    # print_edge(graph, r_e, "after reverse: ")
            else:
                u = v_neg
                pick_dict[graph.vp.id[u]] = '-'
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
                    fifo_queue.append([node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]])

    # verify sorted graph
    print("final verifying graph..")
    assert len(pick_dict) == len(node_dict)
    for key, item in pick_dict.items():
        v_pos, v_neg = node_dict[key]
        if item == '+':
            if v_neg.in_degree() + v_neg.out_degree() > 0:
                print_vertex(graph, v_neg, "pick ambiguous found")
                print("force selection, remove opposite edges")
                for e in set(v_neg.all_edges()):
                    print_edge(graph, e, "force edge deletion")
                    tmp_s = e.source()
                    tmp_t = e.target()
                    print(graph.vp.ori[tmp_s], graph.vp.ori[tmp_t])
                    edge_dict.pop((graph.vp.id[tmp_s], graph.vp.ori[tmp_s], 
                        graph.vp.id[tmp_t], graph.vp.ori[tmp_t]))
                    graph.remove_edge(e)
        else:
            if v_pos.in_degree() + v_pos.out_degree() > 0:
                print_vertex(graph, v_pos, "pick ambiguous found")
                for e in set(v_pos.all_edges()):
                    print_edge(graph, e, "force edge deletion")
                    tmp_s = e.source()
                    tmp_t = e.target()
                    print(graph.vp.ori[tmp_s], graph.vp.ori[tmp_t])
                    edge_dict.pop((graph.vp.id[tmp_s], graph.vp.ori[tmp_s], 
                        graph.vp.id[tmp_t], graph.vp.ori[tmp_t]))
                    graph.remove_edge(e)
    print("Graph is verified")

    simp_node_dict = {}
    for seg_no, pick in pick_dict.items():
        if pick == '+':
            picked = node_dict[seg_no][0]
        else:
            picked = node_dict[seg_no][1]
        graph.vp.ori[picked] = 1
        simp_node_dict[seg_no] = picked

    simp_edge_dict = {}
    for (u, _, v, _), e in edge_dict.items():
        simp_edge_dict[(u,v)] = e
    print("done")
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
    gfa = gfapy.Gfa().from_file(filename=gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))

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
                gfa.write("S\t{0}\t{1}\tDP:f:{2}\n".format(name, graph.vp.seq[v], graph.vp.dp[v]))

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

def flye_info2path(info_file, tempdir):
    print("translating flye {0} into .paths format..".format(info_file))
    contig_file = tempdir + "/tmp/contigs.paths"
    with open(info_file, 'r') as cfile:
        with open(contig_file, 'w') as ofile:
            cfile.readline() # header
            lines = cfile.readlines()
            for line in lines:
                seg_name, length, cov, circ, repeat, mult, alt_group, graph_path = (line[:-1]).split('\t')
                header = str(seg_name) + "_length_" + str(length) + "_cov_" + str(cov)
                path = ""
                ps = graph_path.split(",")
                ps = list(dict.fromkeys(ps))
                for n in ps:
                    if n != '*' and n != '??':
                        # consider the info segment name start with - if neg orientation
                        if n[0] == '-':
                            path += "edge_" + n[1:] + '-' + ","
                        else:
                            path += "edge_" + n[:] + '+' + ","
                path = path[:-1]
                ofile.write(header + "\n")
                ofile.write(path + "\n")

                rev_path = ""
                for n in reversed(ps):
                    if n != '*' and n != '??':
                        # consider the info segment name start with - if neg orientation
                        if n[0] == '-':
                            rev_path += "edge_" + n[1:] + '+' + ","
                        else:
                            rev_path += "edge_" + n[:] + '-' + ","
                rev_path = rev_path[:-1]
                ofile.write(header + "\'\n")
                ofile.write(rev_path + "\n")
            ofile.close()
        cfile.close()
    print("done")
    return contig_file

def get_contig(contig_file, simp_node_dict: dict, simp_edge_dict: dict, min_len=250, at_least_one_edge=False):
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
            split_name = (name[:-1]).split('_')
            cno = str(split_name[1])
            clen = int(split_name[3])
            ccov = float(split_name[5])

            # contig from both orientation
            contigs = [n[:-1] if n[-1] in ['-', '+'] else n for n in seg_nos[:-1].split(',')]
            contigs_rev = [n[:-1] if n[-1] in ['-', '+'] else n for n in seg_nos_r[:-1].split(',')]
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
                        print("pick forward side for contig: ", cno)
                    elif e1_r not in simp_edge_dict:
                        c = contigs[:]
                        pick = True
                        print("pick reverse side for contig: ", cno)
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
                    skip = False
                    for i in range(len(c)-1):
                        c_i = c[i]
                        c_i_1 = c[i+1]
                        if c_i not in simp_node_dict:
                            print("node {0} not in contig {1}, skip".format(c_i, cno))
                            skip = True
                        if c_i_1 not in simp_node_dict:
                            print("node {0} not in contig {1}, skip".format(c_i_1, cno))
                            skip = True
                        if (c_i, c_i_1) not in simp_edge_dict:
                            print("edge {0} not in contig {1}, skip".format((c_i, c_i_1), cno))
                            skip = True
                        if skip:
                            break
                    if not skip:
                        contig_dict[cno] = [c, clen, ccov]
            else:
                c = contigs
                if c[0] in simp_node_dict:
                    if not at_least_one_edge:
                        contig_dict[cno] = [c, clen, ccov]
                    else:
                        # check connectivity
                        if simp_node_dict[c[0]].in_degree() > 0 or simp_node_dict[c[0]].out_degree() > 0:
                            contig_dict[cno] = [c, clen, ccov]
        contigs_file.close()
    print("done, total input contig: ", len(contig_dict))
    return contig_dict

def gfa_to_fasta(gfa_file, fasta_file):
    if not gfa_file:
        print("No gfa file imported")
        return -1
    subprocess.check_call("touch {0}".format(fasta_file), shell=True)
    with open(gfa_file, 'r') as gfa:
        with open(fasta_file, 'w') as fasta:
            for Line in gfa:
                splited = (Line[:-1]).split('\t')
                if splited[0] == 'S':
                    fasta.write(">{0}\n{1}\n".format(splited[1],splited[2]))
            fasta.close()
        gfa.close()

def contig_dict_to_fasta(graph: Graph, simp_node_dict: dict, contig_dict: dict, output_file):
    """
    Store contig dict into fastq file
    """
    subprocess.check_call("touch {0}".format(
    output_file), shell=True)

    with open(output_file, 'w') as fasta:
        for cno, (contig, clen, ccov) in contig_dict.items():
            contig_name = ">" + str(cno) + "_" + str(clen) + "_" + str(round(ccov, 2)) + "\n"
            seq = path_ids_to_seq(graph, contig, contig_name, simp_node_dict) + "\n"
            fasta.write(contig_name)
            fasta.write(seq)
        fasta.close()

def contig_dict_to_path(contig_dict: dict, output_file, keep_original=False):
    """
    Store contig dict into paths file
    """
    subprocess.check_call("touch {0}".format(output_file), shell=True)
    with open(output_file, 'w') as paths:
        for cno, (contig, clen, ccov) in sorted(contig_dict.items(), key=lambda x: x[1][2]):
            contig_name = "NODE_" + str(cno) + "_" + str(clen) + "_" + str(ccov) + "\n"
            path_ids = ""
            for id in contig:
                if keep_original:
                    for iid in str(id).split('_'):
                        if iid.find('*') != -1:
                            path_ids += iid[:iid.find('*')] + ","
                        else:
                            path_ids += iid + ","
                else:
                    path_ids += str(id) + ","
            path_ids = path_ids[:-1]  + "\n"
            paths.write(contig_name)
            paths.write(path_ids)
        paths.close()