#!/usr/bin/env python3

from graph_tool.all import Graph
import gfapy
import subprocess

import numpy


def gfa_to_graph(gfa_file, init_ori):
    """
    Convert assembly graph gfa file to graph
    Nodes: segment with corresponding 
    """

    print("Parsing GFA format graph")
    gfa = gfapy.Gfa(version='gfa2').from_file(filename=gfa_file)
    print("Parsed gfa file length: {0}, version: {1}".format(len(gfa.lines), gfa.version))

    graph = Graph(directed=True)

    vprop_seq = graph.new_vertex_property("string", val="")
    vprop_dp = graph.new_vertex_property("double")
    vprop_kc = graph.new_vertex_property("int32_t")
    vprop_id = graph.new_vertex_property("string", val="UD")
    vprop_visited = graph.new_vertex_property("int16_t", val=0)
    vprop_ori = graph.new_vertex_property("int16_t") # 1 = +, -1 = -
    vprop_group = graph.new_vertex_property("int16_t", val=-1)
    vprop_partition = graph.new_vertex_property("int16_t", val=0)
    vprop_color = graph.new_vertex_property("string")

    graph.vp.seq = vprop_seq
    graph.vp.dp = vprop_dp
    graph.vp.kc = vprop_kc
    graph.vp.id = vprop_id
    graph.vp.visited = vprop_visited
    graph.vp.ori = vprop_ori
    graph.vp.group = vprop_group
    graph.vp.partition = vprop_partition
    graph.vp.color = vprop_color

    eprop_overlap = graph.new_edge_property("int", val=0)
    eprop_visited = graph.new_edge_property("int", val=0)
    eprop_flow = graph.new_edge_property("float", val=0.0)
    eprop_color = graph.new_edge_property("string")

    graph.ep.overlap = eprop_overlap
    graph.ep.visited = eprop_visited
    graph.ep.flow = eprop_flow
    graph.ep.color = eprop_color

    # graph.list_properties()
    # S
    node_dict = {}
    dp_dict = {}
    edge_dict = {}
    for line in gfa.segments:
        # segment, convert into Node^- and Node^+
        [line_type, seg_no, seg, dp, kc] = str(line).split("\t")
        dp_float = float(dp.split(":")[2])
        kc_float = float(kc.split(":")[2])
        v_pos = graph.add_vertex()
        v_neg = graph.add_vertex()
        graph.vp.seq[v_pos] = seg
        graph.vp.dp[v_pos] = dp_float
        graph.vp.kc[v_pos] = kc_float
        graph.vp.id[v_pos] = seg_no
        graph.vp.ori[v_pos] = 1
        graph.vp.group[v_pos] = -1
        graph.vp.visited[v_pos] = -1
        graph.vp.partition[v_pos] = -1
        graph.vp.color[v_pos] = 'black'
        
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
        [line_type, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
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
        
    # # P
    # for path in gfa.paths:
    #     [line_type, path_no, seg_names, seg_overlap] = str(path).split("\t")
    graph, simp_node_dict, simp_edge_dict = flip_graph_bfs(graph, node_dict, edge_dict, dp_dict, init_ori)
    return graph, simp_node_dict, simp_edge_dict

def graph_to_gfa(graph: Graph, simp_node_dict: dict, edge_dict: dict, min_cov, filename):
    """
    store the swapped graph in simplifed_graph.
    """
    subprocess.check_call("echo "" > {0}".format(
    filename), shell=True)
    
    def to_ori(i):
        return '+' if i == 1 else '-'
    
    id_record = []
    dp_record = []

    with open(filename, 'w') as gfa:
        for v in simp_node_dict.values():
            if graph.vp.color[v] == 'black' and graph.vp.dp[v] >= min_cov:
                
                id_record.append(int(graph.vp.id[v]))
                dp_record.append((int(graph.vp.id[v]), int(graph.vp.dp[v])))
                
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
            if graph.vp.dp[node_u] < min_cov or graph.vp.dp[node_v] < min_cov or graph.ep.flow[e] < min_cov:
                continue

            gfa.write("L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format
            (u, to_ori(graph.vp.ori[node_u]), 
            v, to_ori(graph.vp.ori[node_v]), 
            graph.ep.overlap[e]))
        gfa.close()
    print("--------------------", filename, "is stored--------------------")
    # print(id_record)
    # print(dp_record)
    return 0

def flip_graph_bfs(graph: Graph, node_dict: dict, edge_dict: dict, dp_dict: dict, init_ori):
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
        print("source node id: ", seg_no, ", depth: ", dp_dict[seg_no])
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
    print("-------------------------flip graph orientation----------------------")
    pick_dict = {}
    while set(dp_dict):
        seg_no = source_node_via_dp(dp_dict)
        source_pos, source_neg = node_dict[seg_no]
        graph.vp.visited[source_pos] = 0
        graph.vp.visited[source_neg] = 0
        queue = []
        queue.append([node_dict[seg_no], init_ori]) 

        while queue:
            (v_pos, v_neg), ori = queue.pop()
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
            # add further nodes into the queue TODO, all or out only
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1:
                    graph.vp.visited[adj_node] = 0
                    queue.append([node_dict[graph.vp.id[adj_node]], graph.vp.ori[adj_node]])

    # verify sorted graph
    print("-------------------------verify graph----------------------")
    check = True

    check = check and (len(pick_dict) == len(node_dict))
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

    check = check and len(node_dict) == len(pick_dict)
    if check: print("Graph is verified")
    print("-------------------------end verify------------------------")

    simp_node_dict = {}
    for seg_no, pick in pick_dict.items():
        if pick == '+':
            picked = node_dict[seg_no][0]
        else:
            picked = node_dict[seg_no][1]
        graph.vp.ori[picked] = 1
        simp_node_dict[seg_no] = picked

    simp_edge_dict = simplify_edge_dict(edge_dict)
    print("-------------------------flip graph orientation end------------------")
    return graph, simp_node_dict, simp_edge_dict

# FIXME fix the path
def map_ref_to_graph(ref_file, simp_node_dict: dict, simp_file):
    """
    map reference strain to the graph, debug only
    assumption: graph is stored in acc/simplifed_graph, 
    """
    if not ref_file:
        print("No ref file imported")
        return -1
    with open(simp_file, 'r') as gfa:
        with open("acc/gfa_to_fq.fq", 'w') as fq:
            for Line in gfa:
                splited = Line.split('\t')
                if splited[0] == 'S':
                    quality = "B"*len(splited[2])
                    fq.write("@{0}\n{1}\n+\n{2}\n".format(splited[1],splited[2],quality))
            fq.close()
        gfa.close()
    
    # minimap2. you may need to replace the exec path to minimap to fit your case
    subprocess.check_call("/Users/luorunpeng/opt/miniconda3/envs/vg-flow-env/bin/minimap2 {0} {1} > {2}".format(
    ref_file, "acc/gfa_to_fq.fq", "acc/overlap.paf"), shell=True)

    strain_dict = {}
    with open("acc/overlap.paf", 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            seg_no_int = int(splited[0])
            seg_no = splited[0]
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = splited[5]
            if seg_no not in simp_node_dict:
                continue
            if ((seg_f - seg_s) / seg_l) >= 0.8:
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = []
                strain_dict[ref_no].append(seg_no_int)
    subprocess.check_call("rm {0}".format("acc/gfa_to_fq.fq"), shell=True)
    
    print("strain dict mapping")
    for seg_no, strains in strain_dict.items():
        print("strains: ", seg_no, " Path: ", strains)
        print("-------------------")
    return strain_dict

def contig_overlap(contig_dict: dict, contig_file):
    """
    Construct an overlap matrix among all the non-full length contig
    """
    subprocess.check_call("/Users/luorunpeng/opt/miniconda3/envs/vg-flow-env/bin/minimap2 {0} {1} > {2}".format(
        contig_file, contig_file, "acc/contig_overlap.paf"), shell=True)
    return

def contig_dict_to_fq(graph: Graph, contig_dict: dict, simp_node_dict: dict, overlap_len, output_file, min_len=0):
    """
    Store contig dict into fastq file
    """
    subprocess.check_call("echo "" > {0}".format(
    output_file), shell=True)

    with open(output_file, 'w') as fq:
        for cno, (contig, clen, ccov) in contig_dict.items():
            contig_name = ">" + str(cno) + "_" + str(clen) + "_" + str(ccov) + "\n"
            seq = contig_to_seq(graph, contig, contig_name, simp_node_dict, overlap_len) + "\n"
            fq.write(contig_name)
            fq.write(seq)
        fq.close()


def get_contig(graph: Graph, contig_file, simp_node_dict: dict, simp_edge_dict: dict, min_cov, min_len, overlap, min_node=2):
    """
    Map SPAdes's contig to the graph, return all the contigs.
    if nodes of contig have coverage less than min_cov, try to shrink the contig.
    """
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
                break

            if not name or not seg_nos or not name_r or not seg_nos_r: break
            split_name = name.split('_')
            cno = str(split_name[1])
            clen = int(split_name[3])
            ccov = float(split_name[5])

            # contig from both orientation
            contigs = [n[:-1] if n[-1] in ['-', '+'] else n[:-2] for n in seg_nos.split(',')]
            contigs_rev = [n[:-1] if n[-1] in ['-', '+'] else n[:-2] for n in seg_nos_r.split(',')]

            if len(contigs) < min_node or clen < min_len / 10 or ccov < min_cov:
                continue

            contig_len = len(contigs)
            if contig_len > 1:
                i = 0
                c = []
                pick = False
                while not pick:
                    e1 = (contigs[i],contigs[i+1])
                    e1_r = (contigs_rev[i], contigs_rev[i+1])
                    i = i + 1
                    if e1 not in simp_edge_dict and e1_r not in simp_edge_dict:
                        print("edge is not exist in both direction, continue to select, contig: ", cno)
                    elif e1 not in simp_edge_dict:
                        c = contigs_rev[:]
                        pick = True
                        print("pick forward side for contig: ", cno)
                    elif e1_r not in simp_edge_dict:
                        c = contigs[:]
                        pick = True
                        print("pick reverse side for contig: ", cno)
                    else:
                        print("both direction edge, error edge case, contig: ", cno)
                    
                    if not pick and i == contig_len - 1:
                        # still not pick until last edge
                        print("all the edge is removed, no pick until last point")
                        break
                if not pick:
                    # whole contig is reduced already, no split chance
                    continue
                else:
                    contig_list = contig_split(graph, cno, c, simp_node_dict, simp_edge_dict, overlap)
                    for (sub_cno, sub_contig, sub_clen, sub_ccov) in contig_list:
                        contig_dict[sub_cno] = (sub_contig, sub_clen, sub_ccov)
            else:
                c = contigs
                if c[0] not in simp_node_dict:
                    print("only left up node is removed already: ", cno)
                else:
                    contig_dict[cno] = (c,clen,ccov)
                    print("only left up node is picked for contig: ", cno)
        contigs_file.close()

    return contig_dict

def contig_split(graph: Graph, cno, contig: list, simp_node_dict: dict, simp_edge_dict: dict, overlap, min_node=2):
    """
    Split the contig by removing all the removed node, and form segments
    list of contig tuple: (cno, contig, ccov, len)
    """
    contig_list = []
    s = 0
    idx = 0
    while s < len(contig):
        x = s
        for i in range(s, len(contig)):
            v = contig[i]
            if v in simp_node_dict:
                x = x + 1
            else:
                break
        # x will end up to removed node idx for the contig
        sub = contig[s:x]
        if len(sub) >= min_node:
            cflow = contig_flow(graph, simp_edge_dict, sub)
            ccov = numpy.mean(cflow) if len(cflow) != 0 else 0
            clen = path_len(graph, [simp_node_dict[node] for node in sub], overlap)
            print(clen)
            contig_list.append((cno+"^"+str(idx), sub, clen, ccov))
            idx = idx + 1
        s = x
        for i in range(x, len(contig)):
            v = contig[i]
            if v not in simp_node_dict:
                s = s + 1
            else:
                break
        # s will end up to not-removed node
    return contig_list

def path_len(graph: Graph, path, overlap):
    """
    Find length of the linear path.
    """
    lens = [len(graph.vp.seq[u]) for u in path]
    return sum(lens) - overlap * (len(lens) - 1) if len(lens) > 0 else 0

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

def contig_to_seq(graph: Graph, contig: list, contig_name, simp_node_dict: dict, overlap_len):
    seq = ""
    for i in range(len(contig)):
        c = simp_node_dict[contig[i]]
        if i == len(contig) - 1:
            seq = seq + graph.vp.seq[c]
        else:
            seq = seq + (graph.vp.seq[c])[:-overlap_len]
    # print(contig_name, " ", seq)
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

def print_edge(graph, e, s=""):
    print(s, " edge: ", graph.vp.id[e.source()],graph.vp.ori[e.source()], "->", graph.vp.id[e.target()], graph.vp.ori[e.target()], graph.ep.flow[e], graph.ep.color[e])

def print_vertex(graph, v, s=""):
    print(s, " vertex: ", graph.vp.id[v], ", dp: ", graph.vp.dp[v], ", ori: ", graph.vp.ori[v], ", in_degree: ", v.in_degree(), ", out_degree: ", v.out_degree(), graph.vp.color[v])

def print_contig(cno, clen, ccov, contig, s=""):
    print(s, " Contig: ", cno, ", length: ", clen, ", cov: ", ccov, "Path: ", [int(v) for v in contig])