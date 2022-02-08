#!/usr/bin/env python3

from graph_tool.all import Graph
import gfapy
import subprocess


def gfa_to_graph(gfa_file):
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
    return graph, node_dict, edge_dict, dp_dict

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
            if graph.vp.color[v] == 'black' and graph.vp.dp[v] > min_cov:
                
                id_record.append(int(graph.vp.id[v]))
                dp_record.append(int(graph.vp.dp[v]))
                
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
            if graph.vp.dp[node_u] <= min_cov or graph.vp.dp[node_v] <= min_cov or graph.ep.flow[e] <= min_cov:
                continue

            gfa.write("L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format
            (u, to_ori(graph.vp.ori[node_u]), 
            v, to_ori(graph.vp.ori[node_v]), 
            graph.ep.overlap[e]))
        gfa.close()
    print("--------------------", filename, "is stored--------------------")
    print(id_record)
    print(dp_record)
    return 0

# FIXME fix the path
def map_ref_to_graph(ref, graph: Graph, simp_node_dict: dict):
    """
    map reference strain to the graph, debug only
    assumption: graph is stored in acc/simplifed_graph, 
    """
    if not ref:
        return -1
    with open("acc/simplifed_graph", 'r') as gfa:
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
    ref, "acc/gfa_to_fq.fq", "acc/overlap.paf"), shell=True)

    strain_dict = {}
    with open("acc/overlap.paf", 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            seg_no = splited[0]
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = splited[5]
            if ((seg_f - seg_s) / seg_l) >= 0.8:
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = []
                strain_dict[ref_no].append(seg_no)
    subprocess.check_call("rm {0}".format("acc/gfa_to_fq.fq"), shell=True)
    
    # print("strain dict mapping")
    # for seg_no, strains in strain_dict.items():
    #     if "seq_1" in strains:
    #         print("node no: ", seg_no)
    #         print("strains: ", strains)
    #         print("-------------------")
    return strain_dict


def get_contig(graph: Graph, contig_file, simp_node_dict: dict, edge_dict: dict, min_cov, min_node=5):
    """
    Map SPAdes's contig to the graph, return all the contigs, with dict
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

            if len(contigs) <= min_node:
                continue

            # select contig from single orientation
            e1 = (contigs[0],contigs[1])
            e1_r = (contigs_rev[0], contigs_rev[1])
            c = []
            if e1 not in edge_dict and e1_r not in edge_dict:
                print("edge is not exist in both direction, error handling, TODO")
            elif e1 not in edge_dict:
                c = contigs_rev
            elif e1_r not in edge_dict:
                c = contigs
            else:
                print("both direction edge, error edge case, TODO")

            edge_flow = []
            for i in range(len(c)-1):
                e = edge_dict[(c[i],c[i+1])]
                f = graph.ep.flow[e]
                edge_flow.append(f)

            contig_dict[(cno, clen, ccov)] = (c,edge_flow)
        contigs_file.close()

    return contig_dict

def contig_to_seq(graph: Graph, contig: list, contig_name, simp_node_dict: dict, overlap_len):
    seq = ""
    for i in range(len(contig)):
        c = simp_node_dict[contig[i]]
        if i == len(contig) - 1:
            seq = seq + graph.vp.seq[c]
        else:
            seq = seq + (graph.vp.seq[c])[:-overlap_len]
    print(contig_name, " ", seq)
    return seq

def reverse_seq(seq: str):
    return ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[x] for x in seq[::-1])

def simplify_edge_dict(edge_dict: dict):
    simp_edge_dict = {}
    for (u, _, v, _), e in edge_dict.items():
        simp_edge_dict[(u,v)] = e
    return simp_edge_dict
