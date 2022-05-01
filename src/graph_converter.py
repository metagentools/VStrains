#!/usr/bin/env python3

# from concurrent.futures import thread
from graph_tool.all import Graph
from graph_tool.draw import graph_draw
import gfapy
import subprocess
import sys

import numpy

from hap_construction import DEBUG_MODE


def gfa_to_graph(gfa_file, init_ori=1):
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
    red_graph, red_node_dict, red_edge_dict = reduce_graph(graph, simp_node_dict, simp_edge_dict)
    return red_graph, red_node_dict, red_edge_dict

def reduce_graph(unsimp_graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    graph = Graph(directed=True)

    vprop_seq = graph.new_vertex_property("string", val="")
    vprop_dp = graph.new_vertex_property("double")
    vprop_kc = graph.new_vertex_property("int32_t")
    vprop_id = graph.new_vertex_property("string", val="UD")
    vprop_color = graph.new_vertex_property("string")

    graph.vp.seq = vprop_seq
    graph.vp.dp = vprop_dp
    graph.vp.kc = vprop_kc
    graph.vp.id = vprop_id
    graph.vp.color = vprop_color

    eprop_overlap = graph.new_edge_property("int", val=0)
    eprop_flow = graph.new_edge_property("float", val=0.0)
    eprop_color = graph.new_edge_property("string")

    graph.ep.overlap = eprop_overlap
    graph.ep.flow = eprop_flow
    graph.ep.color = eprop_color

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
    vprop_seq = graph.new_vertex_property("string", val="")
    vprop_dp = graph.new_vertex_property("double")
    vprop_kc = graph.new_vertex_property("int32_t")
    vprop_id = graph.new_vertex_property("string", val="UD")
    vprop_color = graph.new_vertex_property("string")

    graph.vp.seq = vprop_seq
    graph.vp.dp = vprop_dp
    graph.vp.kc = vprop_kc
    graph.vp.id = vprop_id
    graph.vp.color = vprop_color

    eprop_overlap = graph.new_edge_property("int", val=0)
    eprop_flow = graph.new_edge_property("float", val=0.0)
    eprop_color = graph.new_edge_property("string")

    graph.ep.overlap = eprop_overlap
    graph.ep.flow = eprop_flow
    graph.ep.color = eprop_color

    red_node_dict = {}
    red_edge_dict = {}

    # S
    for line in gfa.segments:
        [line_type, seg_no, seg, dp, kc] = str(line).split("\t")
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
        [line_type, seg_no_l, ori_l, seg_no_r, ori_r, overlap_len] = str(edge).split("\t")
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
    for (u,v),e in simp_edge_dict.items():
        u_node = simp_node_dict[u]
        v_node = simp_node_dict[v]
        flow = 0
        if u_node.out_degree() == 1 and v_node.in_degree() == 1:
            flow = max(graph.vp.dp[u_node], graph.vp.dp[v_node])
        elif u_node.out_degree() > 1 and v_node.in_degree() == 1:
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
            flow = max(graph.vp.dp[v_node], (graph.vp.dp[v_node]/u_out_sum)*graph.vp.dp[u_node])
        elif u_node.out_degree() == 1 and v_node.in_degree() > 1:
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
            flow = max(graph.vp.dp[u_node], (graph.vp.dp[u_node]/v_in_sum)*graph.vp.dp[v_node])
        else:
            u_out_sum = numpy.sum([graph.vp.dp[n] for n in u_node.out_neighbors()])
            v_in_sum = numpy.sum([graph.vp.dp[n] for n in v_node.in_neighbors()])
            flow = max((graph.vp.dp[v_node] / u_out_sum) * graph.vp.dp[u_node], (graph.vp.dp[u_node] / v_in_sum) * graph.vp.dp[v_node])
        graph.ep.flow[e] = round(max(graph.ep.flow[e], flow), 2)
    return

def graph_simplification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, min_cov):
    """
    Directly remove all the vertex with coverage less than minimum coverage and related edge

    Node belongs to any contigs should not be removed
    return:
        removed_node_dict
        removed_edge_dict
    """
    print("-------------------------graph simplification----------------------")
    print("Total nodes: ", len(simp_node_dict), " Total edges: ", len(simp_edge_dict))
    removed_node_dict = {}
    removed_edge_dict = {}
    # iterate until no more node be removed from the graph
    for id, node in list(simp_node_dict.items()):
        if graph.vp.dp[node] < min_cov:
            if id in node_to_contig_dict:
                if DEBUG_MODE:
                    print("node: {0} should not be removed although with ccov: {1}".format(id, graph.vp.dp[node]))
                continue
            # print_vertex(graph, node, "Node removed by graph simplification -")

            # delete the node
            simp_node_dict.pop(id)
            graph.vp.color[node] = 'gray'
            removed_node_dict[id] = node
            total_reduce_depth = graph.vp.dp[node]
            # delete related edges
            for out_node in node.out_neighbors():
                out_id = graph.vp.id[out_node]
                if (id, out_id) in edge_to_contig_dict:
                    if DEBUG_MODE:
                        print("edge: {0} should not be removed".format((id, out_id)))
                    continue
                if (id, out_id) in simp_edge_dict:
                    e = simp_edge_dict[(id, out_id)]
                    graph.ep.color[e] = 'gray'
                    simp_edge_dict.pop((id, out_id))
                    removed_edge_dict[(id, out_id)] = e

            for in_node in node.in_neighbors():
                in_id = graph.vp.id[in_node]
                if (in_id, id) in edge_to_contig_dict:
                    if DEBUG_MODE:
                        print("edge: {0} should not be removed".format((in_id, id)))
                    continue
                if (in_id, id) in simp_edge_dict:
                    e = simp_edge_dict[(in_id, id)]
                    graph.ep.color[e] = 'gray'
                    simp_edge_dict.pop((in_id, id))
                    removed_edge_dict[(in_id, id)] = e
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
            # print_edge(graph, e, "adding edge to graph")
            gfa.write("L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format
            (u, "+", 
            v, "+", 
            graph.ep.overlap[e]))
        gfa.close()
    print("--------------------", filename, "is stored--------------------")
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
        if DEBUG_MODE:
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
        fifo_queue = []
        fifo_queue.append([node_dict[seg_no], init_ori]) 

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
    print("-------------------------verify graph----------------------")
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
def map_ref_to_graph(ref_file, simp_node_dict: dict, graph_file, store_mapping=False, output_file="overlap.paf", fasta_file="temp_gfa_to_fasta.fasta"):
    """
    map reference strain to the graph, debug only
    assumption: graph is stored in acc/simplifed_graph, 
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
    
    # minimap2. you may need to replace the exec path to minimap to fit your case
    subprocess.check_call("/Users/luorunpeng/opt/miniconda3/envs/spades-hapConstruction-env/bin/minimap2 {0} {1} -c > {2}".format(
    ref_file, fasta_file, output_file), shell=True)

    strain_dict = {}
    with open(output_file, 'r') as paf:
        for Line in paf:
            splited = Line.split('\t')
            seg_no = str(splited[0])
            seg_no = splited[0]
            seg_l = int(splited[1])
            seg_s = int(splited[2])
            seg_f = int(splited[3])
            ref_no = splited[5]
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if seg_no not in simp_node_dict:
                continue
            if (((seg_f - seg_s) / seg_l) >= 0.8 and (nmatch/nblock) >= 0.8 and mark > 0):
                if ref_no not in strain_dict:
                    strain_dict[ref_no] = []
                strain_dict[ref_no].append(seg_no)
        paf.close()
        
    subprocess.check_call("rm {0}".format(fasta_file), shell=True)
    if not store_mapping:
        subprocess.check_call("rm {0}".format(output_file), shell=True)
    
    print("strain dict mapping")
    for seg_no, strains in strain_dict.items():
        print("strains: ", seg_no, " Path: ", list_to_string(strains))
        print("-------------------")
    return strain_dict

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
                    quality = "B"*len(splited[2])
                    fasta.write("@{0}\n{1}\n+\n{2}\n".format(splited[1],splited[2],quality))
            fasta.close()
        gfa.close()

def minimap_api(ref_file, fasta_file, output_file):
    subprocess.check_call("/Users/luorunpeng/opt/miniconda3/envs/spades-hapConstruction-env/bin/minimap2 {0} {1} -c > {2}".format(
        ref_file, fasta_file, output_file), shell=True)
    return  

def trim_contig_dict(graph: Graph, simp_node_dict: dict, contig_dict: dict, overlap):
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        involved_node_set = set()
        new_contig = []
        for i in range(len(contig)):
            if contig[i] not in involved_node_set:
                new_contig.append(contig[i])
                involved_node_set.add(contig[i])
        contig_dict[cno] = [new_contig, path_len(graph, [simp_node_dict[no] for no in new_contig], overlap), ccov]
    return contig_dict

def contig_dict_to_fasta(graph: Graph, contig_dict: dict, simp_node_dict: dict, overlap_len, output_file):
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

def get_contig(graph: Graph, contig_file, simp_node_dict: dict, simp_edge_dict: dict, min_len, min_cov):
    """
    Map SPAdes's contig to the graph, return all the contigs.
    if nodes of contig have coverage less than min_cov, try to shrink the contig.
    
    Also fix the graph by adding removed node&edge if is supported by contigs
    """
    print("-----------------get contig-----------------------")
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
            contig_len = len(contigs)

            # contig filter
            # use as less in-confident contigs as possible.
            # if clen < min_len / 10 or ccov < min_cov:
            #     continue
            if clen < min_len / 10 and len(contigs) < 2:
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
                    # whole contig is reduced already, no split chance, potential error
                    continue
                else:
                    contig_dict[cno] = [c, clen, ccov]
                    succ_add = True
                    for i in range(len(c)):
                        c_i = c[i]
                        c_i_1 = c[i+1] if (i < len(c) - 1) else None
                        if c_i not in simp_node_dict:
                            print("node {0} not in contig {1}, error".format(c_i, cno))
                            succ_add = False

                        if c_i_1 != None and c_i_1 not in simp_node_dict:
                            print("node {0} not in contig {1}, error".format(c_i_1, cno))
                            succ_add = False
                        
                        if c_i_1 != None and (c_i, c_i_1) not in simp_edge_dict:
                            print("edge {0} not in contig {1}, error".format((c_i, c_i_1), cno))
                            succ_add = False
                    if succ_add:
                        print("contig {0} is successfully added".format(cno))
                    else:
                        print("contig {0} is failed to be added".format(cno))
            else:
                c = contigs
                if c[0] not in simp_node_dict:
                    print("only left up node is removed already: ", cno)
                else:
                    contig_dict[cno] = [c, clen, ccov]
                    print("only left up node is picked for contig: ", cno)
        contigs_file.close()
    
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    print("-----------------get contig end-----------------------")
    return contig_dict, node_to_contig_dict, edge_to_contig_dict

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

def contig_preprocess(graph:Graph, simp_node_dict: dict, simp_edge_dict: dict, overlap, min_cov, contig_dict: dict):
    """
    Remove and split the unsatisfied contigs. legacy
    """
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        contig_dict.pop(cno)
        contig_list = contig_split(graph, cno, contig, simp_node_dict, simp_edge_dict, overlap, min_cov) #FIXME
        if contig_list == []:
            print("No sub contig be found for original contig: ", cno)
        else:
            if len(contig_list) == 1:
                print("No split for cno: ", cno)
            else:
                print("Update sub contigs to the concat contig dict for cno: ", cno)
            for (sub_cno, sub_contig, sub_clen, sub_ccov) in contig_list:
                if sub_cno in contig_dict:
                    print("sub cno: ", sub_cno, " already exist, error")
                else:
                    contig_dict[sub_cno] = [sub_contig, sub_clen, sub_ccov]
    return contig_dict

def udpate_node_to_contig_dict(graph: Graph, node_to_contig_dict: dict, simp_node_dict: dict):
    for no in node_to_contig_dict.keys():
        node = simp_node_dict[no]
        node_to_contig_dict[no][1] = graph.vp.dp[no]
        node_to_contig_dict[no][2] = node
        
def update_edge_to_contig_dict(graph: Graph, edge_to_contig_dict: dict, simp_edge_dict: dict):
    for (u,v) in edge_to_contig_dict.keys():
        e = simp_edge_dict[(u,v)]
        edge_to_contig_dict[(u,v)][1] = graph.ep.flow[e]
        edge_to_contig_dict[(u,v)][2] = e

def clear_flow(graph: Graph, simp_edge_dict: dict):
    for _, e in simp_edge_dict.items():
        graph.ep.flow[e] = 0.0

def contig_split(graph: Graph, cno, contig: list, simp_node_dict: dict, simp_edge_dict: dict, overlap, min_cov, min_node=2):
    """
    Split the contig by removing all the removed node, and form segments
    list of contig tuple: (cno, contig, ccov, len)
    """
    contig_list = []
    s = 0
    idx = 0
    if DEBUG_MODE:
        print("contig len: ", len(contig))
    while s < len(contig):
        x = s
        keep = False
        for i in range(s, len(contig)):
            u = contig[i]
            if i < len(contig) - 1:
                v = contig[i + 1]
                if u in simp_node_dict:
                    if (u,v) in simp_edge_dict:
                        x = x + 1
                    else:
                        keep = True
                        break
                else:
                    break
            else:
                if u in simp_node_dict:
                    x = x + 1
                else:
                    break
        # x will end up to removed node idx for the contig
        sub = contig[s:x + 1] if keep else contig[s:x]
        if DEBUG_MODE:
            print("sub start: ", contig[s], " sub end: ", contig[x-1])
        if len(sub) >= min_node:
            cflow = contig_flow(graph, simp_edge_dict, sub)
            ccov = numpy.median(cflow) if len(cflow) != 0 else 0
            clen = path_len(graph, [simp_node_dict[node] for node in sub], overlap)
            contig_list.append((cno+"^"+str(idx), sub, clen, ccov))
            idx = idx + 1

        s = x
        for i in range(x, len(contig)):
            u = contig[i]
            if i < len(contig) - 1:
                v = contig[i + 1]
                if u not in simp_node_dict or (u,v) not in simp_edge_dict:
                    s = s + 1
                else:
                    break
            else:
                if u not in simp_node_dict:
                    s = s + 1
                else:
                    break
        # s will end up to not-removed node
    return contig_list
def contig_dict_fix(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap):
    """
    fix the contig dict, reassign the contig coverage to minimum used **edge flow** coverage,
    we not use the minimum used **node** since the min node may still be shared by multiple contigs.
    however, min edge flow may still not be the true contig coverage

    """
    for cno, [contig, _, ccov] in list(contig_dict.items()):
        if all([no in simp_node_dict for no in contig]):
            # contig_dict[cno][2] = round(numpy.min(contig_flow(graph, simp_edge_dict, contig)), 2)
            None
        else:
            subcontigs = []
            curr_contig = []
            addLast = False
            for no in contig:
                if no in simp_node_dict:
                    curr_contig.append(no)
                    addLast = True
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
                # subcov = round(numpy.min(contig_flow(graph, simp_edge_dict, subc)), 2)
                subcov = ccov
                contig_dict[cno + "_" + str(i)] = [subc, sublen, subcov]
    return

def contig_cov_fix(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, printout=False):
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        if printout:
            print("---------------------------------------------------------------")
        if len(contig) > 1:
            contig_dict[cno][2] = numpy.min(contig_flow(graph, simp_edge_dict, contig))
        else:
            if contig[0] not in simp_node_dict:
                contig_dict.pop(cno)
            else:
                contig_dict[cno][2] = graph.vp.dp[simp_node_dict[contig[0]]]
        if cno in contig_dict:
            if printout:
                print_contig(cno, clen, contig_dict[cno][2], contig)
                if len(contig) > 1:
                    print("mean: ", numpy.mean(contig_flow(graph, simp_edge_dict, contig)), 
                    " median: ", numpy.median(contig_flow(graph, simp_edge_dict, contig)), 
                    " min: ", numpy.min(contig_flow(graph, simp_edge_dict, contig)))
    return

def graph_reduction_c(graph: Graph, cand_path, cand_cov, threshold):
    """
    reduce the graph coverage based on given path and cov,
    only applied after udp be deployed in the graph
    """
    for i in range(len(cand_path)):
        graph.vp.udp[cand_path[i]] -= cand_cov

    for i in range(len(cand_path) - 1):
        e = graph.edge(cand_path[i], cand_path[i+1])
        # print(e, graph.vp.id[cand_path[i]], graph.vp.id[cand_path[i+1]])
        graph.ep.flow[e] -= cand_cov

def graph_splitting(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold, strict_mode=True):
    """
    n-n branch splitting
    """
    print("-------------------------graph splitting: {0}----------------------".format(("strict" if strict_mode else "relax")))
    print("Threshold: ", threshold)
    split_branches = []
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    for no, node in list(simp_node_dict.items()):
        ind = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
        outd = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
        if ind == outd and ind > 1 and (no in node_to_contig_dict or not strict_mode):
            support_contigs = node_to_contig_dict[no] if no in node_to_contig_dict else {}
            print_vertex(graph, node, "---------- branch node, support by contig {0}".format(support_contigs))
            cproduct = [(ie, oe, numpy.mean([graph.ep.flow[ie], graph.ep.flow[oe]]), abs(graph.ep.flow[ie] - graph.ep.flow[oe])) for ie in node.in_edges() for oe in node.out_edges() if abs(graph.ep.flow[ie] - graph.ep.flow[oe]) < threshold]
            cproduct = sorted(cproduct, key=lambda element: element[2])

            es = [val for _, _, val, _ in cproduct]
            amb = numpy.max(es) - numpy.min(es) if cproduct else 0
            if cproduct == []:
                print("no matching edges to split")
            elif not strict_mode and amb < threshold and no not in node_to_contig_dict:
                print("ambiguous split under non-strict mode, skip: ", amb)
            else:
                cproduct = sorted(cproduct, key=lambda tuple: tuple[2])

                for i, (ie, oe, eval, delta) in enumerate(cproduct):
                    print("---------------")
                    print_edge(graph, ie, "in")
                    print_edge(graph, oe, "out")
                    if ((graph.vp.id[ie.source()], graph.vp.id[ie.target()]) not in simp_edge_dict 
                    or (graph.vp.id[oe.source()], graph.vp.id[oe.target()]) not in simp_edge_dict):
                        print("edge has been removed already")
                        continue
                    else:
                        print("Delta: ", delta)
                    prev_node = ie.source()
                    prev_no = graph.vp.id[prev_node]
                    next_node = oe.target()
                    next_no = graph.vp.id[next_node]
                    involved_contigs = []
                    cross_talk = False
                    for cno in support_contigs:
                        contig, clen, ccov = contig_dict[cno]
                        if prev_no in contig and next_no not in contig:
                            print("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        elif prev_no not in contig and next_no in contig:
                            print("contig {0}, {1} pass cross edge".format(cno, ccov))
                            cross_talk = True
                            break
                        elif prev_no in contig and next_no in contig:
                            involved_contigs.append(cno)
                            print_contig(cno, clen, ccov, contig, "support contig ")
                        else:
                            None
                    if cross_talk:
                        print("- current branch split forbidden, cross talk", prev_no, no, next_no)
                    elif not involved_contigs and strict_mode:
                        print("- current branch split forbidden, no supporting contig path", prev_no, no, next_no)
                    else:
                        print("- branch split performed")
                        split_branches.append(no)

                        subid = "0" + no + "0" + str(i)
                        subdp = numpy.max([graph.ep.flow[ie], graph.ep.flow[oe]])
                        sub_node = graph_add_vertex(graph, simp_node_dict, subid, subdp, graph.vp.seq[node], graph.vp.kc[node])

                        graph.vp.dp[node] -= subdp

                        graph_remove_edge(graph, simp_edge_dict, prev_no, no)
                        graph_remove_edge(graph, simp_edge_dict, no, next_no)
                        graph_add_edge(graph, simp_edge_dict, prev_node, prev_no, sub_node, subid, graph.ep.overlap[ie],
                        graph.ep.flow[ie])
                        graph_add_edge(graph, simp_edge_dict, sub_node, subid, next_node, next_no, graph.ep.overlap[oe],
                        graph.ep.flow[oe])

                        for icno in involved_contigs:
                            count_occurance = contig_dict[icno][0].count(no)
                            if count_occurance == 0:
                                print("error, node {0} not in contig {1}".format(no, icno))
                            elif count_occurance > 1:
                                print("error, node {0} in contig {1} more than once".format(no, icno))
                            else:
                                contig_dict[icno][0][contig_dict[icno][0].index(no)] = subid
                                node_to_contig_dict[no].remove(icno)

                                if not node_to_contig_dict[no]:
                                    print("no contig is supporting the node: ", no)
                                    node_to_contig_dict.pop(no)
                                
                                if subid not in node_to_contig_dict:
                                    node_to_contig_dict[subid] = {icno}
                                else:
                                    node_to_contig_dict[subid].add(icno)
    
    # remove all the isolated low cov node
    for node in list(graph.vertices()):
        ind = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
        outd = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
        if ind == 0 and outd == 0 and graph.vp.dp[node] < threshold:
            graph_remove_vertex(graph, simp_node_dict, graph.vp.id[node], "remove isolated low cov node")
    
    print("No of branch be removed: ", len(set(split_branches)))
    print("Split branches: ", list_to_string(set(split_branches)))
    return len(set(split_branches))

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
        v = int(p[-1])
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
        name = "0" + str(id) + "0"
        clen = path_len(graph, p, overlap)
        cov = numpy.max([graph.vp.dp[n] for n in p])
        simp_path_dict[name] = [pids, clen, cov]
    return simp_path_dict

def simp_path_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, simp_path_dict: dict, contig_dict: dict, overlap):
    """
    reduce all the contig to a single node, and keep all the potential src/tgt edge.

    1. reduce the coverage for each involving node by the amount of contig cov
    2. reconnect end-to-end nodes to the contig node
    """
    graph_backup = graph.copy()
    simp_node_dict_backup = simp_node_dict.copy()

    node_to_simp_node = {}
    for id in simp_node_dict.keys():
        node_to_simp_node[id] = id

    contig_info = []
    # reduce all the simple path to a single node from the graph
    for cno, (contig, clen, ccov) in list(simp_path_dict.items()):
        src = contig[0]
        tgt = contig[-1]
        id = src + "00" + cno + "00" + tgt
        cseq = path_to_seq(graph_backup, [simp_node_dict_backup[n] for n in contig], cno, overlap)
        kc = numpy.median([graph_backup.vp.kc[simp_node_dict_backup[u]] for u in contig])
        in_edges = list((graph_backup.vp.id[e.source()], src) for e in simp_node_dict_backup[src].in_edges())
        out_edges = list((tgt, graph_backup.vp.id[e.target()],) for e in simp_node_dict_backup[tgt].out_edges())
        

        for i in range(len(contig)):
            no = contig[i]
            node_to_simp_node[no] = id
            popnode = simp_node_dict.pop(no)
            graph.vp.color[popnode] = 'gray'
            if i != len(contig) - 1:
                e = simp_edge_dict.pop((contig[i], contig[i+1]))
                graph.ep.color[e] = 'gray'

        cv = graph.add_vertex()
        graph.vp.seq[cv] = cseq
        graph.vp.dp[cv] = ccov
        graph.vp.kc[cv] = kc
        graph.vp.id[cv] = id
        graph.vp.color[cv] = 'black'
        simp_node_dict[id] = cv
        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])
    
    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u,v) in in_edges:
            # print("Previous concat: ", (u,v))
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                ue = graph.add_edge(simp_node_dict[u], node)
                graph.ep.overlap[ue] = overlap
                graph.ep.color[ue] = 'black'
                simp_edge_dict[(u, graph.vp.id[node])] = ue
                # print_edge(graph, ue, "reappend edge")
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict:
                    ue = graph.add_edge(in_node, node)
                    graph.ep.overlap[ue] = overlap
                    graph.ep.color[ue] = 'black'
                    simp_edge_dict[(graph.vp.id[in_node], graph.vp.id[node])] = ue
                    # print_edge(graph, ue, "reappend edge")             

        for (u,v) in out_edges:
            # print("Previous concat: ", (u,v))
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                ve = graph.add_edge(node, simp_node_dict[v])
                graph.ep.overlap[ve] = overlap
                graph.ep.color[ve] = 'black'
                simp_edge_dict[(graph.vp.id[node], v)] = ve
                # print_edge(graph, ve, "reappend edge")
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict:
                    ve = graph.add_edge(node, out_node)
                    graph.ep.overlap[ve] = overlap
                    graph.ep.color[ve] = 'black'
                    simp_edge_dict[(graph.vp.id[node], graph.vp.id[out_node])] = ve
                    # print_edge(graph, ve, "reappend edge")
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

def check_contig_intersection(cno, contig, cno2, contig2):
    """
    check if contig1 and 2 are overlapped end-to-end or intersect as parallel
    return true if two contigs are parallel, false if overlap end-to-end
    """
    # intersection region is proper subset for both contig1 and contig2
    # check intersection
    intersect = set(contig).intersection(set(contig2))
    print("intersect check: {0} vs {1}, count: {2}".format(cno, cno2, len(intersect)))
    if len(intersect) <= 0:
        print("No intersection")
        return False, intersect, 0

    if len(intersect) == len(contig):
        print("{0} is covered by {1}".format(cno, cno2))
        return True, intersect, -1

    if len(intersect) == len(contig2):
        print("{0} is covered by {1}".format(cno2, cno))
        return True, intersect, -1

    intermediate_nodes_index = [False for _ in contig]
    for i in [contig.index(e) for e in intersect]:
        intermediate_nodes_index[i] = True
    print("cno: {0} intersect at: {1}".format(cno, list(enumerate(intermediate_nodes_index))))
    if not intermediate_nodes_index[0] and not intermediate_nodes_index[-1]:
        print("intersection in the middle")
        return True, intersect, -1
    prev_false_index = intermediate_nodes_index.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index)):
        if not intermediate_nodes_index[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                print("intersection in the middle")
                return True, intersect, -1

    intermediate_nodes_index2 = [False for _ in contig2]
    for i in [contig2.index(e) for e in intersect]:
        intermediate_nodes_index2[i] = True
    print("cno2: {0} intersect at: {1}".format(cno2, list(enumerate(intermediate_nodes_index2))))
    if not intermediate_nodes_index2[0] and not intermediate_nodes_index2[-1]:
        print("intersection in the middle")
        return True, intersect, -1
    prev_false_index = intermediate_nodes_index2.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index2)):
        if not intermediate_nodes_index2[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                print("intersection in the middle")
                return True, intersect, -1
    cend = 0
    cend = cend + 1 if intermediate_nodes_index[0] else cend
    cend = cend + 1 if intermediate_nodes_index[-1] else cend 
    print("overlap end-to-end, tolarent, cend: ", cend)
    return False, intersect, cend
    

def contig_dict_simp(contig_dict: dict, threshold):
    print("---------->Drop duplicated contig")
    dup_contig_ids = set()
    for cno in contig_dict.keys():
        contig, _, ccov = contig_dict[cno]
        if cno not in dup_contig_ids:        
            for cno2 in contig_dict.keys():
                contig2, _, ccov2 = contig_dict[cno2]
                if cno != cno2:
                    if contig == contig2:
                        print("duplicated contig found: {0}, {1}".format(cno, cno2))
                        dup_contig_ids.add(cno2)
                    else:
                        intersect = set(contig).intersection(set(contig2))
                        print("intersect check: {0} vs {1}, count: {2}".format(cno, cno2, len(intersect)))
                        if len(intersect) <= 0:
                            print("No intersection")
                        elif len(intersect) == len(contig):
                            if abs(ccov - ccov2) < threshold:
                                print("{0} is covered by {1}".format(cno, cno2))
                                dup_contig_ids.add(cno)
                        elif len(intersect) == len(contig2):
                            if abs(ccov - ccov2) < threshold:
                                print("{0} is covered by {1}".format(cno2, cno))
                                dup_contig_ids.add(cno2)
    for cno in dup_contig_ids:
        contig_dict.pop(cno)

    return contig_dict

def graph_grouping(graph: Graph, simp_node_dict: dict, forward="", reverse="", partition_length_cut_off=0):
    """
    Maximimize graph connectivity by minimizing node with 0 in-degree or out-degree, detect and remove all the cycles.
    Out-of-date, TBD
    """
    # determine the isolated subgraphs, and assign each node with its group No, which indicates they are belong to same group
    def bfs_grouping(graph: Graph, start_node, group_no, groups):
        """
        Perform a breadth-first search and assign the group no to all the connected nodes.
        """
        groups[group_no] = []

        queue = []
        queue.append(start_node)
        while queue:
            v = queue.pop()
            graph.vp.group[v] = group_no
            groups[group_no].append(v)

            for neighbor in v.all_neighbors():
                if graph.vp.group[neighbor] == -1:
                    queue.append(neighbor)
        return graph, group_no, groups

    graph.vp.group = graph.new_vertex_property("int16_t", val=-1)
    group_no = 1
    groups = {}
    for v in simp_node_dict.values():
        # grouping
        if graph.vp.group[v] == -1:
            graph, group_no, groups = bfs_grouping(graph, v, group_no, groups)
            group_no = group_no + 1

    # connect sub-graphs based on pair-end reads information TODO
    return graph, groups

def path_len(graph: Graph, path, overlap):
    """
    Find length of the linear path.
    """
    lens = [len(graph.vp.seq[u]) for u in path]
    return sum(lens) - overlap * (len(lens) - 1) if len(lens) > 0 else 0

def path_cov(graph: Graph, path):
    """
    Compute the mean cov for the path
    """
    pcovs = [graph.vp.dp[n] for n in path]
    return numpy.mean(pcovs) if len(pcovs) != 0 else 0

def path_usage(graph: Graph, node_usage_dict: dict, path):
    """
    define the usage of the path
    """
    sum = 0
    path_len = len(path)
    for node in path:
        usage = node_usage_dict[graph.vp.id[node]]
        sum += usage[0]/usage[1]
    return round(sum / path_len, 2)

## FIXME also increment via the path depth
def increment_node_usage_dict(node_usage_dict: dict, path_ids, cov):
    """
    update the node usage dict by incrementing all the involving node from the path.
    """
    for id in path_ids:
        node_usage_dict[id][0] += cov

def decrement_node_usage_dict(node_usage_dict: dict, path_ids, cov):
    """
    update the node usage dict by incrementing all the involving node from the path.
    """
    for id in path_ids:
        node_usage_dict[id][0] -= cov

def contig_node_cov_rise(graph: Graph, simp_node_dict: dict, contig_dict: dict, node_to_contig_dict: dict):
    """
    for any node that involved in one or more contigs, rise the depth if less than the sum of related contig cov
    """
    for no, cnos in node_to_contig_dict.items():
        sum_covs = numpy.sum([contig_dict[cno][2] for cno in cnos])
        node = simp_node_dict[no]
        if sum_covs > graph.vp.dp[node]:
            if DEBUG_MODE:
                print("Node: {0} dp is really low: {1} vs {2}, rise it up".format(no, graph.vp.dp[node], sum_covs))
            graph.vp.dp[node] = sum_covs
    return

def contig_edges(graph: Graph, edge_dict: dict, contig):
    """
    edge flow for the contig
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

def min_flow_edge(graph: Graph, edge_dict: dict, contig):
    edge_flow = []
    if len(contig) < 2:
        return None, None
    cand_e = None
    cand_flow = sys.maxsize
    for i in range(len(contig)-1):
        e = edge_dict[(contig[i],contig[i+1])]
        f = graph.ep.flow[e]
        if f < cand_flow:
            cand_e = e
            cand_flow = f
    return (cand_e, cand_flow)

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
    for seg_no, v in simp_node_dict.items():
        print_vertex(graph, v, "stat")
    for (_,_), e in simp_edge_dict.items():
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

def graph_add_vertex(graph: Graph, simp_node_dict: dict, id, dp, seq, kc, s="add vertex", color='black'):
    node = graph.add_vertex()
    graph.vp.id[node] = id
    graph.vp.dp[node] = dp
    graph.vp.seq[node] = seq
    graph.vp.kc[node] = kc
    graph.vp.color[node] = color
    simp_node_dict[id] = node
    print_vertex(graph, node, s)
    return node

def graph_remove_vertex(graph, simp_node_dict, id, s="remove vertex", color='gray', printout=False):
    node = simp_node_dict[id]
    graph.vp.color[node] = color
    simp_node_dict.pop(id)
    if printout:
        print_vertex(graph, node, s)
    return node

def graph_add_edge(graph: Graph, simp_edge_dict: dict, src, src_id, tgt, tgt_id, overlap, flow, s="add edge", color='black'):
    edge = graph.add_edge(src, tgt)
    graph.ep.overlap[edge] = overlap
    graph.ep.color[edge] = color
    graph.ep.flow[edge] = flow
    simp_edge_dict[(src_id, tgt_id)] = edge
    print_edge(graph, edge, s)
    return edge

def graph_recolor_edge(graph: Graph, simp_edge_dict: dict, src, src_id, tgt, tgt_id, s="edge recolor", color='black'):
    """
    recolour edge to black as default
    """
    edge = graph.edge(src, tgt)
    graph.ep.color[edge] = color
    simp_edge_dict[(src_id, tgt_id)] = edge
    print_edge(graph, edge, s)
    return

def graph_remove_edge(graph: Graph, simp_edge_dict: dict, src_id, tgt_id, s="remove edge", color='gray'):
    edge = simp_edge_dict.pop((src_id, tgt_id))
    graph.ep.color[edge] = color
    print_edge(graph, edge, s)
    return ((src_id, tgt_id))            

def cliq_graph_init(graph: Graph):
    """
    remove all the grayed node/edge and return a new graph.
    """
    cliq_graph = Graph(directed=True)
    cliq_graph.vp.cno = cliq_graph.new_vertex_property("string", val="")
    cliq_graph.vp.clen = cliq_graph.new_vertex_property("int32_t")
    cliq_graph.vp.ccov = cliq_graph.new_vertex_property("double")
    cliq_graph.vp.text = cliq_graph.new_vertex_property("string")
    cliq_graph.vp.color = cliq_graph.new_vertex_property("string")
    
    cliq_graph.ep.color = cliq_graph.new_edge_property("string")
    cliq_graph.ep.slen = cliq_graph.new_edge_property("int32_t")
    cliq_graph.ep.text = cliq_graph.new_edge_property("string")

    cliq_node_dict = {}
    cliq_edge_dict = {}

    for contig in graph.vertices():
        if graph.vp.color[contig] == 'black':
            # print("appending node: ", cno, " cov: ", graph.vp.ccov[contig])
            node = cliq_graph.add_vertex()
            cliq_graph.vp.cno[node] = graph.vp.cno[contig]
            cliq_graph.vp.clen[node] = graph.vp.clen[contig]
            cliq_graph.vp.ccov[node] = round(graph.vp.ccov[contig],2)
            cliq_graph.vp.text[node] = cliq_graph.vp.cno[node] + ":" + str(cliq_graph.vp.clen[node]) + ":" + str(cliq_graph.vp.ccov[node])
            cliq_graph.vp.color[node] = 'black'
            cliq_node_dict[cliq_graph.vp.cno[node]] = node
    
    for e in graph.edges():
        i = e.source()
        icno = graph.vp.cno[i]
        j = e.target()
        jcno = graph.vp.cno[j]
        if graph.vp.color[i] != 'black' or graph.vp.color[j] != 'black':
            continue
        if graph.ep.color[e] == 'black':
            edge = cliq_graph.add_edge(cliq_node_dict[icno], cliq_node_dict[jcno])
            cliq_graph.ep.color[edge] = 'black'
            cliq_graph.ep.slen[edge] = graph.ep.slen[e]
            cliq_graph.ep.text[edge] = graph.ep.text[e]
            cliq_edge_dict[(icno, jcno)] = edge

    return cliq_graph, cliq_node_dict, cliq_edge_dict

def cliq_graph_add_node(cliq_graph: Graph, cliq_node_dict: dict, cno, clen, ccov, text, color='black'):
    cnode = cliq_graph.add_vertex()
    cliq_graph.vp.cno[cnode] = cno
    cliq_graph.vp.clen[cnode] = clen
    cliq_graph.vp.ccov[cnode] = ccov
    cliq_graph.vp.text[cnode] = text
    cliq_graph.vp.color[cnode] = color

    cliq_node_dict[cno] = cnode

    return cnode

def cliq_graph_remove_node(cliq_graph: Graph, cliq_node_dict: dict, cno, cnode, color='gray'):
    if cno in cliq_node_dict:
        cliq_graph.vp.color[cnode] = color
        cliq_node_dict.pop(cno)
    return cnode

def cliq_graph_add_edge(cliq_graph: Graph, cliq_edge_dict: dict, cno1, cnode1, cno2, cnode2, slen, text, color='black'):
    edge = cliq_graph.add_edge(cnode1, cnode2)
    cliq_graph.ep.color[edge] = color
    cliq_graph.ep.slen[edge] = slen
    cliq_graph.ep.text[edge] = text
    cliq_edge_dict[(cno1, cno2)] = edge
    print("edge {0} {1} has been added".format(cno1, cno2))
    return edge

def cliq_graph_remove_edge(cliq_graph: Graph, cliq_edge_dict: dict, cno1, cno2, edge, color='gray'):
    cliq_graph.ep.color[edge] = color
    if (cno1, cno2) in cliq_edge_dict:
        cliq_edge_dict.pop((cno1, cno2))
    print("edge {0} {1} has been removed".format(cno1, cno2))
    return edge

def draw_cliq_graph(cliq_graph: Graph, nnodes, nedges, tempdir, output_file):
    output_size = 120 * (nnodes + nedges)
    vsize= 30
    esize = 30
    graph_draw(g=cliq_graph, output="{0}{1}".format(tempdir, output_file), bg_color="white", 
    vertex_text=cliq_graph.vp.text, vertex_size=vsize, vertex_font_size=int(vsize * 0.8), 
    edge_text=cliq_graph.ep.text, edge_font_size= int(esize * 0.8), output_size=(output_size, output_size))
    print("cliq graph has been stored in: {0}{1}".format(tempdir, output_file))

def similarity(cov1, cov2):
    return 1 - (abs(cov1-cov2)/(cov1+cov2))
def similarity_e(e, cliq_graph):
    return similarity(cliq_graph.vp.ccov[e.source()], cliq_graph.vp.ccov[e.target()])

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
    path = ""
    for v in contig:
        path += v + ", "
    path = path[:-2]
    print(s, " Contig: ", cno, ", length: ", clen, ", cov: ", ccov, "Path: ", path)


def list_to_string(ids: list, s=""):
    string = s + " - "
    for id in ids:
        string += str(id) + ", "
    return string[:-2] if len(string) >= 2 else ""

def path_to_id_string(graph: Graph, path, s=""):
    return list_to_string([graph.vp.id[node] for node in path], s)

def get_row(matrix: list, rowId):
    return matrix[rowId]

def get_col(matrix: list, colId):
    return [row[colId] for row in matrix]



#########
# def contig_variation_span(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, maxlen, tempdir):
#     seaborn.set_theme(style="darkgrid")

#     for cno, [contig, clen, ccov] in contig_dict.items():
#         print("-------------------------------------------------------------")
#         src = simp_node_dict[contig[0]]
#         tgt = simp_node_dict[contig[-1]]

#         involved_node, involved_edge = contig_variation_path(graph, simp_node_dict, contig_dict, cno, src, tgt, overlap)
#         # the histogram of the data
#         if len(involved_edge) == 0:
#             print("no edge within the contig variation, skip graph construction")
#             continue
#         xpair = [(key, graph.ep.flow[simp_edge_dict[key]]) for key in involved_edge]
#         print(xpair)
        # x_self = contig_flow(graph, simp_edge_dict, contig)

        # print("contig mean: ", numpy.mean(x_self), 
        # " median: ", numpy.median(x_self), 
        # " min: ", numpy.min(x_self),
        # " max: ", numpy.max(x_self))
        # x_all = [graph.ep.flow[simp_edge_dict[key]] for key in involved_edge]
        # print("span mean: ", numpy.mean(x_all), 
        # " median: ", numpy.median(x_all), 
        # " min: ", numpy.min(x_all),
        # " max: ", numpy.max(x_all))

        # x_span = [graph.ep.flow[simp_edge_dict[key]] for key in involved_edge.difference(set(contig_edges(graph, simp_edge_dict, contig)))]        
        # plt.figure(figsize=(32,16))
        
        # if len(x_span) != 0:
        #     x_series = pandas.Series(x_span).value_counts()
        #     x_self_series = pandas.Series(x_self).value_counts()
        #     df=pandas.concat([x_self_series,x_series],axis=1)
        #     df.columns = ["Contig", "Span"]
        #     df=df.reset_index().melt(id_vars=['index'])
        #     df.columns = ["coverage", "type", "count"]
        #     ax = seaborn.barplot(
    #             x='coverage',
    #             y='count',
    #             hue='type',
    #             data=df,
    #             palette=['blue','red'],
    #             alpha=.5,
    #             dodge=False,)
    #     else:
    #         d = pandas.Series({"coverage": x_self})
    #         ax = seaborn.countplot(x="coverage", data=d)
    #     for container in ax.containers:
    #         ax.bar_label(container)
    #     ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    #     plt.title('Histogram of Contig {0}, CLEN {1}, CCOV {2} variation'.format(cno, clen, ccov))
    #     plt.savefig("{0}contig_{1}_variation_hist.png".format(tempdir, cno))
    # return