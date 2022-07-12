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

def trim_contig_dict(graph: Graph, simp_node_dict: dict, contig_dict: dict):
    for cno, [contig, _, ccov] in list(contig_dict.items()):
        involved_node_set = set()
        new_contig = []
        for id in contig:
            if id not in involved_node_set:
                new_contig.append(id)
                involved_node_set.add(id)
        contig_dict[cno] = [new_contig, path_len(graph, [simp_node_dict[no] for no in new_contig]), ccov]
    return contig_dict

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

def contig_dict_fix(graph: Graph, simp_node_dict: dict, contig_dict: dict):
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
                sublen = path_len(graph, [simp_node_dict[c] for c in subc])
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

def contig_dict_remapping(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, id_mapping: dict, prev_ids: list):
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
                contig_dict[cno] = [paths[0], path_len(graph, [simp_node_dict[no] for no in paths[0]]), subcov]
        else:
            contig_dict.pop(cno)
            # print("multi mapping for the current contig: whole contig is ambiguous mapping", cno)
            for i, path in enumerate(paths):
                dupcno = cno + "^" + str(i)
                if dupcno in contig_dict:
                    print("dup cno: ", dupcno, " already exist, error")
                else:
                    subcov = path_cov(graph, simp_node_dict, simp_edge_dict, path)
                    contig_dict[dupcno] = [path, path_len(graph, [simp_node_dict[no] for no in path]), subcov]
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
            if src != target:
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

def simple_paths_to_dict(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    simple_paths = simp_path(graph, simp_node_dict, simp_edge_dict)
    simp_path_dict = {}
    for id, p in enumerate(simple_paths):
        pids = [graph.vp.id[n] for n in p]
        name = str(id)
        clen = path_len(graph, p)
        cov = numpy.max([graph.vp.dp[n] for n in p])
        simp_path_dict[name] = [pids, clen, cov]
        # print("Simple PATH: ", list_to_string(pids))
    return simp_path_dict

def simp_path_compactification(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict):
    """
    reduce all the contig to a single node, and keep all the potential src/tgt edge.

    1. reduce the coverage for each involving node by the amount of contig cov
    2. reconnect end-to-end nodes to the contig node
    """

    simp_path_dict = simple_paths_to_dict(graph, simp_node_dict, simp_edge_dict)

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
        cseq = path_to_seq(graph_backup, [simp_node_dict_backup[n] for n in contig], cno)
        in_edges = list((graph_backup.vp.id[e.source()], src, graph_backup.ep.overlap[e]) for e in simp_node_dict_backup[src].in_edges())
        out_edges = list((tgt, graph_backup.vp.id[e.target()], graph_backup.ep.overlap[e]) for e in simp_node_dict_backup[tgt].out_edges())
        

        for i in range(len(contig)):
            no = contig[i]
            node_to_simp_node[no] = id
            graph_remove_vertex(graph, simp_node_dict, no, printout=False)
            if i != len(contig) - 1:
                graph_remove_edge(graph, simp_edge_dict, contig[i], contig[i+1], printout=False)
        cv = graph_add_vertex(graph, simp_node_dict, id, ccov, cseq, printout=False)
        contig_info.append([src, tgt, cno, cv, in_edges, out_edges])
    
    # recover all the in-out edges surrounding the contigs
    for [_, _, _, node, in_edges, out_edges] in contig_info:
        for (u, v, o) in in_edges:
            if u in simp_node_dict and (u, graph.vp.id[node]) not in simp_edge_dict:
                graph_add_edge(graph, simp_edge_dict, simp_node_dict[u], u, node, graph.vp.id[node], o, printout=False)
            
            for [_, tgt, _, in_node, _, _] in contig_info:
                if tgt == u and (graph.vp.id[in_node], graph.vp.id[node]) not in simp_edge_dict:
                    graph_add_edge(graph, simp_edge_dict, in_node, graph.vp.id[in_node], node, graph.vp.id[node], o, printout=False)            

        for (u, v, o) in out_edges:
            if v in simp_node_dict and (graph.vp.id[node], v) not in simp_edge_dict:
                graph_add_edge(graph, simp_edge_dict, node, graph.vp.id[node], simp_node_dict[v], v, o, printout=False)
            
            for [src, _, _, out_node, _, _] in contig_info:
                if src == v and (graph.vp.id[node], graph.vp.id[out_node]) not in simp_edge_dict:
                    graph_add_edge(graph, simp_edge_dict, node, graph.vp.id[node], out_node, graph.vp.id[out_node], o, printout=False)
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
        contig_dict[cno] = [new_contig, path_len(graph, [simp_node_dict[no] for no in new_contig]), ccov]
    
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

def concat_overlap_contig(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict):
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
                    if path_len(graph, overlaps[cno][0]) < path_len(graph, [simp_node_dict[k] for k in intersects]):
                        overlaps[cno] = (intersects, cno2)
                elif status == 'd':
                    if cno not in [c for _, c in overlaps[cno2]]:
                        if path_len(graph, overlaps[cno][0]) < path_len(graph, [simp_node_dict[k] for k in intersects]):
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
            contig_dict[cnos] = [concat_contig, path_len(graph, [simp_node_dict[k] for k in concat_contig]), path_cov(graph, simp_node_dict, simp_edge_dict, concat_contig)]

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

def path_len(graph: Graph, path):
    """
    Find length of the linear path.
    """
    lens = sum([len(graph.vp.seq[u]) for u in path])
    for i in range(len(path) - 1):
        u = path[i]
        v = path[i+1]
        lens -= graph.ep.overlap[graph.edge(u, v)]
    return lens

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

def path_ids_to_seq(graph: Graph, path_ids: list, path_name, simp_node_dict: dict):
    seq = ""
    for i in range(len(path_ids)):
        u = simp_node_dict[path_ids[i]]
        if i == len(path_ids) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            overlap_len = graph.ep.overlap[graph.edge(u, simp_node_dict[path_ids[i+1]])]
            if overlap_len == 0:
                seq = seq + graph.vp.seq[u]
            else:
                seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq

def path_to_seq(graph: Graph, path: list, path_name):
    seq = ""
    for i in range(len(path)):
        u = path[i]
        if i == len(path) - 1:
            seq = seq + graph.vp.seq[u]
        else:
            overlap_len = graph.ep.overlap[graph.edge(u, path[i+1])]
            if overlap_len == 0:
                seq = seq + graph.vp.seq[u]
            else:
                seq = seq + (graph.vp.seq[u])[:-overlap_len]
    return seq

def reverse_seq(seq: str):
    return ''.join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[x] for x in seq[::-1])

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

def graph_add_vertex(graph: Graph, simp_node_dict: dict, id, dp, seq, s="add vertex", color='black', printout=False):
    node = graph.add_vertex()
    graph.vp.id[node] = id
    graph.vp.dp[node] = dp
    graph.vp.seq[node] = seq
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
    print(s, " vertex: ", graph.vp.id[v], ", dp: ", graph.vp.dp[v], ", in_degree: ", v.in_degree(), ", out_degree: ", v.out_degree(), graph.vp.color[v])

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

def add_global_source_sink(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, store_dict=False):
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
        graph.ep.overlap[e] = 0
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
        graph.ep.overlap[e] = 0
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

def paths_from_src(graph: Graph, simp_node_dict: dict, self_node, src, maxlen):
    """
    retrieve all the path from src node to any node 
    within maxlen restriction, in straight direction
    """
    def dfs_rev(graph: Graph, u, curr_path: list, maxlen, visited, all_path):
        visited[u] = True
        curr_path.append(u)
        curr_len = path_len(graph, curr_path)
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

def paths_to_tgt(graph: Graph, simp_node_dict: dict, self_node, tgt, maxlen):
    """
    retrieve all the path from any node to tgt node
    within maxlen restriction, in reverse direction
    """
    def dfs_rev(graph: Graph, v, curr_path: list, maxlen, visited, all_path):
        visited[v] = True
        curr_path.insert(0, v)
        curr_len = path_len(graph, curr_path)
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