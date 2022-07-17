#!/usr/bin/env python3

from graph_tool.topology import shortest_path
from graph_tool.all import Graph

import numpy

from utils.ns_Utilities import *

def dict_to_hist(graph: Graph, contig_dict: dict, b0, b1):
    print("contig histogram generating..")
    contig_hist = {}
    for (cno, [_, clen, ccov]) in contig_dict.items():
        delta = 3*abs(b0 + b1*ccov)
        lb = ccov - delta
        ub = ccov + delta
        min_eset = {}
        min_vset = set()
        for e in graph.edges():
            if graph.ep.flow[e] >= lb and graph.ep.flow[e] <= ub:
                    min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
                    min_vset.add(e.source())
                    min_vset.add(e.target())
        x = [graph.ep.flow[e] for e in min_eset.values()]
        regions, bins = numpy.histogram(x)
        contig_hist[cno] = [clen, ccov, regions, bins]
    print("done")
    return contig_hist

def set_edge_weight(graph: Graph, ccov, b0, b1, relax=False):
    """
    set edge weight for path finding
    with relax is set, remove possible negative cycle by inversing the negative weight assignment
    """
    for e in graph.edges():
        diff = graph.ep.flow[e] - ccov
        delta = 3*abs(b0 + b1 * graph.ep.flow[e])
        if diff < -delta:
            #P4 worst
            graph.ep.eval[e] = (-diff)/delta
            if relax:
                graph.ep.eval[e] += 1
        elif diff >= -delta and diff <= delta:
            #P1 best
            alen = len(str(graph.vp.id[e.source()]).split('_')) + len(str(graph.vp.id[e.target()]).split('_'))
            if relax:
                graph.ep.eval[e] = 0
            else:
                # negative weight, guarantee selection
                graph.ep.eval[e] = -(alen/(abs(diff) + 1))*delta
        elif diff > delta and diff <= 2*delta:
            #P3
            graph.ep.eval[e] = ((diff - delta) / delta) 
            if relax:
                graph.ep.eval[e] += 1
        else:
            #P2
            graph.ep.eval[e] = 1 if relax else 0
    return

def st_path(graph: Graph, ccov, s, t, b0, b1, is_dag=False, ):
    sp_vlist = []
    try:
        sp_vlist, _ = shortest_path(graph, s, t, graph.ep.eval, negative_weights=True, dag=is_dag)
    except ValueError as ve:
        print(ve)
        # remove negative cycle
        set_edge_weight(graph, ccov, b0, b1, relax=True)
        sp_vlist, _ = shortest_path(graph, s, t, graph.ep.eval, negative_weights=False, dag=is_dag)
    return sp_vlist

def eval_score(flow, ccov, threshold):
    diff = flow - ccov
    if diff < -threshold:
        return "P4"
    elif diff >= -threshold and diff <= threshold:
        return "P1"
    elif diff > threshold and diff <= 2*threshold:
        return "P3"
    elif diff > 2*threshold:
        return "P2"
    return None

def extract_cand_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, b0, b1, threshold, min_clen=600):
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict)
    strain_dict = {}
    contig_hist = dict_to_hist(graph, contig_dict, b0, b1)
    is_dag = graph_is_DAG(graph, simp_node_dict)
    graph.vp.keep = graph.new_vertex_property("boolean")
    graph.ep.eval = graph.new_edge_property("double")
    graph.ep.keep = graph.new_edge_property("boolean")
    usage_dict = {}
    for node in simp_node_dict.keys():
        usage_dict[node] = 0

    # primary combination, self cycle and s-t path contig first
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        if ((('global_src', contig[0]) in simp_edge_dict and (contig[-1], 'global_sink') in simp_edge_dict) 
            or (contig[-1], contig[0]) in simp_edge_dict):
            # s-t path contig
            contig_dict.pop(cno)
            print("st path/self cycle contig, retrieve first: ", cno)
            pcov = path_cov(graph, simp_node_dict, simp_edge_dict, contig)
            graph_reduction_c(graph, [simp_node_dict[n] for n in contig], usage_dict, pcov)
            #filter low coverage strain
            if pcov >= threshold:
                print("cand strain found")
                strain_dict['A' + cno] = [contig, clen, pcov]
            else:
                print("low cov strain, removed")

    while len(contig_dict.keys()) > 0:
        print("----------------------------------------------------------------------------")
        cno = max(contig_dict.keys(), key=lambda k: sum(contig_hist[k][2]))
        contig, clen, ccov = contig_dict.pop(cno)
        delta = 3*abs(b0 + b1 * ccov)
        print("REGION: ", contig_hist[cno][2])
        print("BIN: ", contig_hist[cno][3])
        print("Contig: ", cno, "len: ", clen, "->bound: [", ccov - delta, ccov, ccov + delta, "], delta: ", delta, "***" ,list_to_string(contig))

        if ccov < threshold or (all([usage_dict[k] != 0 for k in contig])):
            print("current contig {0} is used previously {1}".format(ccov, threshold))
            continue
        if clen < min_clen:
            print("current contig {0} is too short: {1}bp (unreliable) for s-t extension".format(cno, clen))
            strain_dict['A' + cno] = [contig, clen, ccov]
            continue

        set_edge_weight(graph, ccov, b0, b1)
        strain = []
        # find self cycle first if exist
        if reachable(graph, simp_node_dict, simp_node_dict[contig[-1]], simp_node_dict[contig[0]]):
            print("concat self cycle")
            s = simp_node_dict[contig[-1]]
            t = simp_node_dict[contig[0]]
            ns = None
            prev_edges = []
            if len(contig) == 1:
                # split the single node contig, zip it again later
                ns = graph.add_vertex()
                nt = simp_node_dict[contig[0]]
                for e in list(nt.in_edges()):
                    prev_edges.append((e.source(), graph.ep.overlap[e], graph.ep.color[e], graph.ep.flow[e], graph.ep.eval[e]))
                    ie = graph.add_edge(e.source(), ns)
                    graph.ep.eval[ie] = graph.ep.eval[e]
                    graph.remove_edge(e)
                s = nt
                t = ns
            for e in graph.edges():
                graph.ep.keep[e] = True
            for (u,v) in contig_edges(contig):
                e = simp_edge_dict[(u,v)]
                graph.ep.keep[e] = False
            for v in graph.vertices():
                graph.vp.keep[v] = True
            for v in [simp_node_dict[n] for n in contig[1:-1]]:
                graph.vp.keep[v] = False
            
            graph.set_filters(graph.ep.keep, graph.vp.keep)
            sp_vlist = st_path(graph, ccov, s, t, b0, b1, is_dag)
            strain.extend([simp_node_dict[n] for n in contig])
            print("found path: ", list_to_string([graph.vp.id[n] for n in sp_vlist]))
            strain.extend(sp_vlist[1:-1])

            if ns != None:
                # re-zip split node
                for e in list(ns.in_edges()):
                    graph.remove_edge(e)
                graph.remove_vertex(ns)
                for (src, o, c, f, w) in prev_edges:
                    ie = graph.add_edge(src, s)
                    graph.ep.overlap[ie] = o
                    graph.ep.color[ie] = c
                    graph.ep.flow[ie] = f
                    graph.ep.eval[ie] = w
                    simp_edge_dict[(graph.vp.id[src], graph.vp.id[s])] = ie
            graph.clear_filters()
        else:
            print("concat st path")
            sphead_vlist = st_path(graph, ccov, global_src, simp_node_dict[contig[0]], b0, b1, is_dag)
            sptail_vlist = st_path(graph, ccov, simp_node_dict[contig[-1]], global_sink, b0, b1, is_dag)
            strain.extend(sphead_vlist[1:-1])
            strain.extend([simp_node_dict[n] for n in contig])
            strain.extend(sptail_vlist[1:-1])
        score = []
        for flow in contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in strain]):
            delta = 3*abs(b0 + b1 * flow)
            s = eval_score(flow, ccov, delta)
            score.append(s)
        print("related edges score: ", score)
        ccov = path_cov(graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in strain])
        plen = path_len(graph, strain)
        print(path_to_id_string(graph, strain, "strain, len: {0}, ccov: {1}".format(plen, ccov)))
        #filter low coverage strain
        if ccov >= threshold:
            print("cand strain found")
            strain_dict['A' + cno] = [[graph.vp.id[n] for n in strain], plen, ccov]
        else:
            print("low cov strain, removed")
        graph_reduction_c(graph, strain, usage_dict, ccov)
        contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict)
        contig_hist = dict_to_hist(graph, contig_dict, b0, b1)
    
    for v in [global_sink, global_src]:
        graph.remove_vertex(v)

    return strain_dict