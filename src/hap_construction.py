#!/usr/bin/env python3

import sys, os
import subprocess
import argparse

from graph_tool.topology import shortest_path
from graph_tool.all import Graph

import numpy
import matplotlib.pyplot as plt
import seaborn
import pandas as pd

from graph_converter import *


# k-means clustering
from sklearn.cluster import KMeans
# linear regression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import silhouette_score

usage = "Construct viral strains under deno vo approach"
author = "Runpeng Luo"

DEBUG_MODE = False
TEMP_DIR = "acc/"

def main():
    parser = argparse.ArgumentParser(prog='hap_construction.py', description=usage)
    parser.add_argument('-gfa', '--gfa_file', dest='gfa_file', type=str, required=True, help='assembly graph under gfa format')
    parser.add_argument('-c', '--contig', dest='contig_file', type=str, required=True, help='contig file from SPAdes, paths format')
    parser.add_argument('-mincov' '--minimum_coverage', dest='min_cov', type=int, help=("minimum coverage for strains"))
    parser.add_argument('-overlap', '--vertex_overlap', dest='overlap', default=127, type=int, help=("adjacent vertex overlap in the graph"))
    parser.add_argument('-ref', "--reference_fa", dest='ref_file', type=str, help='reference strain, fasta format, DEBUG_MODE only')
    parser.add_argument('-o', '--output_dir', dest='output_dir', default='acc/', type=str, help='output directory (default: acc/)')
    # parser.add_argument('-f', '--forward', dest='forward', type=str, required=True, help='Forward reads, fastq format')
    # parser.add_argument('-r', '--reverse', dest='reverse', type=str, required=True, help='Reverse reads, fastq format')
    # parser.add_argument('-l', "--insert_size", dest='insert_size', type=int, required=True, help='Pair-end read distance')
    args = parser.parse_args()
    global TEMP_DIR
    TEMP_DIR = args.output_dir
    
    subprocess.check_call("rm -rf {0} && mkdir {0}".format(TEMP_DIR), shell=True)

    print("----------------------------------INPUT---------------------------------------")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, args.overlap)
    contig_dict = get_contig(args.contig_file, simp_node_dict, simp_edge_dict, 250)
    
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}graph_L0.gfa".format(TEMP_DIR))
    print("----------------------------------TIP REMOVAL---------------------------------------")
    # apply tip removal for cyclic graph.
    if not graph_is_DAG(graph0, simp_node_dict0):
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(graph0, simp_node_dict0, args.overlap)

        contig_dict_fix(graph0, simp_node_dict0, contig_dict, args.overlap)
    else:
        print("Graph is DAG, tip removal skipped.")
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1.gfa".format(TEMP_DIR))

    print("-----------------------------DELTA ESTIMATION-----------------------------")
    b0, b1 = delta_estimation(graph1)
    print("-------------------------------GRAPH SIMPLIFICATION & REBALANCE-----------------------------------")
    #FIXME
    mediandp = numpy.median([graph1.vp.dp[node] for node in simp_node_dict1.values()])
    if args.min_cov:
        THRESHOLD = args.min_cov
    else:
        THRESHOLD = 0.05 * mediandp
    print("MEDIAN NODE DEPTH: ", mediandp, "threshold: ", THRESHOLD)
    graph_simplification(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, THRESHOLD)

    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}st_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}st_graph_L2.gfa".format(TEMP_DIR))

    coverage_rebalance_s(graph2, simp_node_dict2, simp_edge_dict2, True)

    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}sdt_graph_L3.gfa".format(TEMP_DIR))
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)
    
    print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict)
    
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict3, "{0}graph_L0.gfa".format(TEMP_DIR), True, "{0}node_to_ref.paf".format(TEMP_DIR), 
            "{0}temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph3, simp_node_dict3, contig_dict, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    # end stat

    print("-----------------------GRAPH BRANCH SPLIT & COMPACTIFICATION-------------------------------")
    grapha = graph3
    simp_node_dicta = simp_node_dict3
    simp_edge_dicta = simp_edge_dict3
    total_removed_branch_nt = 0
    total_removed_branch_t = 0
    iterCount = 'A'
    num_split = 0
    trivial_split_count = 0
    while True:
        # trivial branch split
        prev_ids = list(simp_node_dicta.keys())
        trivial_split_count, id_mapping = graph_split_trivial(grapha, simp_node_dicta, simp_edge_dicta)

        graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)

        # fix trivial splitted contig
        contig_dict_remapping(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, id_mapping, prev_ids, args.overlap)
        # non-trivial branch split
        num_split = graph_splitting(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, b0, b1, THRESHOLD, args.overlap)
    
        graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))

        simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, args.overlap)

        graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        grapha, simp_node_dicta, simp_edge_dicta = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(grapha, simp_node_dicta, simp_edge_dicta)

        contig_dup_removed(grapha, simp_edge_dicta, contig_dict)
        trim_contig_dict(grapha, simp_node_dicta, contig_dict, args.overlap)
        contig_cov_fix(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            iterCount = chr(ord(iterCount) + 1)
        else:
            coverage_rebalance_s(grapha, simp_node_dicta, simp_edge_dicta, True)
            graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
            assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
            break
    print("Total non-trivial branches removed: ", total_removed_branch_nt, " total trivial branches removed: ", total_removed_branch_t)

    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    contig_dict_to_path(contig_dict, "{0}post_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(graph5, simp_node_dict5, contig_dict, args.overlap, "{0}post_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}post_contigs.fasta".format(TEMP_DIR), "{0}post_contigs_to_strain.paf".format(TEMP_DIR))
    map_ref_to_contig(contig_dict, "{0}post_contigs_to_strain.paf".format(TEMP_DIR))
    # end stat

    print("-----------------------CONTIG PATH EXTENSION-------------------------------")
    strain_dict = extract_cand_path(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, args.overlap, b0, b1, THRESHOLD)
    
    print("-----------------------FINAL CLEAN UP-------------------------------")
    contig_dup_removed(graph5, simp_edge_dict5, strain_dict)  
    trim_contig_dict(graph5, simp_node_dict5, strain_dict, args.overlap)  
    contig_dict_to_fasta(graph5, simp_node_dict5, strain_dict, args.overlap, "{0}final_contigs.fasta".format(TEMP_DIR))
    contig_dict_to_path(strain_dict, "{0}final_contigs.paths".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}final_contigs.fasta".format(TEMP_DIR), "{0}final_contigs_to_strain.paf".format(TEMP_DIR))
    return 0

def delta_estimation(graph: Graph):
    global TEMP_DIR
    print("Start delta estimation")
    xs = []
    ys = []
    cxs = []
    cys = []
    for node in graph.vertices():
        if (sum([x.out_degree() for x in node.in_neighbors()]) == node.in_degree() 
            and sum([y.in_degree() for y in node.out_neighbors()]) == node.out_degree()):
            lv = sum([graph.vp.dp[n] for n in node.in_neighbors()])
            rv = sum([graph.vp.dp[n] for n in node.out_neighbors()])
            m = graph.vp.dp[node]
            xs.extend([lv, rv])
            ys.extend([m, m])
            cxs.append(m)
            cys.append(m)

    plt.figure(figsize=(12,8))
    plt.hist([b-a for b, a in zip(ys,xs)], bins=len(ys))
    plt.title("delta_hist_plot")
    plt.savefig("{0}{1}".format(TEMP_DIR, "delta_hist_plot.png"))

    print("sample size: ", len(xs))
    if len(xs) < 100:
        return 32.23620072586657, 0.009936800927088535

    df = pd.DataFrame({'x': xs, 'y': ys})
    # find n_clusters
    ncs = [i for i in range(2, len(xs)//10)]
    scores = []
    for n_c in ncs:
        clusterer = KMeans(n_clusters=n_c)
        preds = clusterer.fit_predict(df)
        centers = clusterer.cluster_centers_
        score = silhouette_score(df, preds)
        scores.append(score)

    plt.figure(figsize=(12,8))
    plt.plot(ncs, scores, c="black")
    plt.title("K-Means Silhouette Score",size=20)
    plt.xlabel("number of cluster", size=16)
    plt.ylabel("score", size=16)
    # plt.show()
    plt.savefig("{0}{1}".format(TEMP_DIR, "silhouette_score_delta.png"))

    optim_nc = max(ncs, key=lambda x: scores[ncs.index(x)])
    print("Optim n_cluster: ", optim_nc)
    # clustering
    kmc = KMeans(n_clusters=optim_nc)
    kmc_model = kmc.fit(df)
    clust_medians = []
    clust_x = []
    # plot
    colors=["red","blue","green","purple","orange"]
    plt.figure(figsize=(12,8))
    for i in range(numpy.max(kmc_model.labels_)+1):
        cxs = df[kmc_model.labels_==i].iloc[:,0]
        cys = df[kmc_model.labels_==i].iloc[:,1]
        if len(cxs) < 10:
            print("skip cluster: ", kmc_model.cluster_centers_[i])
            continue
        plt.scatter(cxs, cys, label=i, c=colors[i%len(colors)], alpha=0.5)

        vals = []
        for (x1,y1) in zip(cxs, cys):
            for (x2,y2) in zip(cxs, cys):
                if x1 != x2 and y1 != y2:
                    vals.append(max(abs(y1-y2), abs(x1-x2)))
        clust_x.append(kmc_model.cluster_centers_[i][0])
        clust_medians.append(numpy.median(vals))

    plt.scatter(kmc_model.cluster_centers_[:,0], kmc_model.cluster_centers_[:,1], label='Cluster Centers', c="black", s=200)
    plt.title("K-Means Clustering",size=20)
    plt.xlabel(df.columns[0], size=16)
    plt.ylabel(df.columns[1], size=16)
    plt.legend()

    plt.figure(figsize=(12,8))
    plt.scatter(clust_x, clust_medians, label='Cluster Centers', c="black", s=100)
    lr = LinearRegression()
    # fit the data using linear regression
    # the reshaping is because the function expects more columns in x
    model = lr.fit(numpy.array(clust_x).reshape(-1,1),clust_medians)
    b0 = model.intercept_
    b1 = model.coef_[0]
    print("The intercept of this model is:", model.intercept_)
    print("The slope coefficient of this model is:", model.coef_[0])
    print("")
    print("Thus the equation is: Coverage_2 = |", b0, "+", b1, "* Coverage |")

    x_range = [0, max(clust_x)]                      # get the bounds for x
    y_range = [b0, b0+b1*x_range[1]]    # get the bounds for y
    plt.plot(x_range, y_range, c="red")
    plt.title("Regression",size=20)
    plt.xlabel("coverage", size=16)
    plt.ylabel("delta", size=16)
    # plt.show()
    plt.savefig("{0}{1}".format(TEMP_DIR, "cluster_delta.png"))

    return b0, b1

curr_path = []
def node_partition(graph: Graph, simp_node_dict: dict):
    """
    partition the graph into cyclic node and linear intersect cyclic + linear node

    store the cyc node into one graph object, where noncyc node into another object
    """
    def noncyc_path(graph: Graph, src, tgt):
        all_noncyc_paths = []
        visited = {}
        for v in graph.vertices():
            visited[v] = False

        def dfs(u, v):
            global curr_path
            visited[u] = True
            curr_path.append(u)
            # print(path_to_id_string(graph, curr_path, graph.vp.id[u]))
            if u == v:
                print(path_to_id_string(graph, curr_path, "path found"))
                all_noncyc_paths.append(curr_path[:])
            else:
                for next in u.out_neighbors():
                    if not visited[next]:
                        dfs(next, v)
            curr_path = curr_path[:-1]
            visited[u] = False
            return
        dfs(src, tgt)
        return all_noncyc_paths

    # init path acc variable
    global curr_path, TEMP_DIR
    curr_path = []
    # for any source and target node pair in the graph, the tranversing non-cycle/simple path
    # would only include the linear&cycle intersection nodes.

    noncyc_nodes = set()
    simple_paths = []

    srcs = []
    tgts = []
    for node in simp_node_dict.values():
        if node.in_degree() == 0 and node.out_degree() == 0:
            print_vertex(graph, node, "isolated node")
        elif node.in_degree() == 0:
            srcs.append(node)
        elif node.out_degree() == 0:
            tgts.append(node)
        else:
            None
    if srcs == [] or tgts == []:
        print("no src or tgt be found on the graph")
        return None, None

    for src in srcs:
        for tgt in tgts:
            paths = noncyc_path(graph, src, tgt)
            for p in paths:
                [noncyc_nodes.add(graph.vp.id[n]) for n in p]
                simple_paths.append([graph.vp.id[n] for n in p])

    print(list_to_string(list(noncyc_nodes), "non-cyclic+intersection ids"))
    [print(p) for p in simple_paths]

    # partitate into linear graph
    noncyc_graph = graph.copy()
    simp_node_dict_noncyc, simp_edge_dict_noncyc = graph_to_dict(noncyc_graph)
    graph_color_other_to_gray(noncyc_graph, simp_node_dict_noncyc, noncyc_nodes)
    graph_to_gfa(noncyc_graph, simp_node_dict_noncyc, simp_edge_dict_noncyc, "{0}nc_graph_L2p.gfa".format(TEMP_DIR))

    return noncyc_nodes, simple_paths

def tip_removal(graph: Graph, simp_node_dict: dict, overlap):
    """
    retrieve all the source/tail simple path, and merge them into adjacent neighbor path if possible
    
    the collapse step can be done before node depeth rebalance, since it only regards to
    matching score within node seq len

    if is the case, then spades contig may also be modified.
    """
    def remove_tip(graph: Graph, simp_node_dict: dict, from_node, to_path):
        """
        collapse the node with the given path, increment given path depth, remove related information
        about the node.
        """
        graph.vp.color[from_node] = 'gray'
        pending_dp = graph.vp.dp[from_node]
        for node in to_path:
            graph.vp.dp[node] += pending_dp
        simp_node_dict.pop(graph.vp.id[from_node])
        for e in from_node.all_edges():
            graph.ep.color[e] = 'gray'
        print(path_to_id_string(graph, to_path, "Tip Node {0} collapsed to path".format(graph.vp.id[from_node])))
        return

    def cand_collapse_path(graph: Graph, from_node, to_paths, temp_dir, overlap, accept_rate = 0.98):
        """
        use minimap2 -c to evaluation the node-path similarity, sort based on matching score in DESC order
        
        return: the most similar path if there exist a path with score >= accept rate, else return None
        """
        ref_loc = "{0}ref.fa".format(temp_dir)
        query_loc = "{0}query.fa".format(temp_dir)
        overlap_loc = "{0}overlap.paf".format(temp_dir)
        subprocess.check_call('touch {0}; touch {1}'.format(ref_loc, query_loc), shell=True)
        
        id_path_dict = {}
        for id, path in list(enumerate(to_paths)):
            id_path_dict[id] = path

        # retrieve all the path information and save into ref.fa
        with open(ref_loc, 'w') as ref_file:
            for id, path in id_path_dict.items():
                name = ">" + str(id) + "\n"
                seq = path_to_seq(graph, path, id, overlap) + "\n"
                ref_file.write(name)
                ref_file.write(seq)
            ref_file.close()

        # save from node info to query.fa
        with open(query_loc, 'w') as query_file:
            name = ">" + graph.vp.id[from_node] + "\n"
            seq = path_to_seq(graph, [from_node], name, overlap) + "\n"
            query_file.write(name)
            query_file.write(seq)
            query_file.close()

        # minimap to obtain matching score for all node-path
        id_evalscore = {}
        minimap_api(ref_loc, query_loc, overlap_loc)
        with open(overlap_loc, 'r') as overlap_file:
            for Line in overlap_file:
                splited = Line.split('\t')
                path_no = int(splited[5])
                nmatch = int(splited[9])
                nblock = int(splited[10])
                if path_no not in id_evalscore:
                    id_evalscore[path_no] = [nmatch/nblock]
                else:
                    id_evalscore[path_no].append(nmatch/nblock)
            overlap_file.close()
        
        # remove temp file
        subprocess.check_call('rm {0}; rm {1}; rm {2}'.format(ref_loc, query_loc, overlap_loc), shell=True)
        
        id_evalscore_sum = []
        for id, scores in id_evalscore.items():
            mean_score = numpy.mean(scores) if len(scores) != 0 else 0
            id_evalscore_sum.append((id, mean_score))
        
        best_match = sorted(id_evalscore_sum, key=lambda t: t[1], reverse=True)
        print("Tip Node: ", graph.vp.id[from_node], best_match)
        if len(best_match) == 0:
            return None
        elif best_match[0][1] >= accept_rate:
            return id_path_dict[best_match[0][0]]
        else:
            return None
    global TEMP_DIR
    is_removed = True
    # get all the source simple path
    src_nodes = []
    tgt_nodes = []
    isolated_node = []
    for node in simp_node_dict.values():
        if node.in_degree() + node.out_degree() == 0:
            print_vertex(graph, node, "isolated node")
            isolated_node.append(node)
        elif node.in_degree() == 0:
            src_nodes.append(node)
        elif node.out_degree() == 0:
            tgt_nodes.append(node) 
        else:
            None
    
    # src node collapse
    src_nodes = sorted(src_nodes, key=lambda x: graph.vp.dp[x])
    for src in src_nodes:
        print("--------------------------src: {0} --------------".format(graph.vp.id[src]))
        src_len = path_len(graph, [src], overlap)
        potential_paths = []
        # path retrieve
        for out_branch in src.out_neighbors():
            if graph.vp.id[out_branch] not in simp_node_dict:
                continue
            print("current out branch: ", graph.vp.id[out_branch])
            for in_tgt in out_branch.in_neighbors():
                if graph.vp.id[in_tgt] == graph.vp.id[src]:
                    # coincidence path
                    continue
                if graph.vp.id[in_tgt] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                print("current in tgt: ", graph.vp.id[in_tgt])
                potential_paths.extend(paths_to_tgt(graph, simp_node_dict, src, in_tgt, overlap, src_len))
        cand_path = cand_collapse_path(graph, src, potential_paths, TEMP_DIR, overlap)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, src, cand_path)
            is_removed = False
        else:
            print_vertex(graph, src, "Tip cannot be removed, no matching path")

    # target node collapse
    tgt_nodes = sorted(tgt_nodes, key=lambda x: graph.vp.dp[x])
    for tgt in tgt_nodes:
        print("--------------------------tgt: {0} --------------".format(graph.vp.id[tgt]))
        tgt_len = path_len(graph, [tgt], overlap)
        potential_paths = []
        # path retrieve
        for in_branch in tgt.in_neighbors():
            if graph.vp.id[in_branch] not in simp_node_dict:
                continue
            print("current in branch: ", graph.vp.id[in_branch])
            for out_src in in_branch.out_neighbors():
                if graph.vp.id[out_src] == graph.vp.id[tgt]:
                    # coincidence path
                    continue
                if graph.vp.id[out_src] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                print("current out src: ", graph.vp.id[out_src])
                potential_paths.extend(paths_from_src(graph, simp_node_dict, tgt, out_src, overlap, tgt_len))
        cand_path = cand_collapse_path(graph, tgt, potential_paths, TEMP_DIR, overlap)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, tgt, cand_path)
            is_removed = False
        else:
            print_vertex(graph, tgt, "Tip cannot be removed, no matching path")

    return is_removed

def coverage_rebalance_s(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strict=False):
    global TEMP_DIR
    print("------------------------------NODE PARTITION-----------------------------------")
    # store the no-cycle nodes in nc_graph_L3p.gfa
    noncyc_nodes = None
    simple_paths = None
    # graphtool is_DAG() may not work if the graph is not connected as several parts
    if not graph_is_DAG(graph, simp_node_dict):
        noncyc_nodes, simple_paths = node_partition(graph, simp_node_dict)

    print("-------------------------------COVERAGE REBALANCE-----------------------------------")
    # all the previous depth has been stored in the dict
    # ratio: normalised balanced node depth / previous node depth
    _, _, _, ratio = coverage_rebalance(graph, simp_node_dict, simp_edge_dict, strict)
    
    if noncyc_nodes != None and simple_paths != None:
        print("rebalance linear subgraph now..")
        graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}nc_graph_L2p.gfa".format(TEMP_DIR))
        coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc, strict)
        print("Done, start coverage merge")

        for no, node in simp_node_dictnc.items():
            cnode = simp_node_dict[no]
            merge_dp = graphnc.vp.dp[node] + graph.vp.dp[cnode]
            graph.vp.dp[cnode] = merge_dp   
    else:
        print("no linear subgraph available..")

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, strict=False):
    def expectation_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
        for (u,v),e in simp_edge_dict.items():

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

    def maximization_node_depth(graph: Graph, simp_node_dict: dict):
        """
        return True if any updates happened
        """
        is_update = False
        for no, node in simp_node_dict.items():
            us = [e.source() for e in node.in_edges() if graph.ep.color[e] == 'black']
            u_out_degrees = numpy.sum([(len([e for e in u_node.out_edges() if graph.ep.color[e] == 'black'])) for u_node in us]) 
            node_in_degree = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
            in_neighbor_dp_sum = -1
            if u_out_degrees == node_in_degree:
                in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])
            
            vs = [e.target() for e in node.out_edges() if graph.ep.color[e] == 'black']
            v_in_degrees = numpy.sum([(len([e for e in v_node.in_edges() if graph.ep.color[e] == 'black'])) for v_node in vs]) 
            node_out_degree = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
            out_neighbor_dp_sum = -1
            if v_in_degrees == node_out_degree:
                out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])

            curr_dp = graph.vp.dp[node]
            inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
            outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
            graph.vp.dp[node] = numpy.max([curr_dp, inflow, outflow, in_neighbor_dp_sum, out_neighbor_dp_sum])
            # print_vertex(graph, node, "prev dp: {0}".format(curr_dp))
            # print("inflow: ", inflow, "outflow: ", outflow, "in n sum: ", in_neighbor_dp_sum, " out n sum: ", out_neighbor_dp_sum)
            if curr_dp != graph.vp.dp[node]:
                is_update = True
        return is_update
    
    # set cutoff delta
    cutoff = 0
    if graph_is_DAG(graph, simp_node_dict):
        cutoff = 0.00001 * len(simp_node_dict) if strict else 0.0001 * len(simp_node_dict)
    else:
        cutoff = 0.0001 * len(simp_node_dict) if strict else 0.01 * len(simp_node_dict)
    # store previous node depth
    prev_dp_dict = {}
    for no, v in simp_node_dict.items():
        prev_dp_dict[no] = graph.vp.dp[v]
    print("cutoff delta: ", cutoff)

    # EM optimization
    sum_delta = 0
    is_update = True
    while is_update:
        # E Step
        expectation_edge_flow(graph, simp_node_dict, simp_edge_dict)
        # Validation
        sum_delta = 0.0
        sum_dp = numpy.sum([graph.vp.dp[node] for node in graph.vertices()])
        dom_delta = 0.0
        deltas = []
        for no, node in simp_node_dict.items():
            node_in_degree = len([e for e in node.in_edges() if graph.ep.color[e] == 'black'])
            node_out_degree = len([e for e in node.out_edges() if graph.ep.color[e] == 'black'])
            if node_in_degree == 0:
                inflow = graph.vp.dp[node]
            else:
                inflow = numpy.sum([graph.ep.flow[e] for e in node.in_edges() if graph.ep.color[e] == 'black'])
            
            if node_out_degree == 0:
                outflow = graph.vp.dp[node]
            else:
                outflow = numpy.sum([graph.ep.flow[e] for e in node.out_edges() if graph.ep.color[e] == 'black'])

            dominator = (inflow + outflow)/ 2
            if dominator != 0.0:
                sum_delta += abs(inflow - outflow)
                dom_delta += dominator
                deltas.append((no, (graph.vp.dp[node] * abs(inflow - outflow))/dominator))

        if strict:
            sum_delta = sum([k[1] for k in deltas]) / sum_dp
        else:
            sum_delta = sum_delta / dom_delta
        print("sum delta: ", sum_delta, "worst delta: ", sorted(deltas, key=lambda p: p[1], reverse=True)[:10])
        if sum_delta < cutoff:
            break
        # M Step
        is_update = maximization_node_depth(graph, simp_node_dict)

    # final evaluation
    ratios = [(graph.vp.dp[u] / prev_dp_dict[no]) for no, u in simp_node_dict.items()]
    sum_depth_before = numpy.sum([dp for dp in prev_dp_dict.values()])
    sum_ratio = (numpy.sum([graph.vp.dp[u] for u in simp_node_dict.values()]) / sum_depth_before)
    print("Sum Ratio: ", sum_ratio, "Ave Ratio: ", numpy.mean(ratios), "Max Ratio: ", numpy.max(ratios), "Min Ratio: ", numpy.min(ratios), "Delta: ", sum_delta)

    ratio_div = sum_ratio
    print("selected ratio: ", ratio_div)
    for node in simp_node_dict.values():
        graph.vp.dp[node] = graph.vp.dp[node] / ratio_div

    node_ratio_dict = {}
    for no in prev_dp_dict.keys():
        if prev_dp_dict[no] != 0:
            node_ratio_dict[no] = (graph.vp.dp[simp_node_dict[no]] / prev_dp_dict[no]) if prev_dp_dict[no] != 0 else 0
    curr_dp_dict = {}
    for no, v in simp_node_dict.items():
        curr_dp_dict[no] = graph.vp.dp[v]

    return prev_dp_dict, curr_dp_dict, node_ratio_dict, ratio_div

def extract_cand_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, b0, b1, threshold, strain_prefix='A'):
    def dict_to_hist(graph: Graph, contig_dict: dict):
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
                    # if (reachable(graph, simp_node_dict, e.target(), simp_node_dict[contig[0]]) or
                    #     reachable(graph, simp_node_dict, simp_node_dict[contig[-1]], e.source())):
                        min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
                        min_vset.add(e.source())
                        min_vset.add(e.target())
            x = [graph.ep.flow[e] for e in min_eset.values()]
            regions, bins = numpy.histogram(x)
            contig_hist[cno] = [clen, ccov, regions, bins]
        print("done..")
        return contig_hist
    
    def set_edge_weight(graph: Graph, ccov, b0, b1, relax=False):
        for e in graph.edges():
            diff = graph.ep.flow[e] - ccov
            delta = 3*abs(b0 + b1 * graph.ep.flow[e])
            if diff < -delta:
                #P4 worst
                graph.ep.eval[e] = (-diff)/delta
                graph.ep.eval[e] += 1 if relax else 0
            elif diff >= -delta and diff <= delta:
                #P1 best
                alen = len(str(graph.vp.id[e.source()]).split('_')) + len(str(graph.vp.id[e.target()]).split('_'))
                # negative weight, guarantee selection
                graph.ep.eval[e] = -(alen/(abs(diff) + 1))*delta
                # when relax, remove possible negative cycle
                graph.ep.eval[e] = 0 if relax else graph.ep.eval[e]
            elif diff > delta and diff <= 2*delta:
                #P3
                graph.ep.eval[e] = ((diff - delta) / delta)
                graph.ep.eval[e] += 1 if relax else 0
            else:
                #P2
                graph.ep.eval[e] = 0
                graph.ep.eval[e] += 1 if relax else 0
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
    global TEMP_DIR
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, True)
    strain_dict = {}
    contig_hist = dict_to_hist(graph, contig_dict)
    is_dag = graph_is_DAG(graph, simp_node_dict)

    # stat
    x = [cov for [_, _, cov] in contig_dict.values()]
    regions, bins = numpy.histogram(x)
    print(regions)
    print(bins)
    plt.figure(figsize=(64,32))
    plt.hist(x)
    plt.title("contig_cov")
    plt.savefig("{0}{1}".format(TEMP_DIR, "contig_cov.png"))
    # end stat

    graph.ep.eval = graph.new_edge_property("double")
    while len(contig_dict.keys()) > 0:
        print("----------------------------------------------------------------------------")
        cno = max(contig_dict.keys(), key=lambda k: sum(contig_hist[k][2]))
        contig, _, ccov = contig_dict.pop(cno)
        delta = 3*abs(b0 + b1 * ccov)
        print("REGION: ", contig_hist[cno][2])
        print("BIN: ", contig_hist[cno][3])
        print(cno, "->bound: [", ccov - delta, ccov, ccov + delta, "], delta: ", delta, "***" ,list_to_string(contig))

        if ccov < threshold:
            print("current contig {0} is used previously {1}".format(ccov, threshold))
            continue

        set_edge_weight(graph, ccov, b0, b1)
        strain = []
        # find self cycle first if exist
        if reachable(graph, simp_node_dict, simp_node_dict[contig[-1]], simp_node_dict[contig[0]]):
            print("concat self cycle")
            sp_vlist = st_path(graph, ccov, simp_node_dict[contig[-1]], simp_node_dict[contig[0]], b0, b1, is_dag)
            strain.extend([simp_node_dict[n] for n in contig])
            strain.extend(sp_vlist[1:-1])
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
            if s == 'P4':
                print("P4: ", flow)
            score.append(s)
        ccov = path_cov(graph, simp_node_dict, simp_edge_dict, [graph.vp.id[n] for n in strain])
        plen = path_len(graph, strain, overlap)
        print(path_to_id_string(graph, strain, "strain, len: {0}, ccov: {1}".format(plen, ccov)))
        print("related edges score: ", score)
        graph_reduction_c(graph, strain, ccov)
        
        # strain_dict[strain_prefix + cno] = [[graph.vp.id[n] for n in strain], plen, ccov]
        #filter low coverage strain
        if ccov >= threshold:
            print("cand strain found")
            strain_dict[strain_prefix + cno] = [[graph.vp.id[n] for n in strain], plen, ccov]
        else:
            print("low cov strain, removed")
        contig_cov_fix(graph, simp_node_dict, simp_edge_dict, contig_dict)
        contig_hist = dict_to_hist(graph, contig_dict)

    return strain_dict

if __name__ == "__main__":
    main()