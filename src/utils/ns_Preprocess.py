#!/usr/bin/env python3

import sys, os
import subprocess
from graph_tool.all import Graph

import numpy
import matplotlib.pyplot as plt
import pandas as pd

from utils.ns_Utilities import *

from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.metrics import silhouette_score

def tip_removal_s(graph: Graph, simp_node_dict: dict, contig_dict: dict, tempdir, accept_rate = 0.99):
    if not graph_is_DAG(graph, simp_node_dict):
        print("Graph is Cyclic, tip removal start..")
        tip_removed = False
        while not tip_removed:
            tip_removed = tip_removal(graph, simp_node_dict, tempdir, accept_rate)
        print("done")
        contig_dict_fix(graph, simp_node_dict, contig_dict)
    else:
        print("Graph is DAG, tip removal skipped.")

def tip_removal(graph: Graph, simp_node_dict: dict, tempdir, accept_rate):
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

    def cand_collapse_path(graph: Graph, from_node, to_paths, temp_dir):
        """
        use minimap2 -c to evaluation the node-path similarity, sort based on matching score in DESC order
        
        return: the most similar path if there exist a path with score >= accept rate, else return None
        """
        ref_loc = "{0}/ref.fa".format(temp_dir)
        query_loc = "{0}/query.fa".format(temp_dir)
        overlap_loc = "{0}/overlap.paf".format(temp_dir)
        subprocess.check_call('touch {0}; touch {1}'.format(ref_loc, query_loc), shell=True)
        
        id_path_dict = {}
        for id, path in list(enumerate(to_paths)):
            id_path_dict[id] = path

        # retrieve all the path information and save into ref.fa
        with open(ref_loc, 'w') as ref_file:
            for id, path in id_path_dict.items():
                name = ">" + str(id) + "\n"
                seq = path_to_seq(graph, path, id) + "\n"
                ref_file.write(name)
                ref_file.write(seq)
            ref_file.close()

        # save from node info to query.fa
        with open(query_loc, 'w') as query_file:
            name = ">" + graph.vp.id[from_node] + "\n"
            seq = path_to_seq(graph, [from_node], name) + "\n"
            query_file.write(name)
            query_file.write(seq)
            query_file.close()

        # minimap to obtain matching score for all node-path
        id_evalscore = {}
        minimap_api(ref_loc, query_loc, overlap_loc)
        with open(overlap_loc, 'r') as overlap_file:
            for Line in overlap_file:
                splited = (Line[:-1]).split('\t')
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
        # print("Tip Node: ", graph.vp.id[from_node], best_match)
        if len(best_match) == 0:
            return None
        elif best_match[0][1] >= accept_rate:
            return id_path_dict[best_match[0][0]]
        else:
            return None
    is_removed = True
    # get all the source simple path
    src_nodes = []
    tgt_nodes = []
    isolated_node = []
    for node in simp_node_dict.values():
        if node.in_degree() + node.out_degree() == 0:
            # print_vertex(graph, node, "isolated node")
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
        # print("--------------------------src: {0} --------------".format(graph.vp.id[src]))
        src_len = path_len(graph, [src])
        potential_paths = []
        # path retrieve
        for out_branch in src.out_neighbors():
            if graph.vp.id[out_branch] not in simp_node_dict:
                continue
            # print("current out branch: ", graph.vp.id[out_branch])
            for in_tgt in out_branch.in_neighbors():
                if graph.vp.id[in_tgt] == graph.vp.id[src]:
                    # coincidence path
                    continue
                if graph.vp.id[in_tgt] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                # print("current in tgt: ", graph.vp.id[in_tgt])
                potential_paths.extend(paths_to_tgt(graph, simp_node_dict, src, in_tgt, src_len))
        cand_path = cand_collapse_path(graph, src, potential_paths, tempdir)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, src, cand_path)
            is_removed = False
        else:
            print_vertex(graph, src, "Tip cannot be removed, no matching path")

    # target node collapse
    tgt_nodes = sorted(tgt_nodes, key=lambda x: graph.vp.dp[x])
    for tgt in tgt_nodes:
        print("--------------------------tgt: {0} --------------".format(graph.vp.id[tgt]))
        tgt_len = path_len(graph, [tgt])
        potential_paths = []
        # path retrieve
        for in_branch in tgt.in_neighbors():
            if graph.vp.id[in_branch] not in simp_node_dict:
                continue
            # print("current in branch: ", graph.vp.id[in_branch])
            for out_src in in_branch.out_neighbors():
                if graph.vp.id[out_src] == graph.vp.id[tgt]:
                    # coincidence path
                    continue
                if graph.vp.id[out_src] not in simp_node_dict:
                    # collapsed path in previous iteration
                    continue
                # print("current out src: ", graph.vp.id[out_src])
                potential_paths.extend(paths_from_src(graph, simp_node_dict, tgt, out_src, tgt_len))
        cand_path = cand_collapse_path(graph, tgt, potential_paths, tempdir)
        if cand_path != None:
            remove_tip(graph, simp_node_dict, tgt, cand_path)
            is_removed = False
        else:
            print_vertex(graph, tgt, "Tip cannot be removed, no matching path")

    return is_removed

def delta_estimation(graph: Graph, tempdir, cutoff_size=200):
    print("Start delta estimation")
    xs = []
    ys = []
    sample_size = 0
    for node in graph.vertices():
        if (sum([x.out_degree() for x in node.in_neighbors()]) == node.in_degree() 
            and sum([y.in_degree() for y in node.out_neighbors()]) == node.out_degree()):
            if node.in_degree() > 1:
                sample_size += 1
            lv = sum([graph.vp.dp[n] for n in node.in_neighbors()])
            rv = sum([graph.vp.dp[n] for n in node.out_neighbors()])
            m = graph.vp.dp[node]
            xs.extend([lv, rv])
            ys.extend([m, m])
    print("sample size: ", sample_size)
    if sample_size < cutoff_size:
        b0 = 32.23620072586657
        b1 = 0.009936800927088535
        print("use default delta estimation function: Delta = |", b0, "+", b1, "* Coverage |")
        return b0, b1

    plt.figure(figsize=(12,8))
    plt.hist([b-a for b, a in zip(ys,xs)], bins=len(ys))
    plt.title("delta_hist_plot")
    plt.savefig("{0}{1}".format(tempdir, "tmp/delta_hist_plot.png"))

    df = pd.DataFrame({'x': xs, 'y': ys})
    # find n_clusters
    ncs = [i for i in range(2, len(xs)//10)]
    scores = []
    for n_c in ncs:
        clusterer = KMeans(n_clusters=n_c)
        preds = clusterer.fit_predict(df)
        score = silhouette_score(df, preds)
        scores.append(score)

    plt.figure(figsize=(12,8))
    plt.plot(ncs, scores, c="black")
    plt.title("K-Means Silhouette Score",size=20)
    plt.xlabel("number of cluster", size=16)
    plt.ylabel("score", size=16)
    # plt.show()
    plt.savefig("{0}{1}".format(tempdir, "tmp/silhouette_score_delta.png"))

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
    print("Thus the equation is: Delta = |", b0, "+", b1, "* Coverage |")

    x_range = [0, max(clust_x)]                      # get the bounds for x
    y_range = [b0, b0+b1*x_range[1]]    # get the bounds for y
    plt.plot(x_range, y_range, c="red")
    plt.title("Regression",size=20)
    plt.xlabel("coverage", size=16)
    plt.ylabel("delta", size=16)
    # plt.show()
    plt.savefig("{0}{1}".format(tempdir, "tmp/cluster_delta.png"))

    return b0, b1

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

            for e in set(node.all_edges()):
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