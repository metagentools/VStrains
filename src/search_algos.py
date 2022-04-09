#!/usr/bin/env python3

from graph_tool import Graph
import sys

import graph_converter

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
        for neighbour in v.out_neighbors():
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
    for id, node in simp_node_dict.items():
        visited[node] = False
        recStack[node] = False
    for id, node in simp_node_dict.items():
        if not visited[node]:
            if isCyclicUtil(graph, node, visited, recStack):
                print("graph is cyclic")
                return False
    print("graph is not cyclic")
    return True

def reachable(graph: Graph, simp_node_dict: dict, src, src_cno, tgt, tgt_cno, contig_dict: dict):
    """
    determine whether src can possibly reach the tgt
    """
    print("reachable check: {0} - {1}".format(graph.vp.id[src], graph.vp.id[tgt]))
    visited = {}
    for no in simp_node_dict.keys():
        visited[no] = False
    for c in contig_dict[src_cno][0]:
        visited[c] = True
    for c in contig_dict[tgt_cno][0]:
        visited[c] = True

    visited[graph.vp.id[src]] = False
    visited[graph.vp.id[tgt]] = False

    queue = [src]
    reached = False
    while queue:
        curr = queue.pop()
        visited[graph.vp.id[curr]] = True
        if curr == tgt:
            reached = True
            break
        for out in curr.out_neighbors():
            if not visited[graph.vp.id[out]]:
                graph.vp.prev[out] = graph.vp.id[curr]
                queue.append(out)
    
    rec_path = None

    if not reached:
        print("not reachable")
    else:
        rec_path = [graph.vp.id[tgt]]
        node = tgt
        while graph.vp.prev[node] != graph.vp.id[src]:
            prev_id = graph.vp.prev[node]
            rec_path.insert(0, prev_id)
            node = simp_node_dict[prev_id]
        rec_path.insert(0, graph.vp.id[src])
        print("PATH: ", graph_converter.list_to_string(rec_path))

        for cno, (contig, _, _) in contig_dict.items():
            if cno in [src_cno, tgt_cno]:
                # self involved, tolerant
                continue
            if all(c in rec_path for c in contig):
                reached = False
                print("cno: {0} is contained in the path".format(cno))

    # clean prev
    for node in simp_node_dict.values():
        graph.vp.prev[node] = ""

    if reached:
        print("reachable")
    return reached, rec_path


def dijkstra_sp(graph: Graph, simp_node_dict: dict, source, sink, closest_cov, overlap: int):
    """
    Use basic dijkstra algorithm to compute the shortest path between source and sink
    """
    print("Dijkastra path finding: closest cov guide: ", closest_cov)
    dist = {}
    prev = {}
    Q = set()
    for node in simp_node_dict.values():
        dist[node] = sys.maxsize
        prev[node] = None
        Q.add(node)
    dist[source] = 0

    while Q:
        u = None
        # retrieve min
        for n, d in sorted(dist.items(), key=lambda x: x[1]):
            if n in Q:
                u = n
                break

        Q.remove(u)
        if u == sink:
            break

        for v in u.out_neighbors():
            if v in Q:
                # relax
                alt = dist[u] + pow(graph.ep.flow[graph.edge(u, v)] - closest_cov, 2)
                # graph_converter.path_len(graph, [u, v], overlap)
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = u
    
    if dist[sink] == sys.maxsize:
        # not reachable
        return None, None

    node = sink
    sp = []
    while prev[node]:
        sp.insert(0, node)
        node = prev[node]
    if not sp:
        print("SP not found")
        return None, None
    else:
        sp.insert(0, source)
        print(graph_converter.path_to_id_string(graph, sp, "SP be found: "))
        print("plen: ", graph_converter.path_len(graph, sp[1:-1], overlap))

        return sp[1:-1], graph_converter.path_len(graph, sp[1:-1], overlap)

def get_concat_len(head_cno, head_clen, tail_cno, tail_clen, plen, overlap):
    total_len = 0
    if head_cno == tail_cno:
        if plen == 0:
            total_len = head_clen - overlap
        else:
            total_len = head_clen + plen - 2*overlap
        print("total cyclic shortest length: ", total_len)
    else:
        if plen == 0:
            total_len = head_clen + tail_clen - overlap
        else:
            total_len = head_clen + tail_clen + plen - 2*overlap
    return total_len

def transitive_graph_reduction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    for k in simp_node_dict.keys():
        for i in simp_node_dict.keys():
            for j in simp_node_dict.keys():
                if i != j and i != k and j != k:
                    if (i, k) in simp_edge_dict and (k, j) in simp_edge_dict and (i, j) in simp_edge_dict:
                        graph.ep.color[simp_edge_dict[(i, j)]] = 'gray'
                        simp_edge_dict.pop((i, j))                  
    return

# curr_simp_path = []
# def retrieve_simple_paths(graph: Graph, simp_node_dict: dict, src, tgt, previsited: list, max_len, overlap):
#     print("src: {0}, tgt: {1}, max_len: {2}".format(graph.vp.id[src], graph.vp.id[tgt], max_len))
#     simple_path = []
#     visited = {}
#     for v in graph.vertices():
#         visited[v] = False
#     for p in previsited:
#         visited[simp_node_dict[p]] = True
#     visited[src] = False
#     visited[tgt] = False

#     def dfs(u, v):
#         global curr_simp_path
#         visited[u] = True
#         curr_simp_path.append(u)
#         # print(path_to_id_string(graph, curr_simp_path, "temp path"))
#         if path_len(graph, curr_simp_path, overlap) >= max_len:
#             print("max len reached, drop the path")
#         elif u == v:
#             print(path_to_id_string(graph, curr_simp_path, "path found"))
#             simple_path.append(curr_simp_path[:])
#         else:
#             for next in u.out_neighbors():
#                 if not visited[next]:
#                     dfs(next, v)
#         curr_simp_path = curr_simp_path[:-1]
#         visited[u] = False
#         return
#     dfs(src, tgt)

#     return simple_path

# def predict_sp_max_len(graph: Graph, simp_node_dict: dict, contig_dict: dict, tail_cno, head_cno, max_len, overlap):
#     tail_contig, tail_clen, _ = contig_dict[tail_cno]
#     src_len = len(graph.vp.seq[simp_node_dict[tail_contig[-1]]])
#     head_contig, head_clen, _ = contig_dict[head_cno]
#     tgt_len = len(graph.vp.seq[simp_node_dict[head_contig[0]]])
#     if head_cno == tail_cno:
#         return max_len - tail_clen + src_len + tgt_len
#     else:
#         return max_len - head_clen - tail_clen + src_len + tgt_len + overlap 