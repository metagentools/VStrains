#!/usr/bin/env python3

from graph_tool import GraphView, _in_degree
# import graph_tool
from graph_tool.all import Graph
from graph_tool.search import dfs_iterator
from graph_tool.topology import is_DAG, topological_sort, all_circuits
from graph_tool.draw import graph_draw

from utils.ns_Utilities import *
from search_algos import *
from tmp.legacy import *
"""
This file is used to store out of date scripts to keep a record for further usage needed.
"""

def graph_grouping(graph: Graph, simp_node_dict: dict, forward, reverse, partition_length_cut_off):
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

    group_no = 1
    groups = {}
    for v in simp_node_dict.values():
        # grouping
        if graph.vp.group[v] == -1:
            graph, group_no, groups = bfs_grouping(graph, v, group_no, groups)
            group_no = group_no + 1

    print("number of groups:", len(groups))
    # for key, item in groups.items():
    #     print("group number: ", key, " member: ", [(graph.vp.id[v],graph.vp.partition[v]) for v in item])

    # connect sub-graphs based on pair-end reads information
    # TODO
    return graph, groups

def graph_partition(graph: Graph, simp_node_dict: dict, edge_dict: dict, dp_dict: dict):
    """
    start from source node, run BFS to determine the partition for each nodes until reaching the sink node.
    partition is determined based on seg-length
    """
    None
    # queue = []
    # queue.append(src)
    # while queue:
    #     u = queue.pop()
    #     if graph.vp.id[u] == 'src':

    #         for v in u.out_neighbors():
    #             queue.append(v)
    #     elif graph.vp.id[u] == 'sink':
    #         None
    #     else:
    #         None

    # # pick a source node based on depth.
    # seg_no = source_node_via_dp(dp_dict)
    # src = simp_node_dict[seg_no]
    # graph.vp.visited[src] = 0
    # max_len_per_partition = 1000
    # max_node_per_partition = len(simp_node_dict) // 10
    # # estimate #partitions in the graph.
    # print("max len per partition: ", max_len_per_partition, ", max node per partition: ", max_node_per_partition)

    # queue = []
    # queue.append(src)
    # partition_dict = {}
    # partition_len_dict = {}
    # partition = 1
    # length = 0
    # partition_dict[partition] = []
    # while queue:
    #     u = queue.pop()
    #     # #print_vertex(graph, u, "current vertex")
    #     # label with partition tag
    #     # consider all the 
    #     curr_seq = graph.vp.seq[u]
    #     curr_len = len(curr_seq)

    #     if len(partition_dict[partition]) >= max_node_per_partition:
    #         partition = partition + 1
    #         length = 0
        
    #     if partition not in partition_dict:
    #         partition_dict[partition] = []

    #     if length + curr_len > max_len_per_partition and length + curr_len > (max_len_per_partition * 1.5):
    #         # overfit to put into the current partition 
    #         next_partition = partition + 1
    #         while next_partition in partition_dict:
    #             if len(partition_dict[next_partition]) >= max_node_per_partition:
    #                 next_partition = next_partition + 1
    #             else:
    #                 break
    #         if next_partition not in partition_dict:
    #             partition_dict[next_partition] = []
    #         # print("overfit in current partition ", partition, ", jump to next partition ", next_partition)
    #         partition_dict[next_partition].append(u)
    #         graph.vp.partition[u] = next_partition
    #     else:
    #         # best fit
    #         # print("best fit in partition ", partition)
    #         partition_dict[partition].append(u)
    #         graph.vp.partition[u] = partition
        
    #     for v in u.out_neighbors():
    #         if graph.vp.partition[v] == -1:
    #             queue.append(v)
    # # for s, v in simp_node_dict.items():
    # #     # #print_vertex(graph, v, "partition check")
    # #     if graph.vp.partition[v] == -1:
    # #         print("current node with seg no: ", s, " is in partition ", graph.vp.partition[v])
    return graph

def src_sink(graph: Graph, simp_node_dict: dict, edge_dict: dict, contig_dict: dict):
    """
    add manufactured source/sink node for finding maximum flow over graph
    """
    src = graph.add_vertex()
    sink = graph.add_vertex()
    graph.vp.id[src] = 'src'
    graph.vp.id[sink] = 'sink'
    graph.vp.ori[src] = 1
    graph.vp.ori[sink] = 1
    graph.vp.color[src] = 'white'
    graph.vp.color[sink] = 'white'

    src_is_connect = False
    sink_is_connect = False
    inf_value = 2**20 # may replace with numpy inf TODO
    for u in simp_node_dict.values():
        if u.in_degree() == 0 and u.out_degree() == 0:
            #print_vertex(graph, u, "isolated vertex")
            continue

        if u.in_degree() == 0:
            #print_vertex(graph, u, "vertex with 0 in degree, connect to src")
            e = graph.add_edge(src, u)
            graph.ep.flow[e] = inf_value
            src_is_connect = True
        
        if u.out_degree() == 0:
            #print_vertex(graph, u, "vertex with 0 out degree, connect to sink")
            e = graph.add_edge(u, sink)
            graph.ep.flow[e] = inf_value
            sink_is_connect = True
    # TODO create a src and sink node if is cyclic
    # attempt to break the cycle by removing one edge contains the maximum depth among all the edges
    if not src_is_connect or not sink_is_connect:
        # obtain all the involving node from contigs
        contig_node_set = set()
        for (contig, contig_rev) in contig_dict.values():
            for n in contig:
                contig_node_set.add(n)

        largest_circuit = max(all_circuits(graph), key=lambda c: len(c))
        circuit_nodes = [simp_node_dict[graph.vp.id[v]] for v in largest_circuit]
        for v in sorted(circuit_nodes, key=lambda n: graph.vp.dp[n], reverse=True):
            if src_is_connect and sink_is_connect:
                break
            if graph.vp.id[v] in contig_node_set:
                print("already in the contig set, never consider to break the relevant edge")
                continue
            if v.out_degree() == 1 and v.in_degree() != 0 and not sink_is_connect:
                # break the edge, connect to sink
                u = list(v.out_neighbors())[0]
                #print_edge(graph, edge_dict[(v, u)], "edge be removed")
                graph.remove_edge(edge_dict[(v, u)])
                edge_dict.pop((v, u))
                e = graph.add_edge(v, sink)
                graph.ep.flow[e] = inf_value
                sink_is_connect = True

            if v.in_degree() == 1 and v.out_degree() != 0 and not src_is_connect:
                # break the edge, connect to src        
                u = list(v.in_neighbors())[0]
                #print_edge(graph, edge_dict[(u, v)], "edge be removed")
                graph.remove_edge(edge_dict[(u, v)])
                edge_dict.pop((u, v))
                e = graph.add_edge(src, v)
                graph.ep.flow[e] = inf_value
                src_is_connect = True

    if not src_is_connect:
        print("source is not connected still")
    if not sink_is_connect:
        print("sink is not connected still")
    if src_is_connect and sink_is_connect:
        print("src and sink is initialised")
    if not is_DAG(graph):
        print("The graph still have cycles")
    return graph, simp_node_dict, edge_dict, src, sink

def max_flow(graph: Graph, simp_node_dict: dict, edge_dict: dict, contig_dict: dict):
    print("find max flow")
    c1 = contig_dict[("1", 9557, 1897.189396)]
    src = simp_node_dict[c1[0]]
    
    sink = simp_node_dict[c1[-1]]
    cap = graph.ep.flow
    # res = flow.boykov_kolmogorov_max_flow(graph, src, sink, cap)
    # res.a = cap.a - res.a  # the actual flow
    # max_flow = sum(res[e] for e in sink.in_edges())
    # print(max_flow)

    # for i in range(len(c1)-1):
    #     u = simp_node_dict[c1[i]]
    #     v = simp_node_dict[c1[i+1]]
    #     e = edge_dict[(u,v)]
    #     f = graph.ep.flow[e]
        #print_edge(graph, e, "edge flow: {0}, res flow: {1}".format(f, res[e]))
    return

def graph_dfs(graph: Graph, source):
    """
    Count the maximum depth path the source node can reach via directed edge in the graph
    """
    visited = {}
    for u in graph.vertices():
        visited[u] = False
    
    def dfs_helper(graph: Graph, u, visited):
        if u.out_degree() == 0:
            return [u]
        else:
            rtn = []
            for v in u.out_neighbors():
                # print_vertex(graph, v)
                if not visited[v]:
                    visited[v] = True
                    cmp = dfs_helper(graph, v, visited)
                    cmp.insert(0, (v, len(graph.vp.seq[v])))
                    rtn = rtn if len(cmp) < len(rtn) else cmp
            return rtn

    return dfs_helper(graph, source, visited)

def graph_reduction(graph: Graph, contig_dict: dict, simp_node_dict: dict, simp_edge_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, output_file, min_cov, min_len, overlap):
    """
    reduce the node/edge weight based on existing contig found by SPAdes.
    only contig with minimum strain length satisfied be removed
    return non-removed contig dict
    """

    print("-------------------------graph reduction----------------------")
    udpate_node_to_contig_dict(node_to_contig_dict, simp_node_dict)
    update_edge_to_contig_dict(edge_to_contig_dict, simp_edge_dict)
    cand_strains_dict = {}
    temp_contigs_dict = {}
    for cno, (contig, clen, ccov) in contig_dict.items():
        if clen >= min_len:
            cand_strains_dict[cno] = (contig, clen, ccov)
            print("full-length contig found: ", cno, clen, ccov)
            contig_reduction(graph, contig, cno, clen, ccov, simp_node_dict, simp_edge_dict, min_cov)  
        else:
            temp_contigs_dict[cno] = [contig, clen, ccov]
            print("imcomplete contig found: ", cno, clen, ccov) 

    # reappend nodes
    for no, [cnos, dp, node] in list(node_to_contig_dict.items()):

        for strain_cno in cand_strains_dict.keys():
            if strain_cno in cnos:
                cnos.remove(strain_cno)

        # update cnos
        node_to_contig_dict[no][0] = cnos
        if len(cnos) == 0:
            node_to_contig_dict.pop(no)
            continue
        
        if no not in simp_node_dict and no in node_to_contig_dict:
            print("contig node {0} be removed, append it back".format(no))
            graph.vp.dp[node] = dp
            graph.vp.color[node] = "black"
            simp_node_dict[no] = node
            print_vertex(graph, node, "from graph reduction: ")

    # reappends edges
    for (u,v), [cnos, flow, edge] in list(edge_to_contig_dict.items()):

        for strain_cno in cand_strains_dict.keys():
            if strain_cno in cnos:
                cnos.remove(strain_cno)

        # update cnos
        edge_to_contig_dict[(u,v)][0] = cnos
        if len(cnos) == 0:
            edge_to_contig_dict.pop((u,v))
            continue
        
        if (u,v) not in simp_edge_dict and (u,v) in edge_to_contig_dict:
            print("contig edge {0} be removed, append it back".format((u,v)))
            graph.ep.flow[edge] = flow
            graph.ep.color[edge] = "black"
            simp_edge_dict[(u,v)] = edge

    # store output graph 
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, output_file)

    print("-----------------------graph reduction end--------------------")
    return cand_strains_dict, temp_contigs_dict

def assign_edge_flow(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, consistent_node: dict, inconsistent_node: dict):
    """
    Assign the edge flow based on node weight and contig alignment.
    """
    #TODO assign non-branch path node dp first, keep the start and end node dp.
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for no, node in simp_node_dict.items():
        if no not in consistent_node:
            continue
        w = graph.vp.dp[node]
        in_d = len([n for n in node.in_neighbors() if graph.vp.color[n] == 'black'])
        node_in_edges = [e for e in node.in_edges() if graph.ep.color[e] == 'black']
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if graph.ep.flow[e] == 0.0 and (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = round((w + graph.vp.dp[src])/2, 2)
                    else:
                        graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1

        out_d = len([n for n in node.out_neighbors() if graph.vp.color[n] == 'black'])
        node_out_edges = [e for e in node.out_edges() if graph.ep.color[e] == 'black']
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = round((w + graph.vp.dp[src])/2, 2)
                    else:
                        graph.ep.flow[e] = w
                    un_assigned_edge = un_assigned_edge - 1
    print("un-assigned edges after init iteration : ", un_assigned_edge)

    # coverage iteration
    converage_flag = 0
    while True:          
        for no, node in simp_node_dict.items():
            if no not in consistent_node:
                continue
            in_d = len([n for n in node.in_neighbors() if graph.vp.color[n] == 'black'])
            node_in_edges = [e for e in node.in_edges() if graph.ep.color[e] == 'black']
            in_w = graph.vp.dp[node]
            in_e = []
            for e in node_in_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    in_d = in_d - 1
                    in_w = in_w - f
                else:
                    in_e.append(e)
            if in_d == 1:
                for e in in_e:
                    if graph.ep.flow[e] == 0.0 and (graph.vp.id[e.source()], no) in simp_edge_dict:
                        if in_w <= 0:
                            print("in low edge flow error: ", graph.vp.id[e.source()], " -> ", graph.vp.id[e.target()], in_w)
                        graph.ep.flow[e] = in_w
                        un_assigned_edge = un_assigned_edge - 1

            out_d = len([n for n in node.out_neighbors() if graph.vp.color[n] == 'black'])
            node_out_edges = [e for e in node.out_edges() if graph.ep.color[e] == 'black']
            out_w = graph.vp.dp[node]
            out_e = []
            for e in node_out_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    out_d = out_d - 1
                    out_w = out_w - f
                else:
                    out_e.append(e)
            if out_d == 1:
                for e in out_e:
                    if graph.ep.flow[e] == 0.0 and (no, graph.vp.id[e.target()]) in simp_edge_dict:
                        if out_w <= 0:
                            print_vertex(graph, node, "Current node")
                            print("out low edge flow error: ", graph.vp.id[e.source()], " -> ", graph.vp.id[e.target()], out_w)
                        graph.ep.flow[e] = out_w
                        un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  
    print("un-assigned edges after node-weight coverage iteration : ", un_assigned_edge)

    # FIXME for inconsistent node part

    for (u,v), e in simp_edge_dict.items():
        if graph.ep.flow[e] == 0.0:
            if u in inconsistent_node and v in inconsistent_node:
                continue

            u = e.source()
            v = e.target()
            u_flow_remain = graph.vp.dp[u]
            u_degree = len([n for n in u.out_neighbors() if graph.vp.color[n] == 'black'])
            u_out_edges = [e for e in u.out_edges() if graph.ep.color[e] == 'black']

            for ue in u_out_edges:
                if graph.ep.flow[ue] != 0.0:
                    # assigned edge
                    u_degree = u_degree - 1
                    u_flow_remain = u_flow_remain - graph.ep.flow[ue]

            v_flow_remain = graph.vp.dp[v]
            v_degree = len([n for n in v.in_neighbors() if graph.vp.color[n] == 'black'])
            v_in_edges = [e for e in v.in_edges() if graph.ep.color[e] == 'black']

            for ve in v_in_edges:
                if graph.ep.flow[ve] != 0.0:
                    # assigned edge
                    v_degree = v_degree - 1
                    v_flow_remain = v_flow_remain - graph.ep.flow[ve]

            u_flow_remain = u_flow_remain if u_flow_remain > 0 else 0
            v_flow_remain = v_flow_remain if v_flow_remain > 0 else 0
            if u_flow_remain == 0 or v_flow_remain == 0:
                assign_flow = (u_flow_remain / u_degree) + (v_flow_remain / v_degree)
            else:
                assign_flow = ((u_flow_remain / u_degree) + (v_flow_remain / v_degree)) / 2
            if assign_flow <= 0:
                print("low edge flow error: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                print("u_flow_remain: ", u_flow_remain, " u_degree: ", u_degree)
                print("v_flow_remain: ", v_flow_remain, " v_degree: ", v_degree)
            else:
                print("manual assigned edge: ", graph.vp.id[u], " -> ", graph.vp.id[v], assign_flow)
                if graph.ep.flow[e] == 0.0:
                    un_assigned_edge = un_assigned_edge - 1
                    graph.ep.flow[e] = assign_flow
    print("un-assigned edges after manual assign iteration : ", un_assigned_edge)
    u_nodes = []
    for (u,v), e in simp_edge_dict.items():
        if graph.ep.flow[e] == 0.0:
            u_nodes.append(int(u))
            u_nodes.append(int(v))
        if u == '49355':
            print_edge(graph, e, "")
    print("unassigned related node: ", u_nodes)
    print("-----------------------assign edge flow end--------------------")

def node_rebalance_via_simple_path(graph: Graph, simp_edge_dict: dict):
    # get all the simple path (no branch path)
    simple_paths = simp_path(graph, simp_edge_dict)
    for sp in simple_paths:
        if DEBUG_MODE:
            print("path: ", [int(graph.vp.id[u]) for u in sp])
        dp_mean_p = numpy.mean([graph.vp.dp[u] for u in sp])
        # rebalance the node depth
        for u in sp:
            graph.vp.dp[u] = dp_mean_p

def node_dp_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    """
    Re-balance all the node depth, such that for any given node u, sum(dp(in_neighbors)) == dp(u) == sum(dp(out_neigbors))
    where dp(u) = max(sum(dp(in_neighbors)), dp(u), sum(dp(out_neigbors)))
    plan to run BFS from the source node (how to pick?) until all the node is rebalanced.
    source node may pick from the node with 0 indegree.
    """
    def source_node_via_dp(dp_dict: dict):
        """
        return the pos-neg node with minimum depth
        """
        seg_no = max(dp_dict, key=dp_dict.get)
        print("source node id: ", seg_no, ", depth: ", dp_dict[seg_no])
        return seg_no
    print("-------------------------node depth rebalance----------------------")
    graph.vp.visited = graph.new_vertex_property("int16_t", val=-1)
    dp_dict = {}
    for no, u in simp_node_dict.items():
        dp_dict[no] = graph.vp.dp[u]
        graph.vp.visited[u] = -1

    while set(dp_dict):
        for no in dp_dict.keys():
            dp_dict[no] = graph.vp.dp[simp_node_dict[no]]
       
        seg_no = source_node_via_dp(dp_dict)
        src = simp_node_dict[seg_no]
        graph.vp.visited[src] = 0
        fifo_queue = [src]
        # BFS
        while fifo_queue:
            u = fifo_queue.pop()
            dp_dict.pop(graph.vp.id[u])

            u_dp = graph.vp.dp[u]

            in_es = [e for e in u.in_edges() if graph.ep.color[e] == 'black']
            in_flow_sum = numpy.sum([graph.ep.flow[e] for e in in_es])

            out_es = [e for e in u.out_edges() if graph.ep.color[e] == 'black']
            out_flow_sum = numpy.sum([graph.ep.flow[e] for e in out_es])

            max_dp = max(in_flow_sum, u_dp, out_flow_sum)
            graph.vp.dp[u] = max_dp
            if DEBUG_MODE:
                print("Node: {0} dp from {1} to {2}".format(graph.vp.id[u], u_dp, round(max_dp)))
                print("in {0} dp {1} out {2}".format(in_flow_sum, u_dp, out_flow_sum))

            if len(in_es) == 1:
                graph.ep.flow[in_es[0]] = max_dp
            elif len(in_es) > 1:
                in_flow_load = max_dp
                for i in range(len(in_es) - 1):
                    e = in_es[i]
                    flow = round((graph.ep.flow[e]/in_flow_sum) * max_dp)
                    graph.ep.flow[e] = flow
                    in_flow_load -= flow
                graph.ep.flow[in_es[-1]] = round(in_flow_load)
            else:
                print("No in edges, skip")

            if len(out_es) == 1:
                graph.ep.flow[out_es[0]] = max_dp
            elif len(out_es) > 1:
                out_flow_load = max_dp
                for i in range(len(out_es) - 1):
                    e = out_es[i]
                    flow = round((graph.ep.flow[e]/out_flow_sum) * max_dp)
                    graph.ep.flow[e] = flow
                    out_flow_load -= flow
                graph.ep.flow[out_es[-1]] = round(out_flow_load)
            else:
                print("No out edges, skip")

            graph.vp.visited[u] = 1

            # append neighbors into fifo queue
            for adj_node in u.all_neighbors():
                if graph.vp.visited[adj_node] == -1 and graph.vp.color[adj_node] == 'black':
                    graph.vp.visited[adj_node] = 0
                    fifo_queue.append(adj_node)
    print("Done...")
    print("Check node dp balance now...")
    # sum check
    false_count = 0
    for no, u in simp_node_dict.items():
        dp = graph.vp.dp[u]
        in_es = [n for n in u.in_edges() if graph.ep.color[n] == 'black']
        in_flow_sum = numpy.sum([graph.ep.flow[v] for v in in_es])

        out_es = [n for n in u.out_edges() if graph.ep.color[n] == 'black']
        out_flow_sum = numpy.sum([graph.ep.flow[v] for v in out_es])

        if len(in_es) != 0:
            if abs(in_flow_sum - dp) > 100:
                if len(out_es) == 0:
                    # edge case, manual fix
                    print("mannual fix, sink node", no)
                    graph.vp.dp[no] = max(in_flow_sum, dp)
                else:
                    print("Current node: ", no, dp)
                    false_count += 1
                    for e in u.in_edges():
                        if graph.ep.color[e] == 'black':
                            print_edge(graph, e, "in edge")  
                    print("Dp diff: ", in_flow_sum - dp)

        if len(out_es) != 0:
            if abs(out_flow_sum - dp) > 100:
                if len(in_es) == 0:
                    # edge case, manual fix
                    print("mannual fix, source node", no)
                    graph.vp.dp[no] = max(out_flow_sum, dp)
                else:
                    print("Current node: ", no, dp)
                    false_count += 1
                    for e in u.out_edges():
                        if graph.ep.color[e] == 'black':
                            print_edge(graph, e, "out edge")  
                    print("Dp diff: ", out_flow_sum - dp)             
    print("False case: ", false_count)
    print("-------------------------node depth rebalance end----------------------")
    return


# refine all the node coverage by reassign all the edge flow
    for (u,v), e in simp_edge_dict.items():
        graph.ep.flow[e] = 0.0
    
    print("Final refine")
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict, consistent_node, inconsistent_node)

    # re assign the node depth = max(curr_dp, out_neighbor_dp_sum if no branch, in_neighbor_dp_sum if no branch, out edge flow sum, in edge flow sum)
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        node = queue.pop()
        visited.append(graph.vp.id[node])
        curr_dp = graph.vp.dp[node]

        us = list(node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        in_neighbor_dp_sum = -1
        if u_out_degrees == node.in_degree():
            in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

        vs = list(node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        out_neighbor_dp_sum = -1
        if v_in_degrees == node.out_degree():
            out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in node.out_edges()])
        graph.vp.dp[node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[node]))
        assign_edge_flow2(graph, simp_node_dict, simp_edge_dict)
        
        for n in node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)

def assign_edge_flow2(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    un_assigned_edge = len(simp_edge_dict)
    print("-------------------------assign edge flow----------------------")
    print("Assign edge flow: Total unassigned edges: ", un_assigned_edge)
    # it is necessary to distinguish two phase to avoid assembly graph mistake, or do we ignore the mistake?
    # init iteration
    for no, node in simp_node_dict.items():
        w = graph.vp.dp[node]
        in_d = len([n for n in node.in_neighbors()])
        node_in_edges = [e for e in node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if src.out_degree() == 1:
                    graph.ep.flow[e] = numpy.max([w, graph.ep.flow[e], graph.vp.dp[src]])
                    un_assigned_edge = un_assigned_edge - 1
        out_d = len([n for n in node.out_neighbors()])
        node_out_edges = [e for e in node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if tgt.in_degree() == 1:
                    graph.ep.flow[e] = numpy.max([w, graph.ep.flow[e], graph.vp.dp[tgt]])
                    un_assigned_edge = un_assigned_edge - 1

    # coverage iteration
    converage_flag = 0
    while True:          
        for no, node in simp_node_dict.items():
            in_d = len([n for n in node.in_neighbors()])
            node_in_edges = [e for e in node.in_edges()]
            in_w = graph.vp.dp[node]
            in_e = []
            for e in node_in_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    in_d = in_d - 1
                    in_w = in_w - f
                else:
                    in_e.append(e)
            if in_d == 1:
                for e in in_e:
                    graph.ep.flow[e] = max(in_w,graph.ep.flow[e])
                    un_assigned_edge = un_assigned_edge - 1

            out_d = len([n for n in node.out_neighbors()])
            node_out_edges = [e for e in node.out_edges()]
            out_w = graph.vp.dp[node]
            out_e = []
            for e in node_out_edges:
                f = graph.ep.flow[e]
                if f != 0.0:
                    out_d = out_d - 1
                    out_w = out_w - f
                else:
                    out_e.append(e)
            if out_d == 1:
                for e in out_e:
                        graph.ep.flow[e] = max(out_w,graph.ep.flow[e])
                        un_assigned_edge = un_assigned_edge - 1
        if converage_flag == un_assigned_edge:
            break
        else:
            converage_flag = un_assigned_edge  

    # double cross consistent node assign
    for (u,v), e in simp_edge_dict.items():

        src = e.source()
        tgt = e.target()

        u_flow_remain = graph.vp.dp[src]
        u_degree = len([n for n in src.out_neighbors()])
        u_out_edges = [e for e in src.out_edges()]
        for ue in u_out_edges:
            if graph.ep.flow[ue] != 0.0:
                # assigned edge
                u_degree = u_degree - 1
                u_flow_remain = u_flow_remain - graph.ep.flow[ue]

        v_flow_remain = graph.vp.dp[tgt]
        v_degree = len([n for n in tgt.in_neighbors()])
        v_in_edges = [e for e in tgt.in_edges()]
        for ve in v_in_edges:
            if graph.ep.flow[ve] != 0.0:
                # assigned edge
                v_degree = v_degree - 1
                v_flow_remain = v_flow_remain - graph.ep.flow[ve]

        u_flow_remain = u_flow_remain if u_flow_remain > 0 else 0
        v_flow_remain = v_flow_remain if v_flow_remain > 0 else 0
        udiv = (u_flow_remain / u_degree) if u_degree != 0 else 0
        vdiv = (v_flow_remain / v_degree) if v_degree != 0 else 0
        if udiv != 0 and vdiv != 0:
            assign_flow = (udiv + vdiv) / 2
        elif udiv != 0:
            assign_flow = udiv
        elif vdiv != 0:
            assign_flow = vdiv
        else:
            assign_flow = 0
        un_assigned_edge = un_assigned_edge - 1
        graph.ep.flow[e] = max(assign_flow,graph.ep.flow[e])
    print("-----------------------assign edge flow v2--------------------")
    return

def graph_split_by_breadth(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict):
    srcs = []
    for no, node in simp_node_dict.items():
        if node.in_degree() == 0:
            srcs.append(node)
    
    queue = []
    visited = []
    srcs = sorted(srcs, key=lambda x: graph.vp.dp[x], reverse=True)
    queue.append(srcs[0])
    split_by_level = {}
    split_by_level[0] = [graph.vp.id[srcs[0]]]
    curr_level = 1
    while queue:
        node = queue.pop()
        visited.append(graph.vp.id[node])
        split_by_level[curr_level] = []
        for out_n in node.out_neighbors():
            if not graph.vp.id[out_n] in visited:
                split_by_level[curr_level].append(graph.vp.id[out_n])
                queue.append(out_n)
        curr_level += 1
    sorted_level_pair = sorted(split_by_level.items(), key = lambda x: x[0])
    allNodes = []
    for level, nodes in sorted_level_pair:
        print("Current level: ", level)
        print("Nodes: ", [int(v) for v in nodes])
        [allNodes.append((int(v))) for v in nodes]
    print("All nodes: ", allNodes)

    for no, node in simp_node_dict.items():
        curr_dp = graph.vp.dp[node]
        in_edges = list(node.in_edges())
        in_sum = numpy.sum([graph.ep.flow[e] for e in in_edges])
        out_edges = list(node.out_edges())
        out_sum = numpy.sum([graph.ep.flow[e] for e in out_edges])
        if in_edges != [] and curr_dp > in_sum:
            for e in in_edges:
                graph.ep.flow[e] = round(curr_dp * (graph.ep.flow[e] / in_sum), 2)
        
        if out_edges != [] and curr_dp > out_sum:
            for e in out_edges:
                graph.ep.flow[e] = round(curr_dp * (graph.ep.flow[e] / out_sum), 2)

    for edge in simp_edge_dict.values():
        print_edge(graph, edge, "After refine")

    # re assign the node depth = max(curr_dp, out_neighbor_dp_sum if no branch, in_neighbor_dp_sum if no branch, out edge flow sum, in edge flow sum)
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        curr_node = queue.pop()
        visited.append(graph.vp.id[curr_node])
        curr_dp = graph.vp.dp[curr_node]

        us = list(curr_node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        in_neighbor_dp_sum = -1
        if u_out_degrees == curr_node.in_degree():
            in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

        vs = list(curr_node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        out_neighbor_dp_sum = -1
        if v_in_degrees == curr_node.out_degree():
            out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
        graph.vp.dp[curr_node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[curr_node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))

        in_d = len([n for n in curr_node.in_neighbors()])
        node_in_edges = [e for e in curr_node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[src])
                    else:
                        graph.ep.flow[e] = w

        out_d = len([n for n in curr_node.out_neighbors()])
        node_out_edges = [e for e in curr_node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[tgt] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[tgt])
                    else:
                        graph.ep.flow[e] = w

        for n in curr_node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)
    print("Next")
    queue = []
    queue.append((sorted(simp_node_dict.values(),key=lambda n: graph.vp.dp[n], reverse=True))[0])
    visited = []
    while queue:
        curr_node = queue.pop()
        visited.append(graph.vp.id[curr_node])
        curr_dp = graph.vp.dp[curr_node]

        us = list(curr_node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        in_neighbor_dp_sum = -1
        if u_out_degrees == curr_node.in_degree():
            in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

        vs = list(curr_node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        out_neighbor_dp_sum = -1
        if v_in_degrees == curr_node.out_degree():
            out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
        
        in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
        out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
        graph.vp.dp[curr_node] = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
        print(graph.vp.id[curr_node], " node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))

        in_d = len([n for n in curr_node.in_neighbors()])
        node_in_edges = [e for e in curr_node.in_edges()]
        if in_d == 1:
            for e in node_in_edges:
                src = e.source()
                if (graph.vp.id[src], no) in simp_edge_dict:
                    if graph.vp.id[src] in consistent_node and src.out_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[src])
                    else:
                        graph.ep.flow[e] = w

        out_d = len([n for n in curr_node.out_neighbors()])
        node_out_edges = [e for e in curr_node.out_edges()]
        if out_d == 1:
            for e in node_out_edges:
                tgt = e.target()
                if (no, graph.vp.id[tgt]) in simp_edge_dict:
                    if graph.vp.id[tgt] in consistent_node and tgt.in_degree() == 1:
                        graph.ep.flow[e] = max(w, graph.vp.dp[tgt])
                    else:
                        graph.ep.flow[e] = w

        for n in curr_node.all_neighbors():
            if graph.vp.id[n] not in visited:
                queue.append(n)

def coverage_rebalance(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, min_cov, cutoff=125):
    """
    rebalance all the node coverage to ensure flow consistency
    """
    # helper functions
    def node_expansion_both(graph: Graph, con_nodes, incon_nodes, curr_node):
        """
        from given node, BFS search all the in/out branches, if the curr node is covered by c node layer
        then return true, if any inconsistent node be found then return false
        """
        visited = []
        queue = [curr_node]
        while queue:
            node = queue.pop()
            id = graph.vp.id[node]
            visited.append(id)
            if id in con_nodes:
                #skip
                None
            elif id in incon_nodes:
                # print("From src node {0} reach in consistent node {1} return false".format(graph.vp.id[curr_node], id))
                return False
            else:
                for u in node.all_neighbors():
                    if graph.vp.id[u] not in visited:
                        queue.append(u)
        return True
    
    def subfix(graph: Graph, curr_node, no_changes_in, consistent_node: dict, inconsistent_node: dict, simp_node_dict: dict):
        no_changes = no_changes_in
        for edge in curr_node.all_edges():
            u = graph.vp.id[edge.source()]
            v = graph.vp.id[edge.target()]

            if u in inconsistent_node and v in inconsistent_node:
                None
            elif v in inconsistent_node:
                # u is consistent
                u_node = simp_node_dict[u]
                all_incon_dp = 0
                incon_edges = []
                out_flow_remain = graph.vp.dp[u_node]
                for out_edge in u_node.out_edges():
                    tgt = out_edge.target()
                    if graph.ep.flow[out_edge] != 0.0:
                        out_flow_remain -= graph.ep.flow[out_edge]
                    else:
                        all_incon_dp += graph.vp.dp[tgt]
                        incon_edges.append(out_edge)
                for incon_edge in incon_edges:
                    tgt = incon_edge.target()
                    maxflow =  max(graph.ep.flow[incon_edge], round((graph.vp.dp[tgt]/all_incon_dp) * out_flow_remain, 2))
                    if maxflow != graph.ep.flow[incon_edge]:
                        graph.ep.flow[incon_edge] = maxflow
                        no_changes = False
            elif u in inconsistent_node:
                # v is consistent
                v_node = simp_node_dict[v]
                all_incon_dp = 0
                incon_edges = []
                in_flow_remain = graph.vp.dp[v_node]
                for in_edge in v_node.in_edges():
                    src = in_edge.source()
                    if graph.ep.flow[in_edge] != 0.0:
                        in_flow_remain -= graph.ep.flow[in_edge]
                    else:
                        all_incon_dp += graph.vp.dp[src]
                        incon_edges.append(in_edge)
                for incon_edge in incon_edges:
                    src = incon_edge.source()
                    maxflow = numpy.max([graph.ep.flow[incon_edge], round((graph.vp.dp[src]/all_incon_dp) * in_flow_remain, 2)])
                    if maxflow != graph.ep.flow[incon_edge]:
                        graph.ep.flow[incon_edge] = maxflow
                        no_changes = False     
            else:
                None

        if curr_id in consistent_node:
            w = graph.vp.dp[curr_node]
            node_in_edges = list(curr_node.in_edges())
            if len(node_in_edges) < 1:
                None
            elif len(node_in_edges) == 1:
                e = node_in_edges[0]
                src = e.source()
                if src.out_degree() == 1:
                    if graph.vp.id[src] in consistent_node:
                        maxflow = numpy.max([graph.ep.flow[e], w, graph.vp.dp[src]])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
                    else:
                        maxflow = numpy.max([graph.ep.flow[e], w])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
            else:
                # mutliple in edges
                total_dp = numpy.sum([graph.vp.dp[u] for u in curr_node.in_neighbors()])
                for e in node_in_edges:
                    src = e.source()
                    maxflow = numpy.max([graph.ep.flow[e], round((graph.vp.dp[src]/total_dp) * w)])
                    if graph.ep.flow[e] != maxflow:
                        graph.ep.flow[e] = maxflow
                        no_changes = False                            

            node_out_edges = list(curr_node.out_edges())
            if len(node_out_edges) < 1:
                None
            elif len(node_out_edges) == 1:
                e = node_out_edges[0]
                tgt = e.target()
                if tgt.in_degree() == 1:
                    if graph.vp.id[tgt] in consistent_node:
                        maxflow = numpy.max([graph.ep.flow[e], w, graph.vp.dp[tgt]])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
                    else:
                        maxflow = numpy.max([graph.ep.flow[e], w])
                        if graph.ep.flow[e] != maxflow:
                            graph.ep.flow[e] = maxflow
                            no_changes = False
            else:
                # mutliple in edges
                total_dp = numpy.sum([graph.vp.dp[u] for u in curr_node.out_neighbors()])
                for e in node_out_edges:
                    tgt = e.target()
                    maxflow = numpy.max([graph.ep.flow[e], round((graph.vp.dp[tgt]/total_dp) * w)])
                    if graph.ep.flow[e] != maxflow:
                        graph.ep.flow[e] = maxflow
                        no_changes = False
        return no_changes
    #### Helper functions
    print("-------------------------coverage rebalance----------------------")
    print("Total node: ", len(simp_node_dict))
    consistent_count = 0

    consistent_node = {}
    inconsistent_node = {}
    untracked_node = {}

    # select the most confident nodes
    for no, node in simp_node_dict.items():
        dp = graph.vp.dp[node]
        untracked = False

        us = list(node.in_neighbors())
        u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
        out_consistent = False
        if u_out_degrees == 0:
            out_consistent = False
        else:
            if u_out_degrees == node.in_degree():
                u_sum = numpy.sum([graph.vp.dp[u] for u in us])
                out_consistent = abs(u_sum - dp) < cutoff
            else:
                out_consistent = True
                untracked = True

        vs = list(node.out_neighbors())
        v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
        in_consistent = False
        if v_in_degrees == 0:
            in_consistent = False
        else:
            if v_in_degrees == node.out_degree():
                v_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                in_consistent = abs(v_sum - dp) < cutoff
            else:
                in_consistent = True
                untracked = True

        if in_consistent and out_consistent:
            if untracked:
                untracked_node[graph.vp.id[node]] = node
            else:
                consistent_node[graph.vp.id[node]] = node
                consistent_count += 1
        else:
            inconsistent_node[graph.vp.id[node]] = node
    
    print("Total consistent node after most confident extraction: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node.keys()])
    print([int(n) for n in inconsistent_node.keys()]) 
    print([int(n) for n in untracked_node.keys()])

    # deal with untracked node
    # from current untracked node, BFS to both end, for every branch, stop until non-untracked
    # node be found, if found a in-consistent node during the search, then keep it untracked/mark it
    # inconsistent, otherwise mark it consistent.
    for no, node in list(untracked_node.items()):
        untracked_node.pop(no)
        if node_expansion_both(graph, consistent_node, inconsistent_node, node):
            dp = graph.vp.dp[node]

            us = list(node.in_neighbors())
            u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
            out_consistent = False
            out_amb = False
            if u_out_degrees == 0:
                out_consistent = False
            else:
                if u_out_degrees == node.in_degree():
                    u_sum = numpy.sum([graph.vp.dp[u] for u in us])
                    out_consistent = abs(u_sum - dp) < cutoff
                else:
                    out_amb = True

            vs = list(node.out_neighbors())
            v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
            in_consistent = False
            in_amb = False
            if v_in_degrees == 0:
                in_consistent = False
            else:
                if v_in_degrees == node.out_degree():
                    v_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                    in_consistent = abs(v_sum - dp) < cutoff
                else:
                    in_amb = True
            
            if out_amb and in_amb:
                # no way to determine its consistency
                inconsistent_node[no] = node
            else:
                if in_consistent and out_consistent:
                    consistent_node[no] = node
                    consistent_count += 1
                elif in_consistent and out_amb:
                    consistent_node[no] = node
                    consistent_count += 1
                elif out_consistent and in_amb:
                    consistent_node[no] = node
                    consistent_count += 1
                else:
                    inconsistent_node[no] = node
        else:
            inconsistent_node[no] = node
    print("Total consistent node after untracked classification: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node])
    print([int(n) for n in inconsistent_node]) 
    print([int(n) for n in untracked_node])

    # assign edge flow only to consistent node related edge
    assign_edge_flow(graph, simp_node_dict, simp_edge_dict, consistent_node, inconsistent_node)
    
    sorted_consistent_node = sorted(consistent_node.items(), key=lambda x: graph.vp.dp[x[1]], reverse=True)
    # fix the inconsistent node depth, and also try to 
    # increment all the imbalance consistent node
    i = 0
    break_point = consistent_count
    no_changes = False
    cyclic = not is_DAG(graph)

    while len(inconsistent_node) != 0 or not no_changes:
        no_changes = True
        no, node = sorted_consistent_node[i]
        print("Current iteration: ", i, "Total consistent node: ", consistent_count)
        print([int(n) for n in inconsistent_node])
        print("Current solid node: ", no)

        queue = [node]
        visited = []
        while queue:
            curr_node = queue.pop()
            curr_id = graph.vp.id[curr_node]
            curr_dp = graph.vp.dp[curr_node]
            visited.append(curr_id)
            # fix the curr node if is in consistent
            if curr_id in inconsistent_node:
                fix = True
                vtotal = 0
                for e in curr_node.in_edges():
                    if graph.ep.flow[e] != 0.0:
                        vtotal += graph.ep.flow[e]
                    else:
                        # there exist an unassigned edge around the curr node
                        fix = False
                        break
                # perform fix
                if fix:
                    graph.vp.dp[curr_node] = vtotal
                    if DEBUG_MODE:
                        print(curr_id, "inconsistent node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))
                    inconsistent_node.pop(curr_id)
                    consistent_node[curr_id] = curr_node
                    consistent_count += 1
                    no_changes = False
            else:
                us = list(curr_node.in_neighbors())
                u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
                in_neighbor_dp_sum = -1
                if u_out_degrees == curr_node.in_degree():
                    in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

                vs = list(curr_node.out_neighbors())
                v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
                out_neighbor_dp_sum = -1
                if v_in_degrees == curr_node.out_degree():
                    out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                
                in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
                out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
                maxdp = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
                if curr_dp != maxdp:
                    graph.vp.dp[curr_node] = maxdp
                    no_changes = False
            
            no_changes = subfix(graph, curr_node, no_changes, consistent_node, inconsistent_node, simp_node_dict)
            # append all the out neighbors
            for n in curr_node.out_neighbors():
                if graph.vp.id[n] not in visited:
                    queue.append(n)

        queue = [node]
        visited = []
        while queue:
            curr_node = queue.pop()
            curr_id = graph.vp.id[curr_node]
            curr_dp = graph.vp.dp[curr_node]
            visited.append(curr_id)
            # fix the curr node if is in consistent
            if curr_id in inconsistent_node:
                fix = True
                vtotal = 0
                for e in curr_node.out_edges():
                    if graph.ep.flow[e] != 0.0:
                        vtotal += graph.ep.flow[e]
                    else:
                        # there exist an unassigned edge around the curr node
                        fix = False
                        break
                # perform fix
                if fix:
                    graph.vp.dp[curr_node] = vtotal
                    if DEBUG_MODE:
                        print(curr_id, "inconsistent node depth from {0} to {1}".format(curr_dp, graph.vp.dp[curr_node]))
                    inconsistent_node.pop(curr_id)
                    consistent_node[curr_id] = curr_node
                    consistent_count += 1
                    no_changes = False
            else:
                us = list(curr_node.in_neighbors())
                u_out_degrees = numpy.sum([u.out_degree() for u in us]) 
                in_neighbor_dp_sum = -1
                if u_out_degrees == curr_node.in_degree():
                    in_neighbor_dp_sum = numpy.sum([graph.vp.dp[u] for u in us])

                vs = list(curr_node.out_neighbors())
                v_in_degrees = numpy.sum([v.in_degree() for v in vs]) 
                out_neighbor_dp_sum = -1
                if v_in_degrees == curr_node.out_degree():
                    out_neighbor_dp_sum = numpy.sum([graph.vp.dp[v] for v in vs])
                
                in_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.in_edges()])
                out_eflows = numpy.sum([graph.ep.flow[e] for e in curr_node.out_edges()])
                maxdp = numpy.max([curr_dp, in_neighbor_dp_sum, out_neighbor_dp_sum, in_eflows, out_eflows])
                if curr_dp != maxdp:
                    graph.vp.dp[curr_node] = maxdp
                    no_changes = False
            
            no_changes = subfix(graph, curr_node, no_changes, consistent_node, inconsistent_node, simp_node_dict)
            # append all the in neighbors
            for n in curr_node.in_neighbors():
                if graph.vp.id[n] not in visited:
                    queue.append(n)        
        i += 1
        if i >= len(sorted_consistent_node):
            if break_point == consistent_count and (no_changes or cyclic):
                # no more changes
                break
            else:
                i = 0
                break_point = consistent_count

    print("Total consistent node after solid fix: ", consistent_count, len(simp_node_dict))
    print([int(n) for n in consistent_node])
    print([int(n) for n in inconsistent_node]) 
    print([int(n) for n in untracked_node])
    if DEBUG_MODE:
        for edge in simp_edge_dict.values():
            print_edge(graph, edge, "")
    return

###### main


    # pairwise_contig_concatenation(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, contig_node_dict, contig_concat_plans)
    # strain_dict = {}
    # iter_condition = set()
    # concat_contig_dict_iter = contig_dict.copy()
    # concat_strain_dict_iter = {}
    # contig pair-wise concatenation iteration
    # with iteration stopped when no further concatenation occurred
    # while True:
    #     print("Current iteration: ", level_no)
    #     graph_iter, simp_node_dict_iter, simp_edge_dict_iter = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no)))
    #     graph_simplification(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, {}, {}, args.min_cov/2)
    #     assign_edge_flow(graph_iter, simp_node_dict_iter, simp_edge_dict_iter)
    #     output_file = "{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no+1))
    #     # contig pair-wise concatenation
    #     concat_strain_dict_iter, concat_contig_dict_iter = contig_merge(graph_iter, simp_node_dict_iter, simp_edge_dict_iter, concat_contig_dict_iter, node_to_contig_dict, edge_to_contig_dict, output_file, args.min_cov, args.min_len, args.max_len, args.overlap)

    #     contig_dict_to_fasta(graph_iter, concat_contig_dict_iter, simp_node_dict_iter, args.overlap, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(level_no)))
    #     minimap_api(args.ref_file, "{0}L{1}_contigs.fasta".format(TEMP_DIR, str(level_no)), "{0}L{1}_contigs_to_strain.paf".format(TEMP_DIR, str(level_no)))

    #     strain_dict.update(concat_strain_dict_iter.copy())
    #     level_no += 1

    #     iter_compare = set(concat_contig_dict_iter.keys())
    #     if iter_condition == iter_compare:
    #         print("end of contig pair-wise concatenation iteration: ", level_no)
    #         break
    #     else:
    #         iter_condition = iter_compare
    
    # # further step to recover the rest of strain from the reduced graph

    # graph_red, simp_node_dict_red, simp_edge_dict_red = flipped_gfa_to_graph("{0}graph_L{1}.gfa".format(TEMP_DIR, str(level_no)))
    # # graph_simplification(graph_red, simp_node_dict_red, simp_edge_dict_red, node_to_contig_dict, edge_to_contig_dict, args.min_cov)
    # assign_edge_flow(graph_red, simp_node_dict_red, simp_edge_dict_red)
    
    # # extract the last mile paths from the graph
    # final_strain_dict = path_extraction(graph_red, simp_node_dict_red, simp_edge_dict_red, node_usage_dict, args.overlap, args.min_cov, args.min_len)   
    
    # graph_to_gfa(graph_red, simp_node_dict_red, simp_edge_dict_red, "{0}final_graph.gfa".format(TEMP_DIR))

    # strain_dict.update(final_strain_dict)

    # # print strain and minimap overlaps
    # contig_dict_to_fasta(pre_graph_v2, strain_dict, simp_node_dict_pre_v2, args.overlap, "{0}cand_strains.fasta".format(TEMP_DIR))
    # contig_dict_to_path(strain_dict, "{0}cand_strains.paths".format(TEMP_DIR))
    # minimap_api(args.ref_file, "{0}cand_strains.fasta".format(TEMP_DIR), "{0}cand_strain_to_strain.paf".format(TEMP_DIR))
    
    # # print out node usage stat
    # node_usage_pair = sorted(node_usage_dict.items(), key=lambda x: x[1][1] - x[1][0], reverse=True)
    # for id, used in node_usage_pair:
    #     print("id: {0} has been used {1} times".format(id, used))

def get_concat_plan(contig_dict: dict, max_len):
    contig_impossible_dict = {}
    all_contig_ids = contig_dict.keys()
    for no in contig_dict.keys():
        contig_impossible_dict[no] = set()
    for tail_cno, [tail_contig, tail_clen, _] in contig_dict.items():
        for head_cno, [head_contig, head_clen, _] in contig_dict.items():
            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len:
                    contig_impossible_dict[tail_cno].add(head_cno)
                if list(set(tail_contig) & set(head_contig)) != []:
                    contig_impossible_dict[tail_cno].add(head_cno)
    contig_concat_plans = {}
    for key, item in contig_impossible_dict.items():
        ps = all_contig_ids - item
        print("cno: ", key, " can concat with following: ", ps)
        contig_concat_plans[key] = ps
    return contig_concat_plans

def contig_merge(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, temp_contigs_dict: dict, node_to_contig_dict: dict, edge_to_contig_dict: dict, node_usage_dict: dict, output_file, min_cov, min_len, max_len, overlap):

    print("-----------------------contig pair-wise concatenaton--------------------")

    is_linear_strain = is_DAG(graph)
    if not is_linear_strain:
        print("graph is cyclic, cyclic strian may exist")
    else:
        print("graph is not cyclic, lienar strain may exist")

    # find pair-wise shortest path among contigs and construct the pairwise contig path matrix
    # TODO if sp found, do we also reduce its coverage?
    dist_matrix = {}
    concat_dicision = {}
    for tail_cno, [tail_contig, tail_clen, tail_ccov] in temp_contigs_dict.items():
        for head_cno, [head_contig, head_clen, head_ccov] in temp_contigs_dict.items():
            print("------------------------------------------------------")
            print("Tail Contig: ", tail_cno, " -> Head Contig: ", head_cno)
            if tail_cno != head_cno:
                if tail_clen + head_clen > max_len:
                    print("Contig ", tail_cno, "-", "Contig ", head_cno, " with len over maxlen")
                    continue
                intersect = list(set(tail_contig) & set(head_contig))
                if intersect != []:
                    print("Contig ", tail_cno, "-", "Contig ", head_cno, " Intersection with ", len(intersect), " nodes")
                    if DEBUG_MODE:
                        print(intersect)
                    continue
            else:
                # definitely no self-to-self path.
                if is_linear_strain:
                    print("No cyclic strain possible, skip")
                    continue
            print("start path construction")
            s_path, s_len = distance_search(graph, simp_node_dict, node_usage_dict, tail_contig[:-1], tail_contig[-1], head_contig[1:], head_contig[0], overlap)
            if s_path != None:
                s_path_ids = [graph.vp.id[v] for v in s_path]
                print("shortest path length: ", s_len, "path: ", [int(v) for v in s_path_ids])
                s_path_edge_flow = contig_flow(graph, simp_edge_dict,  s_path_ids)
                s_path_ccov = numpy.mean(s_path_edge_flow) if len(s_path_edge_flow) != 0 else 0
                dist_matrix[(tail_cno, head_cno)] = (s_path_ids, s_len, s_path_ccov)

                concat_ccov = min(head_ccov, tail_ccov, s_path_ccov) if s_path_ccov != 0.0 else min(head_ccov, tail_ccov)
                concat_len = head_clen + tail_clen - overlap if s_len == 0 else head_clen + s_len + tail_clen - overlap * 2
                print("coverage: head contig eflow: ", head_ccov, " s path eflow: ", s_path_ccov, " tail contig eflow: ", tail_ccov)
                print("potential concatenated length: ", concat_len)

                if concat_len > max_len:
                    print("exceed maximum strain length")
                else:
                    print("length satisfied, concat_ccov: ", concat_ccov)
                    concat_dicision[(tail_cno, head_cno)] = concat_ccov
            else:
                print("shortest path not found")
    print("------------------------------------------------------")
    
    # TODO re-evaluate the concat dicision, sort the concat dicision via concat_ccov, with reverse order or not, TBD
    sorted_pair = sorted(concat_dicision.items(), key=lambda x:x[1])
    concat_dicision_s = [t[0] for t in sorted_pair]

    #start concatenation
    # FIXME fix the concat dicision, test for used contig only concat n times.
    concat_strain_dict = {}
    concat_contig_dict = {}
    skip_key = set()
    used_contig = set()
    #FIXME not only used_contig, but also the nodes been used throughout the path.
    for (tail_cno, head_cno) in concat_dicision_s:
        print("------------------------------------------------------")
        print("Concatentate contigs: ", tail_cno, " <-> ", head_cno)
        [tail_contig, tail_clen, tail_ccov] = temp_contigs_dict[tail_cno]
        [head_contig, head_clen, head_ccov] = temp_contigs_dict[head_cno]

        if head_cno in used_contig:
            print("contig {0} be used from previous concatenation, skip".format(head_cno))
            continue
        if tail_cno in used_contig:
            print("contig {0} be used from previous concatenation, skip".format(tail_cno))
            continue

        if (head_cno, tail_cno) in skip_key:
            continue
        if (tail_cno, head_cno) in skip_key:
            continue
        skip_key.add((head_cno, tail_cno))
        skip_key.add((tail_cno, head_cno))

        is_linear = False
        if (head_cno, tail_cno) not in concat_dicision_s:
            print("no reverse concatenation exists, potential lienar strain")
            is_linear = True
        else: 
            print("access contig in both direction, potential cyclic strain")
            is_linear = False
        
        (s_path_ids_l, s_len_l, s_path_ccov_l) = dist_matrix[(tail_cno, head_cno)]

        concat_cno = tail_cno + "_" + head_cno
        overlap_count = 1 if s_len_l == 0 else 2
        if is_linear:
            #linear strain concatenation
            concat_clen = tail_clen + s_len_l + head_clen - overlap_count * overlap
            concat_ccov = min(v for v in [head_ccov, tail_ccov, s_path_ccov_l] if v > 0.0)
            concat_c = tail_contig + s_path_ids_l + head_contig
        else:
            #cyclic strain concatentation
            (s_path_ids_r, s_len_r, s_path_ccov_r) = dist_matrix[(head_cno, tail_cno)]
            overlap_count = overlap_count + 1 if s_len_r == 0 else overlap_count + 2

            concat_clen = tail_clen + s_len_l + head_clen + s_len_r - overlap_count * overlap
            concat_ccov = min(v for v in [head_ccov, tail_ccov, s_path_ccov_l, s_path_ccov_r] if v > 0.0)
            concat_c = tail_contig + s_path_ids_l + head_contig + s_path_ids_r
        
        # reduce the pre-defined contig's ccov by concat_ccov
        if head_ccov - concat_ccov < min_cov:
            used_contig.add(head_cno)
        else:
            temp_contigs_dict[head_cno][2] = temp_contigs_dict[head_cno][2] - concat_ccov

        if tail_ccov - concat_ccov < min_cov:
            used_contig.add(tail_cno)
        else:
            temp_contigs_dict[tail_cno][2] = temp_contigs_dict[tail_cno][2] - concat_ccov

        if concat_clen >= min_len:
            concat_strain_dict[concat_cno] = (concat_c, concat_clen, concat_ccov)
        else:
            concat_contig_dict[concat_cno] = [concat_c, concat_clen, concat_ccov]
    
    # re-append rest of contig
    for cno, [contig, clen, ccov] in temp_contigs_dict.items():
        if cno not in used_contig:
            if clen >= min_len:
                print("Full len contig directly store into strain dict", cno, clen, ccov)
                used_contig.add(cno)
                concat_strain_dict[cno] = (contig, clen, ccov)
            else:
                print("Re-appending contig to temp_contig_dict: ", cno, clen, ccov)
                concat_contig_dict[cno] = [contig, clen, ccov]
        else:
            print("contig {0} is used".format(cno))
    
    udpate_node_to_contig_dict(graph, node_to_contig_dict, simp_node_dict)
    update_edge_to_contig_dict(graph, edge_to_contig_dict, simp_edge_dict)
    # add concat contig cno to dict
    for cno, [c, _, _] in concat_contig_dict.items():
        for n in c:
            if n not in node_to_contig_dict:
                node_to_contig_dict[n] = [{cno},graph.vp.dp[simp_node_dict[n]], simp_node_dict[n]]
            else:
                node_to_contig_dict[n][0].add(cno)     
           
        for i in range(len(c)):
            c_i = c[i]
            c_i_1 = c[i+1] if (i < len(c) - 1) else None
            if c_i_1 != None:
                if (c_i, c_i_1) not in edge_to_contig_dict:
                    edge_to_contig_dict[(c_i, c_i_1)] = [{cno}, graph.ep.flow[simp_edge_dict[(c_i, c_i_1)]], simp_edge_dict[(c_i, c_i_1)]]
                else:
                    edge_to_contig_dict[(c_i, c_i_1)][0].add(cno)

    # update edge_to_contig_dict
    for (u,v), [cnos, _, _] in list(edge_to_contig_dict.items()):
        for cno in used_contig:
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove used contig cno {0} from edge: {1} cno list".format(cno, (u,v)))
        for cno in concat_strain_dict.keys():
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remov cand strain cno {0} from edge: {1} cno list".format(cno, (u,v)))

        if len(cnos) == 0:
            edge_to_contig_dict.pop((u,v))
        else:
            edge_to_contig_dict[(u,v)][0] = cnos

    # update node_to_contig_dict
    for no, [cnos, _, _] in list(node_to_contig_dict.items()):
        for cno in used_contig:
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove used contig cno {0} from node: {1} cno list".format(cno, no))
        for cno in concat_strain_dict.keys():
            if cno in cnos:
                cnos.remove(cno)
                if DEBUG_MODE:
                    print("Remove cand strain cno {0} from node: {1} cno list".format(cno, no))

        if len(cnos) == 0:
            node_to_contig_dict.pop(no)
        else:
            node_to_contig_dict[no][0] = cnos

    # reduce all the full length concated contigs as cand strain
    for cno, (c, clen, ccov) in concat_strain_dict.items():
        contig_reduction(graph, c, cno, clen, ccov, simp_node_dict, simp_edge_dict, min_cov)
        increment_node_usage_dict(node_usage_dict, c, ccov)

    for no, [cnos, dp, node] in node_to_contig_dict.items():
        if no not in simp_node_dict:
            simp_node_dict[no] = node
            sumcov = numpy.sum([concat_contig_dict[c][2] for c in cnos])
            graph.vp.dp[node] = sumcov
            graph.vp.color[node] = "black"
    # recover edge
    for (u,v), [cnos, flow, edge] in edge_to_contig_dict.items():
        if (u,v) not in simp_edge_dict:
            simp_edge_dict[(u,v)] = edge
            sumcov = numpy.sum([concat_contig_dict[c][2] for c in cnos])
            graph.ep.flow[edge] = sumcov
            graph.ep.color[edge] = "black"

    for cno, [contig, clen, ccov] in concat_contig_dict.items():
        print("------------------------------------------------------")
        print_contig(cno, clen, ccov, contig, "partial length residue contig found after concatenation")
        u = simp_node_dict[contig[0]]
        v = simp_node_dict[contig[-1]]
        if u.in_degree() == 0 and v.out_degree() == 0:
            print("no further concatentation exist")

    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, output_file)

    print("--------------------contig pair-wise concatenaton end--------------------")
    return concat_strain_dict, concat_contig_dict

def path_extraction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, node_usage_dict: dict, overlap, min_cov, min_len):
    """
    extract the last mile path from the graph, with support from residue contigs
    """
    print("--------------------path extraction--------------------")

    paths_per_group_dict = {}

    graph, groups = graph_grouping(graph, simp_node_dict)
    print("number of groups:", len(groups))
    final_strain_dict = {}
    
    for gno, group in groups.items():
        print("group number: ", gno, " member: ", [graph.vp.id[v] for v in group])
        if len(group) < 2:
            # single node group
            if len(group) == 1:
                node = group[0]
                pcov = graph.vp.dp[node]
                id = graph.vp.id[node]
                plen = len(graph.vp.seq[node])
                print("Single node Path: ", int(id), "path len: ", plen, "cov: ", pcov)
                # ratio%
                usage_ratio = round((node_usage_dict[id][0]/node_usage_dict[id][1])*100, 2)
                if pcov >= min_cov and usage_ratio > 0:
                    print("accept")
                    final_strain_dict[id] = ([id], plen, pcov)
            else:
                continue
        srcs = []
        sinks = []
        isolations = []
        middles = []
        paths = []
        # classify the nodes
        for u in group:
            no = graph.vp.id[u]
            if u.in_degree() == 0 and u.out_degree() == 0:
                isolations.append(no)
            elif u.in_degree() == 0:
                srcs.append(no)
            elif u.out_degree() == 0:
                sinks.append(no)
            else:
                middles.append(no)
        #FIXME potential issue: no src/sink node
        for src in srcs:
            for sink in sinks:
                print("Path extraction")
                p, plen = distance_search(graph, simp_node_dict, node_usage_dict, [], src, [], sink, overlap)
                if p != None:
                    s_path_ids = [graph.vp.id[v] for v in p]
                    s_path_ids.append(sink)
                    s_path_ids.insert(0, src)             
                    plen = plen + len(graph.vp.seq[simp_node_dict[src]]) + len(graph.vp.seq[simp_node_dict[sink]]) - 2 * overlap
                    paths.append((s_path_ids, plen))
        paths_per_group_dict[gno] = paths

    for gno, paths in paths_per_group_dict.items():
        print("Current group: ", gno)
        for p_ids, plen in paths:
            print("------------------------------------------------------")
            print("Path: ", [int(u) for u in p_ids])
            pcov = numpy.mean([graph.vp.dp[simp_node_dict[u]] for u in p_ids if u in simp_node_dict])
            pno = str(p_ids[0]) + "_" + str(p_ids[-1])
            print("path no: {0} path len: {1} path cov: {2}".format(pno, plen, pcov))

            if pcov >= min_cov:
                print("path accepted")
                final_strain_dict[pno] = (p_ids, plen, pcov)
                contig_reduction(graph, p_ids, pno, plen, pcov, simp_node_dict, simp_edge_dict, min_cov)
                increment_node_usage_dict(node_usage_dict, p_ids, pcov)

    print("------------------------------------------------------")
    print("--------------------path extraction end--------------------")
    return final_strain_dict

def contig_reduction(graph: Graph, contig, cno, clen, ccov, simp_node_dict: dict, simp_edge_dict: dict, min_cov):
    """
    reduce and update the graph by given contig
    """
    print("*---Contig reduction: ", cno, clen, ccov)
    next_node_index = 1
    for node in contig:
        if next_node_index >= len(contig):
            break
        adj_node = contig[next_node_index]

        u = simp_node_dict[node] if node in simp_node_dict else None
        v = simp_node_dict[adj_node] if adj_node in simp_node_dict else None
        e = simp_edge_dict[(node,adj_node)] if (node,adj_node) in simp_edge_dict else None 

        # edge may be eliminated from previous execution already
        if u == None or v == None or e == None:
            if e != None:
                graph.ep.flow[e] = 0
                graph.ep.color[e] = 'gray'
                simp_edge_dict.pop((node,adj_node))
                print_edge(graph, e, "edge removed 1st")
            continue
        if DEBUG_MODE:
            print_edge(graph, e, "current edge eval")
        # reduce the depth for involving edge
        if graph.ep.flow[e] - ccov <= min_cov:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            print_edge(graph, e, "edge removed 2nd")
        else:
            graph.ep.flow[e] = round(graph.ep.flow[e] - ccov, 2)

        # reduce the depth for involving node, gray color for forbiddened node
        if graph.vp.dp[u] - ccov <= min_cov:
            graph.vp.dp[u] = 0
            graph.vp.color[u] = 'gray'
            simp_node_dict.pop(node)
        else:
            graph.vp.dp[u] = round(graph.vp.dp[u] - ccov, 2)

        # last node in the contig, reduce its depth
        if next_node_index == len(contig) - 1:
            if graph.vp.dp[v] - ccov <= min_cov: 
                graph.vp.dp[v] = 0
                graph.vp.color[v] = 'gray'
                simp_node_dict.pop(adj_node)
            else:
                graph.vp.dp[v] = round(graph.vp.dp[v] - ccov, 2)
        
        # update edges
        if (graph.vp.color[u] == 'gray' or graph.vp.color[v] == 'gray') and (node,adj_node) in simp_edge_dict:
            graph.ep.flow[e] = 0
            graph.ep.color[e] = 'gray'
            simp_edge_dict.pop((node,adj_node))
            print_edge(graph, e, "edge removed 3rd")

        next_node_index = next_node_index + 1
    return

ss_path = None
cmp_path = None
def distance_search(graph: Graph, simp_node_dict: dict, contig_covered_node_ids: set, source, sink, overlap: int):
    """
    Compute minimal distance and its path between source node and sink node
    optimise the function with contig overlap check TODO
    FIXME re-define shortest path
    """
    def dfs_helper(graph: graph, u, sink, visited):
        """
        Return None if no shortest path is founded
        """
        if graph.vp.id[u] == sink:
            print("found sink")
            return [u]
        elif u.out_degree() == 0:
            return None
        else:
            global ss_path
            global cmp_path
            for v in u.out_neighbors():
                if not visited[v] and graph.vp.id[v] in simp_node_dict:
                    # print_vertex(graph, v, "curr node: ")
                    visited[v] = True
                    cmp_path = dfs_helper(graph, v, sink, visited)
                    if cmp_path != None:
                        # path lead to sink
                        cmp_path.insert(0, u)
                        cmp_len = path_len(graph, cmp_path, overlap)
                        if ss_path == None:
                            ss_path = cmp_path
                        else:
                            ss_len = path_len(graph, ss_path, overlap)
                            ss_path = ss_path if ss_len < cmp_len else cmp_path
                        # print(path_to_id_string(graph, ss_path, "curr sp: "))
                    visited[v] = False
        return ss_path

    print("source: ", source, "sink: ", sink)
    global ss_path
    ss_path = None
    global cmp_path
    cmp_path = None

    s_path = None
    s_len = 0
    print("start ssp")
    visited = {}
    for u in graph.vertices():
        visited[u] = False
    # avoid double contig path
    for s in contig_covered_node_ids:
        visited[simp_node_dict[s]] = True
    # avoid cycle
    visited[simp_node_dict[source]] = True
    visited[simp_node_dict[sink]] = False

    u = simp_node_dict[source]

    s_path = dfs_helper(graph, u, sink, visited)
    print(path_to_id_string(graph, s_path, "path: ") if s_path != None else "path: ")
    if s_path == None:
        print("Path not found")
    elif len(s_path) >= 2:
        s_path = s_path[1:-1]
        # compute directed path between ith tail to tth head (distance = len * #nodes in the path - (#nodes in the path - 1) * overlap)
        s_len = path_len(graph, s_path, overlap)
        print("shortest path found, len: ", s_len)
    else:
        s_path = None
        print("error path found")
    return s_path, s_len

def optim_path_linear(graph: Graph, simp_node_dict: dict, contig_dict: dict, srccno, tgtcno, src, tgt, cand_cov, bound_length, overlap, threshold, max_len, is_circular=False):
        """
        find the optimal path from src tgt, where intermediate nodes with cov ~ cand_cov would be prefered
        """
        print("Start path finding: {0} -> {1}, cand cov: {2}, bound: {3}".format(srccno, tgtcno, cand_cov, bound_length))
        path = [src]
        pathqueue = deque()
        pathqueue.append(path)
        rtn_paths = []

        visited = {}
        for node in simp_node_dict.values():
            visited[node] = False
        if not is_circular:
            for no in contig_dict[srccno][0][:-1]:
                visited[simp_node_dict[no]] = True
            for no in contig_dict[tgtcno][0][1:]:
                visited[simp_node_dict[no]] = True
        else:
            for no in contig_dict[srccno][0][1:-1]:
                visited[simp_node_dict[no]] = True

        while pathqueue:
            curr_path = pathqueue.popleft()
            if curr_path[-1] == tgt:
                print(path_to_id_string(graph, curr_path[1:-1], "curr path"))
                rtn_paths.append(curr_path[1:-1])
                continue
            if path_len(graph, curr_path[1:-1], overlap) > bound_length:
                print(path_to_id_string(graph, curr_path[1:-1], "curr path len overbound, skip impossible path: "))
                continue

            for next in curr_path[-1].out_neighbors():
                if next not in curr_path and not visited[next]:
                    # print_vertex(graph, next, "next: ")
                    split_path = curr_path[:]
                    split_path.append(next)

                    pathqueue.append(split_path)
        
        # rtn_paths = sorted(rtn_paths, key=lambda path: numpy.sum([pow(f-cand_cov,2) for f in contig_uflow(graph, simp_edge_dict, [graph.vp.id[n] for n in path])]))
        for p in rtn_paths:
            plen = path_len(graph, p, overlap)
            if not is_circular:
                total_len = get_concat_len(srccno, contig_dict[srccno][1], tgtcno, contig_dict[tgtcno][1], plen, overlap)
            else:
                total_len = contig_dict[srccno][1] + plen - 2*overlap
            num_similarity = 0
            for node in p:
                num_similarity = num_similarity + 1 if abs(graph.vp.udp[node] - cand_cov) < threshold else num_similarity
            print("total len: ", total_len, " similarity: ", num_similarity, path_to_id_string(graph, p, "--->path: "))
            if total_len <= max_len:
                return rtn_paths, p, total_len, num_similarity
        return None, None, None, None

def contig_pairwise_concatenation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
cliq_graph: Graph, cliq_node_dict: dict, cliq_edge_dict: dict, sp_path_dict: dict, 
min_cov, min_len, max_len, overlap, threshold, tempdir):
    """
    1. concat the self cycle and reduce some self cycles if self len < minlen
    2. concat the most confident contig pair in order
    3. if cycle exist within the graph, break the cycle by removing the most weak edge.
    4. gradually concat from source nodes until all the contig be concated.
    """
    # helper functions
    def graph_reduction_c(cand_path, cand_cov):
        for i in range(len(cand_path) - 1):
            u = cand_path[i]
            v = cand_path[i + 1]
            e = graph.edge(u, v)
            graph.vp.udp[u] -= cand_cov
            graph.vp.udp[v] -= cand_cov
            graph.ep.flow[e] -= cand_cov

    def contig_pair_reduction(cno1, cno2, cand_cov, cand_path, cand_len, cno_mapping: dict):
        """
        1. reduce the graph and cliq graph via founded path and cand_cov, 
        for cliq graph, then duplicate/remove
        the cno1/2 node and merge into a single node with cand_cov, keep 
        all the connections other than 1-2
        """
        # original graph udp/edge flow reduction
        graph_reduction_c(cand_path, cand_cov)

        # cliq graph reduction
        cnode1 = cliq_node_dict[cno1]
        cliq_graph.vp.ccov[cnode1] -= cand_cov
        cnode2 = cliq_node_dict[cno2]
        cliq_graph.vp.ccov[cnode2] -= cand_cov

        print("L1 node cov after deduction: ", cliq_graph.vp.ccov[cnode1])
        print("L2 node cov after deduction: ", cliq_graph.vp.ccov[cnode2])

        if cliq_graph.vp.ccov[cnode1] <= threshold and cliq_graph.vp.ccov[cnode2] <= threshold:
            print("both node be used up, merge L1 to L2, keep the L1 in edges and L2 out edges only")
            cno_merged = cno1 + "->" + cno2

            if cno1 in cno_mapping:
                if cno1 in cno_mapping[cno1]:
                    cno_mapping[cno1].remove(cno1)
                cno_mapping[cno1].add(cno_merged)
            else:
                cno_mapping[cno1] = {cno_merged}

            if cno2 in cno_mapping:
                if cno2 in cno_mapping[cno2]:
                    cno_mapping[cno2].remove(cno2)
                cno_mapping[cno2].add(cno_merged)
            else:
                cno_mapping[cno2] = {cno_merged}

            if cno_merged not in cno_mapping:
                cno_mapping[cno_merged] = {cno_merged}

            prev1 = contig_dict.pop(cno1)
            prev2 = contig_dict.pop(cno2)
            contig_dict[cno_merged] = [prev1[0]+[graph.vp.id[v] for v in cand_path]+prev2[0], cand_len, cand_cov]
            
            cnode_merged = cliq_graph_add_node(cliq_graph, cliq_node_dict, 
            cno_merged, cand_len, cand_cov, 
            cno_merged + ":" + str(cand_len) + ":" + str(cand_cov))

            # removed the used up contig node
            cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno1, cnode1)
            cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno2, cnode2)

            # removed the related L1 edges
            for edge1out in cnode1.out_edges():
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict,
                cliq_graph.vp.cno[edge1out.source()], cliq_graph.vp.cno[edge1out.target()], edge1out)

            # replace the related L1 in edges
            for edge1in in cnode1.in_edges():
                if cliq_graph.ep.color[edge1in] != 'black':
                    continue
                src1 = edge1in.source()
                src1cno = cliq_graph.vp.cno[src1]
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                src1cno, cliq_graph.vp.cno[edge1in.target()], edge1in)
                
                # do not add self cycle edges
                if src1cno != cno_merged:
                    cliq_graph_add_edge(cliq_graph, cliq_edge_dict, src1cno, src1, 
                    cno_merged, cnode_merged, cliq_graph.ep.slen[edge1in], cliq_graph.ep.text[edge1in])
            
            # removed the related L2 in edges
            for edge2in in cnode2.in_edges():
                if cliq_graph.ep.color[edge2in] != 'black':
                    continue
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                cliq_graph.vp.cno[edge2in.source()], cliq_graph.vp.cno[edge2in.target()], edge2in)
            
            # replace the related L2 out edges
            for edge2out in cnode2.out_edges():
                if cliq_graph.ep.color[edge2out] != 'black':
                    continue
                tgt2 = edge2out.target()
                tgt2cno = cliq_graph.vp.cno[tgt2]
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                cliq_graph.vp.cno[edge2out.source()], tgt2cno, edge2out)
                if cno_merged != tgt2cno:
                    cliq_graph_add_edge(cliq_graph, cliq_edge_dict, cno_merged, cnode_merged,
                    tgt2cno, tgt2, cliq_graph.ep.slen[edge2out], cliq_graph.ep.text[edge2out])

        elif cliq_graph.vp.ccov[cnode1] <= threshold:
            print("L1 node be used up, split L2 node, and merge L1 to L2")
            # split L2
            cno_dup = cno1 + "->" + cno2

            if cno1 in cno_mapping:
                if cno1 in cno_mapping[cno1]:
                    cno_mapping[cno1].remove(cno1)
                cno_mapping[cno1].add(cno_dup)
            else:
                cno_mapping[cno1] = {cno_dup}
                
            if cno2 in cno_mapping:
                cno_mapping[cno2].add(cno_dup)
            else:
                cno_mapping[cno2] = {cno_dup}

            if cno_dup not in cno_mapping:
                cno_mapping[cno_dup] = {cno_dup}

            prev1 = contig_dict.pop(cno1)
            contig_dict[cno_dup] = [prev1[0]+[graph.vp.id[v] for v in cand_path]+contig_dict[cno2][0], cand_len, cand_cov]
            contig_dict[cno2][2] -= cand_cov

            cnode_dup = cliq_graph_add_node(cliq_graph, cliq_node_dict, 
            cno_dup, cand_len, cand_cov, 
            cno_dup + ":" + str(cand_len) + ":" + str(cand_cov))

            # remove the L1 node
            cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno1, cnode1)

            # remove the related L1 edges
            for edge1 in cnode1.all_edges():
                if cliq_graph.ep.color[edge1] != 'black':
                    continue
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                cliq_graph.vp.cno[edge1.source()], cliq_graph.vp.cno[edge1.target()], edge1)
            
            # append the related L2 out edges to dup node
            for edge2out in cnode2.out_edges():
                if cliq_graph.ep.color[edge2out] != 'black':
                    continue
                tgt2 = edge2out.target()
                tgt2cno = cliq_graph.vp.cno[tgt2]
                if tgt2cno != cno_dup:
                    cliq_graph_add_edge(cliq_graph, cliq_edge_dict, cno_dup, cnode_dup,
                    tgt2cno, tgt2, cliq_graph.ep.slen[edge2out], cliq_graph.ep.text[edge2out])

        elif cliq_graph.vp.ccov[cnode2] <= threshold:
            print("L2 node be used up, split L1 node, and merge L1 to L2")  
            # split L1
            cno_dup = cno1 + "->" + cno2

            if cno1 in cno_mapping:
                cno_mapping[cno1].add(cno_dup)
            else:
                cno_mapping[cno1] = {cno_dup}

            if cno2 in cno_mapping:
                if cno2 in cno_mapping[cno2]:
                    cno_mapping[cno2].remove(cno2)
                cno_mapping[cno2].add(cno_dup)
            else:
                cno_mapping[cno2] = {cno_dup}

            if cno_dup not in cno_mapping:
                cno_mapping[cno_dup] = {cno_dup}

            prev2 = contig_dict.pop(cno2)
            contig_dict[cno_dup] = [contig_dict[cno1][0]+[graph.vp.id[v] for v in cand_path]+prev2[0], cand_len, cand_cov]
            contig_dict[cno1][2] -= cand_cov

            cnode_dup = cliq_graph_add_node(cliq_graph, cliq_node_dict, 
            cno_dup, cand_len, cand_cov, 
            cno_dup + ":" + str(cand_len) + ":" + str(cand_cov))

            cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno2, cnode2)

            # remove the related L2 inedges
            for edge2 in cnode2.in_edges():
                if cliq_graph.ep.color[edge2] != 'black':
                    continue
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                cliq_graph.vp.cno[edge2.source()], cliq_graph.vp.cno[edge2.target()], edge2)
            
            # replace the related L2 out edges to dup node
            for edge2out in cnode2.out_edges():
                if cliq_graph.ep.color[edge2out] != 'black':
                    continue
                tgt2 = edge2out.target()
                tgt2cno = cliq_graph.vp.cno[tgt2]
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
                cliq_graph.vp.cno[edge2out.source()], tgt2cno, edge2out)
                if tgt2cno != cno_dup:
                    cliq_graph_add_edge(cliq_graph, cliq_edge_dict, cno_dup, cnode_dup,
                    tgt2cno, tgt2, cliq_graph.ep.slen[edge2out], cliq_graph.ep.text[edge2out])
        else:
            print("error: L1 node: ", cliq_graph.vp.ccov[cnode1], ", L2 node: ", cliq_graph.vp.ccov[cnode2])

    def buffer_concatenation(concat_buffer: list, cno_mapping: dict):
        for (cno1, cno2, cov, delta) in concat_buffer:
            print("------------------------------------------------------------------------")
            print("-->Before mapping: {0}: {1} - {2}: {3}".format(cno1, cno_mapping[cno1], cno2, cno_mapping[cno2]))
            for pcno in list(cno_mapping[cno1]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno1].remove(pcno)
                    print("pcno removed: ", pcno)
            for pcno in list(cno_mapping[cno2]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno2].remove(pcno)  
                    print("pcno removed: ", pcno)

            # cno1m = sorted(cno_mapping[cno1], key=lambda x: abs(cov - cliq_graph.vp.ccov[cliq_node_dict[x]]))[0]
            # cno2m = sorted(cno_mapping[cno2], key=lambda x: abs(cov - cliq_graph.vp.ccov[cliq_node_dict[x]]))[0]
            
            pairs = []
            for cno1m in cno_mapping[cno1]: 
                for cno2m in cno_mapping[cno2]:
                    if (cno1m, cno2m) in cliq_edge_dict:
                        pairs.append((cno1m, cno2m))
            if not pairs:
                print("contig has been used from previous step: {0} {1}".format(cno1, cno2))
                continue

            if cno1 in contig_usage:
                contig_usage[cno1] += 1
            if cno2 in contig_usage:
                contig_usage[cno2] += 1

            cno1m, cno2m = min(pairs, key=lambda p: 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[0]]] - cov, 2) + 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[1]]] - cov, 2))
            
            print("-->PAIR UP {0} - {1}, cov: {2}, diff: {3}".format(cno1m, cno2m, cov, delta))
            src = simp_node_dict[contig_dict[cno1m][0][-1]]
            tgt = simp_node_dict[contig_dict[cno2m][0][0]]
            cand_path, plen = dijkstra_sp(graph, simp_node_dict, src, tgt, cov, threshold, overlap)
            cand_len = get_concat_len(cno1m, contig_dict[cno1m][1], cno2m, contig_dict[cno2m][1], plen, overlap)

            contig_pair_reduction(cno1m, cno2m, cov, cand_path, cand_len, cno_mapping)

        return
    ###################################################################################################
    # udp, node depth with related contig cov be deducted
    graph.vp.udp = graph.new_vertex_property("double")
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    for no, node in simp_node_dict.items():
        graph.vp.udp[node] = graph.vp.dp[node]
        if no in node_to_contig_dict:
            for cno in node_to_contig_dict[no]:
                graph.vp.udp[node] -= contig_dict[cno][2]
    contig_usage = {}
    for cno in contig_dict.keys():
        contig_usage[cno] = 0

    # retrieve all the self cycle first
    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                if cliq_graph.vp.clen[contig_node] < min_len:
                    # remove self cycle edge with self cycle + outer connection feature when clen < minlen
                    cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_graph.edge(contig_node, contig_node))
                    print("remove self edge+outer connection {0} -> {0} with node length < minlen".format(cno))
            else:
                print("Circular PATH: ", cno)
                cov = cliq_graph.vp.ccov[contig_node]
                src = simp_node_dict[contig_dict[cno][0][-1]]
                tgt = simp_node_dict[contig_dict[cno][0][0]]
                cand_path, plen = dijkstra_sp(graph, simp_node_dict, src, tgt, cov, threshold, overlap)
                cand_len = get_concat_len(cno, contig_dict[cno][1], cno, contig_dict[cno][1], plen, overlap)

                if cand_path != None and cand_len != None:
                    graph_reduction_c(cand_path, cov)
                    print(path_to_id_string(graph, [simp_node_dict[cno] for cno in contig_dict[cno][0]] + cand_path, "cov: {0}".format(cov)))
                    contig_dict[cno][0].extend([graph.vp.id[n] for n in cand_path])
                    contig_dict[cno][1] = cand_len
                    contig_dict[cno][2] = cov
                    cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_edge_dict[(cno, cno)])
                else:
                    print("Path not found, error")

    cno_mapping = {}
    for id, node in cliq_node_dict.items():
        cno_mapping[id] = {id}

    concat_buffer = []
    for (cno1, cno2) in cliq_edge_dict.keys():
        print("---------------------------------------------------------------")
        delta = abs(cliq_graph.vp.ccov[cliq_node_dict[cno1]] - cliq_graph.vp.ccov[cliq_node_dict[cno2]])
        cov = min(cliq_graph.vp.ccov[cliq_node_dict[cno1]], cliq_graph.vp.ccov[cliq_node_dict[cno2]])
        print("{0} - {1}, cov: {2}, diff: {3}".format(cno1, cno2, cov, delta))
        if delta < threshold:
            concat_buffer.append((cno1, cno2, cov, delta))
    
    concat_buffer = sorted(concat_buffer, key=lambda tuple: tuple[3])
    print("all most confident sorted concats are: ", concat_buffer)
    buffer_concatenation(concat_buffer, cno_mapping)

    for cno, contig in cliq_node_dict.items():
        contig, clen, ccov = contig_dict[cno]
        print_contig(cno, clen, ccov, contig)
    for (u, v), e in cliq_edge_dict.items():
        print("EDGE: ", u, v)

    if not graph_is_DAG(cliq_graph, cliq_node_dict):
        print("------>graph contains cycles, cyclic concatenation processing")
        cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
        cycles = all_circuits(cliq_graph, True)
        for cyc in cycles:
            print("cycle: ", [cliq_graph.vp.cno[n] for n in cyc])
            skip_cycle = False
            all_edges = {}
            for i in range(len(cyc)):
                u = cyc[i]
                v = cyc[(i+1) % len(cyc)]
                e = cliq_graph.edge(u, v)
                if graph.ep.color[e] != 'black':
                    skip_cycle = True
                    break
                all_edges[(cliq_graph.vp.cno[u],cliq_graph.vp.cno[v])] = cliq_graph.edge(u, v)
            if skip_cycle:
                continue
            removed_edge_tuple = min(all_edges.items(), key=lambda e_tuple: cliq_graph.ep.slen[e_tuple[1]])
            print("Removed edge: ", removed_edge_tuple[0][0], "->", removed_edge_tuple[0][1])
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, removed_edge_tuple[0][0], removed_edge_tuple[0][1], removed_edge_tuple[1])

    print("--------------Start graduate concatentation------------------")
    for cno, cnos in cno_mapping.items():
        print("cno: ", cno, "maps to: ", cnos)
    # concat 
    # process all the non self-cycle contig until no isolated node.
    # direct pair two adjacent contig only if coverage difference is within pairwise threshold
    while True:
        # clean up the cno_mapping
        cno_mapping = {}
        for id in cliq_node_dict.keys():
            cno_mapping[id] = {id}
        # all the src contig
        L1_contigs = {}
        for no, contig_node in list(cliq_node_dict.items()):
            if cliq_graph.vp.color[contig_node] != 'black':
                continue
            ind = len([e for e in contig_node.in_edges() if cliq_graph.ep.color[e] == 'black'])
            outd = len([e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black'])
            if ind == 0 and outd != 0:
                # src node
                L1_contigs[no] = contig_node

        L2_contigs = {}
        st_pairs = {}
        for no, contig_node in list(L1_contigs.items()):
            out_edges = [e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black']
            for out_e in out_edges:
                out_contig_node = out_e.target()
                if cliq_graph.vp.color[out_contig_node] != 'black':
                    continue
                outid = cliq_graph.vp.cno[out_contig_node]
                if not outid in L2_contigs:
                    L2_contigs[outid] = out_contig_node
                st_pairs[(no, outid)] = (contig_node, out_contig_node)
        
        if not L1_contigs or not L2_contigs:
            print("no more linear concatentation left, break")
            break

        print("L1: ", [no for no in L1_contigs.keys()])
        print("L2: ", [no for no in L2_contigs.keys()])
        
        # detect the most confident contig pair first
        concat_buffer = []
        for (cno1, cno2), (node1, node2) in st_pairs.items():
            delta = abs(cliq_graph.vp.ccov[node1] - cliq_graph.vp.ccov[node2])
            concat_buffer.append((cno1, cno2, min(cliq_graph.vp.ccov[node1], cliq_graph.vp.ccov[node2]), delta))
        
        concat_buffer = sorted(concat_buffer, key=lambda tuple: tuple[3])
        print("all most confident sorted concats are: ", concat_buffer)
        buffer_concatenation(concat_buffer, cno_mapping)

    print("--------------> contig after concatentations...")
    for cno, contig_node in list(cliq_node_dict.items()):
        print("--------------------------------------------------------")
        if cno in contig_usage:
            count = contig_usage[cno]
            print("cno: ", cno, "original usage count: ", count)
            if count > 0:
                print("original contig has been used, remove it to reduce the potential duplication")
                cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno, contig_node)
                contig_dict.pop(cno)
                continue
        contig, clen, ccov = contig_dict[cno]
        cflows = contig_flow(graph, simp_edge_dict, contig)
        print_contig(cno, clen, ccov, contig)
        if len(cflows) > 0:
            print("mean: ", numpy.mean(cflows), 
            " median: ", numpy.median(cflows), 
            " min: ", numpy.min(cflows))

        if simp_node_dict[contig[0]].in_degree() == 0 and simp_node_dict[contig[-1]].out_degree() == 0:
            print("not extensible")
        elif simp_node_dict[contig[0]] in list(simp_node_dict[contig[-1]].out_neighbors()):
            print("self cycle, not entensible")
        else:
            print("extensible")

    for (u, v), e in cliq_edge_dict.items():
        print("potential uncaught error: EDGE: ", u, v)

    for no, node in simp_node_dict.items():
        print("Node: {0}, full: {1}, left: {2}, usage: {3}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), round((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2))
    
    # simplify the graph
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), tempdir, "cliq_graphL2.png")

    return contig_dict


    {
    # graph transitive reduction FIXME select the option the reduce.
    # print("total edges: ", len(cliq_edge_dict))
    # transitive_graph_reduction(graph, simp_node_dict, simp_edge_dict, cliq_graph, contig_dict, cliq_node_dict, cliq_edge_dict, sp_path_dict)
    # print("total edges after reduction: ", len([e for e in cliq_edge_dict.values() if cliq_graph.ep.color[e] == 'black']))
    # cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    # edge reduction
    # ensure every contig have at most 1 in connection and 1 out connection
    # if graph_is_DAG(graph, simp_node_dict):
    #     max_edges = []
    #     for no, node in cliq_node_dict.items():
    #         print("current contig evaluating: ", no)
    #         # in connection
    #         ine = [e for e in node.in_edges() if cliq_graph.ep.color[e] == 'black']
    #         if len(ine) > 1:
    #             ine = sorted(ine, 
    #             key=lambda e: similarity_e(e, cliq_graph), reverse=True)
    #             max_sim = similarity_e(ine[0], cliq_graph)
    #             print("MAX IN EDGE: {0}->{1}, plen: {2}".
    #             format(cliq_graph.vp.cno[ine[0].source()], cliq_graph.vp.cno[ine[0].target()], cliq_graph.ep.slen[ine[0]]))
    #             print("MAX SIM: ", max_sim)
    #             for e in ine:
    #                 cs = similarity_e(e, cliq_graph)
    #                 if cs < max_sim:
    #                     cliq_graph.ep.color[e] = 'gray'
    #                     print("drop edge {0} -> {1}, sim: {2}, slen: {3}".format(cliq_graph.vp.cno[e.source()], cliq_graph.vp.cno[e.target()], cs, cliq_graph.ep.slen[e]))
    #                 else:
    #                     max_edges.append(e)
    #         # out connection
    #         oute = [e for e in node.out_edges() if cliq_graph.ep.color[e] == 'black']
    #         if len(oute) > 1:
    #             oute = sorted(oute, 
    #             key=lambda e: similarity_e(e, cliq_graph), reverse=True)
    #             max_sim = similarity_e(oute[0], cliq_graph)
    #             print("MAX OUT EDGE: {0}->{1}, plen: {2}".format(cliq_graph.vp.cno[oute[0].source()], cliq_graph.vp.cno[oute[0].target()], cliq_graph.ep.slen[oute[0]]))
    #             print("MAX SIM: ", max_sim)
    #             for e in oute:
    #                 cs = similarity_e(e, cliq_graph)
    #                 if cs < max_sim:
    #                     cliq_graph.ep.color[e] = 'gray'    
    #                     print("drop edge {0} -> {1}, sim: {2}, slen: {3}".format(cliq_graph.vp.cno[e.source()], cliq_graph.vp.cno[e.target()], cs, cliq_graph.ep.slen[e]))     
    #                 else:
    #                     max_edges.append(e)

    #     # recover node in/out connection if edges be removed lead to no connection
    #     sorted_sim_edges_f = sorted([e for e in cliq_graph.edges() if cliq_graph.ep.color[e] != 'black'], key=lambda e: similarity_e(e, cliq_graph), reverse=True)
    #     for fe in sorted_sim_edges_f:
    #         u = fe.source()
    #         v = fe.target()
    #         # print("current check edge {0} -> {1}".format(cliq_graph.vp.cno[u], cliq_graph.vp.cno[v]))
    #         uoute = [e for e in u.out_edges() if cliq_graph.ep.color[e] == 'black']
    #         vine = [e for e in v.in_edges() if cliq_graph.ep.color[e] == 'black']
    #         if len(uoute) == 0 and len(vine) == 0:
    #             cliq_graph.ep.color[fe] = 'black'
    #             print("reappend edge {0} -> {1}, {2}".format(cliq_graph.vp.cno[u], cliq_graph.vp.cno[v], similarity_e(fe, cliq_graph)))
    }


# def contig_variation_tree(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold):
#     contig, clen, ccov = contig_dict['6']
#     # head_queue = [simp_node_dict[contig[0]]]
#     # head_nodes = []
#     # visited = {}
#     # for v in graph.vertices():
#     #     visited[v] = False
#     # while head_queue:
#     #     curr = head_queue.pop()
#     #     head_nodes.append(curr)
#     #     visited[curr] = True
#     #     for prev in curr.in_neighbors():
#     #         if not visited[prev]:
#     #             head_queue.append(prev)
#     mine, minflow = min_flow_edge(graph, simp_edge_dict, contig)
#     if mine == None:
#         minflow = graph.vp.dp[simp_node_dict[contig[-1]]]
#     qe = [contig[-1], minflow]
#     tail_queue = [qe]
#     tail_nodes = set()
#     visited = {}
#     outd = 0
#     for v in graph.vertices():
#         visited[graph.vp.id[v]] = False
#     while tail_queue:
#         currid, cflow = tail_queue.pop()
#         currnode = simp_node_dict[currid]
        
#         outflow = 0
#         for oute in currnode.out_edges():
#             next = oute.target()
#             if not visited[graph.vp.id[next]]:
#                 tail_nodes.add(graph.vp.id[next])
#                 outd += 1
#                 outflow = cflow * (graph.ep.flow[oute] / graph.vp.dp[currnode]) 
#                 tail_queue.append([graph.vp.id[next], outflow])

#         graph.vp.dp[currnode] -= cflow
#         if graph.vp.dp[currnode] < threshold:
#             visited[currid] = True

#     # print(path_to_id_string(graph, head_nodes, "c6 head tree: "))
#     print(path_to_id_string(graph, list(tail_nodes), "c6 tail tree: "))
#     print("outd: ", outd)


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


    # TODO
    # print("--------------------------------------GRAPH TRIVIAL SPLIT----------------------------------")
    # prev_ids = list(simp_node_dict5.keys())
    # trivial_split_count, id_mapping = graph_split_trivial(graph5, simp_node_dict5, simp_edge_dict5)
    # graph_to_gfa(graph5, simp_node_dict5, simp_edge_dict5, "{0}fcbsdt_graph_L6.gfa".format(TEMP_DIR))
    # graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}fcbsdt_graph_L6.gfa".format(TEMP_DIR))
    # assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)
    # contig_dict_resol(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, id_mapping, prev_ids, args.overlap)

    # simp_path_dict = simple_paths_to_dict(graph5, simp_node_dict5, simp_edge_dict5, args.overlap)
    # simp_path_compactification(graph5, simp_node_dict5, simp_edge_dict5, simp_path_dict, contig_dict, args.overlap)
    # contig_dict = contig_dict_simp(contig_dict, THRESHOLD)

    # graph5, simp_node_dict5, simp_edge_dict5, contig_dict = reindexing(graph5, simp_node_dict5, simp_edge_dict5, contig_dict)

    # contig_dict_to_path(contig_dict, "{0}post_contig.paths".format(TEMP_DIR))

    # graph_to_gfa(graph5, simp_node_dict5, simp_edge_dict5, "{0}sfcbsdt_graph_L7.gfa".format(TEMP_DIR))
    # graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}sfcbsdt_graph_L7.gfa".format(TEMP_DIR))
    # assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)

    # num_split, branch_id_mapping = graph_splitting(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, THRESHOLD, strict_mode=False, oddBranch=True)

    # graph_to_gfa(graph5, simp_node_dict5, simp_edge_dict5, "{0}ssfcbsdt_graph_L8.gfa".format(TEMP_DIR))



    # while bubbles
    # gen_bubble_detection(graph, simp_node_dict, simp_edge_dict, overlap)
    branches, bubbles = minimal_bubble_detection(graph)
    for (a, b), bubble in bubbles.items():
        print("------------------------------------------------")
        print("Bubble {0} <-> {1}".format(a, b))
        bin = {}
        strains = {}
        prev_strain_loc = {}
        for var in bubble:
            varid = graph.vp.id[var]
            bin[varid] = graph.vp.dp[var]
            if varid in node_to_contig_dict:
                for cno in node_to_contig_dict[varid]:
                    if cno not in strains:
                        strains[cno] = contig_dict[cno][2]
                    if cno not in prev_strain_loc:
                        prev_strain_loc[cno] = varid
                print("    bin: {0}, bin size: {1}, curr usage: {2}%\n    strain: {3}\n".format(varid, graph.vp.dp[var], pre_usages[varid],
                [(cno, contig_dict[cno][2]) for cno in node_to_contig_dict[varid]] if varid in node_to_contig_dict else ""))
            else:
                print("    bin: {0}, bin size: {1}, curr usage: {2}%\n    strain: None\n".format(varid, graph.vp.dp[var], pre_usages[varid]))
        if len(strains) > 0:
            bin_assignment, strain_placement = LPSolveBubbleLocal(bin, strains, threshold)
            for bin, contained_strain in bin_assignment.items():
                for cno in contained_strain:
                    curr_contig = contig_dict[cno][0]
                    rep_contig = contig_replacement(curr_contig, prev_strain_loc[cno], bin)
                    print("cno: {0} {1} --> {2}".format(cno, list_to_string(curr_contig), list_to_string(rep_contig)))
                    contig_dict[cno][0] = rep_contig
    print("------------------------------------------------")
    # handle all the local minimal bubble, do local swap
    # post-processing, based on swapped information, split the bubble edges evenly
    # store graph, combine simple path, store graph, map contig node
    # end loop

    # # stat on current usage
    # # udp, node depth with related contig cov be deducted
    # graph.vp.udp = graph.new_vertex_property("double")
    # for no, node in simp_node_dict.items():
    #     graph.vp.udp[node] = graph.vp.dp[node]
    # print("------------------> all the extended contigs: ")
    # for cno, [contig, clen, ccov] in list(contig_dict.items()):
    #     print_contig(cno, clen, ccov, contig, "-----> extended contig ")
    #     call = [simp_node_dict[n] for n in contig]
    #     if global_src in simp_node_dict[contig[0]].in_neighbors():
    #         call.insert(0, global_src)
    #     if global_sink in simp_node_dict[contig[-1]].out_neighbors():
    #         call.append(global_sink)
    #     graph_reduction_c(graph, call, ccov, threshold)

    # for v in graph.vertices():
    #     if graph.vp.udp[v] < threshold:
    #         graph.vp.color[v] = 'gray'
    # for e in graph.edges():
    #     if graph.ep.flow[e] < threshold:
    #         graph.ep.color[e] = 'gray'
    # partial_used_nodes = []
    # free_nodes = []
    # overused_nodes = []
    # pre_usages = {}
    # node_to_contig_dict, _ = contig_map_node(contig_dict)
    # for no, node in simp_node_dict.items():
    #     print("----------------------------------------------------")
    #     if no == 'global_src' or no == 'global_sink':
    #         continue
    #     ratio = round(((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2)
    #     print("Node: {0}, full: {1}, left: {2}, usage: {3}, color: {4}, CNOs: {5}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), ratio, graph.vp.color[node], node_to_contig_dict[no] if no in node_to_contig_dict else ""))
    #     if ratio < 100 and graph.vp.udp[node] > threshold:
    #         partial_used_nodes.append(no)
    #     if ratio <= 0:
    #         free_nodes.append(no)
    #     if ratio > 100 and graph.vp.udp[node] < -threshold:
    #         overused_nodes.append(no)
    #     pre_usages[no] = ratio
    # overall_pre_usage = numpy.mean(list(pre_usages.values()))
    # print("Free nodes: ", list_to_string(free_nodes))
    # print("Partial used nodes: ", list_to_string(partial_used_nodes))
    # print("Over used nodes: ", list_to_string(overused_nodes))
    # print("Usage: ", overall_pre_usage)

    # # analytics part
    # global TEMP_DIR
    # seaborn.set_theme(style="darkgrid")
    # plt.figure(figsize=(128,64))
    # df = pandas.DataFrame(
    #     {'Id': pre_usages.keys(),
    #     'pre': pre_usages.values(),
    #     })
    # ax = seaborn.barplot(x='Id', y='pre', data=df)
    # for container in ax.containers:
    #     ax.bar_label(container)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    # plt.title('Bar plot of Node Usage')
    # plt.savefig("{0}barplot_usage_pre.png".format(TEMP_DIR))
    # # end stat

    # if noncyc_nodes != None and simple_paths != None:
    #     print("rebalance linear subgraph now..")
    #     graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}nc_graph_L2p.gfa".format(TEMP_DIR))
    #     coverage_rebalance(graphnc, simp_node_dictnc, simp_edge_dictnc)
    #     print("Done, start coverage merge")

    #     for no, node in simp_node_dictnc.items():
    #         cnode = simp_node_dict1[no]
    #         merge_dp = graphnc.vp.dp[node] + graph1.vp.dp[cnode]
    #         graph1.vp.dp[cnode] = merge_dp   
    # else:
    #     print("no linear subgraph available..")


    # print("------------------------------NODE PARTITION-----------------------------------")
    # # store the no-cycle nodes in nc_graph_L3p.gfa
    # noncyc_nodes = None
    # simple_paths = None
    # # graphtool is_DAG() may not work if the graph is not connected as several parts
    # if not graph_is_DAG(graph1, simp_node_dict1):
    #     noncyc_nodes, simple_paths = node_partition(graph1, simp_node_dict1, simp_edge_dict1, TEMP_DIR)

def a(contig_dict, args, THRESHOLD):
    while True:
        # trivial branch split
        prev_ids = list(simp_node_dicta.keys())
        trivial_split_count, id_mapping = graph_split_trivial(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

        graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)

        contig_dict_remapping(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, id_mapping, prev_ids, args.overlap)

        # non-trivial branch split
        num_split, branch_id_mapping = graph_splitting(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, THRESHOLD)
        # clean up graph
        graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphc, simp_node_dictc, simp_edge_dictc)

        contig_cov_fix(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, branch_id_mapping)

        # simp path compactification
        simp_path_dict = simple_paths_to_dict(graphc, simp_node_dictc, simp_edge_dictc, args.overlap)
        simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, simp_path_dict, contig_dict, args.overlap)

        graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)

        # contig_dict_simp(contig_dict, THRESHOLD)
        trim_contig_dict(graph5, simp_node_dict5, contig_dict, args.overlap)

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            grapha = graph5
            simp_node_dicta = simp_node_dict5
            simp_edge_dicta = simp_edge_dict5
            iterCount = chr(ord(iterCount) + 1)
        else:
            break



#backup
    # while True:
    #     # branch split
    #     num_split, branch_id_mapping = graph_splitting(grapha, simp_node_dicta, simp_edge_dicta, contig_dict, THRESHOLD)
    #     # clean up graph
    #     graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
    #     graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
    #     assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)
        
    #     contig_cov_fix(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, branch_id_mapping)

    #     simp_path_dict = simple_paths_to_dict(graphb, simp_node_dictb, simp_edge_dictb, args.overlap)
    #     simp_path_compactification(graphb, simp_node_dictb, simp_edge_dictb, simp_path_dict, contig_dict, args.overlap)
        
    #     contig_dict = contig_dict_simp(contig_dict, THRESHOLD)

    #     graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}1s.gfa".format(TEMP_DIR, iterCount))
    #     graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1s.gfa".format(TEMP_DIR, iterCount))
    #     assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)
        
    #     # normal split FIXME
    #     prev_ids = list(simp_node_dictb.keys())
    #     trivial_split_count, id_mapping = graph_split_trivial(graphb, simp_node_dictb, simp_edge_dictb, contig_dict)

    #     graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
    #     graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
    #     assign_edge_flow(graphc, simp_node_dictc, simp_edge_dictc)

    #     contig_dict_resol(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, id_mapping, prev_ids, args.overlap)

    #     simp_path_dict = simple_paths_to_dict(graphc, simp_node_dictc, simp_edge_dictc, args.overlap)
    #     simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, simp_path_dict, contig_dict, args.overlap)

    #     graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
    #     graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
    #     assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)

    #     contig_dict = contig_dict_simp(contig_dict, THRESHOLD)
    #     contig_dict = trim_contig_dict(graph5, simp_node_dict5, contig_dict, args.overlap)
        
    #     if num_split != 0 or trivial_split_count != 0:
    #         total_removed_branch_nt += num_split
    #         total_removed_branch_t += trivial_split_count
    #         grapha = graph5
    #         simp_node_dicta = simp_node_dict5
    #         simp_edge_dicta = simp_edge_dict5
    #         iterCount = chr(ord(iterCount) + 1)
    #     else:
    #         break

    # print("-------------------------------REAPPEND CYC EDGES-----------------------------------")
    # for (i, o) in removed_edges:
    #     if i in simp_node_dict3 and o in simp_node_dict3:
    #         # node still exist
    #         graph_add_edge(graph3, simp_edge_dict3, simp_node_dict3[i], i, simp_node_dict3[o], o, args.overlap, 0)
    # assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)


def contig_branch_splitting(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold):
    """
    pre-request: ccov stored in contig dict should already be the mincov among the contig path.
    """
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)

    splitters_per_contig = {}
    for cno, [contig, clen, ccov] in contig_dict.items():
        print_contig(cno, clen, ccov, contig, "currnet split contig")
        splitters = []
        trivials = []
        # within range [1, len_contig - 1], guarantee to split the branch if any
        currIndex = 1
        while currIndex < len(contig) - 1:
            nid = contig[currIndex]
            nnode = simp_node_dict[nid]
            nnode = simp_node_dict[nid]
            if graph.vp.color[nnode] == 'black':
                if is_splitter_branch(nnode):
                    print("Found a splitter branch: ", nid)
                    splitters.append(nid)
                    # split the branch
                    prev_nid = contig[currIndex - 1]
                    next_nid = contig[currIndex + 1]
                    assert (prev_nid, nid) in simp_edge_dict and (nid, next_nid) in simp_edge_dict
                    prev_sp_contigs = node_to_contig_dict[prev_nid]
                    curr_sp_contigs = node_to_contig_dict[nid]
                    next_sp_contigs = node_to_contig_dict[next_nid]
                    involved_contig = (prev_sp_contigs.intersection(curr_sp_contigs)).intersection(next_sp_contigs)
                    left_edge = simp_edge_dict[(prev_nid, nid)]
                    right_edge = simp_edge_dict[(nid, next_nid)]
                    print("prev: {0}, edge flow: {1}, curr: {2}, edge flow: {3}, next: {4}".
                    format(prev_nid, graph.ep.flow[left_edge], nid, graph.ep.flow[right_edge], next_nid))
                    print("involved contig: ", involved_contig)

                    # subid = "X" + nid + "X"
                    # subdp = numpy.mean([graph.ep.flow[ie], graph.ep.flow[oe]])
                    # sub_node = graph_add_vertex(graph, simp_node_dict, subid, subdp, graph.vp.seq[node], graph.vp.kc[node])

                    # if no not in no_mapping:
                    #     no_mapping[no] = []
                    # no_mapping[no].append(str(subid))

                    # graph.vp.dp[node] -= subdp

                    # graph_remove_edge(graph, simp_edge_dict, prev_no, no)
                    # graph_remove_edge(graph, simp_edge_dict, no, next_no)
                    # graph_add_edge(graph, simp_edge_dict, prev_node, prev_no, sub_node, subid, graph.ep.overlap[ie],
                    # graph.ep.flow[ie])
                    # graph_add_edge(graph, simp_edge_dict, sub_node, subid, next_node, next_no, graph.ep.overlap[oe],
                    # graph.ep.flow[oe])
                    # map all involved contig
                elif is_trivial_branch(nnode):
                    print("Found a trivial branch: ", nid)
                    trivials.append(nid)

            currIndex += 1        
        splitters_per_contig[cno] = [splitters,trivials]
    
    for cno, [splitters, trivials] in splitters_per_contig.items():
        print("Contig: ", cno, list_to_string(contig_dict[cno][0], "c: "))
        print(list_to_string(splitters, "splitters: "))
        print(list_to_string(trivials, "trivials: "))
    return

def split_contig(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, threshold):
    """
    split contig out
    """
    dup_record = {}
    for no in simp_node_dict.keys():
        dup_record[no] = 0
    for cno in contig_dict.keys():
        contig, _, ccov = contig_dict[cno]
        contig = list(contig)
        print("--->split contig: ", cno)
        for i in range(1,len(contig) - 1):
            no = contig_dict[cno][0][i]
            prev_no = contig_dict[cno][0][i-1]
            prev_old_no = contig[i-1]
            ie = simp_edge_dict[(prev_old_no, no)]
            graph.ep.flow[ie] -= ccov
            # if graph.ep.flow[ie] <= 0:
            #     graph_remove_edge(graph, simp_edge_dict, graph.vp.id[ie.source()], graph.vp.id[ie.target()])

            old_vertex = simp_node_dict[no]

            graph.vp.dp[old_vertex] -= ccov
            new_vertex = graph_add_vertex(graph, simp_node_dict, str(no) + "X" + str(dup_record[no]), 
                ccov, graph.vp.seq[old_vertex], graph.vp.kc[old_vertex])
            graph_add_edge(graph, simp_edge_dict, simp_node_dict[prev_no], 
                prev_no, new_vertex, graph.vp.id[new_vertex], graph.ep.overlap[ie], ccov)
            if i == len(contig) - 2:
                next_no = contig_dict[cno][0][i+1]
                oe = simp_edge_dict[(no, next_no)]
                graph.ep.flow[oe] -= ccov
                # if graph.ep.flow[oe] <= 0:
                #     graph_remove_edge(graph, simp_edge_dict, graph.vp.id[oe.source()], graph.vp.id[oe.target()])
                graph_add_edge(graph, simp_edge_dict, new_vertex, graph.vp.id[new_vertex], 
                    simp_node_dict[next_no], next_no, graph.ep.overlap[ie], ccov) 
            
            # if graph.vp.dp[old_vertex] <= 0:
            #     graph_remove_vertex(graph, simp_node_dict, no)
            contig_dict[cno][0][i] = graph.vp.id[new_vertex]
            dup_record[no] += 1
    return


def extract_cand_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, overlap, threshold):
    # top sort the nodes
    lb = min([graph.ep.flow[e] for e in graph.edges()])
    ub = lb*2
    print("[{0},{1})".format(lb, ub))
    gs, gt = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, True)
    # graph.ep.isMin = graph.new_edge_property("bool")
    min_vset = set()
    min_eset = {}
    for e in graph.edges():
        if graph.ep.flow[e] > lb and graph.ep.flow[e] <= ub:
            # graph.ep.isMin[e] = True
            min_vset.add(e.source())
            min_vset.add(e.target())
            min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e
    print(path_to_id_string(graph, list(min_vset)))
    global TEMP_DIR
    x = [graph.ep.flow[e] for e in min_eset.values()]
    plt.figure(figsize=(128,64))
    # bins=round(graph.ep.flow[max(graph.edges(), key=lambda e: graph.ep.flow[e])]/threshold)
    ax = plt.hist(x)
    plt.axvline(x=lb, color='r')
    plt.axvline(x=ub, color='r')
    plt.title("min_eset")
    plt.savefig("{0}{1}".format(TEMP_DIR, "min_eset.png"))

    sorted_eset = sorted(min_eset.values(), key=lambda e: graph.ep.flow[e])
    for e in sorted_eset:
        print_edge(graph, e)
    # for (au, av), e_tail in min_eset.items():
    #     for (bu, bv), e_head in min_eset.items():

    return




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
        trivial_split_count, id_mapping = graph_split_trivial(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

        graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        graphb, simp_node_dictb, simp_edge_dictb = flipped_gfa_to_graph("{0}split_graph_L{1}1.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphb, simp_node_dictb, simp_edge_dictb)

        contig_dict_resol(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, id_mapping, prev_ids, args.overlap)

        # non-trivial branch split
        num_split, branch_id_mapping = graph_splitting(graphb, simp_node_dictb, simp_edge_dictb, contig_dict, THRESHOLD)

        graph_to_gfa(graphb, simp_node_dictb, simp_edge_dictb, "{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        graphc, simp_node_dictc, simp_edge_dictc = flipped_gfa_to_graph("{0}split_graph_L{1}2.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(graphc, simp_node_dictc, simp_edge_dictc)

        contig_cov_fix(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, branch_id_mapping)

        # simp path compactification
        simp_path_compactification(graphc, simp_node_dictc, simp_edge_dictc, contig_dict, args.overlap)

        graph_to_gfa(graphc, simp_node_dictc, simp_edge_dictc, "{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        grapha, simp_node_dicta, simp_edge_dicta = flipped_gfa_to_graph("{0}split_graph_L{1}3.gfa".format(TEMP_DIR, iterCount))
        assign_edge_flow(grapha, simp_node_dicta, simp_edge_dicta)

        contig_dict_simp(grapha, simp_edge_dicta, contig_dict, THRESHOLD)
        trim_contig_dict(grapha, simp_node_dicta, contig_dict, args.overlap)

        if num_split != 0 or trivial_split_count != 0:
            total_removed_branch_nt += num_split
            total_removed_branch_t += trivial_split_count
            iterCount = chr(ord(iterCount) + 1)
        else:
            break
    print("Total non-trivial branches removed: ", total_removed_branch_nt, " total trivial branches removed: ", total_removed_branch_t)

    contig_dict_to_path(contig_dict, "{0}pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(grapha, contig_dict, simp_node_dicta, args.overlap, "{0}pre_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}pre_contigs.fasta".format(TEMP_DIR), "{0}pre_contigs_to_strain.paf".format(TEMP_DIR))
    
    grapha, simp_node_dicta, simp_edge_dicta, contig_dict = reindexing(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

    contig_dict_to_path(contig_dict, "{0}post_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(grapha, contig_dict, simp_node_dicta, args.overlap, "{0}post_contigs.fasta".format(TEMP_DIR))
    minimap_api(args.ref_file, "{0}post_contigs.fasta".format(TEMP_DIR), "{0}post_contigs_to_strain.paf".format(TEMP_DIR))

    graph_to_gfa(grapha, simp_node_dicta, simp_edge_dicta, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
    graph5, simp_node_dict5, simp_edge_dict5 = flipped_gfa_to_graph("{0}rbsdt_graph_L5.gfa".format(TEMP_DIR))
    assign_edge_flow(graph5, simp_node_dict5, simp_edge_dict5)

    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}rbsdt_graph_L5.gfa".format(TEMP_DIR), True, "{0}node_to_ref_red.paf".format(TEMP_DIR), "{0}temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    gfa_to_fasta("{0}rbsdt_graph_L5.gfa".format(TEMP_DIR), "{0}graph_L5_node.fasta".format(TEMP_DIR))



    print("-------------------------------ACYCLIC TRANSLATION-----------------------------------")
    removed_edges = cyclic_to_dag(graph1, simp_node_dict1, contig_dict, args.overlap)
    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}t_graph_L1_dag.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1_dag.gfa".format(TEMP_DIR))

def branch_split_formal(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, threshold):
    def is_partial(e):
        return abs(cap(e)) < threshold
    def cap(e):
        return graph.ep.flow[e] - graph.ep.usage[e]
    def closest_edge(e, es):
        return min(es, key=lambda e2: abs(cap(e2) - cap(e)))

    graph.ep.usage = graph.new_edge_property("int", val=0)
    dup_records = {}
    split_plan = {}
    for no in simp_node_dict.keys():
        dup_records[no] = 0
        split_plan[no] = []

    for no, node in simp_node_dict.items():
        Ein = list(filter(is_partial, node.in_edges()))
        Eout = list(filter(is_partial, node.out_edges()))
        if len(Ein) <= 1 and len(Eout) <= 1:
            # elim non-branch node
            continue
        print("Current eval branch: ", no)
        while len(Ein) > 0 and len(Eout) > 0:
            min_edge = min(set(Ein).union(set(Eout)), key=cap)
            min_cap = cap(min_edge)
            graph.ep.usage[min_edge] += min_cap
            c_edge = None
            if min_edge in Ein:
                # min in edge
                c_edge = closest_edge(min_edge, Eout)
                graph.ep.usage[c_edge] += min_cap
                if not is_partial(min_edge):
                    Ein.remove(min_edge)
                if not is_partial(c_edge):
                    Eout.remove(c_edge)
                split_plan[no].append([min_edge, c_edge, min_cap])
            else:
                # min out edge
                c_edge = closest_edge(min_edge, Ein)
                graph.ep.usage[c_edge] += min_cap
                if not is_partial(min_edge):
                    Eout.remove(min_edge)
                if not is_partial(c_edge):
                    Ein.remove(c_edge)
                split_plan[no].append([c_edge, min_edge, min_cap])
        print("Split plan for ", no)
        for i, o, u in split_plan[no]:
            print("------------------------------>")
            print_edge(graph, i, "left")
            print_edge(graph, o, "right")
            print("Capacity: ", u)
    return 0


            # for i, path in enumerate(paths):
            #     dupcno = cno + "^" + str(i)
            #     if dupcno in contig_dict:
            #         print("dup cno: ", dupcno, " already exist, error")
            #     else:
            #         eflow = contig_flow(graph, simp_edge_dict, path)
            #         subcov = 0
            #         if len(eflow) < 1:
            #             # only one vertex
            #             subcov = graph.vp.dp[simp_node_dict[path[0]]]
            #         else:
            #             subcov = min(eflow)
            #         contig_dict[dupcno] = [path, path_len(graph, [simp_node_dict[no] for no in path], overlap), subcov]
            #     print("duplicated mapped contig: ", dupcno, list_to_string(path))

        # contig_list = []
        # k = -1
        # for i in range(len(contig)):
        #     currnodes = [path[i] for path in paths]
        #     if all(e==currnodes[0] for e in currnodes):
        #         if k == i - 1:
        #             # just concat to the last sublist from contig_list
        #             if len(contig_list) == 0:
        #                 contig_list = [[currnodes[0]]]
        #             else:
        #                 contig_list[-1].append(currnodes[0])
        #             k += 1
        #         else:
        #             contig_list.append([currnodes[0]])
        #             k = i
        # contig_dict.pop(cno)
        # if len(contig_list) < 1:
        #     print("no mapping for the current contig: whole contig is ambiguous mapping", cno)
        #     for i, path in enumerate(paths):
        #         dupcno = cno + "^" + str(i)
        #         if dupcno in contig_dict:
        #             print("dup cno: ", dupcno, " already exist, error")
        #         else:
        #             eflow = contig_flow(graph, simp_edge_dict, path)
        #             subcov = 0
        #             if len(eflow) < 1:
        #                 # only one vertex
        #                 subcov = graph.vp.dp[simp_node_dict[path[0]]]
        #             else:
        #                 subcov = min(eflow)
        #             contig_dict[dupcno] = [path, path_len(graph, [simp_node_dict[no] for no in path], overlap), subcov]
        #         print("duplicated mapped contig: ", dupcno, list_to_string(path))
        # elif len(contig_list) == 1:
        #     print("single mapping")
        #     print(list_to_string(contig_list[0]))
        #     contig_dict[cno] = [contig_list[0], path_len(graph, [simp_node_dict[no] for no in contig_list[0]], overlap), ccov]
        # else:
        #     for i, subcontig in enumerate(contig_list):
        #         subcno = cno + "^" + str(i)
        #         if subcno in contig_dict:
        #             print("sub cno: ", subcno, " already exist, error")
        #         else:
        #             contig_dict[subcno] = [subcontig, path_len(graph, [simp_node_dict[no] for no in subcontig], overlap), ccov]
        #         print("sub mapped contig: ", subcno, list_to_string(subcontig))


    
   # print("-------------------------------CONTIG CLIQUE GRAPH BUILD-----------------------------------")
    
    # draw_edgeflow(graph5, simp_edge_dict5, TEMP_DIR, 'Bar plot of Edge flow', 'barplot_edge_flow_sp.png')

    # cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict, adj_matrix, cno_to_index = contig_clique_graph_build(graph5, simp_node_dict5, simp_edge_dict5,
    #                                                 contig_dict, args.max_len, THRESHOLD, args.overlap)
    # # draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), TEMP_DIR, "cliq_graphL1.png")
    # cliq_graph, cliq_node_dict, cliq_edge_dict = clique_graph_clean(cliq_graph, cliq_node_dict, cliq_edge_dict, adj_matrix, cno_to_index, THRESHOLD)
    # # draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), TEMP_DIR, "cliq_graphL2.png")
    # print("-------------------------------CONTIG PAIRWISE CONCATENATION-----------------------------------")
    # concat_contig_dict = contig_pairwise_concatenation(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, 
    # cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict, 
    # args.min_cov, args.min_len, args.max_len, args.overlap, THRESHOLD, TEMP_DIR)

    # # partial result
    # contig_dict_to_fasta(graph5, concat_contig_dict, simp_node_dict5, args.overlap, "{0}concat_contig.fasta".format(TEMP_DIR))
    # contig_dict_to_path(concat_contig_dict, "{0}concat_contig.paths".format(TEMP_DIR))
    # minimap_api(args.ref_file, "{0}concat_contig.fasta".format(TEMP_DIR), "{0}concat_contig_to_strain.paf".format(TEMP_DIR))

    # print("-------------------------------------STRAIN EXTENSION--------------------------------------")
    # extended_contig_dict = strain_extension(graph5, simp_node_dict5, simp_edge_dict5, concat_contig_dict, args.min_cov, args.max_len, THRESHOLD, args.overlap)

    # print("-------------------------------------LOCAL OPTIMISATION--------------------------------------")
    # local_search_optimisation(graph5, simp_node_dict5, simp_edge_dict5, extended_contig_dict, 
    # args.overlap, THRESHOLD)
    # # print("----------------------------------FINAL STRAIN EXTRACTION----------------------------------")
    # # final_strain_extraction(graph5, simp_node_dict5, simp_edge_dict5, extended_contig_dict, pre_usages, THRESHOLD, args.overlap)
    # # print("-------------------------------LOCAL SEARCH OPTIMISATION-----------------------------------")
    # # local_search_optimisation(graph5, simp_node_dict5, simp_edge_dict5, concat_contig_dict, args.min_cov, args.min_len, args.max_len, args.overlap, THRESHOLD)
    # # stat
    # extended_contig_dict = trim_contig_dict(graph5, simp_node_dict5,  extended_contig_dict, args.overlap)


    
    # grapha, simp_node_dicta, simp_edge_dicta, contig_dict = reindexing(grapha, simp_node_dicta, simp_edge_dicta, contig_dict)

    # contig_dict_to_path(contig_dict, "{0}post_contigs.paths".format(TEMP_DIR))
    # contig_dict_to_fasta(grapha, contig_dict, simp_node_dicta, args.overlap, "{0}post_contigs.fasta".format(TEMP_DIR))
    # minimap_api(args.ref_file, "{0}post_contigs.fasta".format(TEMP_DIR), "{0}post_contigs_to_strain.paf".format(TEMP_DIR))


    # print("-------------------------------ACYCLIC TRANSLATION-----------------------------------")
    # #FIXME
    # # removed_edges = cyclic_to_dag(graph1, simp_node_dict1, contig_dict, args.overlap)
    # graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}t_graph_L1_dag.gfa".format(TEMP_DIR))
    # graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}t_graph_L1_dag.gfa".format(TEMP_DIR))



# cliq graph version
def allowed_concat_init(graph: Graph, contig_dict: dict, simp_node_dict: dict, max_len, threshold, overlap):
    """
    Decide whether any contig pair should be connected
    1. two contig cannot concat if:
        1. two contig have intermediate parallel sharing nodes
        2. no path exist between two contig
    """
    self_concat_off = graph_is_DAG(graph, simp_node_dict)
    impossible_concat_dict = {}
    sp_path_dict = {}

    for no in contig_dict.keys():
        impossible_concat_dict[no] = set()
    for tail_cno, [tail_contig, tail_clen, tail_ccov] in contig_dict.items():
        for head_cno, [head_contig, head_clen, head_ccov] in contig_dict.items():
            print("---------> tail cno: ", tail_cno, " vs head cno: ", head_cno)
            print("tail cov: {0} - head cov: {1}".format(tail_ccov, head_ccov))

            tail_node = simp_node_dict[contig_dict[tail_cno][0][-1]]
            head_node = simp_node_dict[contig_dict[head_cno][0][0]]

            # in-out degree reachable filter
            if tail_node.out_degree() == 0 or head_node.in_degree() == 0:
                impossible_concat_dict[tail_cno].add(head_cno)
            
            # contig intersection filter
            intersects = None
            cend = None
            status = None
            if tail_cno != head_cno:
                isParallel, intersects, cend, status = check_contig_intersection(tail_cno, tail_contig, head_cno, head_contig)
                if isParallel:
                    impossible_concat_dict[tail_cno].add(head_cno)
                if status != None:                 
                    if status != 'n' and status != 'o':
                        if status != 'b':
                            # end-to-end overlap case, forward direction
                            sp_path_dict[(tail_cno, head_cno)] = None
                            continue
                        else:
                            impossible_concat_dict[tail_cno].add(head_cno)
            else:
                if self_concat_off:
                    impossible_concat_dict[tail_cno].add(head_cno)
            if head_cno not in impossible_concat_dict[tail_cno]:
                total_len = 0
                if intersects != None and cend != None:
                    # approx overlap length
                    intersect_len = sum([len(graph.vp.seq[simp_node_dict[id]]) for id in intersects]) - overlap * (len(intersects) + cend)
                    total_len = -intersect_len
                sp, plen, pmark = dijkstra_sp(graph, tail_node, head_node, min(tail_ccov, head_ccov), threshold, overlap)
                if sp != None:
                    if head_cno == tail_cno:
                        if plen == 0:
                            total_len += head_clen - overlap
                        else:
                            total_len += head_clen + plen - 2*overlap
                        print("total cyclic shortest length: ", total_len)
                    else:
                        if plen == 0:
                            total_len += head_clen + tail_clen - overlap
                        else:
                            total_len += head_clen + tail_clen + plen - 2*overlap
                        print("total linear shortest length: ", total_len)
                    if total_len >= max_len:
                        print("even sp exceed the maxlen: ", max_len)
                        impossible_concat_dict[tail_cno].add(head_cno)
                    else:
                        print("SP length within upper bound max len")
                        sp_path_dict[(tail_cno, head_cno)] = (sp, plen, pmark)
                else:
                    print("SP not found, impossible")
                    impossible_concat_dict[tail_cno].add(head_cno)

    all_contig_ids = contig_dict.keys()
    contig_concat_plans = {}
    for key, item in impossible_concat_dict.items():
        contig_concat_plans[key] = all_contig_ids - item
    return contig_concat_plans, sp_path_dict

def contig_clique_graph_build(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, max_len, threshold, overlap):
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
    
    for cno, (contig, clen, ccov) in contig_dict.items():
        contig_node = cliq_graph.add_vertex()
        cliq_graph.vp.cno[contig_node] = cno
        cliq_graph.vp.clen[contig_node] = clen
        cliq_graph.vp.ccov[contig_node] = ccov
        cliq_graph.vp.color[contig_node] = 'black'
        cliq_graph.vp.text[contig_node] = cno + ":" + str(clen) + ":" + str(round(ccov, 2))
        cliq_node_dict[cno] = contig_node

    concat_plan, sp_path_dict = allowed_concat_init(graph, contig_dict, simp_node_dict, max_len, threshold, overlap)
    
    cno2cno_adjMtx_csv = []
    cno2cno_adjMtx_csv.append([cno+"/"+str(round(contig_dict[cno][2])) for cno in concat_plan.keys()])
    cno2cno_adjMtx_csv[0].insert(0, "X")
    cno2cno_adjMtx = []

    cno_to_index = {}
    cno_to_index_csv = {}
    for i, cno in enumerate(cno2cno_adjMtx_csv[0]):
        if i > 0:
            cno_to_index_csv[cno.split("/")[0]] = i
            cno_to_index[cno.split("/")[0]] = i - 1

    for tail_cno, head_cnos in concat_plan.items():
        print("------> tail cno: ", tail_cno, " can concat with following head cnos: ", head_cnos)
        curr_row_csv = [" " for _ in concat_plan.keys()]
        curr_row_csv.insert(0, tail_cno+"/"+str(round(contig_dict[tail_cno][2])))
        curr_row = [[sys.maxsize,'X'] for _ in concat_plan.keys()]

        src_node = cliq_node_dict[tail_cno]
        for head_cno in head_cnos:
            plen = int(sp_path_dict[(tail_cno, head_cno)][1]) if sp_path_dict[(tail_cno, head_cno)] != None else -1
            abs_cdif = abs(contig_dict[tail_cno][2]-contig_dict[head_cno][2])
            curr_row_csv[cno_to_index_csv[head_cno]] = str(round(abs_cdif)) + ":" + str(plen)
            curr_row[cno_to_index[head_cno]] = [abs_cdif, 'W']
            tgt_node = cliq_node_dict[head_cno]
            contig_edge = cliq_graph.add_edge(src_node, tgt_node)
            cliq_graph.ep.slen[contig_edge] = plen
            cliq_graph.ep.color[contig_edge] = 'black'
            cliq_graph.ep.text[contig_edge] = str(cliq_graph.ep.slen[contig_edge])
            cliq_edge_dict[(tail_cno, head_cno)] = contig_edge
        cno2cno_adjMtx.append(curr_row)
        cno2cno_adjMtx_csv.append(curr_row_csv)
    return cliq_graph, cliq_node_dict, cliq_edge_dict, sp_path_dict, cno2cno_adjMtx, cno_to_index

def clique_graph_clean(cliq_graph: Graph, cliq_node_dict: dict, cliq_edge_dict: dict, 
adj_matrix, cno_to_index: dict, threshold):
    """
    adj matrix, elem in 5 colors: X, gray(invalid), white(has connection), red(candidate connection), blue(fixed connection)
    FIXME not all contig should be concated, check the diff < threshold, overlap contig pair have highest priority
    """

    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                remove_self_cycle = False
                if not remove_self_cycle:
                    for inn in contig_node.in_neighbors():
                        if inn != contig_node and contig_node in (inn.in_neighbors()):
                            print("remove sc set 01, ")
                            remove_self_cycle = True
                            break
                if not remove_self_cycle:
                    for onn in contig_node.out_neighbors():
                        if onn != contig_node and contig_node in (onn.out_neighbors()):
                            print("remove sc set 02, ")
                            remove_self_cycle = True
                            break    
                if remove_self_cycle:
                    # remove self cycle edge with self cycle + outer connection feature
                    adj_matrix[cno_to_index[cno]][cno_to_index[cno]][1] = 'G'
                    print("remove self edge+outer connection {0} -> {0}".format(cno))
                else:
                    print("Keep self cycle since no feedback from neighbor ", cno)
    index_to_cno = {}
    for cno, i in cno_to_index.items():
        index_to_cno[i] = cno

    for rowId, row in enumerate(adj_matrix):
        for colId, [flow, color] in enumerate(row):
            if color == 'W':
                if flow >= threshold:
                    print("manual forbidden concat: ", index_to_cno[rowId], index_to_cno[colId], flow)
                    adj_matrix[rowId][colId][1] = 'G'
                else:
                    print("Potential concat: ", index_to_cno[rowId], index_to_cno[colId], flow)
    has_changes = True
    dim = len(adj_matrix)
    # iterate the adjacent matrix in rowwise and columnwise
    while has_changes:
        print("looping")
        has_changes = False
        # row wise
        for rowId in range(dim):
            row = get_row(adj_matrix, rowId)
            hasFixed = False
            colId = None
            minFlow = sys.maxsize
            for cid, [flow, color] in enumerate(row):
                if color == 'B':
                    # current row fixed
                    hasFixed = True
                    break
                elif color == 'W':
                    if flow < minFlow:
                        colId = cid
                        minFlow = flow
                elif color == 'X' or color == 'G':
                    # invalid, skip
                    None
                else:
                    print("Error: ", rowId, index_to_cno[rowId], cid, index_to_cno[cid])
                    assert color != 'R'
            if not hasFixed:
                if colId != None:
                    # found a minimum block to assign
                    adj_matrix[rowId][colId][1] = 'R'
                    has_changes = True
                else:
                    # no more spot for the row
                    None
        for colId in range(dim):
            col = get_col(adj_matrix, colId)
            # hasFixed
            rowId = None
            minFlow = sys.maxsize
            # if only one red+blue among the column, then assign it to blue/nochange, otherwise select the cand flow red -> blue, 
            # (also challenge the blue) and recolor the other red to False
            cands = []
            for rid, [flow, color] in enumerate(col):
                if color == 'R' or color == 'B':
                    cands.append((rid, [flow, color]))
            if len(cands) == 0:
                # relax
                None
            elif len(cands) == 1:
                rid, [flow, color] = cands[0]
                if color == 'R':
                    adj_matrix[rid][colId] = [flow, 'B']
                    has_changes = True
            else:
                mrid, [mflow, _] = min(cands, key=lambda p: p[1][0])
                for rid, [flow, color] in cands:
                    adj_matrix[rid][colId] = [flow, 'G']
                adj_matrix[mrid][colId] = [mflow, 'B']
                has_changes = True
    
    for rowId, row in enumerate(adj_matrix):
        print([c for _, c in row])
        for colId, [abs_diff, color] in enumerate(row):
            if color != 'X':
                e = cliq_edge_dict[(index_to_cno[rowId], index_to_cno[colId])]
                if color != 'B':
                    cliq_graph.ep.color[e] = 'gray'
                else:
                    if abs_diff < threshold:
                        cliq_graph.ep.color[e] = 'black'
                    else:
                        cliq_graph.ep.color[e] = 'gray'
                    print("Decided Edge: {0} -> {1}, asb_diff: {2}".format(index_to_cno[rowId], index_to_cno[colId], abs_diff))
    
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    return cliq_graph, cliq_node_dict, cliq_edge_dict



def contig_pairwise_concatenation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
cliq_graph: Graph, cliq_node_dict: dict, cliq_edge_dict: dict, sp_path_dict: dict, 
min_cov, min_len, max_len, overlap, threshold, tempdir):
    """
    1. concat the self cycle and reduce some self cycles if self len < minlen
    2. concat the most confident contig pair in order
    3. if cycle exist within the graph, break the cycle by removing the most weak edge.
    4. gradually concat from source nodes until all the contig be concated.
    """
    # helper functions
    def contig_pair_reduction(cno1, cno2, cand_cov, cand_path, cand_len, cno_mapping: dict):
        """
        1. reduce the graph and cliq graph via founded path and cand_cov, 
        for cliq graph, then duplicate/remove
        the cno1/2 node and merge into a single node with cand_cov, keep 
        all the connections other than 1-2
        """
        # cliq graph reduction
        cnode1 = cliq_node_dict[cno1]
        cliq_graph.vp.ccov[cnode1] -= cand_cov
        cnode2 = cliq_node_dict[cno2]
        cliq_graph.vp.ccov[cnode2] -= cand_cov

        print("L1 node cov after deduction: ", cliq_graph.vp.ccov[cnode1])
        print("L2 node cov after deduction: ", cliq_graph.vp.ccov[cnode2])

        print("merge L1 to L2, keep the L1 in edges and L2 out edges only")
        cno_merged = cno1 + "->" + cno2

        if cno1 in cno_mapping:
            if cno1 in cno_mapping[cno1]:
                cno_mapping[cno1].remove(cno1)
            cno_mapping[cno1].add(cno_merged)
        else:
            cno_mapping[cno1] = {cno_merged}

        if cno2 in cno_mapping:
            if cno2 in cno_mapping[cno2]:
                cno_mapping[cno2].remove(cno2)
            cno_mapping[cno2].add(cno_merged)
        else:
            cno_mapping[cno2] = {cno_merged}

        if cno_merged not in cno_mapping:
            cno_mapping[cno_merged] = {cno_merged}

        prev1 = contig_dict.pop(cno1)
        prev2 = contig_dict.pop(cno2)
        contig_dict[cno_merged] = [prev1[0]+[graph.vp.id[v] for v in cand_path]+prev2[0], cand_len, cand_cov]
        
        cnode_merged = cliq_graph_add_node(cliq_graph, cliq_node_dict, 
        cno_merged, cand_len, cand_cov, 
        cno_merged + ":" + str(cand_len) + ":" + str(cand_cov))

        # removed the used up contig node
        cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno1, cnode1)
        cliq_graph_remove_node(cliq_graph, cliq_node_dict, cno2, cnode2)

        # removed the related L1 edges
        for edge1out in cnode1.out_edges():
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict,
            cliq_graph.vp.cno[edge1out.source()], cliq_graph.vp.cno[edge1out.target()], edge1out)

        # replace the related L1 in edges
        for edge1in in cnode1.in_edges():
            if cliq_graph.ep.color[edge1in] != 'black':
                continue
            src1 = edge1in.source()
            src1cno = cliq_graph.vp.cno[src1]
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            src1cno, cliq_graph.vp.cno[edge1in.target()], edge1in)
            
            # do not add self cycle edges
            if src1cno != cno_merged:
                cliq_graph_add_edge(cliq_graph, cliq_edge_dict, src1cno, src1, 
                cno_merged, cnode_merged, cliq_graph.ep.slen[edge1in], cliq_graph.ep.text[edge1in])
        
        # removed the related L2 in edges
        for edge2in in cnode2.in_edges():
            if cliq_graph.ep.color[edge2in] != 'black':
                continue
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            cliq_graph.vp.cno[edge2in.source()], cliq_graph.vp.cno[edge2in.target()], edge2in)
        
        # replace the related L2 out edges
        for edge2out in cnode2.out_edges():
            if cliq_graph.ep.color[edge2out] != 'black':
                continue
            tgt2 = edge2out.target()
            tgt2cno = cliq_graph.vp.cno[tgt2]
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, 
            cliq_graph.vp.cno[edge2out.source()], tgt2cno, edge2out)
            if cno_merged != tgt2cno:
                cliq_graph_add_edge(cliq_graph, cliq_edge_dict, cno_merged, cnode_merged,
                tgt2cno, tgt2, cliq_graph.ep.slen[edge2out], cliq_graph.ep.text[edge2out])
        return
    def buffer_concatenation(concat_buffer: list, cno_mapping: dict):
        for (cno1, cno2, cov, delta) in concat_buffer:
            print("------------------------------------------------------------------------")
            print("-->Before mapping: {0}: {1} - {2}: {3}".format(cno1, cno_mapping[cno1], cno2, cno_mapping[cno2]))
            for pcno in list(cno_mapping[cno1]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno1].remove(pcno)
                    print("pcno removed: ", pcno)
            for pcno in list(cno_mapping[cno2]):
                if pcno not in cliq_node_dict:
                    cno_mapping[cno2].remove(pcno)  
                    print("pcno removed: ", pcno)
            pairs = []
            for cno1m in cno_mapping[cno1]: 
                for cno2m in cno_mapping[cno2]:
                    if (cno1m, cno2m) in cliq_edge_dict:
                        pairs.append((cno1m, cno2m))
            if not pairs:
                print("contig has been used from previous step: {0} {1}".format(cno1, cno2))
                continue

            cno1m, cno2m = min(pairs, key=lambda p: 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[0]]] - cov, 2) + 
                pow(cliq_graph.vp.ccov[cliq_node_dict[p[1]]] - cov, 2))
            
            print("-->PAIR UP {0} - {1}, cov: {2}, diff: {3}".format(cno1m, cno2m, cov, delta))
            if (cno1m.split('->')[-1], cno2m.split('->')[0]) in sp_path_dict:
                cno1m_rel = cno1m.split('->')[-1]
                cno2m_rel = cno2m.split('->')[0]
                if sp_path_dict[(cno1m_rel, cno2m_rel)] != None:
                    cand_path, plen, pmark = sp_path_dict[(cno1m_rel, cno2m_rel)]
                else:
                    # cno1m -> end overlap with cno2m
                    intersect = set(contig_dict[cno1m][0]).intersection(set(contig_dict[cno2m][0]))
                    print("{0} -> {1} end to end intersects: {2}".format(cno1m, cno2m, intersect))
                    intermediate_nodes_index = [False for _ in contig_dict[cno1m][0]]
                    for i in [contig_dict[cno1m][0].index(e) for e in intersect]:
                        intermediate_nodes_index[i] = True
                    # in order
                    cand_path = [simp_node_dict[node_id] for i, node_id in enumerate(contig_dict[cno1m][0]) if intermediate_nodes_index[i]]
                    plen = path_len(graph, cand_path, overlap)
                    contig1m, clen1m, ccov1m = contig_dict[cno1m]
                    contig2m, clen2m, ccov2m = contig_dict[cno2m]
                    contig1m = contig1m[:-len(cand_path)]
                    clen1m = path_len(graph, [simp_node_dict[node_id] for node_id in contig1m], overlap)
                    contig_dict[cno1m] = [contig1m, clen1m, ccov1m]

                    contig2m = contig2m[len(cand_path):]
                    clen2m = path_len(graph, [simp_node_dict[node_id] for node_id in contig2m], overlap)
                    contig_dict[cno2m] = [contig2m, clen2m, ccov2m]
                    print("overlap contig concat: ", cno1m, contig1m, cno2m, contig2m, intersect)
            else:
                src = simp_node_dict[contig_dict[cno1m][0][-1]]
                tgt = simp_node_dict[contig_dict[cno2m][0][0]]
                cand_path, plen, pmark = dijkstra_sp(graph, src, tgt, cov, threshold, overlap)
            cand_len = get_concat_len(cno1m, contig_dict[cno1m][1], cno2m, contig_dict[cno2m][1], plen, overlap)

            contig_pair_reduction(cno1m, cno2m, cov, cand_path, cand_len, cno_mapping)

        return
    ###################################################################################################        

    # retrieve all the self cycle first
    for cno, contig_node in list(cliq_node_dict.items()):
        if contig_node in list(contig_node.out_neighbors()):
            if len(list(contig_node.all_edges())) > 2:
                # remove self cycle edge with self cycle + outer connection feature
                cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_graph.edge(contig_node, contig_node))
                print("should be removed already ... remove self edge+outer connection {0} -> {0}, legacy".format(cno))
            else:
                print("Circular PATH: ", cno)
                cov = cliq_graph.vp.ccov[contig_node]
                if (cno, cno) in sp_path_dict:
                    cand_path, plen, pmark = sp_path_dict[(cno, cno)]
                else:
                    src = simp_node_dict[contig_dict[cno][0][-1]]
                    tgt = simp_node_dict[contig_dict[cno][0][0]]
                    cand_path, plen, pmark = dijkstra_sp(graph, src, tgt, cov, threshold, overlap)
                cand_len = get_concat_len(cno, contig_dict[cno][1], cno, contig_dict[cno][1], plen, overlap)

                if cand_path != None and cand_len != None:
                    contig_dict[cno][0].extend([graph.vp.id[n] for n in cand_path])
                    contig_dict[cno][1] = cand_len
                    contig_dict[cno][2] = cov
                    cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, cno, cno, cliq_edge_dict[(cno, cno)])
                    print(path_to_id_string(graph, [simp_node_dict[cno] for cno in contig_dict[cno][0]] + cand_path, "cov: {0}".format(cov)))
                else:
                    print("Path not found, error")

    cno_mapping = {}
    for id in cliq_node_dict.keys():
        cno_mapping[id] = {id}

    for cno, contig in cliq_node_dict.items():
        contig, clen, ccov = contig_dict[cno]
        print_contig(cno, clen, ccov, contig)
    for (u, v), e in cliq_edge_dict.items():
        print("EDGE: ", u, v)

    if not graph_is_DAG(cliq_graph, cliq_node_dict):
        print("------>graph contains cycles, cyclic concatenation processing")
        cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
        cycles = all_circuits(cliq_graph, True)
        for cyc in cycles:
            print("cycle: ", [cliq_graph.vp.cno[n] for n in cyc])
            skip_cycle = False
            all_edges = {}
            for i in range(len(cyc)):
                u = cyc[i]
                v = cyc[(i+1) % len(cyc)]
                e = cliq_graph.edge(u, v)
                if cliq_graph.ep.color[e] != 'black':
                    skip_cycle = True
                    break
                all_edges[(cliq_graph.vp.cno[u],cliq_graph.vp.cno[v])] = cliq_graph.edge(u, v)
            if skip_cycle:
                continue
            removed_edge_tuple = min(all_edges.items(), key=lambda e_tuple: cliq_graph.ep.slen[e_tuple[1]])
            print("Removed edge: ", removed_edge_tuple[0][0], "->", removed_edge_tuple[0][1])
            cliq_graph_remove_edge(cliq_graph, cliq_edge_dict, removed_edge_tuple[0][0], removed_edge_tuple[0][1], removed_edge_tuple[1])

    print("--------------Start graduate concatentation------------------")
    for cno, cnos in cno_mapping.items():
        print("cno: ", cno, "maps to: ", cnos)
    # concat 
    # process all the non self-cycle contig until no isolated node.
    # direct pair two adjacent contig only if coverage difference is within pairwise threshold
    while True:
        # clean up the cno_mapping
        cno_mapping = {}
        for id in cliq_node_dict.keys():
            cno_mapping[id] = {id}
        # all the src contig
        L1_contigs = {}
        for no, contig_node in list(cliq_node_dict.items()):
            if cliq_graph.vp.color[contig_node] != 'black':
                continue
            ind = len([e for e in contig_node.in_edges() if cliq_graph.ep.color[e] == 'black'])
            outd = len([e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black'])
            if ind == 0 and outd != 0:
                # src node
                L1_contigs[no] = contig_node

        L2_contigs = {}
        st_pairs = {}
        for no, contig_node in list(L1_contigs.items()):
            out_edges = [e for e in contig_node.out_edges() if cliq_graph.ep.color[e] == 'black']
            for out_e in out_edges:
                out_contig_node = out_e.target()
                if cliq_graph.vp.color[out_contig_node] != 'black':
                    continue
                outid = cliq_graph.vp.cno[out_contig_node]
                if not outid in L2_contigs:
                    L2_contigs[outid] = out_contig_node
                st_pairs[(no, outid)] = (contig_node, out_contig_node)
        
        if not L1_contigs or not L2_contigs:
            print("no more linear concatentation left, break")
            break

        print("L1: ", [no for no in L1_contigs.keys()])
        print("L2: ", [no for no in L2_contigs.keys()])
        
        # detect the most confident contig pair first
        concat_buffer = []
        for (cno1, cno2), (node1, node2) in st_pairs.items():
            delta = abs(cliq_graph.vp.ccov[node1] - cliq_graph.vp.ccov[node2])
            concat_buffer.append((cno1, cno2, min(cliq_graph.vp.ccov[node1], cliq_graph.vp.ccov[node2]), delta))
        
        concat_buffer = sorted(concat_buffer, key=lambda tuple: tuple[3])
        print("all most confident sorted concats are: ", concat_buffer)
        buffer_concatenation(concat_buffer, cno_mapping)


    for (u, v), e in cliq_edge_dict.items():
        print("potential uncaught error: EDGE: ", u, v)
    
    # simplify the graph
    cliq_graph, cliq_node_dict, cliq_edge_dict = cliq_graph_init(cliq_graph)
    # draw_cliq_graph(cliq_graph, len(cliq_node_dict), len(cliq_edge_dict), tempdir, "cliq_graphL3.png")

    return contig_dict
    
    
    # store the no-cycle nodes in nc_graph_L3p.gfa
    # noncyc_nodes = None
    # simple_paths = None
    # # graphtool is_DAG() may not work if the graph is not connected as several parts
    # if not graph_is_DAG(graph, simp_node_dict):
    #     noncyc_nodes, simple_paths = node_partition(graph, simp_node_dict, tempdir)


    # if noncyc_nodes != None and simple_paths != None:
    #     print("rebalance linear subgraph now..")
    #     graphnc, simp_node_dictnc, simp_edge_dictnc = flipped_gfa_to_graph("{0}/gfa/nc_graph_L2p.gfa".format(tempdir))
    #     coverage_rebalance_ave(graphnc, simp_node_dictnc, simp_edge_dictnc)
    #     print("Done, start coverage merge")

    #     for no, node in simp_node_dictnc.items():
    #         cnode = simp_node_dict[no]
    #         merge_dp = graphnc.vp.dp[node] + graph.vp.dp[cnode]
    #         graph.vp.dp[cnode] = merge_dp   
    # else:
    #     print("no linear subgraph available..")
