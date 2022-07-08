#!/usr/bin/env python3

from graph_tool import Graph
import sys
import heapq
from collections import deque
from graph_tool.topology import all_circuits, all_shortest_paths, min_spanning_tree, topological_sort
from numpy import sort
import graph_converter
import gurobipy as gp

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
                print("graph is cyclic")
                return False
    print("graph is not cyclic")
    return True

def retrieve_cycle(graph: Graph, n=1):
    """
    retrieve a cycle, if any, else return None, sequential
    """
    cycles = []
    def processDFSTree(graph: Graph, stack: list, visited: list, n):
        for out_e in stack[-1].out_edges():
            if graph.ep.color[out_e] != 'black':
                continue
            if n == 0:
                return n
            next = out_e.target()
            if visited[next] == 'instack':
                # revisit a visited node, cycle
                n -= 1
                store_cycle(stack, next)
            elif visited[next] == 'unvisited':
                visited[next] = 'instack'
                stack.append(next)
                n = processDFSTree(graph, stack, visited, n)
        visited[stack[-1]] = 'done'
        stack.pop()
        return n

    def store_cycle(stack: list, next):
        stack = stack[stack.index(next):]
        print("Cycle: ", graph_converter.list_to_string([graph.vp.id[node] for node in stack]))
        cycles.append(stack)

    visited = {}
    for node in graph.vertices():
        visited[node] = 'unvisited'
    
    for v in graph.vertices():
        if visited[v] == 'unvisited':
            stack = [v]
            visited[v] = 'instack'
            n = processDFSTree(graph, stack, visited, n)
    
    return cycles if len(cycles) > 0 else None

def reachable(graph: Graph, simp_node_dict: dict, src, tgt):
    """
    determine whether src can possibly reach the tgt
    """
    visited = {}
    for no in simp_node_dict.keys():
        visited[no] = False

    queue = [src]
    while queue:
        curr = queue.pop()
        visited[graph.vp.id[curr]] = True
        if curr == tgt:
            return True
        for oute in curr.out_edges():
            out = oute.target()
            if not visited[graph.vp.id[out]]:
                queue.append(out)
    return False
    
def reachable_with_path(graph: Graph, simp_node_dict: dict, src, src_cno, tgt, tgt_cno, contig_dict: dict):
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


def dijkstra_sp(graph: Graph, source, sink, closest_cov, threshold, overlap: int):
    """
    Use basic dijkstra algorithm to compute the shortest path between source and sink
    """
    print("Dijkastra path finding: closest cov guide: ", closest_cov)
    dist = {}
    prev = {}
    Q = set()
    for node in graph.vertices():
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
                # obtain current edge cost
                edge_flow = graph.ep.flow[graph.edge(u, v)]
                diff = edge_flow - closest_cov
                # alt = dist[u] + pow(diff, 2)
                if abs(diff) > threshold and diff < 0:
                    alt = dist[u] + pow(diff, 2)
                else:
                    alt = dist[u] + abs(diff)
                # relax
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = u
    
    if dist[sink] == sys.maxsize:
        # not reachable
        return None, None, None

    node = sink
    sp = []
    mark = 0
    while prev[node]:
        sp.insert(0, node)
        mark += (graph.ep.flow[graph.edge(prev[node], node)] - closest_cov)**2
        node = prev[node]
    if not sp:
        print("SP not found")
        return None, None, None
    else:
        sp.insert(0, source)
        print(graph_converter.path_to_id_string(graph, sp, "SP be found: "))
        mark = 1/(1+mark**0.5)
        print("plen: ", graph_converter.path_len(graph, sp[1:-1], overlap), "pmark: ", mark)

        return sp[1:-1], graph_converter.path_len(graph, sp[1:-1], overlap), mark

def dijkstra_sp_v2(graph: Graph, simp_node_dict: dict, source, sink, overlap: int):
    """
    Use basic dijkstra algorithm to compute the max coverage between source and sink
    """
    print("Dijkastra path finding: Max cov possible guide")
    if graph.vp.color[source] != 'black' or graph.vp.color[sink] != 'black':
        print("s/t unavailable")
        return None, None
    dist = {}
    prev = {}
    Q = set()
    for node in simp_node_dict.values():
        if graph.vp.color[node] == 'black':
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
        # found sink
        if u == sink:
            break
        for v in u.out_neighbors():
            if v in Q:
                # obtain current edge cost
                e = graph.edge(u, v)
                if graph.ep.color[e] != 'black':
                    continue
                alt = dist[u] + pow(graph.ep.flow[e],2)
                # relax
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = u
    
    if dist[sink] == sys.maxsize:
        print("not reachable")
        return None, None

    node = sink
    sp = []
    while prev[node]:
        sp.insert(0, node)
        node = prev[node]
    if not sp:
        print("SP not found")
        return None, None
    elif node != source:
        print("SP src not found")
        return None, None
    else:
        sp.insert(0, source)
        print(graph_converter.path_to_id_string(graph, sp, "SP be found: "))
        print("plen: ", graph_converter.path_len(graph, sp[1:-1], overlap))

        return sp[1:-1], graph_converter.path_len(graph, sp[1:-1], overlap)

def eval_score(flow, ccov, threshold, min_e):
    diff = flow - ccov
    if diff < -threshold:
        return "P4"
    elif diff >= -threshold and diff <= threshold:
        return "P1"
    elif diff > threshold and diff <= min_e - threshold:
        return "P3"
    elif diff > min_e - threshold:
        return "P2"
    return None

def dijkstra_sp_v3(graph: Graph, source, sink, closest_cov, threshold, overlap: int, min_e):
    """
    Use basic dijkstra algorithm to compute the shortest path between source and sink
    """
    print("Dijkastra path finding: closest cov guide: ", closest_cov)
    dist = {}
    prev = {}
    Q = set()
    for node in graph.vertices():
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
                # obtain current edge cost
                edge_flow = graph.ep.flow[graph.edge(u, v)]
                diff = edge_flow - closest_cov
                if diff < -threshold:
                    #P4
                    alt = dist[u] + (- diff)/threshold
                elif diff >= -threshold and diff <= threshold:
                    #P1
                    alt = dist[u] - (len(str(graph.vp.id[v]).split('_'))/(abs(diff) + 1))*threshold
                elif diff > threshold and diff <= 2*threshold:
                    #P3
                    # alt = dist[u] + 1
                    alt = dist[u] + ((diff - threshold) / threshold)
                elif diff > 2*threshold:
                    #P2
                    alt = dist[u] + 0
                # relax
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = u
    
    if dist[sink] == sys.maxsize:
        # not reachable
        return None, None, None

    node = sink
    sp = []
    mark = 0
    while prev[node]:
        sp.insert(0, node)
        mark += (graph.ep.flow[graph.edge(prev[node], node)] - closest_cov)**2
        node = prev[node]
    if not sp:
        print("SP not found")
        return None, None, None
    else:
        sp.insert(0, source)
        print(graph_converter.path_to_id_string(graph, sp, "SP be found: "))
        mark = 1/(1+mark**0.5)
        print("plen: ", graph_converter.path_len(graph, sp[1:-1], overlap), "pmark: ", mark)

        return sp[1:-1], graph_converter.path_len(graph, sp[1:-1], overlap), mark

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

def transitive_graph_reduction(cliq_node_dict: dict, cliq_edge_dict: dict):
    for k in cliq_node_dict.keys():
        for i in cliq_node_dict.keys():
            for j in cliq_node_dict.keys():
                if i != j and i != k and j != k:
                    if (i, k) in cliq_edge_dict and (k, j) in cliq_edge_dict and (i, j) in cliq_edge_dict:
                        ik = cliq_edge_dict[(i,k)]
                        kj = cliq_edge_dict[(k,j)]
                        ij = cliq_edge_dict[(i,j)]
                        cliq_edge_dict[(i, j)] = 'gray'
                        cliq_edge_dict.pop((i, j))                                            

    return

def paths_from_src(graph: Graph, simp_node_dict: dict, self_node, src, overlap, maxlen):
    """
    retrieve all the path from src node to any node 
    within maxlen restriction, in straight direction
    """
    def dfs_rev(graph: Graph, u, curr_path: list, maxlen, visited, all_path):
        visited[u] = True
        curr_path.append(u)
        curr_len = graph_converter.path_len(graph, curr_path, overlap)
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

def paths_to_tgt(graph: Graph, simp_node_dict: dict, self_node, tgt, overlap, maxlen):
    """
    retrieve all the path from any node to tgt node
    within maxlen restriction, in reverse direction
    """
    def dfs_rev(graph: Graph, v, curr_path: list, maxlen, visited, all_path):
        visited[v] = True
        curr_path.insert(0, v)
        curr_len = graph_converter.path_len(graph, curr_path, overlap)
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

def st_variation_path(graph: Graph, src, src_contig, tgt, tgt_contig, closest_cov, overlap):
    """
    find the optimal path from src tgt, where intermediate nodes with cov ~ cand_cov would be prefered
    """
    print("Start st variation finding: {0} -> {1}".format(graph.vp.id[src], graph.vp.id[tgt]))

    rtn_paths = []
    local_visited = {}
    global_visited = {}
    for node in graph.vertices():
        global_visited[graph.vp.id[node]] = False
        local_visited[graph.vp.id[node]] = False
    # avoid searching back to contig itself
    for id in src_contig[:-1]:
        global_visited[id] = True
    for id in tgt_contig[1:]:
        global_visited[id] = True
    
    # mark source node as visited
    global_visited[graph.vp.id[src]] = True
    local_visited[graph.vp.id[src]] = True 

    path = [src]
    pathqueue = deque()
    pathqueue.append([path, local_visited, len(graph.vp.seq[src]), 0]) 
    
    while pathqueue:
        curr_path, local_visited, curr_len, curr_mark = pathqueue.popleft()
        if curr_path[-1] == tgt and curr_path[0] == src:
            rtn_paths.append([curr_path, curr_len, 1/(1 + curr_mark**0.5)])
            continue

        for e in curr_path[-1].out_edges():
            next = e.target()
            if graph.ep.color[e] != 'black' and graph.vp.color[next] != 'black':
                continue
            if global_visited[graph.vp.id[next]] or local_visited[graph.vp.id[next]]:
                continue
            if next not in curr_path:
                split_path = curr_path[:]
                split_local_visited = dict(local_visited)
                split_path.append(next)
                split_local_visited[graph.vp.id[next]] = True
                split_len = curr_len + len(graph.vp.seq[next]) - overlap
                split_mark = curr_mark + (graph.ep.flow[e] - closest_cov)**2
                pathqueue.append([split_path, split_local_visited, split_len, split_mark])

    # for p, plen in rtn_paths:
    #     if p[0] == src and p[-1] == tgt:
    #         print("Length: ", plen, graph_converter.path_to_id_string(graph, p, "-Path be found: "))

    return rtn_paths

def path_replacement_account(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, usage_dict: dict, s, t):
    """
    obtain any complete partial used path and overused path
    """
    def dfs_rev(graph: Graph, v, curr_path: list, visited: dict, all_path: list, usage_dict: dict, cond):
        if len(curr_path) > 0 and curr_path[-1] == v:
            all_path.append(list(curr_path[1:-1]))
        else:
            for next in curr_path[-1].out_neighbors():
                if not visited[next] and (usage_dict[graph.vp.id[next]][2] == cond or next == v):
                    visited[next] = True
                    curr_path.append(next)
                    dfs_rev(graph, v, curr_path, visited, all_path, usage_dict, cond)
                    curr_path.pop()
                    visited[next] = False
        return
    all_path_o = []
    visited = {}
    for node in graph.vertices():
        visited[node] = False
    dfs_rev(graph, t, [s], visited, all_path_o, usage_dict, "over")

    all_path_p = []
    for node in graph.vertices():
        visited[node] = False
    dfs_rev(graph, t, [s], visited, all_path_p, usage_dict, "partial")

    return all_path_o, all_path_p


def retrieve_bubble(graph: Graph, simp_edge_dict: dict, s, t):
    """
    since the bubble is at most length =1, simply return 
    all the connected intermediate node.
    """
    bubbles = []
    for child in t.in_neighbors():
        if graph.vp.color[child] == 'black':
            prev_conn = (graph.vp.id[s], graph.vp.id[child])
            if prev_conn in simp_edge_dict:
                if graph.ep.color[simp_edge_dict[prev_conn]] == 'black':
                    bubbles.append(child)
    print("Bubbles: ", graph_converter.path_to_id_string(graph, bubbles))
    return bubbles

def gen_bubble_detection(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, overlap):
    def print_branch_queue(branch_queue):
        for close, ind, outd in branch_queue:
            print("|{0}, in: {1}, out:{2}|".format(graph.vp.id[close], ind, outd))
    global_src, global_sink = graph_converter.add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap)
    # use topological sort to figure out the order of the branch node
    if not graph_is_DAG(graph, simp_node_dict):
        print("Graph is not DAG, obtain spanning tree")
        tree = min_spanning_tree(graph)
        graph.set_edge_filter(tree)
    sort = topological_sort(graph)
    node_order = {}
    print("TOP SORT")
    for i, n in enumerate(sort):
        node = graph.vertex(n)
        if node.in_degree() > 1 or node.out_degree() > 1:
            node_order[node] = i
            print("branch: ", graph.vp.id[n])
        else:
            node_order[node] = sys.maxsize
    print("TOP SORT END")

    visited = {}
    prev = {}
    for node in graph.vertices():
        visited[node] = False
        prev[node] = None
    visited[global_src] = True
    branch_queue = []
    normal_queue = [[global_src, 0, global_src.out_degree(), True]]
    loopCount = 30
    while normal_queue:
        loopCount -= 1
        if loopCount <= 0:
            break
        currnode, ind, outd, isBranch = normal_queue.pop()
        print("----->Curr pop: {0}, {1}, {2}, isBranch: {3}".format(graph.vp.id[currnode], ind, outd, isBranch))
        for v in currnode.out_neighbors():
            if not visited[v]:
                visited[v] = True
                print("append brand new nodes: ", graph.vp.id[v])
                b = not (v.in_degree() <= 1 and v.out_degree() <= 1)
                normal_queue.insert(0, [v, v.in_degree(), v.out_degree(), b])
                if v.in_degree() == 1:
                    prev[v] = [currnode, isBranch]
        
        if isBranch:
            # check if it is the minimum branch within normal queue
            min_branch_inQ, _, _, isBranchQ= min(normal_queue, key=lambda t: node_order[t[0]])
            if node_order[min_branch_inQ] < node_order[currnode] and isBranchQ:
                print("there are lower order branch dig in the queue, push the current back")
                normal_queue.insert(0, [currnode, ind, outd, isBranch])
            else:
                print("curr branch is the lowest one, push it to branch queue")
                branch_queue.insert(0, [currnode, ind, outd])

            if len(branch_queue) > 1:
                print("BEFORE BUBBLE CHECK: ")
                print_branch_queue(branch_queue)
                # last inserted, highest order
                leftClose, lind, loutd = branch_queue[0]
                print("Leftmost: ", graph.vp.id[leftClose], lind, loutd)
                for i in range(1, len(branch_queue)):
                    rightClose, rind, routd = branch_queue[i]
                    # bubble found, deal with it.
                    print("Bubble found: ", graph.vp.id[rightClose], rind, routd, "<->", lind)
                    # do some work
                    retrieve_bubble(graph, simp_edge_dict, rightClose, leftClose)
                    #...
                    # clean up
                    midd = min(lind, routd)
                    lind -= midd
                    routd -= midd
                    branch_queue[i] = [rightClose, rind, routd]
                    print("end: ", rind, routd, "<->", lind)
                branch_queue[0] = [leftClose, lind, loutd]
                branch_queue_copy = []
                for close, ind, outd in branch_queue:
                    if ind != 0 or outd != 0:
                        branch_queue_copy.append([close, ind, outd])
                branch_queue = branch_queue_copy
                print("AFTER BUBBLE CHECK: ")
                print_branch_queue(branch_queue)

def minimal_bubble_detection(graph: Graph):
    # retrieve all the minimal bubble
    # every bubble is structure as (start branch -> single node child -> end branch)
    branches = graph_converter.retrieve_branch(graph)
    bubbles = {}
    for no, branch in branches.items():
        bubble_dict = {}
        for child in branch.out_neighbors():
            if graph.vp.id[child] not in branches:
                if child.out_degree() == 0:
                    # iff the child  == global sink
                    print("Child is sink node: ", graph.vp.id[child])
                elif child.out_degree() == 1:
                    accBranch = list(child.out_neighbors())[0]
                    if graph.vp.id[accBranch] in branches:
                        if accBranch not in bubble_dict:
                            bubble_dict[accBranch] = [child]
                        else:
                            bubble_dict[accBranch].append(child)
                    else:
                        print("NEXT BRANCH ERROR, NOT A LISTED BRANCH", graph.vp.id[branch], graph.vp.id[child], graph.vp.id[accBranch])
                else:
                    print("CHILD ERROR, SHOULD BE MARKED AS A BRANCH", graph.vp.id[branch], graph.vp.id[child])
        for nextBranch, bubble in bubble_dict.items():
            if len(bubble) > 1:
                bubbles[(no, graph.vp.id[nextBranch])] = bubble
    return branches, bubbles

def LPSolveBubbleLocal(bin_dict: dict, strain_dict: dict, threshold):
    """
    solve the bin-assignment problem by assigning strain only once to best match bin
    such that maximize the bin usage if possible
    """
    print("Start Linear Programming solver..")

    bin_list = list(bin_dict.items())
    strain_list = list(strain_dict.items())
    bin = [v for _, v in bin_list]
    strain = [v for _, v in strain_list]

    threshold = (sum(strain) - sum(bin)) + threshold if sum(strain) > sum(bin) else threshold
    
    print("bin sum: ", sum(bin), " strain sum: ", sum(strain))
    print("bin: ", bin_dict)
    print("strains: ", strain_dict)
    print("alpha: ", threshold)

    # FIXME lock step should fix, to avoid the order influence
    locked = {}
    for i in range(len(strain)):
        locked[i] = False

    m = gp.Model("BubbleSwap")
    xs = []
    obj = gp.LinExpr()
    for b in bin:
        le = gp.LinExpr()
        x = m.addVars([i for i in range(len(strain))] ,vtype=gp.GRB.BINARY)
        xs.append(x)
        # check how many lockable item
        lockable = []
        # FIXME better lock method, may use double-match
        for i, v in x.items():
            if abs(strain[i] - b) < threshold and not locked[i]:
                lockable.append(i)
            le += v*strain[i]
        if len(lockable) != 0:
            lockable_s = sorted(lockable, key=lambda lk : abs(strain[lk] - b))
            print("sorted lockable: ", lockable_s, "bin: ", b)
            m.addConstr(x[lockable_s[0]] == 1)
            locked[lockable_s[0]] = True
            print("pre-locked: bin: {0}, strain: {1}".format(b, strain[lockable_s[0]]))
        # m.addConstr(le <= (b+threshold))
        obj += (b - le)**2

    # add constraint, such that each strain only used once
    for i in range(len(strain)):
        le = gp.LinExpr()
        for item in xs:
            le += item[i]
        m.addConstr(le == 1)

    m.setObjective(obj, gp.GRB.MINIMIZE)

    try:
        m.optimize()
    except gp.GurobiError():
        print("Optimize failed due to non-convexity")
    #  obtain the decision variable result
    #  which strain is contained in which bin
    bin_assignment = {}
    strain_placement = {}
    for i, x in enumerate(xs):
        # each x is telling the current ith bin assignment
        bin_id = bin_list[i][0]
        bin_assignment[bin_id] = []
        for j, var in x.items():
            if int(var.x) == 1:
                # i.e., current jth strain is assigned to ith bin
                strain_id = strain_list[j][0]
                bin_assignment[bin_id].append(strain_id)
                if strain_id not in strain_placement:
                    strain_placement[strain_id] = [bin_id]
                else:
                    strain_placement[strain_id].append(bin_id)
                
    print("Bin assignment: ", bin_assignment)
    print("Strain placement: ", strain_placement)
    print([int(var.x) for var in m.getVars()])
    print(m.objVal)
    print(m.display())
    return bin_assignment, strain_placement