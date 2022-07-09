#!/usr/bin/env python3
from graph_converter import *
from search_algos import *
from graph_tool import Graph
import heapq
from collections import deque
from graph_tool.topology import all_circuits, all_shortest_paths, min_spanning_tree, topological_sort
from numpy import sort
import gurobipy as gp

def node_status(used, capacity, threshold):
    if used - capacity > threshold:
        return 'over'
    else:
        if abs(used - capacity) < threshold:
            return 'full'
        else:
            return 'partial'

def test_comb(graph: Graph, comb: list, usage_dict: dict, all_cnos: list, all_paths: dict, prev_cno_mapping: dict, contig_dict: dict, threshold):
    """
    legacy
    """
    usage_dict_c = copy.deepcopy(usage_dict) #copy
    new_state_dict = {}
    for i, key in enumerate(comb):
        # update contig dict
        cno = all_cnos[i]
        prev_p = [graph.vp.id[n] for n in all_paths[prev_cno_mapping[cno]][0]]
        new_p = [graph.vp.id[n] for n in all_paths[key][0]]

        for prev_id in prev_p:
            usage_dict_c[prev_id][0] -= contig_dict[cno][2]
            usage_dict_c[prev_id][2] = node_status(usage_dict_c[prev_id][0], usage_dict_c[prev_id][1], threshold)
        for new_id in new_p:
            usage_dict_c[new_id][0] += contig_dict[cno][2]
            usage_dict_c[new_id][2] = node_status(usage_dict_c[new_id][0], usage_dict_c[new_id][1], threshold)
    for [path, _] in all_paths.values():
        for n in path:
            if graph.vp.id[n] not in new_state_dict:
                new_state_dict[graph.vp.id[n]] = usage_dict_c[graph.vp.id[n]]
    return new_state_dict

def draw_edgeflow(graph: Graph, edge_dict: dict, tempdir, title, filename):
    """
    draw count plot for graph edge coverage, legacy
    """
    seaborn.set_theme(style="darkgrid")
    plt.figure(figsize=(128,64))
    drawdict = {}
    for id, e in edge_dict.items():
        drawdict[id] = graph.ep.flow[e]
    sorted_draw_dict = sorted(drawdict.items(), key=lambda x: x[1])
    df = pandas.DataFrame(
        {'Id': [id for id, _ in sorted_draw_dict],
        'Flow': [f for _, f in sorted_draw_dict]
        })
    ax = seaborn.barplot(x='Id', y='Flow', data=df)
    for container in ax.containers:
        ax.bar_label(container)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    plt.title(title)
    plt.savefig("{0}{1}".format(tempdir, filename))

def cyclic_to_dag(graph: Graph, simp_node_dict: dict, contig_dict: dict, overlap):
    """
    convert graph to dag by delete minimum coverage not-in-contig edge in cycle until reduced to dag, legacy
    """
    def remove_edge(fst, snd):
        print("removing edge: {0} -> {1} to reduce a cycle".format(graph.vp.id[fst], graph.vp.id[snd]))
        graph.ep.color[graph.edge(fst, snd)] = 'gray'
        removed_edges.append((graph.vp.id[fst], graph.vp.id[snd]))
        return fst, snd
    removed_edges = []
    while not graph_is_DAG(graph, simp_node_dict):
        cycle = retrieve_cycle(graph)[0]
        left_node = None
        right_node = None
        removed_cycle = False
        _, edge_to_contig_dict = contig_map_node(contig_dict)
        for min_node in sorted(cycle, key=lambda v: graph.vp.dp[v]):
            # remove min_node prev edge
            prev_node = cycle[(cycle.index(min_node) - 1) % len(cycle)]
            next_node = cycle[(cycle.index(min_node) + 1) % len(cycle)]

            if graph.vp.dp[prev_node] < graph.vp.dp[next_node]:
                if (graph.vp.id[prev_node], graph.vp.id[min_node]) not in edge_to_contig_dict:
                    left_node, right_node = remove_edge(prev_node, min_node)
                    removed_cycle = True
            else:
                if (graph.vp.id[min_node], graph.vp.id[next_node]) not in edge_to_contig_dict:
                    left_node, right_node =  remove_edge(min_node, next_node)
                    removed_cycle = True
            if removed_cycle:
                break
        # final round, if all edge is involved in contig, split the minimal one, worst case
        if not removed_cycle:
            min_node = min(cycle, key=lambda v: graph.vp.dp[v])
            prev_node = cycle[(cycle.index(min_node) - 1) % len(cycle)]
            next_node = cycle[(cycle.index(min_node) + 1) % len(cycle)]
            if graph.vp.dp[prev_node] < graph.vp.dp[next_node]:
                left_node, right_node =  remove_edge(prev_node, min_node)
                print("involved contig: ", edge_to_contig_dict[(graph.vp.id[left_node], graph.vp.id[right_node])])
            else:
                left_node, right_node =  remove_edge(min_node, right_node)
                print("involved contig: ", edge_to_contig_dict[(graph.vp.id[left_node], graph.vp.id[right_node])])

        if (graph.vp.id[left_node], graph.vp.id[right_node]) in edge_to_contig_dict:
            for cno in edge_to_contig_dict[(graph.vp.id[left_node], graph.vp.id[right_node])]:
                contig, clen, ccov = contig_dict.pop(cno)
                left_contig = contig[:contig.index(graph.vp.id[left_node]) + 1]
                right_contig = contig[contig.index(graph.vp.id[right_node]):]
                contig_dict[cno+"l"] = [left_contig, path_len(graph, [simp_node_dict[id] for id in left_contig], overlap), ccov]
                contig_dict[cno+"r"] = [right_contig, path_len(graph, [simp_node_dict[id] for id in right_contig], overlap), ccov]
    return removed_edges

def align_reads_to_contig(contig_dict: dict, contig_file, forward_file, reverse_file, temp_dir, accept_rate = 0.95):
    """
    Use minimap2 to align pair-end reads information to contig, legacy
    """
    c2c_mat = {}
    for cno in contig_dict.keys():
        c2c_mat[cno] = {}
        for cno2 in contig_dict.keys():
            c2c_mat[cno][cno2] = 0
    subprocess.check_call("touch {0}aln.paf".format(temp_dir), shell=True)
    # align reads to contig to strong the contig pairwise connectivity.
    subprocess.check_call("/Users/luorunpeng/opt/miniconda3/envs/spades-hapConstruction-env/bin/minimap2 -x sr {0} {1} {2} > {3}aln.paf".format(contig_file, forward_file, reverse_file, temp_dir), shell=True)
    with open("{0}aln.paf".format(temp_dir), 'r') as aln:
        curr_read = None
        curr_contigs = None
        for i, line in enumerate(aln):
            splited = line.split('\t')
            seg_no = str(splited[0])
            ref_no = str(splited[5])
            nmatch = int(splited[9])
            nblock = int(splited[10])
            mark = int(splited[11])
            if i == 0:
                curr_read = seg_no
                curr_contigs = set()
            elif curr_read != seg_no:
                # process prev reads
                for c1, c2 in list(product(list(curr_contigs), repeat=2)):
                    if c1 != c2:
                        c2c_mat[c1][c2] += 1             
                curr_read = seg_no
                curr_contigs = set()

            if nmatch/nblock >= accept_rate:
                curr_contigs.add(ref_no.split('_')[0])
        aln.close()
    for cno, row in c2c_mat.items():
        print(cno, [(k,v) for k,v in row.items()])
    return c2c_mat

# extract candidate path version 1
def extract_cand_path(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, strain_dict: dict, overlap, threshold, itercount='A'):
    global TEMP_DIR
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap, True)
    for e in sorted(graph.edges(), key=lambda a: graph.ep.flow[a]):
        print_edge(graph, e)
    lb = min([graph.ep.flow[e] for e in graph.edges()])
    ub = lb*2
    print("[{0},{1})".format(lb, ub))
    threshold = ub/20
    print("local threshold: ", threshold)

    min_eset = {}
    for e in graph.edges():
        if graph.ep.flow[e] > lb and graph.ep.flow[e] <= ub:
            min_eset[(graph.vp.id[e.source()], graph.vp.id[e.target()])] = e

    x = [graph.ep.flow[e] for e in min_eset.values()]
    regions, bins = numpy.histogram(x)
    # find peak region
    peak_region_index = list(regions).index(max(regions))
    peak_lb = bins[peak_region_index+1] - threshold
    # peak_lb = bins[peak_region_index]
    # peak_ub = bins[peak_region_index+1]
    peak_ub = bins[peak_region_index] + threshold
    #################################STAT###################################
    print(regions)
    print(bins)
    print(peak_lb, bins[peak_region_index], bins[peak_region_index+1], peak_ub)
    plt.figure(figsize=(128,64))
    plt.hist(x, bins)
    plt.axvline(x=peak_lb, color='r')
    plt.axvline(x=bins[peak_region_index], color='b')
    plt.axvline(x=bins[peak_region_index+1], color='b')
    plt.axvline(x=peak_ub, color='r')
    plt.title("min_eset")
    plt.savefig("{0}{1}".format(TEMP_DIR, "min_eset.png"))
    ########################################################################
    peak_eset = dict()
    peak_vset = set()
    connect_set = dict()
    for k, e in min_eset.items():
        if graph.ep.flow[e] >= peak_lb and graph.ep.flow[e] < peak_ub:
            peak_eset[k] = e
    for (u,v), e in peak_eset.items():
        connect_set[e] = []
        peak_vset.add(u)
        peak_vset.add(v)    
    print(list_to_string(list(peak_vset), "Peak related nodes: "))

    # connectivity
    for (au, av), e_tail in peak_eset.items():
        for (bu, bv), e_head in peak_eset.items():
            if au == bu and av == bv:
                continue
            if reachable(graph, simp_node_dict, simp_node_dict[av], simp_node_dict[bu]):
                connect_set[e_tail].append(e_head)
    
    # transitive reduction
    for k in connect_set:
        for i in connect_set:
            for j in connect_set:
                if k != i and i != j and k != j:
                    if k in connect_set[i] and j in connect_set[k] and j in connect_set[i]:
                        connect_set[i].remove(j)
    
    graphm = Graph(directed=True)
    graphm.vp.id = graphm.new_vertex_property("string")
    graphm.vp.flow = graphm.new_vertex_property("int32_t")
    graphm.vp.text = graphm.new_vertex_property("string")
    graphm.ep.color = graphm.new_edge_property("string")
    graphm_vdict = dict()
    graphm_edict = dict()
    for u in connect_set.keys():
        vertex = graphm.add_vertex()
        graphm.vp.id[vertex] = str(graph.vp.id[u.source()]) + "->" + str(graph.vp.id[u.target()])
        graphm.vp.flow[vertex] = graph.ep.flow[u]
        graphm.vp.text[vertex] = graphm.vp.id[vertex] + ":" + str(graphm.vp.flow[vertex])
        graphm_vdict[u] = vertex
    for u, vs in connect_set.items():
        for v in vs:
            e = graphm.add_edge(graphm_vdict[u], graphm_vdict[v])
            graphm_edict[(u,v)] = e
    
    #################################STAT###################################
    output_size = 120 * (len(list(graphm.vertices()) )+ len(list(graphm.edges())))
    vsize= 30
    graph_draw(g=graphm, output="{0}{1}".format(TEMP_DIR, "graphm.png"), bg_color="white", 
    vertex_text=graphm.vp.text, vertex_size=vsize, vertex_font_size=int(vsize * 0.8), 
    output_size=(output_size, output_size))
    print("min peak graph has been stored in: {0}{1}".format(TEMP_DIR, "graphm.png"))
    ########################################################################
    adjMtx = [[] for _ in connect_set.keys()]
    e2i = {}
    i2e = {}
    for i, e in enumerate(connect_set.keys()):
        e2i[e] = i
        i2e[i] = e
    for u, vs in connect_set.items():
        row = [[sys.maxsize,'X'] for _ in connect_set.keys()]
        for v in vs:
            abs_dif = abs(graph.ep.flow[u] - graph.ep.flow[v])
            row[e2i[v]] = [abs_dif, 'W']
        adjMtx[e2i[u]] = row
    
    has_changes = True
    dim = len(adjMtx)
    # iterate the adjacent matrix in rowwise and columnwise
    while has_changes:
        print("looping")
        has_changes = False
        # row wise
        for rowId in range(dim):
            row = get_row(adjMtx, rowId)
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
                    print("Error: ", rowId, i2e[rowId], cid, i2e[cid])
                    assert color != 'R'
            if not hasFixed:
                if colId != None:
                    # found a minimum block to assign
                    adjMtx[rowId][colId][1] = 'R'
                    has_changes = True
                else:
                    # no more spot for the row
                    None
        for colId in range(dim):
            col = get_col(adjMtx, colId)
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
                    adjMtx[rid][colId] = [flow, 'B']
                    has_changes = True
            else:
                mrid, [mflow, _] = min(cands, key=lambda p: p[1][0])
                for rid, [flow, color] in cands:
                    adjMtx[rid][colId] = [flow, 'G']
                adjMtx[mrid][colId] = [mflow, 'B']
                has_changes = True
    adj_list = {}
    for e in peak_eset.values():
        adj_list[e] = None
    for rowId, row in enumerate(adjMtx):
        for colId, [_, color] in enumerate(row):
            if color != 'X':
                e = graphm_edict[(i2e[rowId], i2e[colId])]
                if color != 'B':
                    graphm.ep.color[e] = 'gray'
                else:
                    graphm.ep.color[e] = 'black'
                    if adj_list[i2e[rowId]] != None:
                        print("Error, branch detected")
                    else:
                        adj_list[i2e[rowId]] = i2e[colId]
    
    for e in sorted(graphm.edges()):
        if graphm.ep.color[e] != 'black':
            graphm.remove_edge(e)
    #################################STAT###################################
    output_size = 120 * (len(list(graphm.vertices()) )+ len(list(graphm.edges())))
    vsize = 30
    graph_draw(g=graphm, output="{0}{1}".format(TEMP_DIR, "graphm_v2.png"), bg_color="white", 
    vertex_text=graphm.vp.text, vertex_size=vsize, vertex_font_size=int(vsize * 0.8), 
    output_size=(output_size, output_size))
    print("min peak graph has been stored in: {0}{1}".format(TEMP_DIR, "graphm_v2.png"))
    ########################################################################
    epaths = []
    for u in adj_list.keys():
        if graphm_vdict[u].in_degree() != 0:
            # intermediate node
            continue
        epath = []
        acc = u
        while acc != None:
            epath.append(acc)
            acc = adj_list[acc]
        epaths.append(epath)
    for i, epath in enumerate(sorted(epaths, key=len, reverse=True)):
        strain = []
        print("--->{0}{1} Strain finding.. length: {2}".format(itercount, i, len(epath)))
        #FIXME
        ccov = max([graph.ep.flow[e] for e in epath])
        if epath[0].source() != global_src:
            # generate path from global src to epath 1st node
            s = global_src
            t = epath[0].source()
            sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
            strain.extend(sp)
            strain.append(t)
        for i in range(len(epath) - 1):
            s = epath[i].target()
            t = epath[i+1].source()
            if s == t:
                strain.append(s)
            else:
                # find intermediate path between s_e to t_e
                sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
                strain.append(s)
                strain.extend(sp)
                strain.append(t)
        if epath[-1].target() != global_sink:
            # generate path from epath last node to global sink
            s = epath[-1].target()
            t = global_sink
            sp, plen, pmark = dijkstra_sp_v3(graph, s, t, ccov, threshold, overlap, lb)
            strain.append(s)
            strain.extend(sp)
        print("Strain - ccov: {0}, {1}".format(ccov, path_to_id_string(graph, strain)))
        # strain_dict[itercount + str(i)] = [[graph.vp.id[n] for n in strain], path_len(graph, strain, overlap), ccov]
        # for node in strain:
        #     eval_score(graph.vp.id[node], graph.vp.dp[node] - ccov, threshold, lb, True)
        # graph_reduction_c(graph, strain, ccov)
        for e in epath:
            strain_dict[graph.vp.id[e.source()]] = [[graph.vp.id[e.source()]], path_len(graph, [e.source()], overlap), graph.vp.dp[e.source()]]
            strain_dict[graph.vp.id[e.target()]] = [[graph.vp.id[e.target()]], path_len(graph, [e.target()], overlap), graph.vp.dp[e.target()]]
        break

    return

# split contig
def split_contig(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, threshold):
    """
    split contig out
    """
    dup_record = {}
    reduced_edges = []
    reduced_vertices = []
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
            reduced_edges.append(ie)           
            if graph.ep.flow[ie] <= threshold:
                graph_remove_edge(graph, simp_edge_dict, graph.vp.id[ie.source()], graph.vp.id[ie.target()])

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
                reduced_edges.append(oe)
                if graph.ep.flow[oe] <= threshold:
                    graph_remove_edge(graph, simp_edge_dict, graph.vp.id[oe.source()], graph.vp.id[oe.target()])
                graph_add_edge(graph, simp_edge_dict, new_vertex, graph.vp.id[new_vertex], 
                    simp_node_dict[next_no], next_no, graph.ep.overlap[ie], ccov) 
            
            reduced_vertices.append(old_vertex)
            if graph.vp.dp[old_vertex] <= threshold:
                graph_remove_vertex(graph, simp_node_dict, no)
            contig_dict[cno][0][i] = graph.vp.id[new_vertex]
            dup_record[no] += 1
    # print("DELTA: ", threshold)
    # for e in sorted(graph.edges(), key=lambda ee: graph.ep.flow[ee]):
    #     if e in reduced_edges:
    #         print("contig split")
    #     else:
    #         print("original")
    #     print_edge(graph, e, "lower than delta" if graph.ep.flow[e] < threshold 
    #         else "greater than delta")
    
    # for v in sorted(graph.vertices(), key=lambda vv: graph.vp.dp[vv]):
    #     if v in reduced_vertices:
    #         print("contig split")
    #     else:
    #         print("original")
    #     print_vertex(graph, v, "lower than delta" if graph.vp.dp[v] < threshold 
    #         else "greater than delta")
    return

def local_search_optimisation(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
overlap, threshold):   
    """
    need to check whether current path flow below maximum flow 
    """
    
    usage_dict = {}
    # init status to all the node
    node_to_contig_dict, _ = contig_map_node(contig_dict)

    full_count = 0
    partial_count = 0
    over_count = 0
    for no, node in simp_node_dict.items():
        if no not in node_to_contig_dict:
            usage_dict[no] = [0, graph.vp.dp[node], 'partial']
            partial_count += 1
        else:
            sum_used = 0
            for cno in node_to_contig_dict[no]:
                sum_used += contig_dict[cno][2]
            if sum_used - graph.vp.dp[node] > threshold:
                # overused
                usage_dict[no] = [sum_used, graph.vp.dp[node], 'over']
                over_count += 1
            else:
                if abs(sum_used - graph.vp.dp[node]) < threshold:
                    usage_dict[no] = [sum_used, graph.vp.dp[node], 'full']
                    full_count += 1
                else:
                    usage_dict[no] = [sum_used, graph.vp.dp[node], 'partial']
                    partial_count += 1
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(usage_dict[no][0]), round(usage_dict[no][1]), usage_dict[no][2]))

    print("overflow nodes: ", over_count)
    print("full nodes: ", full_count)
    print("partial used nodes: ", partial_count)
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    swap = True
    itercount = 0
    while swap:
        swap = False
        itercount += 1
        print("iteration: ", itercount)
        for sno, snode in simp_node_dict.items():
            for eno, enode in simp_node_dict.items():
                # check if sno & eno adjacent
                if sno == eno or snode in enode.all_neighbors():
                    continue
                # not adjacent nodes, check any bounded over-path and partial-path exist
                over_paths, partial_paths = path_replacement_account(graph, simp_node_dict, simp_edge_dict, usage_dict, snode, enode)
                if len(over_paths) <= 0 or len(partial_paths) <= 0:
                    continue
                print("-----*curr s: {0} t: {1}".format(sno, eno))
                all_cnos = set()
                all_paths = {}
                prev_cno_mapping = {}
                for i, op in enumerate(over_paths):
                    min_bound = min([graph.vp.dp[n] for n in op])
                    print(path_to_id_string(graph, op, "over cap: {0}".format(min_bound)))

                    cno_set = None
                    for n in op:
                        if graph.vp.id[n] in node_to_contig_dict:
                            acc_set = set(node_to_contig_dict[graph.vp.id[n]])
                            cno_set = acc_set if cno_set == None else cno_set.intersection(acc_set)
                        else:
                            cno_set = set()
                            break
                    print("involved strain: ", [(cno, contig_dict[cno][2]) for cno in cno_set])
                    all_cnos = all_cnos.union(cno_set)
                    all_paths[('o', i)] = (op, min_bound)
                    for cno in cno_set:
                        prev_cno_mapping[cno] = ('o', i)

                for i, pp in enumerate(partial_paths):
                    min_bound = min([graph.vp.dp[n] for n in pp])
                    print(path_to_id_string(graph, pp, "partial cap: {0}".format(min_bound)))
                    cno_set = None
                    for n in pp:
                        if graph.vp.id[n] in node_to_contig_dict:
                            acc_set = set(node_to_contig_dict[graph.vp.id[n]])
                            cno_set = acc_set if cno_set == None else cno_set.intersection(acc_set)
                        else:
                            cno_set = set()
                            break
                    print("involved strain: ", [(cno, contig_dict[cno][2]) for cno in cno_set])
                    all_cnos = all_cnos.union(cno_set)
                    all_paths[('p', i)] = (pp, min_bound)
                    for cno in cno_set:
                        prev_cno_mapping[cno] = ('p', i)
                
                all_comb = list(product(list(all_paths.keys()), repeat=len(all_cnos)))
                all_cnos = list(all_cnos)

                prev_state_dict = {}
                for [path, _] in all_paths.values():
                    for n in path:
                        if graph.vp.id[n] not in prev_state_dict:
                            prev_state_dict[graph.vp.id[n]] = copy.deepcopy(usage_dict[graph.vp.id[n]])
                prev_overcount = len([k for k, [_, _, s] in prev_state_dict.items() if s == 'over'])
                print([(cno, contig_dict[cno][2]) for cno in all_cnos])
                optim_score = None
                optim_overcount = None
                optim_comb = []
                for i, comb in enumerate(all_comb):
                    # select the best combination that provides optimal capacity usage
                    comb_subscores = {}
                    for key, [_, cap] in all_paths.items():
                        comb_subscores[key] = [0, cap]
                    for j in range(len(all_cnos)):
                        flow = contig_dict[all_cnos[j]][2]
                        comb_subscores[comb[j]][0] += flow
                    #FIXME how to assign the score for one comb, include the check on over-node reduction
                    score_flow = pow(sum([flow/cap for [flow, cap] in comb_subscores.values()]) - len(comb_subscores), 2)
                    print("comb: ", [(path_to_id_string(graph, all_paths[key][0], str(all_paths[key][1])), "cno: " + all_cnos[i], contig_dict[all_cnos[i]][2]) for i, key in enumerate(comb)])
                    new_state_dict = test_comb(graph, comb, usage_dict, all_cnos, all_paths, prev_cno_mapping, contig_dict, threshold)
                    score_overcount = len([k for k, [_, _, s] in new_state_dict.items() if s == 'over'])
                    print("prev states: ", prev_state_dict)
                    print("new states: ", new_state_dict)
                    print("previous over node count: ", len([k for k, [_, _, s] in prev_state_dict.items() if s == 'over']))
                    print("new over node count: ", score_overcount)

                    if optim_score == None:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount
                    elif score_overcount < optim_overcount:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount 
                    elif score_overcount == optim_overcount and score_flow < optim_score:
                        optim_score = score_flow
                        optim_comb = comb
                        optim_overcount = score_overcount                     

                print("Optim comb: ", [(path_to_id_string(graph, all_paths[key][0], str(all_paths[key][1])), "cno: " + all_cnos[i], contig_dict[all_cnos[i]][2]) for i, key in enumerate(optim_comb)])
                print("Optim score: ", optim_score)
                print("Optim overnode: ", optim_overcount)
                if prev_overcount <= optim_overcount:
                    print("no more over usage node be reduced")
                    continue
                swap = True
                # replacement
                for i, key in enumerate(optim_comb):
                    # update contig dict
                    cno = all_cnos[i]
                    prev_p = [graph.vp.id[n] for n in all_paths[prev_cno_mapping[cno]][0]]
                    new_p = [graph.vp.id[n] for n in all_paths[key][0]]
                    print("prev: {0}, {1}, {2}".format(cno, list_to_string(contig_dict[cno][0]), list_to_string(prev_p)))
                    contig = contig_replacement_c(contig_dict[cno][0], prev_p, new_p)
                    print("after: {0}, {1}, {2}".format(cno, list_to_string(contig), list_to_string(new_p)))
                    clen = path_len(graph, contig, overlap)
                    contig_dict[cno] = [contig, clen, contig_dict[cno][2]]
                    # update usage dict
                    #TODO
                    for prev_id in prev_p:
                        # if prev_id in new_p:
                        #     continue
                        usage_dict[prev_id][0] -= contig_dict[cno][2]
                        prev_state = usage_dict[prev_id][2]
                        usage_dict[prev_id][2] = node_status(usage_dict[prev_id][0], usage_dict[prev_id][1], threshold)
                        print("prev-state changed: {0} {1}->{2}".format(prev_id, prev_state, usage_dict[prev_id][2]))
                    for new_id in new_p:
                        # if new_id in prev_p:
                        #     continue
                        usage_dict[new_id][0] += contig_dict[cno][2]
                        prev_state = usage_dict[new_id][2]
                        usage_dict[new_id][2] = node_status(usage_dict[new_id][0], usage_dict[new_id][1], threshold)
                        print("new-state changed: {0} {1}->{2}".format(new_id, prev_state, usage_dict[new_id][2]))
                # update the node to contig mapping
                node_to_contig_dict, _ = contig_map_node(contig_dict)
        break
    over_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'over'])
    full_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'full'])
    partial_count = len([n for n in usage_dict.keys() if usage_dict[n][2] == 'partial'])
    print("overflow nodes: ", over_count)
    print("full nodes: ", full_count)
    print("partial used nodes: ", partial_count)
    for no, [used, cap, state] in usage_dict.items():
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(used), round(cap), state))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    return usage_dict

def final_strain_extraction(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, usage_dict: dict, threshold, overlap):
    """
    extract rest of partial used node path from global source to sink if any
    """
    print("--------------Start strain final extraction------------------")
    global_src = simp_node_dict['global_src']
    global_sink = simp_node_dict['global_sink']

    for no, [_, _, state] in usage_dict.items():
        if state == 'full' or state == 'over':
            node = simp_node_dict[no]
            graph.vp.color[node] = 'gray'
            for edge in node.all_edges():
                graph.ep.color[edge] = 'gray'
    
    path = []
    curr_id = 0
    while path != None:
        path, plen = dijkstra_sp_v2(graph, simp_node_dict, global_src, global_sink, overlap)
        if path != None:
            # FIXME
            cflows = contig_flow(graph, simp_edge_dict, [graph.vp.id[n] for n in path])
            redcov = numpy.min(cflows)
            for node in path:
                node_id = graph.vp.id[node]
                usage_dict[node_id][0] += redcov
                usage_dict[node_id][2] = node_status(usage_dict[node_id][0], 
                    usage_dict[node_id][1], threshold)
                if usage_dict[node_id][2] == 'full' or usage_dict[node_id][2] == 'over':
                    graph.vp.color[node] = 'gray'
                    for edge in node.all_edges():
                        graph.ep.color[edge] = 'gray'
            contig_dict["st" + str(curr_id)] = [[graph.vp.id[n] for n in path], plen, redcov]
            curr_id += 1
        else:
            print("Path not found")
    for no, [used, cap, state] in usage_dict.items():
        print("node: {0}, curr usage: {1}, capacity: {2}, status: {3}".format(
            no, round(used), round(cap), state))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'over'], "OVER: "))
    print(list_to_string([k for k, [_, _, s] in usage_dict.items() if s == 'partial'], "PAR: "))
    # # double check usage
    # post_usages = {}
    # post_free_nodes = []
    # post_partial_used_nodes = []
    # for no, node in simp_node_dict.items():
    #     print("----------------------------------------------------")
    #     if no == 'global_src' or no == 'global_sink':
    #         continue
    #     ratio = round(((graph.vp.dp[node] - graph.vp.udp[node]) * 100 / graph.vp.dp[node]), 2)
    #     print("Node: {0}, full: {1}, left: {2}, usage: {3}, {4}".format(no, round(graph.vp.dp[node], 2), round(graph.vp.udp[node], 2), ratio, graph.vp.color[node]))
    #     if ratio < 100 and graph.vp.udp[node] > threshold:
    #         post_partial_used_nodes.append(no)
    #     if ratio <= 0:
    #         post_free_nodes.append(no)
    #     post_usages[no] = ratio
    # overall_post_usage = numpy.mean(list(post_usages.values()))
    # print("Free nodes: ", list_to_string(post_free_nodes))
    # print("Partial used nodes: ", list_to_string(post_partial_used_nodes))
    # print("Overall Usage Post: {0}".format(overall_post_usage))
    # # do some plotting
    # df = pandas.DataFrame(
    #     {'Id': [i for i in range(len(pre_usages.keys()))],
    #     'pre': pre_usages.values(),
    #     'post': post_usages.values()
    #     })
    # tidy = df.melt(id_vars='Id').rename(columns=str.title)
    # ax = seaborn.barplot(x='Id', y='Value', hue='Variable', data=tidy)
    # for container in ax.containers:
    #     ax.bar_label(container)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right")
    # plt.title('Bar plot of Node Usage')
    # plt.savefig("{0}barplot_usage_post.png".format(TEMP_DIR))
    # # post process
    # for node in graph.vertices():
    #     graph.vp.dp[node] = graph.vp.udp[node]
    # simp_node_dict.pop(graph.vp.id[global_src])
    # graph.vp.color[global_src] = 'gray'
    # simp_node_dict.pop(graph.vp.id[global_sink])
    # graph.vp.color[global_sink] = 'gray'
    return

def strain_extension(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, contig_dict: dict, 
min_cov, max_len, threshold, overlap):
    """
    Extend the strain length
    1. map all the strain back into the graph
    """
    print("--------------Start strain extension------------------")
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap)

    # extend all the contig from both end
    for cno, [contig, clen, ccov] in list(contig_dict.items()):
        print_contig(cno, clen, ccov, contig, "-----> current extending contig: ")
        if simp_node_dict[contig[0]].in_degree() == 0 and simp_node_dict[contig[-1]].out_degree() == 0:
            print("not extensible")
        elif simp_node_dict[contig[0]] in list(simp_node_dict[contig[-1]].out_neighbors()):
            print("self cycle, not entensible")
        else:
            print("extensible")
            contig_head = simp_node_dict[contig[0]]
            if contig_head.in_degree() != 0 and global_src not in contig_head.in_neighbors():
                print("---> {0} ~~> contig head: {1}".format(graph.vp.id[global_src], graph.vp.id[contig_head]))
                sp, _, _ = dijkstra_sp(graph, global_src, contig_head, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [[graph.vp.id[n] for n in sp] + contig_dict[cno][0], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")

            contig_tail = simp_node_dict[contig[-1]]  
            if contig_tail.out_degree() != 0 and global_sink not in contig_tail.out_neighbors():
                print("---> contig tail: {0} ~~> {1}".format(graph.vp.id[contig_tail], graph.vp.id[global_sink]))
                sp, _, _ = dijkstra_sp(graph, contig_tail, global_sink, ccov, threshold, overlap)
                if sp != None:
                    plen = path_len(graph, sp, overlap)
                    print("path len: ", plen)
                    contig_dict[cno] = [contig_dict[cno][0] + [graph.vp.id[n] for n in sp], plen + contig_dict[cno][1] - overlap, ccov]
                else:
                    print("path not found")
            
            # # re-assign the strain cov to min flow among the contig
            # assert len(contig_dict[cno][0]) >= 1
            # if len(contig_dict[cno][0]) == 1:
            #     redcov = graph.vp.dp[contig_dict[cno][0][0]]
            # else:
            #     redcov = numpy.min(contig_flow(graph, simp_edge_dict, contig_dict[cno][0]))
            contig_dict[cno][2] = ccov
    remove_global_source_sink(graph, global_src, global_sink)
    return contig_dict

def contig_replacement_c(contig: list, a, b):
    rtn = []
    rtn.extend(contig[0:contig.index(a[0])])
    rtn.extend(b)
    next = contig.index(a[-1]) + 1
    if next < len(contig):
        rtn.extend(contig[next:])
    return rtn

def minimal_bubble_detection(graph: Graph):
    # retrieve all the minimal bubble
    # every bubble is structure as (start branch -> single node child -> end branch)
    branches = retrieve_branch(graph)
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

def retrieve_branch(graph: Graph):
    branches = {}
    for node in graph.vertices():
        if node.in_degree() > 1 or node.out_degree() > 1:
            branches[graph.vp.id[node]] = node
    return branches

def contig_replacement(contig: list, a, b):
    """
    replace all occurance of a from contig to b
    """
    return [(n if n != a else b) for n in contig]

def get_row(matrix: list, rowId):
    return matrix[rowId]

def get_col(matrix: list, colId):
    return [row[colId] for row in matrix]

def copy_dict_hard(d: dict):
    """
    hard copy a dict with list element, key is immutable. legacy
    """
    rd = {}
    for k, v in d.items():
        rd[k] = v.copy()
    return rd

def similarity(cov1, cov2):
    return 1 - (abs(cov1-cov2)/(cov1+cov2))

def similarity_e(e, cliq_graph):
    return similarity(cliq_graph.vp.ccov[e.source()], cliq_graph.vp.ccov[e.target()])

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

def graph_recolor_edge(graph: Graph, simp_edge_dict: dict, src, src_id, tgt, tgt_id, s="edge recolor", color='black'):
    """
    recolour edge to black as default
    """
    edge = graph.edge(src, tgt)
    graph.ep.color[edge] = color
    simp_edge_dict[(src_id, tgt_id)] = edge
    print_edge(graph, edge, s)
    return

def min_flow_edge(graph: Graph, edge_dict: dict, contig):
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
        # print("No intersection")
        return False, intersect, 0, 'n'

    if len(intersect) == len(contig):
        print("{0} is covered by {1}".format(cno, cno2))
        return True, intersect, -1, 'o'

    if len(intersect) == len(contig2):
        print("{0} is covered by {1}".format(cno2, cno))
        return True, intersect, -1, 'o'

    intermediate_nodes_index = [False for _ in contig]
    for i in [contig.index(e) for e in intersect]:
        intermediate_nodes_index[i] = True
    print("cno: {0} intersect at: {1}".format(cno, list(enumerate(intermediate_nodes_index))))
    if not intermediate_nodes_index[0] and not intermediate_nodes_index[-1]:
        print("intersection in the middle")
        return True, intersect, -1, 'o'
    prev_false_index = intermediate_nodes_index.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index)):
        if not intermediate_nodes_index[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                print("intersection in the middle")
                return True, intersect, -1, 'o'

    intermediate_nodes_index2 = [False for _ in contig2]
    for i in [contig2.index(e) for e in intersect]:
        intermediate_nodes_index2[i] = True
    print("cno2: {0} intersect at: {1}".format(cno2, list(enumerate(intermediate_nodes_index2))))
    if not intermediate_nodes_index2[0] and not intermediate_nodes_index2[-1]:
        print("intersection in the middle")
        return True, intersect, -1, 'o'
    prev_false_index = intermediate_nodes_index2.index(False)
    for j in range(prev_false_index + 1, len(intermediate_nodes_index2)):
        if not intermediate_nodes_index2[j]:
            if prev_false_index + 1 == j:
                prev_false_index = j
            else:
                print("intersection in the middle")
                return True, intersect, -1, 'o'
    cend = 0
    direction = None
    if intermediate_nodes_index[0]:
        cend += 1
        direction = 'b'
    if intermediate_nodes_index[-1]:
        cend += 1
        direction = 'f' if direction == None else 'd'
    print("overlap end-to-end, tolarent, cend: ", cend)
    return False, intersect, cend, direction

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

def contig_split(graph: Graph, cno, contig: list, simp_node_dict: dict, simp_edge_dict: dict, overlap, min_node=2):
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
    print("Bubbles: ", path_to_id_string(graph, bubbles))
    return bubbles

def gen_bubble_detection(graph: Graph, simp_node_dict: dict, simp_edge_dict: dict, overlap):
    def print_branch_queue(branch_queue):
        for close, ind, outd in branch_queue:
            print("|{0}, in: {1}, out:{2}|".format(graph.vp.id[close], ind, outd))
    global_src, global_sink = add_global_source_sink(graph, simp_node_dict, simp_edge_dict, overlap)
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

def path_replacement_account(graph: Graph, usage_dict: dict, s, t):
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
    #         print("Length: ", plen, path_to_id_string(graph, p, "-Path be found: "))

    return rtn_paths

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
        print(path_to_id_string(graph, sp, "SP be found: "))
        mark = 1/(1+mark**0.5)
        print("plen: ", path_len(graph, sp[1:-1], overlap), "pmark: ", mark)

        return sp[1:-1], path_len(graph, sp[1:-1], overlap), mark

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
        print(path_to_id_string(graph, sp, "SP be found: "))
        print("plen: ", path_len(graph, sp[1:-1], overlap))

        return sp[1:-1], path_len(graph, sp[1:-1], overlap)

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

def dijkstra_sp_v3(graph: Graph, source, sink, closest_cov, b0, b1, overlap: int):
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
                threshold = abs(b0 + b1 * edge_flow)
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