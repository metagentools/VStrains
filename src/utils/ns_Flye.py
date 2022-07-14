#!/usr/bin/env python3

from utils.ns_Preprocess import delta_estimation, graph_simplification, tip_removal_s
from utils.ns_CovBalance import coverage_rebalance_s, assign_edge_flow
from utils.ns_Path import extract_cand_path
from utils.ns_Split import iterated_graph_split
from utils.ns_IO import (
    flye_info2path, 
    graph_to_gfa, 
    flipped_gfa_to_graph, 
    get_contig, gfa_to_graph, 
    contig_dict_to_path, 
    contig_dict_to_fasta,
)

from utils.ns_Utilities import *

def run(args):
    TEMP_DIR = args.output_dir

    print("----------------------------------INPUT---------------------------------------")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file)
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "{0}/gfa/graph_L0.gfa".format(TEMP_DIR))
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph("{0}/gfa/graph_L0.gfa".format(TEMP_DIR))
    path_file = flye_info2path(args.info_file, TEMP_DIR)
    contig_dict = get_contig(path_file, simp_node_dict0, simp_edge_dict0, args.min_len, True)
    print("----------------------------------TIP REMOVAL---------------------------------------")
    tip_removal_s(graph0, simp_node_dict0, contig_dict, TEMP_DIR)
    graph_to_gfa(graph0, simp_node_dict0, simp_edge_dict0, "{0}/gfa/t_graph_L1.gfa".format(TEMP_DIR))
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph("{0}/gfa/t_graph_L1.gfa".format(TEMP_DIR))

    print("-----------------------------DELTA ESTIMATION-----------------------------")
    b0, b1 = delta_estimation(graph1, TEMP_DIR)

    print("-------------------------------GRAPH SIMPLIFICATION & REBALANCE-----------------------------------")
    if args.min_cov > 0:
        THRESHOLD = args.min_cov
        print("user-defined node minimum coverage: ", THRESHOLD)
    else:
        THRESHOLD = 0.05 * numpy.median([graph1.vp.dp[node] for node in graph1.vertices()])
        print("computed node minimum coverage: ", THRESHOLD)

    graph_simplification(graph1, simp_node_dict1, simp_edge_dict1, contig_dict, THRESHOLD)

    graph_to_gfa(graph1, simp_node_dict1, simp_edge_dict1, "{0}/gfa/st_graph_L2.gfa".format(TEMP_DIR))
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph("{0}/gfa/st_graph_L2.gfa".format(TEMP_DIR))

    coverage_rebalance_s(graph2, simp_node_dict2, simp_edge_dict2, TEMP_DIR, True)

    graph_to_gfa(graph2, simp_node_dict2, simp_edge_dict2, "{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR))
    # graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph("{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR))
    # expectation_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)
    
    # print("-------------------------------CONTIG COVERAGE REBALANCE-----------------------------------")
    # contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict)
    
    # # stat evaluation
    # contig_dict_to_path(contig_dict, "{0}/tmp/pre_contigs.paths".format(TEMP_DIR))
    # contig_dict_to_fasta(graph3, simp_node_dict3, contig_dict, "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR))
    # if args.ref_file:
    #     map_ref_to_graph(args.ref_file, simp_node_dict3, "{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR), False, "{0}/paf/node_to_ref.paf".format(TEMP_DIR), 
    #         "{0}tmp/temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR))
    #     minimap_api(args.ref_file, "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR), "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR))
    #     map_ref_to_contig(contig_dict, "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR))
    # # end stat

    # print("-----------------------GRAPH BRANCH SPLIT & COMPACTIFICATION-------------------------------")
    # graph5, simp_node_dict5, simp_edge_dict5 = iterated_graph_split(graph3, simp_node_dict3, simp_edge_dict3, contig_dict, TEMP_DIR, b0, b1, THRESHOLD)

    # # stat evaluation
    # contig_dict_to_path(contig_dict, "{0}/tmp/post_contigs.paths".format(TEMP_DIR))
    # contig_dict_to_fasta(graph5, simp_node_dict5, contig_dict, "{0}/tmp/post_contigs.fasta".format(TEMP_DIR))
    # if args.ref_file:
    #     map_ref_to_graph(args.ref_file, simp_node_dict5, "{0}/gfa/rbsdt_graph_L5.gfa".format(TEMP_DIR), False, "{0}/paf/node_to_ref_red.paf".format(TEMP_DIR), 
    #         "{0}/tmp/temp_gfa_to_fasta.fasta".format(TEMP_DIR))
    #     minimap_api(args.ref_file, "{0}/tmp/post_contigs.fasta".format(TEMP_DIR), "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR))
    #     map_ref_to_contig(contig_dict, "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR))
    # # end stat

    # print("-----------------------CONTIG PATH EXTENSION-------------------------------")
    # strain_dict = extract_cand_path(graph5, simp_node_dict5, simp_edge_dict5, contig_dict, b0, b1, THRESHOLD)
    
    # print("-----------------------FINAL CLEAN UP-------------------------------")
    # contig_dup_removed_s(strain_dict)  
    # trim_contig_dict(graph5, simp_node_dict5, strain_dict)  
    # contig_dict_to_fasta(graph5, simp_node_dict5, strain_dict, "{0}/strain.fasta".format(TEMP_DIR))
    # contig_dict_to_path(strain_dict, "{0}/strain.paths".format(TEMP_DIR), True)
    # if args.ref_file:
    #     minimap_api(args.ref_file, "{0}/strain.fasta".format(TEMP_DIR), "{0}/paf/strain_to_ref.paf".format(TEMP_DIR))

    return 0