#!/usr/bin/env python3

from utils.vsAware_Path import extract_cand_path
from utils.vsAware_Split import iterated_graph_split
from utils.vsAware_Utilities import *
from utils.vsAware_Preprocess import (
    delta_estimation,
    graph_simplification,
    tip_removal_s,
    reindexing,
)
from utils.vsAware_CovBalance import coverage_rebalance_s, assign_edge_flow
from utils.vsAware_IO import (
    graph_to_gfa,
    flipped_gfa_to_graph,
    gfa_to_graph,
    contig_dict_to_path,
    contig_dict_to_fasta,
    flye_info_parser,
)


def run(args, logger):
    TEMP_DIR = args.output_dir

    logger.info("vsAware-Flye started")

    logger.info(">>>STAGE: parsing graph and contigs")
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file, logger)
    graph_to_gfa(
        graph,
        simp_node_dict,
        simp_edge_dict,
        logger,
        "{0}/gfa/graph_L0.gfa".format(TEMP_DIR),
    )
    graph0, simp_node_dict0, simp_edge_dict0 = flipped_gfa_to_graph(
        "{0}/gfa/graph_L0.gfa".format(TEMP_DIR), logger
    )
    graph0, simp_node_dict0, simp_edge_dict0, idx_mapping = reindexing(
        graph0, simp_node_dict0, simp_edge_dict0
    )

    contig_dict, contig_info = flye_info_parser(
        graph0,
        simp_node_dict0,
        simp_edge_dict0,
        idx_mapping,
        logger,
        args.info_file,
        args.min_len,
    )
    copy_contig_dict = {}
    for cno, [contig, clen, ccov] in contig_dict.items():
        copy_contig_dict[cno] = [list(contig), clen, ccov]

    logger.info(">>>STAGE: preprocess")
    tip_removal_s(graph0, simp_node_dict0, contig_dict, logger, TEMP_DIR)
    graph_to_gfa(
        graph0,
        simp_node_dict0,
        simp_edge_dict0,
        logger,
        "{0}/gfa/t_graph_L1.gfa".format(TEMP_DIR),
    )
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph(
        "{0}/gfa/t_graph_L1.gfa".format(TEMP_DIR), logger
    )

    b0, b1 = delta_estimation(graph1, logger, TEMP_DIR)

    if args.min_cov > 0:
        THRESHOLD = args.min_cov
        logger.info("user-defined node minimum coverage: {0}", format(THRESHOLD))
    else:
        THRESHOLD = 0.05 * numpy.median(
            [graph1.vp.dp[node] for node in graph1.vertices()]
        )
        logger.info("computed node minimum coverage: {0}".format(THRESHOLD))

    graph_simplification(
        graph1, simp_node_dict1, simp_edge_dict1, contig_dict, logger, THRESHOLD
    )
    graph_to_gfa(
        graph1,
        simp_node_dict1,
        simp_edge_dict1,
        logger,
        "{0}/gfa/st_graph_L2.gfa".format(TEMP_DIR),
    )
    graph2, simp_node_dict2, simp_edge_dict2 = flipped_gfa_to_graph(
        "{0}/gfa/st_graph_L2.gfa".format(TEMP_DIR), logger
    )

    coverage_rebalance_s(graph2, simp_node_dict2, simp_edge_dict2, logger)
    graph_to_gfa(
        graph2,
        simp_node_dict2,
        simp_edge_dict2,
        logger,
        "{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR),
    )
    graph3, simp_node_dict3, simp_edge_dict3 = flipped_gfa_to_graph(
        "{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR), logger
    )
    assign_edge_flow(graph3, simp_node_dict3, simp_edge_dict3)

    contig_cov_fix(graph3, simp_node_dict3, simp_edge_dict3, contig_dict, None)

    contig_dict_to_path(contig_dict, "{0}/tmp/pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(
        graph3,
        simp_node_dict3,
        contig_dict,
        "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR),
    )
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(
            args.ref_file,
            simp_node_dict3,
            "{0}/gfa/cst_graph_L3.gfa".format(TEMP_DIR),
            False,
            "{0}/paf/node_to_ref.paf".format(TEMP_DIR),
            "{0}/tmp/temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR),
        )
        minimap_api(
            args.ref_file,
            "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR),
            "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR),
        )
        map_ref_to_contig(
            contig_dict, "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR)
        )
    # end stat

    logger.info(">>>STAGE: graph branch split & compactification")
    graphf, simp_node_dictf, simp_edge_dictf = iterated_graph_split(
        graph3,
        simp_node_dict3,
        simp_edge_dict3,
        contig_dict,
        logger,
        TEMP_DIR,
        b0,
        b1,
        THRESHOLD,
    )

    contig_dict_to_path(contig_dict, "{0}/tmp/post_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(
        graphf,
        simp_node_dictf,
        contig_dict,
        "{0}/tmp/post_contigs.fasta".format(TEMP_DIR),
    )
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(
            args.ref_file,
            simp_node_dictf,
            "{0}/gfa/rbsdt_graph_L5.gfa".format(TEMP_DIR),
            False,
            "{0}/paf/node_to_ref_red.paf".format(TEMP_DIR),
            "{0}/tmp/temp_gfa_to_fasta.fasta".format(TEMP_DIR),
        )
        minimap_api(
            args.ref_file,
            "{0}/tmp/post_contigs.fasta".format(TEMP_DIR),
            "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR),
        )
        map_ref_to_contig(
            contig_dict, "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR)
        )
    # end stat

    logger.info(">>>STAGE: contig path extension")
    strain_dict = extract_cand_path(
        graphf, simp_node_dictf, simp_edge_dictf, contig_dict, logger, b0, b1, THRESHOLD
    )

    logger.info(">>>STAGE: final process")
    contig_dup_removed_s(strain_dict, logger)
    trim_contig_dict(graphf, simp_node_dictf, strain_dict, logger)
    # recover repeat nodes back to contig
    strain_repeat_resol(
        graph0, simp_node_dict0, strain_dict, contig_info, copy_contig_dict, logger
    )

    logger.info(">>>STAGE: generate result")
    contig_dict_to_fasta(
        graph0, simp_node_dict0, strain_dict, "{0}/strain.fasta".format(TEMP_DIR)
    )
    contig_dict_to_path(strain_dict, "{0}/strain.paths".format(TEMP_DIR), True)
    if args.ref_file:
        minimap_api(
            args.ref_file,
            "{0}/strain.fasta".format(TEMP_DIR),
            "{0}/paf/strain_to_ref.paf".format(TEMP_DIR),
        )
    logger.info("vsAware-Flye finished")
    return 0
