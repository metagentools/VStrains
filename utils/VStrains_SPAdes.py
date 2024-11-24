#!/usr/bin/env python3

from utils.VStrains_Utilities import *
from utils.VStrains_Preprocess import (
    graph_simplification,
    reindexing,
    threshold_estimation,
)
from utils.VStrains_IO import (
    graph_to_gfa,
    flipped_gfa_to_graph,
    gfa_to_graph,
    contig_dict_to_path,
    contig_dict_to_fasta,
    spades_paths_parser,
    process_pe_info,
    store_reinit_graph,
)
from utils.VStrains_Decomposition import *
from utils.VStrains_Extension import path_extension, best_matching
import os
import sys


def run(args, logger):
    TEMP_DIR = args.output_dir

    logger.info("VStrains-SPAdes started")

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
    graph_to_gfa(
        graph0,
        simp_node_dict0,
        simp_edge_dict0,
        logger,
        "{0}/gfa/graph_L0r.gfa".format(TEMP_DIR),
    )

    # cut-off coverage, graph preprocess parameter
    THRESHOLD = 0
    if args.min_cov != None:
        THRESHOLD = args.min_cov
        logger.info("user-defined node minimum coverage: {0}".format(THRESHOLD))
    else:
        THRESHOLD = threshold_estimation(graph0, logger, TEMP_DIR)
        logger.info("computed node minimum coverage: {0}".format(THRESHOLD))

    contig_dict, contig_info = spades_paths_parser(
        graph0,
        simp_node_dict0,
        simp_edge_dict0,
        idx_mapping,
        logger,
        args.path_file,
        args.min_len,
        THRESHOLD,
    )
    copy_contig_dict = {}
    for cno, [contig, clen, ccov] in contig_dict.items():
        copy_contig_dict[cno] = [list(contig), clen, ccov]
    # debug only
    contig_dict_to_path(contig_dict, "{0}/tmp/init_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(
        graph0,
        simp_node_dict0,
        contig_dict,
        "{0}/tmp/init_contigs.fasta".format(TEMP_DIR),
    )
    if args.ref_file:
        minimap_api(
            args.ref_file,
            "{0}/tmp/init_contigs.fasta".format(TEMP_DIR),
            "{0}/paf/init_contigs_to_strain.paf".format(TEMP_DIR),
        )
    # debug only
    logger.info(">>>STAGE: preprocess")
    graph_simplification(
        graph0, simp_node_dict0, simp_edge_dict0, None, logger, THRESHOLD
    )
    graph_to_gfa(
        graph0,
        simp_node_dict0,
        simp_edge_dict0,
        logger,
        "{0}/gfa/s_graph_L1.gfa".format(TEMP_DIR),
    )
    graph1, simp_node_dict1, simp_edge_dict1 = flipped_gfa_to_graph(
        "{0}/gfa/s_graph_L1.gfa".format(TEMP_DIR), logger
    )

    # filter out contig that contains erroroness nodes
    for cno, [contig, _, _] in list(contig_dict.items()):
        if any([c not in simp_node_dict1 for c in contig]):
            contig_dict.pop(cno)
            logger.debug("unreliable contig with low coverage: {0}".format(cno))

    # get graph kmer size
    ksize = graph1.ep.overlap[list(graph1.edges())[0]] if graph1.num_edges() > 0 else 0
    logger.info("graph kmer size: {0}".format(ksize))
    if ksize <= 0:
        logger.error("invalid kmer-size, the graph does not contain any edges, exit..")
        sys.exit(1)

    # obtain paired end information
    script_path = "{0}/VStrains_PE_Inference.py".format(
        os.path.abspath(os.path.dirname(__file__))
    )
    subprocess.check_call(
        "python {0} -g {1} -o {2} -f {3} -r {4} -k {5}".format(
            script_path,
            "{0}/gfa/s_graph_L1.gfa".format(TEMP_DIR),
            "{0}/aln".format(TEMP_DIR),
            args.fwd,
            args.rve,
            ksize,
        ),
        shell=True,
    )
    logger.info("paired end information stored")
    pe_info_file = "{0}/aln/pe_info".format(TEMP_DIR)
    st_info_file = "{0}/aln/st_info".format(TEMP_DIR)
    pe_info, dcpy_pe_info = process_pe_info(
        simp_node_dict1.keys(), pe_info_file, st_info_file
    )

    edge_cleaning(graph1, simp_edge_dict1, contig_dict, pe_info, logger)

    graph2, simp_node_dict2, simp_edge_dict2 = store_reinit_graph(
        graph1,
        simp_node_dict1,
        simp_edge_dict1,
        logger,
        "{0}/gfa/es_graph_L2.gfa".format(TEMP_DIR),
    )

    contig_dict_to_path(contig_dict, "{0}/tmp/pre_contigs.paths".format(TEMP_DIR))
    contig_dict_to_fasta(
        graph2,
        simp_node_dict2,
        contig_dict,
        "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR),
    )
    # stat evaluation
    if args.ref_file:
        map_ref_to_graph(
            args.ref_file,
            simp_node_dict2,
            "{0}/gfa/es_graph_L2.gfa".format(TEMP_DIR),
            logger,
            True,
            "{0}/paf/node_to_ref.paf".format(TEMP_DIR),
            "{0}/tmp/temp_gfa_to_fasta_pre.fasta".format(TEMP_DIR),
        )
        minimap_api(
            args.ref_file,
            "{0}/tmp/pre_contigs.fasta".format(TEMP_DIR),
            "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR),
        )
        map_ref_to_contig(
            contig_dict, logger, "{0}/paf/pre_contigs_to_strain.paf".format(TEMP_DIR)
        )
    # end stat

    # split the branches using link information
    graphf, simp_node_dictf, simp_edge_dictf = iter_graph_disentanglement(
        graph2,
        simp_node_dict2,
        simp_edge_dict2,
        contig_dict,
        pe_info,
        args.ref_file,
        logger,
        0.05 * numpy.median([graph2.vp.dp[node] for node in graph2.vertices()]),
        TEMP_DIR,
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
            "{0}/gfa/split_graph_final.gfa".format(TEMP_DIR),
            logger,
            True,
            "{0}/paf/node_to_ref_final.paf".format(TEMP_DIR),
            "{0}/tmp/temp_gfa_to_fasta_post.fasta".format(TEMP_DIR),
        )
        minimap_api(
            args.ref_file,
            "{0}/tmp/post_contigs.fasta".format(TEMP_DIR),
            "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR),
        )
        map_ref_to_contig(
            contig_dict, logger, "{0}/paf/post_contigs_to_strain.paf".format(TEMP_DIR)
        )
    # end stat
    logger.info(">>>STAGE: contig path extension")

    # refine partial links using best match
    full_link = best_matching(
        graphf, simp_node_dictf, simp_edge_dictf, contig_dict, pe_info, logger
    )

    # update graph coverage on non-trivial branch, maximize
    increment_nt_branch_coverage(graphf, simp_node_dictf, logger)

    graph_to_gfa(
        graphf,
        simp_node_dictf,
        simp_edge_dictf,
        logger,
        "{0}/gfa/split_graph_final.gfa".format(TEMP_DIR),
    )

    # extend the graph
    p_delta = 0.05 * numpy.median([graphf.vp.dp[node] for node in graphf.vertices()])
    strain_dict, usages = path_extension(
        graphf,
        simp_node_dictf,
        simp_edge_dictf,
        contig_dict,
        full_link,
        dcpy_pe_info,
        logger,
        p_delta,
        TEMP_DIR,
    )

    logger.info(">>>STAGE: final process")
    contig_resolve(strain_dict)
    graphl, simp_node_dictl, simp_edge_dictl = flipped_gfa_to_graph(
        "{0}/gfa/es_graph_L2.gfa".format(TEMP_DIR), logger
    )
    trim_contig_dict(graphl, simp_node_dictl, strain_dict, logger)
    contig_dup_removed_s(strain_dict, logger)
    contig_dict_to_path(
        strain_dict, "{0}/tmp/tmp_strain.paths".format(TEMP_DIR), None, False
    )

    # recover repeat nodes back to contig
    strain_repeat_resol(
        graph0, simp_node_dict0, strain_dict, contig_info, copy_contig_dict, logger
    )

    logger.info(">>>STAGE: generate result")
    contig_dict_to_fasta(
        graph0, simp_node_dict0, strain_dict, "{0}/strain.fasta".format(TEMP_DIR)
    )
    contig_dict_to_path(
        strain_dict, "{0}/strain.paths".format(TEMP_DIR), idx_mapping, True
    )
    if args.ref_file:
        minimap_api(
            args.ref_file,
            "{0}/strain.fasta".format(TEMP_DIR),
            "{0}/paf/strain_to_ref.paf".format(TEMP_DIR),
        )
    logger.info("VStrains-SPAdes finished")
    return 0
