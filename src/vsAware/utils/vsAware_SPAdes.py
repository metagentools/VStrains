#!/usr/bin/env python3

from utils.vsAware_Utilities import *
from utils.vsAware_Preprocess import (
    graph_simplification,
    reindexing,
    threshold_estimation,
)
from utils.vsAware_IO import (
    graph_to_gfa,
    flipped_gfa_to_graph,
    gfa_to_graph,
    contig_dict_to_path,
    contig_dict_to_fasta,
    spades_paths_parser,
    process_pe_info,
    store_reinit_graph
)
from utils.vsAware_Decomposition import *
from utils.vsAware_Extension import path_extension, best_matching
from utils.vsAware_IO import strain_dict_to_fasta


__author__ = "Runpeng Luo"
__copyright__ = "Copyright 2022-2025, vsAware Project"
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Runpeng Luo"
__email__ = "John.Luo@anu.edu.au"
__status__ = "Production"


def run(args, logger):
    TEMP_DIR = args.output_dir

    logger.info("vsAware-SPAdes started")

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
    # map_ref_to_graph(
    #     args.ref_file,
    #     simp_node_dict0,
    #     "{0}/gfa/graph_L0r.gfa".format(TEMP_DIR),
    #     logger,
    #     True,
    #     "{0}/paf/node_to_ref_0.paf".format(TEMP_DIR),
    #     "{0}/tmp/temp_gfa_to_fasta_0.fasta".format(TEMP_DIR),
    # )

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
    # FIXME we may need to use paired end information to better filter out erroroness nodes and edges
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

    # obtain paired end information
    pe_info = process_pe_info(simp_node_dict1.keys(), args.pe_info, args.st_info)
    dcpy_pe_info = {}
    for (uid, wid), u in pe_info.items():
        dcpy_pe_info[(uid, wid)] = u

    # TODO may employ paired end information to remove erroroness edges
    # if not clean, the false edge may lead to a NT 
    # branch, which will split wrong. e.g., 5hiv node (58)
    edge_cleaning(graph1, simp_edge_dict1, contig_dict, pe_info, logger)

    graph2, simp_node_dict2, simp_edge_dict2 = store_reinit_graph(graph1, simp_node_dict1, simp_edge_dict1, logger, "{0}/gfa/es_graph_L2.gfa".format(TEMP_DIR))
    # coverage_rebalance_s(graph2, simp_node_dict2, simp_edge_dict2, logger)

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
    graphf, simp_node_dictf, simp_edge_dictf = iter_graph_decomposition(
        graph2, simp_node_dict2, simp_edge_dict2,
        contig_dict, pe_info, None, logger, 0.05 * numpy.median([graph0.vp.dp[node] for node in graph2.vertices()]), TEMP_DIR)

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
    full_link = best_matching(graphf, simp_node_dictf, simp_edge_dictf, contig_dict, pe_info, logger)
    # extend the graph

    # update graph coverage on non-trivial branch, maximize
    increment_nt_branch_coverage(graphf, simp_node_dictf, logger)

    graph_to_gfa(
        graphf,
        simp_node_dictf,
        simp_edge_dictf,
        logger,
        "{0}/gfa/split_graph_final.gfa".format(TEMP_DIR),
    )

    strain_dict, usages = path_extension(graphf, simp_node_dictf, simp_edge_dictf, contig_dict, full_link, dcpy_pe_info, logger, 0.05 * numpy.median([graphf.vp.dp[node] for node in graphf.vertices()]), args.ref_file, TEMP_DIR)

    contig_dup_removed_s(strain_dict, logger)
    # contig_dict_to_path(strain_dict, "{0}/tmp/tmp_strain.paths".format(TEMP_DIR), None, False)
    strain_dict_to_fasta(
        strain_dict, "{0}/tmp/tmp_strain.fasta".format(TEMP_DIR)
    )

    minimap_api(
        args.ref_file,
        "{0}/tmp/tmp_strain.fasta".format(TEMP_DIR),
        "{0}/paf/strain_to_ref.paf".format(TEMP_DIR),
    )




    # contig_scaffording(graphf, simp_node_dictf, strain_dict, usages, logger, TEMP_DIR)
    
    # graphp, simp_node_dictp, simp_edge_dictp = flipped_gfa_to_graph(
    #     "{0}/gfa/cts_graph_L3.gfa".format(TEMP_DIR), logger
    # )
    # usages_p = contig_resolve(simp_node_dictp.keys(), strain_dict, usages)

    # contig_dup_removed_s(strain_dict, usages_p, logger)
    # trim_contig_dict(graphp, simp_node_dictp, strain_dict, logger)

    # # contig_scaffording(graphp, simp_node_dictp, strain_dict, usages_p, logger, TEMP_DIR)

    # contig_dict_to_path(strain_dict, "{0}/tmp/tmp_strain.paths".format(TEMP_DIR))
    # contig_dict_to_fasta(
    #     graphp, simp_node_dictp, strain_dict, "{0}/tmp/tmp_strain.fasta".format(TEMP_DIR)
    # )

    # if args.fwd and args.rve:
    #     minimap_api("{0}/tmp/tmp_strain.fasta".format(TEMP_DIR), args.fwd,  "{0}/tmp/fwd_aln.paf".format(TEMP_DIR))
    #     minimap_api("{0}/tmp/tmp_strain.fasta".format(TEMP_DIR), args.rve,  "{0}/tmp/rve_aln.paf".format(TEMP_DIR))
    #     link_dict = aln2contig(strain_dict, logger, "{0}/tmp/fwd_aln.paf".format(TEMP_DIR), "{0}/tmp/rve_aln.paf".format(TEMP_DIR))


    # concat_overlap_contig(graphf, simp_node_dictf, simp_edge_dictf, strain_dict, logger)
    # logger.info(">>>STAGE: final process")
    # # recover repeat nodes back to contig
    # strain_repeat_resol(
    #     graph0, simp_node_dict0, strain_dict, contig_info, copy_contig_dict, logger
    # )

    # logger.info(">>>STAGE: generate result")
    # contig_dict_to_fasta(
    #     graph0, simp_node_dict0, strain_dict, "{0}/strain.fasta".format(TEMP_DIR)
    # )
    # contig_dict_to_path(
    #     strain_dict, "{0}/strain.paths".format(TEMP_DIR), idx_mapping, True
    # )
    # if args.ref_file:
    #     minimap_api(
    #         args.ref_file,
    #         "{0}/strain.fasta".format(TEMP_DIR),
    #         "{0}/paf/strain_to_ref.paf".format(TEMP_DIR),
    #     )
    # logger.info("vsAware-SPAdes finished")
    return 0

def increment_nt_branch_coverage(graph: Graph, simp_node_dict: dict, logger: Logger):
    nt_branches = get_non_trivial_branches(graph, simp_node_dict)
    for no, node in nt_branches.items():
        prev_dp = graph.vp.dp[node]
        if (
            sum([x.out_degree() for x in node.in_neighbors()]) == node.in_degree()
            and sum([y.in_degree() for y in node.out_neighbors()]) == node.out_degree()
        ):
            sum_in_dp = sum(graph.vp.dp[n] for n in node.in_neighbors())
            sum_out_dp = sum(graph.vp.dp[n] for n in node.out_neighbors())
            graph.vp.dp[node] = max([prev_dp, sum_in_dp, sum_out_dp])
            logger.debug("Simple NT Branch:{0}, cov: {1} -> {2}".format(no, prev_dp, graph.vp.dp[node]))

        else:
            sum_in_flow = sum(graph.ep.flow[e] for e in node.in_edges())
            sum_out_flow = sum(graph.ep.flow[e] for e in node.out_edges())
            graph.vp.dp[node] = max([prev_dp, sum_in_flow, sum_out_flow])
            logger.debug("Non-Simple NT Branch:{0}, cov: {1} -> {2}".format(no, prev_dp, graph.vp.dp[node]))

def contig_resolve(ids, contig_dict: dict, usages: dict):
    rid = ""
    for cno in contig_dict.keys():
        [contig, clen, ccov] = contig_dict[cno]
        rcontig = []
        for id in contig:
            for iid in str(id).split("&"):
                if iid.find("*") != -1:
                    rid = iid[: iid.find("*")]
                else:
                    rid = iid
                rcontig.append(rid)
        contig_dict[cno] = [rcontig, clen, ccov]
    rtn_usages = dict.fromkeys(ids, 0)

    for id, u in usages.items():
        for iid in str(id).split("&"):
            if iid.find("*") != -1:
                rid = iid[: iid.find("*")]
            else:
                rid = iid
            rtn_usages[rid] += u
    return rtn_usages

def contig_scaffording(graph: Graph, simp_node_dict: dict, contig_dict: dict, usages: dict, logger: Logger, temp_dir):
    non_trivial_branches = get_non_trivial_branches(graph, simp_node_dict)
    node_to_contig_dict, edge_to_contig_dict = contig_map_node(contig_dict)
    for no, node in non_trivial_branches.items():
        us = [graph.vp.id[src] for src in node.in_neighbors()]
        ws = [graph.vp.id[tgt] for tgt in node.out_neighbors()]
        logger.debug("---------------------------------------------")
        logger.debug("current non trivial branch: {0}, in-degree: {1}, out-degree: {2}".format(no, len(us), len(ws)))
        if len(us) == len(ws):
            logger.debug("NN Branch")
        elif len(us) < len(ws):
            logger.debug("NM Branch")
            if all([usages[wid] == 1 for wid in ws]):
                logger.debug("M side fully used")
            elif any([usages[wid] > 1 for wid in ws]):
                logger.debug("M side over used")
            else:
                logger.debug("M side partial used")
        else:
            logger.debug("MN Branch")
            if all([usages[uid] == 1 for uid in us]):
                logger.debug("M side fully used")
            elif any([usages[uid] > 1 for uid in us]):
                logger.debug("M side over used")
            else:
                logger.debug("M side partial used")

        logger.debug("-----> U-Side")
        for uid in us:
            logger.debug("{0} - [{1}], used by {2}".format(uid, usages[uid], node_to_contig_dict.get(uid, None)))
        logger.debug("-----> W-Side")
        for wid in ws:
            logger.debug("{0} - [{1}], used by {2}".format(wid, usages[wid], node_to_contig_dict.get(wid, None))) 
    return