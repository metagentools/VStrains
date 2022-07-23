from ns_Preprocess import graph_to_gfa, gfa_to_graph
from ns_Utilities import map_ref_to_graph
import argparse

__author__ = "Runpeng Luo"
__copyright__ = ""
__credits__ = ["Runpeng Luo", "Yu Lin"]
__license__ = ""
__version__ = "0.0.1"
__maintainer__ = "Runpeng Luo"
__email__ = ""
__status__ = ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="minimap_api.py", description="")
    parser.add_argument(
        "-gfa",
        "--gfa_file",
        dest="gfa_file",
        type=str,
        required=True,
        help="assembly graph under gfa format",
    )
    parser.add_argument(
        "-ref",
        "--reference_fa",
        dest="ref_file",
        type=str,
        help="reference strain, fasta format, DEBUG_MODE only",
    )

    args = parser.parse_args()
    graph, simp_node_dict, simp_edge_dict = gfa_to_graph(args.gfa_file)
    graph_to_gfa(graph, simp_node_dict, simp_edge_dict, "flipped_graph.gfa")

    map_ref_to_graph(args.ref_file, simp_node_dict, args.gfa_file)
