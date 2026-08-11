import unittest

from butterfly_mvp import (
    assemble_component_transcripts,
    edge_support_from_reads,
    find_start_nodes,
    greedy_walk,
    path_to_transcript,
    run_butterfly,
)
from chrysalis_mvp import build_de_bruijn_graph, cluster_kmers
from inchworm_mvp import reverse_complement


class EdgeSupportFromReadsTests(unittest.TestCase):
    def test_counts_matching_edge_from_forward_read(self):
        graph = {"ATGC": {"TGCA"}, "TGCA": set()}
        support = edge_support_from_reads(graph, ["ATGCA"], 4)
        self.assertEqual(support, {("ATGC", "TGCA"): 1})

    def test_counts_matching_edge_from_reverse_complement_read(self):
        graph = {"ATGC": {"TGCA"}, "TGCA": set()}
        read = reverse_complement("ATGCA")
        support = edge_support_from_reads(graph, [read], 4)
        self.assertEqual(support, {("ATGC", "TGCA"): 1})

    def test_ignores_reads_that_match_neither_orientation(self):
        graph = {"ATGC": {"TGCA"}, "TGCA": set()}
        support = edge_support_from_reads(graph, ["TTTTT"], 4)
        self.assertEqual(support, {})

    def test_accumulates_multiple_confirming_reads(self):
        graph = {"ATGC": {"TGCA"}, "TGCA": set()}
        reads = ["ATGCA", "ATGCA", reverse_complement("ATGCA")]
        support = edge_support_from_reads(graph, reads, 4)
        self.assertEqual(support, {("ATGC", "TGCA"): 3})


class FindStartNodesTests(unittest.TestCase):
    def test_finds_the_single_node_with_no_incoming_edges(self):
        graph = {"A": {"B"}, "B": {"C"}, "C": set()}
        self.assertEqual(find_start_nodes(graph), ["A"])

    def test_returns_the_lone_node_for_a_no_edge_graph(self):
        graph = {"ATGCA": set()}
        self.assertEqual(find_start_nodes(graph), ["ATGCA"])

    def test_returns_empty_list_for_a_pure_cycle(self):
        graph = {"A": {"B"}, "B": {"C"}, "C": {"A"}}
        self.assertEqual(find_start_nodes(graph), [])


class GreedyWalkTests(unittest.TestCase):
    def test_walks_a_simple_linear_graph_to_completion(self):
        graph = {"A": {"B"}, "B": {"C"}, "C": set()}
        support = {("A", "B"): 5, ("B", "C"): 5}
        self.assertEqual(greedy_walk(graph, support, "A"), ["A", "B", "C"])

    def test_stops_at_a_true_dead_end(self):
        graph = {"A": set()}
        self.assertEqual(greedy_walk(graph, {}, "A"), ["A"])

    def test_prevents_revisiting_a_node_in_a_cycle(self):
        # A -> B -> C, C's best edge loops back to A (high support), but a
        # weaker escape edge C -> D exists. The walk must not loop forever.
        graph = {"A": {"B"}, "B": {"C"}, "C": {"A", "D"}, "D": set()}
        support = {("A", "B"): 5, ("B", "C"): 5, ("C", "A"): 10, ("C", "D"): 1}
        path = greedy_walk(graph, support, "A")
        self.assertEqual(path, ["A", "B", "C", "D"])
        self.assertEqual(len(path), len(set(path)))

    def test_prefers_the_higher_support_neighbor_when_two_candidates_exist(self):
        graph = {"A": {"B", "C"}, "B": set(), "C": set()}
        support = {("A", "B"): 1, ("A", "C"): 5}
        self.assertEqual(greedy_walk(graph, support, "A"), ["A", "C"])

    def test_stops_extending_once_best_candidate_support_is_below_threshold(self):
        graph = {"A": {"B"}, "B": set()}
        support = {("A", "B"): 1}
        self.assertEqual(greedy_walk(graph, support, "A", min_edge_support=2), ["A"])


class PathToTranscriptTests(unittest.TestCase):
    def test_reconstructs_sequence_from_overlapping_kmers(self):
        self.assertEqual(path_to_transcript(["ATGC", "TGCA", "GCAT"]), "ATGCAT")

    def test_returns_the_single_node_unchanged_for_a_length_one_path(self):
        self.assertEqual(path_to_transcript(["ATGC"]), "ATGC")

    def test_returns_empty_string_for_an_empty_path(self):
        self.assertEqual(path_to_transcript([]), "")


class AssembleComponentTranscriptsTests(unittest.TestCase):
    def test_trivial_single_node_component_yields_that_node_as_transcript(self):
        graph = {"ATGCA": set()}
        self.assertEqual(assemble_component_transcripts(graph, [], 5), ["ATGCA"])

    def test_assembles_a_simple_linear_component(self):
        contig = "AAGCCCAATAAACCACTCTGACTG"
        graph = build_de_bruijn_graph([contig], 5)
        # Dense sliding windows so every edge in the contig is covered.
        reads = [contig[i:i + 12] for i in range(len(contig) - 11)]
        transcripts = assemble_component_transcripts(graph, reads, 5)
        self.assertEqual(transcripts, [contig])

    def test_disambiguates_a_fork_using_read_support(self):
        # ATG forks to TGC or TGA; reads mostly confirm ATG -> TGC.
        graph = {"ATG": {"TGC", "TGA"}, "TGC": set(), "TGA": set()}
        reads = ["ATGC", "ATGC", "ATGC", "ATGA"]
        transcripts = assemble_component_transcripts(graph, reads, 3)
        self.assertEqual(transcripts, ["ATGC"])

    def test_drops_a_walk_whose_only_candidate_edge_is_below_min_edge_support(self):
        graph = {"A": {"B"}, "B": set()}
        transcripts = assemble_component_transcripts(graph, [], 5, min_edge_support=1)
        self.assertEqual(transcripts, [])

    def test_handles_a_pure_cycle_component_by_picking_a_deterministic_start_node(self):
        # min_edge_support=0 isolates the fallback-start behavior from the
        # separate low-support-drop rule (tested on its own above), since
        # there's no read support at all for this abstract A/B/C graph.
        graph = {"A": {"B"}, "B": {"C"}, "C": {"A"}}
        transcripts = assemble_component_transcripts(graph, [], 5, min_edge_support=0)
        self.assertEqual(len(transcripts), 1)

    def test_keeps_distinct_transcripts_from_different_start_nodes_even_if_they_share_a_suffix(self):
        # A and D are separate starts that both merge into B - they must
        # NOT be collapsed into one transcript just because they converge,
        # since that would silently drop a real distinct isoform/prefix.
        # min_edge_support=0 since no reads are supplied here.
        graph = {"A": {"B"}, "D": {"B"}, "B": set()}
        transcripts = assemble_component_transcripts(graph, [], 1, min_edge_support=0)
        self.assertEqual(sorted(transcripts), ["AB", "DB"])


class RunButterflyTests(unittest.TestCase):
    def test_assembles_transcripts_across_multiple_components(self):
        contigs = ["AAGCCCAATAAACCACTCTGACTG", "TTTTTGGGGGCCCCCAAAAATTTTTG"]
        graphs = cluster_kmers(contigs, 5)
        reads = [
            contig[i:i + 12]
            for contig in contigs
            for i in range(len(contig) - 11)
        ]
        transcripts = run_butterfly(graphs, reads, 5)
        self.assertEqual(len(transcripts), len(graphs))

    def test_returns_empty_list_for_no_components(self):
        self.assertEqual(run_butterfly([], [], 5), [])

    def test_deduplicates_across_components(self):
        graphs = [{"A": {"B"}, "B": set()}, {"A": {"B"}, "B": set()}]
        transcripts = run_butterfly(graphs, [], 1, min_edge_support=0)
        self.assertEqual(transcripts, ["AB"])


if __name__ == "__main__":
    unittest.main()
