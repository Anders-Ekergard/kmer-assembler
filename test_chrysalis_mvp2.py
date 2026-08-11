import unittest

from chrysalis_mvp import cluster_kmers, orient_component
from inchworm_mvp import kmers_in_contig, reverse_complement


class OrientComponentTests(unittest.TestCase):
    def test_orient_component_flips_contigs_to_match_the_first_one(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        component = [forward, reverse_complement(forward)]
        oriented = orient_component(component, 5)
        self.assertEqual(oriented, [forward, forward])

    def test_orient_component_handles_a_single_contig(self):
        self.assertEqual(orient_component(["ATGCA"], 3), ["ATGCA"])

    def test_orient_component_handles_empty_input(self):
        self.assertEqual(orient_component([], 3), [])


class ClusterKmersTests(unittest.TestCase):
    def test_cluster_kmers_merges_a_contig_and_its_reverse_complement(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        contigs = [forward, reverse_complement(forward)]
        graphs = cluster_kmers(contigs, 5)

        self.assertEqual(len(graphs), 1)
        graph = graphs[0]
        self.assertEqual(set(graph.keys()), kmers_in_contig(forward, 5))

    def test_cluster_kmers_keeps_unrelated_contigs_in_separate_graphs(self):
        graphs = cluster_kmers(["ATGCA", "GCATT", "TTTTT"], 3)
        self.assertEqual(len(graphs), 2)


if __name__ == "__main__":
    unittest.main()
