import unittest

from chrysalis_mvp import cluster_kmers, group_contigs_by_kmer_overlap, orient_component
from inchworm_mvp import kmers_in_contig, reverse_complement


class GroupContigsByKmerOverlapTests(unittest.TestCase):
    def test_groups_transitively_connected_contigs(self):
        contigs = ["ATG", "TGC", "GCA"]
        clusters = group_contigs_by_kmer_overlap(contigs, 2)
        self.assertEqual(clusters, [["ATG", "TGC", "GCA"]])

    def test_keeps_unrelated_contigs_separate(self):
        contigs = ["ATG", "TGC", "TTT"]
        clusters = group_contigs_by_kmer_overlap(contigs, 2)
        self.assertEqual(clusters, [["ATG", "TGC"], ["TTT"]])

    def test_groups_a_contig_with_its_reverse_complement(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        contigs = [forward, reverse_complement(forward), "TTTTTTT"]
        clusters = group_contigs_by_kmer_overlap(contigs, 5)
        self.assertEqual(clusters, [[forward, reverse_complement(forward)], ["TTTTTTT"]])


class MinSharedKmersThresholdTests(unittest.TestCase):
    def test_default_threshold_merges_contigs_that_share_a_single_incidental_kmer(self):
        # Two otherwise unrelated contigs (think: two different organisms'
        # transcripts in a metatranscriptomic sample) that happen to share
        # exactly one 3-mer, "TTT" - representative of a chance collision at
        # low k rather than a genuine overlap.
        contig_x = "AAATTT"
        contig_y = "TTTCCC"
        clusters = group_contigs_by_kmer_overlap([contig_x, contig_y], 3)
        self.assertEqual(clusters, [[contig_x, contig_y]])

    def test_raising_the_threshold_keeps_an_incidental_single_kmer_match_separate(self):
        contig_x = "AAATTT"
        contig_y = "TTTCCC"
        clusters = group_contigs_by_kmer_overlap([contig_x, contig_y], 3, min_shared_kmers=2)
        self.assertEqual(clusters, [[contig_x], [contig_y]])

    def test_a_genuine_multi_kmer_overlap_still_merges_under_a_higher_threshold(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        contigs = [forward, reverse_complement(forward), "TTTTTTT"]
        clusters = group_contigs_by_kmer_overlap(contigs, 5, min_shared_kmers=3)
        self.assertEqual(clusters, [[forward, reverse_complement(forward)], ["TTTTTTT"]])


class OrientComponentTests(unittest.TestCase):
    def test_flips_contigs_to_match_the_first_one(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        component = [forward, reverse_complement(forward)]
        oriented = orient_component(component, 5)
        self.assertEqual(oriented, [forward, forward])

    def test_handles_a_single_contig(self):
        self.assertEqual(orient_component(["ATGCA"], 3), ["ATGCA"])

    def test_handles_empty_input(self):
        self.assertEqual(orient_component([], 3), [])


class ClusterKmersTests(unittest.TestCase):
    def test_merges_a_contig_and_its_reverse_complement(self):
        forward = "AAGCCCAATAAACCACTCTGACTG"
        contigs = [forward, reverse_complement(forward)]
        graphs = cluster_kmers(contigs, 5)

        self.assertEqual(len(graphs), 1)
        self.assertEqual(set(graphs[0].keys()), kmers_in_contig(forward, 5))

    def test_keeps_unrelated_contigs_in_separate_graphs(self):
        graphs = cluster_kmers(["ATGCA", "GCATT", "TTTTT"], 3)
        self.assertEqual(len(graphs), 2)

    def test_forwards_min_shared_kmers_to_the_grouping_step(self):
        # "AAATTT" and "TTTCCC" share exactly one 3-mer - cluster_kmers must
        # keep them apart once a higher threshold is requested, same as
        # group_contigs_by_kmer_overlap does directly.
        graphs = cluster_kmers(["AAATTT", "TTTCCC"], 3, min_shared_kmers=2)
        self.assertEqual(len(graphs), 2)


if __name__ == "__main__":
    unittest.main()
