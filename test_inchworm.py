import unittest

from Inchworm import assemble_transcripts
from inchworm_mvp import share_kmer


class InchwormAssemblyTests(unittest.TestCase):
    def test_greedy_assembly_builds_a_single_transcript(self):
        reads = ["ATG", "TGC", "GCA", "CAT"]
        transcripts = assemble_transcripts(reads, min_overlap=2)
        self.assertEqual(transcripts, ["ATGCAT"])


class InchwormMvpTests(unittest.TestCase):
    def test_share_kmer_returns_true_when_contigs_share_a_kmer(self):
        self.assertTrue(share_kmer("ATGCGT", "GCGTAA", 3))

    def test_share_kmer_returns_false_when_contigs_do_not_share_a_kmer(self):
        self.assertFalse(share_kmer("ATGCGT", "TTTAAA", 3))


def test_curt_kmer():
    kmer = "AAATC"
    k = 4
    current = "TATTC"
 
    return kmer[-(k-1):], current[:k-1], kmer[:k-1], current[-(k-1):]
if __name__ == "__main__":
    #unittest.main()
    results = list(test_curt_kmer())
    for r in results:
        print(r)