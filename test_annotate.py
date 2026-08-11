import tempfile
import unittest
from pathlib import Path

from annotate_mvp import (
    annotate_transcripts,
    best_match,
    kmer_containment,
    load_reference_genes,
)
from inchworm_mvp import reverse_complement


class LoadReferenceGenesTests(unittest.TestCase):
    def test_parses_multiple_records_including_multi_line_sequences(self):
        fasta = ">gene_a\nATGC\nGTAA\n>gene_b\nTTTT\n"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ref.fasta"
            path.write_text(fasta)
            genes = load_reference_genes(str(path))
        self.assertEqual(genes, {"gene_a": "ATGCGTAA", "gene_b": "TTTT"})


class KmerContainmentTests(unittest.TestCase):
    def test_full_containment_scores_one(self):
        self.assertEqual(kmer_containment("ATGCA", "AAATGCATTT", 3), 1.0)

    def test_disjoint_sequences_score_zero(self):
        self.assertEqual(kmer_containment("ATGCA", "TTTTTTT", 3), 0.0)

    def test_partial_overlap_scores_between_zero_and_one(self):
        score = kmer_containment("ATGCATT", "ATGCAGG", 3)
        self.assertTrue(0.0 < score < 1.0)

    def test_query_shorter_than_k_scores_zero(self):
        self.assertEqual(kmer_containment("AT", "ATGCATGC", 3), 0.0)


class BestMatchTests(unittest.TestCase):
    def test_picks_the_reference_with_highest_containment(self):
        reference_genes = {"gene_a": "AAATGCATTT", "gene_b": "GGGGGGGGGG"}
        name, score = best_match("ATGCA", reference_genes, 3)
        self.assertEqual(name, "gene_a")
        self.assertEqual(score, 1.0)

    def test_matches_a_reverse_complement_transcript(self):
        # The transcript is stored as the reverse complement of a sequence
        # that fully matches gene_a - best_match must flip it back to see that.
        reference_genes = {"gene_a": "AAATGCATTT"}
        transcript = reverse_complement("ATGCA")
        name, score = best_match(transcript, reference_genes, 3)
        self.assertEqual(name, "gene_a")
        self.assertEqual(score, 1.0)

    def test_returns_none_when_no_reference_shares_a_kmer(self):
        reference_genes = {"gene_a": "AAATGCATTT"}
        name, score = best_match("CCCCCCC", reference_genes, 3)
        self.assertIsNone(name)
        self.assertEqual(score, 0.0)


class AnnotateTranscriptsTests(unittest.TestCase):
    def test_reports_a_match_above_the_threshold(self):
        reference_genes = {"gene_a": "AAATGCATTT"}
        result = annotate_transcripts(["ATGCA"], reference_genes, 3, min_score=0.5)
        self.assertEqual(result, [{"transcript": "ATGCA", "match": "gene_a", "score": 1.0}])

    def test_suppresses_a_match_below_the_threshold(self):
        reference_genes = {"gene_a": "AAATGCATTT"}
        result = annotate_transcripts(["CCCCCCC"], reference_genes, 3, min_score=0.5)
        self.assertEqual(result, [{"transcript": "CCCCCCC", "match": None, "score": 0.0}])


if __name__ == "__main__":
    unittest.main()
