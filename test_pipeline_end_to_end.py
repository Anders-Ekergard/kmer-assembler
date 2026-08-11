"""
End-to-end integration test: a small synthetic "metatranscriptomic" sample
with a known ground truth, run through the full pipeline (Inchworm ->
Chrysalis -> Butterfly -> Annotate -> Express -> CSV).

The three synthetic genes stand in for three unrelated organisms whose
transcripts all end up pooled together in one sequencing run, same as a
real metatranscriptomic sample. They were picked by random search (offline,
not part of this file) to guarantee zero shared k-mers AND zero shared
(k-1)-mers with each other or their own reverse complements at k=6 - i.e.
no coincidental overlap at all, real or spurious, so any merging Chrysalis
does is a bug, not noise in the fixture. This is what lets the "no
chimera" assertion below be a trustworthy regression check rather than a
coin flip.
"""
import csv
import tempfile
import unittest
from pathlib import Path

from annotate_mvp import annotate_transcripts
from api.pipeline import run_pipeline
from export_mvp import build_report_rows, write_csv
from express_mvp import compare_samples, quantify_transcripts

GENE_ALPHA = "CACCTGGTGATCCTATGCTTGTGA"
GENE_BETA = "GTACCCAGAAAATAGCGACGGACC"
GENE_GAMMA = "GCGGTGTTAAGTGTCGAGCTACAT"
K = 6
READ_LEN = 14


def _tile_reads(sequence: str, read_len: int, depth: int) -> list[str]:
    """Slide a window over sequence and repeat the tiling `depth` times, to simulate a given sequencing depth."""
    tiles = [sequence[i:i + read_len] for i in range(len(sequence) - read_len + 1)]
    return tiles * depth


class MockCommunityEndToEndTests(unittest.TestCase):
    """
    Ground truth: gene_alpha present in both samples, higher in stimulus_b
    (real induction). gene_beta present at an IDENTICAL raw read count in
    both samples (no real change) - but gene_gamma showing up only in
    stimulus_b dilutes stimulus_b's total signal, so gene_beta's TPM share
    still drops there. That's not a bug: it's the same composition effect
    covered when TPM normalization was added to express_mvp.py, kept here
    on purpose rather than tuned away, since it's a real thing a user will
    hit. gene_gamma is absent from stimulus_a entirely (a "species" that
    only appears under one stimulus).
    """

    @classmethod
    def setUpClass(cls):
        cls.reads_a = (
            _tile_reads(GENE_ALPHA, READ_LEN, depth=1)
            + _tile_reads(GENE_BETA, READ_LEN, depth=3)
        )
        cls.reads_b = (
            _tile_reads(GENE_ALPHA, READ_LEN, depth=4)
            + _tile_reads(GENE_BETA, READ_LEN, depth=3)
            + _tile_reads(GENE_GAMMA, READ_LEN, depth=3)
        )
        pooled_reads = cls.reads_a + cls.reads_b

        cls.reference_genes = {"gene_alpha": GENE_ALPHA, "gene_beta": GENE_BETA, "gene_gamma": GENE_GAMMA}

        # Co-assembly: pool both samples' reads for Inchworm/Chrysalis/Butterfly,
        # so there's one consistent transcript set to compare expression against -
        # see the per-sample-vs-co-assembly design question flagged in project memory.
        result = run_pipeline(
            pooled_reads,
            k=K,
            min_kmer_count=2,
            use_reverse_complement=True,
            min_edge_support=1,
            min_annotation_score=0.5,
            min_shared_kmers=3,
        )
        cls.transcripts = result["butterfly"]
        # run_pipeline annotates against the app's own reference_genes.fasta, which
        # knows nothing about our synthetic genes - re-annotate against the true
        # reference so the test's ground truth doesn't depend on that file's contents.
        cls.annotations = annotate_transcripts(cls.transcripts, cls.reference_genes, K, min_score=0.5)

        quantified = quantify_transcripts(cls.transcripts, {"stimulus_a": cls.reads_a, "stimulus_b": cls.reads_b})
        cls.compared = compare_samples(quantified, "stimulus_a", "stimulus_b")

    def _entry_for_gene(self, gene_name: str) -> dict:
        match_by_transcript = {a["transcript"]: a["match"] for a in self.annotations}
        for entry in self.compared:
            if match_by_transcript.get(entry["transcript"]) == gene_name:
                return entry
        raise AssertionError(f"no assembled transcript was annotated as {gene_name}")

    def test_recovers_exactly_three_distinct_transcripts_no_chimera(self):
        # The core regression check this test exists for: three unrelated
        # synthetic genes pooled together must come back as three separate
        # transcripts, not fewer (merged into a chimera) or more (a gene
        # split into fragments).
        self.assertEqual(len(self.transcripts), 3)

    def test_every_transcript_is_correctly_annotated_to_its_source_gene(self):
        matches = {a["match"] for a in self.annotations}
        self.assertEqual(matches, {"gene_alpha", "gene_beta", "gene_gamma"})
        for annotation in self.annotations:
            self.assertEqual(annotation["score"], 1.0)

    def test_gene_alpha_shows_real_induction_in_stimulus_b(self):
        entry = self._entry_for_gene("gene_alpha")
        self.assertGreater(entry["counts"]["stimulus_a"], 0)
        self.assertGreater(entry["counts"]["stimulus_b"], entry["counts"]["stimulus_a"])
        self.assertGreater(entry["log2fc"], 0)

    def test_gene_beta_has_equal_raw_counts_but_a_composition_effect_in_tpm(self):
        entry = self._entry_for_gene("gene_beta")
        self.assertEqual(entry["counts"]["stimulus_a"], entry["counts"]["stimulus_b"])
        # Despite identical raw counts, gene_gamma's arrival in stimulus_b
        # inflates that sample's total signal, so gene_beta's TPM share (and
        # therefore its fold change) still moves - a real TPM composition
        # effect, not a bug in the comparison.
        self.assertLess(entry["log2fc"], 0)

    def test_gene_gamma_is_absent_from_stimulus_a_and_present_in_b(self):
        entry = self._entry_for_gene("gene_gamma")
        self.assertEqual(entry["counts"]["stimulus_a"], 0)
        self.assertGreater(entry["counts"]["stimulus_b"], 0)

    def test_csv_export_has_one_row_per_transcript_with_the_right_matches(self):
        rows = build_report_rows(self.annotations, self.compared)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.csv"
            write_csv(rows, str(path))
            with open(path, newline="") as f:
                reader = list(csv.DictReader(f))

        self.assertEqual(len(reader), 3)
        self.assertEqual(
            set(reader[0].keys()),
            {
                "transcript", "length", "match", "annotation_score",
                "count_stimulus_a", "count_stimulus_b",
                "tpm_stimulus_a", "tpm_stimulus_b", "log2fc",
            },
        )
        self.assertEqual({row["match"] for row in reader}, {"gene_alpha", "gene_beta", "gene_gamma"})


if __name__ == "__main__":
    unittest.main()
