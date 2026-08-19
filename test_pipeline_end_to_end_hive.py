"""
Scaled-up end-to-end validation on REAL gene sequences, separate from
test_pipeline_end_to_end.py's synthetic, engineered-zero-overlap fixture.
Purpose: prove the k-mer pipeline's parameters and runtime hold up at
real-gene scale BEFORE the HIVE demo is wired into
index.html/reference_genes_hive.fasta (see project plan).

STAT1 (_stat1_raw.fasta, NM_007315.4, "transcript variant alpha" - the
commonly cited/canonical STAT1 mRNA) is the real disease-associated gene for
the eventual demo. GAPDH, the chosen housekeeping control gene, hasn't been
sourced yet (a separate step) - FGF12 (_fgf12_raw.fasta, NM_004113.6, already
in the repo) stands in here purely as a second real gene of comparable
length/complexity, to validate parameters/runtime independent of which gene
ends up as the final control. This test's job is to prove the pipeline
mechanics work at real-gene scale, not to make a biological claim about
FGF12 - swapping in GAPDH later should not require different parameters,
since it's the sequence's realism/scale that matters here, not its identity.

Full-length STAT1 (4116nt) / FGF12 (5353nt) were benchmarked as unusably
slow with the current expand_forward/expand_backward implementation: every
contig-extension step rescans the ENTIRE remaining k-mer dict, which is
roughly O(n^2) in sequence length. Measured locally: ~4.7s total at a 900nt
trim, ~45s at 2400nt. A 400nt window keeps the full pipeline well under a
second while still exercising real, non-engineered sequence composition
(repeats, GC content, etc.) that the synthetic fixture deliberately avoids.
This is a known scaling limitation of the current pure-Python approach, not
something this test works around - it's the reason a trimmed window is used
instead of the full transcript.
"""
import time
import unittest

from annotate_mvp import annotate_transcripts
from api.pipeline import run_pipeline
from hive_demo_data import load_gene_sequence, tile_reads

TRIM = 400
K = 21
READ_LEN = 60
MAX_SECONDS = 15  # generous ceiling; observed well under 1s locally at this trim


class RealGeneScaleEndToEndTests(unittest.TestCase):
    """
    Ground truth: STAT1 tiled at higher depth in the "hive" sample than in
    "control" (real induction, matching its role as a HIVE-associated hub
    DEG in the course project). The control-gene stand-in is tiled at an
    IDENTICAL depth in both samples (no real change expected).
    """

    @classmethod
    def setUpClass(cls):
        stat1_full = load_gene_sequence("_stat1_raw.fasta")
        control_gene_full = load_gene_sequence("_fgf12_raw.fasta")  # stand-in, see module docstring
        cls.stat1_seq = stat1_full[:TRIM]
        cls.control_gene_seq = control_gene_full[:TRIM]
        cls.reference_genes = {"STAT1": cls.stat1_seq, "control_gene": cls.control_gene_seq}

        reads_control = (
            tile_reads(cls.stat1_seq, READ_LEN, depth=1)
            + tile_reads(cls.control_gene_seq, READ_LEN, depth=3)
        )
        reads_hive = (
            tile_reads(cls.stat1_seq, READ_LEN, depth=4)
            + tile_reads(cls.control_gene_seq, READ_LEN, depth=3)
        )

        t0 = time.perf_counter()
        cls.result = run_pipeline(
            {"control": reads_control, "hive": reads_hive},
            k=K,
            min_kmer_count=2,
            use_reverse_complement=True,
            min_edge_support=1,
            min_annotation_score=0.5,
            min_shared_kmers=3,
        )
        cls.elapsed = time.perf_counter() - t0

        cls.transcripts = cls.result["butterfly"]
        # run_pipeline annotates against the app's own reference_genes.fasta
        # (still the 3 synthetic placeholders at this point in the rollout),
        # which knows nothing about STAT1/the control gene - re-annotate
        # against the true reference, same approach as
        # test_pipeline_end_to_end.py.
        cls.annotations = annotate_transcripts(cls.transcripts, cls.reference_genes, K, min_score=0.5)
        cls.compared = cls.result["expression"]

    def _entry_for_gene(self, gene_name: str) -> dict:
        match_by_transcript = {a["transcript"]: a["match"] for a in self.annotations}
        for entry in self.compared:
            if match_by_transcript.get(entry["transcript"]) == gene_name:
                return entry
        raise AssertionError(f"no assembled transcript was annotated as {gene_name}")

    def test_runs_within_a_sane_time_budget(self):
        self.assertLess(self.elapsed, MAX_SECONDS, f"pipeline took {self.elapsed:.2f}s, budget was {MAX_SECONDS}s")

    def test_recovers_exactly_two_distinct_transcripts_no_chimera(self):
        self.assertEqual(len(self.transcripts), 2)

    def test_both_transcripts_are_correctly_annotated_and_full_length(self):
        matches = {a["match"] for a in self.annotations}
        self.assertEqual(matches, {"STAT1", "control_gene"})
        for annotation in self.annotations:
            self.assertEqual(annotation["score"], 1.0)
        self.assertEqual({len(t) for t in self.transcripts}, {TRIM})

    def test_run_pipelines_own_live_reference_also_recognizes_stat1(self):
        # Unlike test_pipeline_end_to_end.py's synthetic genes, the live
        # reference_genes_hive.fasta wired into api/pipeline.py now genuinely
        # contains real STAT1 - so, unlike the original test, this checks
        # run_pipeline's own built-in annotations (not just the re-annotation
        # against cls.reference_genes above), to prove the production wiring
        # itself works end-to-end, not just the pipeline mechanics in isolation.
        live_matches = {a["match"] for a in self.result["annotations"] if a["match"]}
        self.assertIn("STAT1", live_matches)

    def test_stat1_shows_real_induction_in_hive_sample(self):
        entry = self._entry_for_gene("STAT1")
        self.assertGreater(entry["counts"]["control"], 0)
        self.assertGreater(entry["counts"]["hive"], entry["counts"]["control"])
        self.assertGreater(entry["log2fc"], 0)

    def test_control_gene_has_equal_raw_counts_but_a_composition_effect_in_tpm(self):
        # Same TPM composition effect documented in test_pipeline_end_to_end.py
        # for gene_beta: identical raw counts in both samples, but STAT1's 4x
        # depth increase in "hive" inflates that sample's total signal, so the
        # control gene's TPM share (and therefore its log2fc) still moves.
        # That's expected behavior of TPM normalization, not a bug - the raw
        # count equality below is the actual "no real change" signal.
        entry = self._entry_for_gene("control_gene")
        self.assertEqual(entry["counts"]["control"], entry["counts"]["hive"])
        self.assertLess(entry["log2fc"], 0)


if __name__ == "__main__":
    unittest.main()
