import unittest

from express_mvp import (
    add_tpm,
    compare_samples,
    count_read_support,
    log2_fold_change,
    quantify_transcripts,
    reads_per_kilobase,
)
from inchworm_mvp import reverse_complement


class CountReadSupportTests(unittest.TestCase):
    def test_counts_reads_that_are_a_substring_of_the_transcript(self):
        transcript = "AAGCCCAATAAACCACTCTGACTG"
        reads = ["AAGCCCAAT", "CCACTCTGACTG", "TTTTTTTT"]
        self.assertEqual(count_read_support(transcript, reads), 2)

    def test_counts_a_read_via_its_reverse_complement(self):
        transcript = "AAGCCCAATAAACCACTCTGACTG"
        read = reverse_complement("AAGCCCAAT")
        self.assertEqual(count_read_support(transcript, [read]), 1)

    def test_does_not_double_count_a_palindromic_style_match(self):
        # A read whose forward orientation already matches must be counted
        # once, not once per orientation checked.
        transcript = "AAGCCCAATAAACCACTCTGACTG"
        reads = ["AAGCCCAAT"] * 3
        self.assertEqual(count_read_support(transcript, reads), 3)

    def test_returns_zero_for_no_matching_reads(self):
        self.assertEqual(count_read_support("AAGCCCAATAAA", ["TTTTTTTT"]), 0)


class ReadsPerKilobaseTests(unittest.TestCase):
    def test_normalizes_count_by_length(self):
        self.assertEqual(reads_per_kilobase(10, 500), 20.0)

    def test_zero_length_transcript_scores_zero(self):
        self.assertEqual(reads_per_kilobase(10, 0), 0.0)


class QuantifyTranscriptsTests(unittest.TestCase):
    def test_reports_counts_and_rpk_per_sample(self):
        transcripts = ["AAGCCCAATAAACCACTCTGACTG"]
        samples = {
            "stimulus_a": ["AAGCCCAAT"],
            "stimulus_b": ["AAGCCCAAT", "CCACTCTGACTG"],
        }
        result = quantify_transcripts(transcripts, samples)
        self.assertEqual(len(result), 1)
        entry = result[0]
        self.assertEqual(entry["length"], 24)
        self.assertEqual(entry["counts"], {"stimulus_a": 1, "stimulus_b": 2})
        self.assertAlmostEqual(entry["rpk"]["stimulus_a"], 1 / (24 / 1000))
        self.assertAlmostEqual(entry["rpk"]["stimulus_b"], 2 / (24 / 1000))


class Log2FoldChangeTests(unittest.TestCase):
    def test_positive_when_sample_b_is_higher(self):
        self.assertGreater(log2_fold_change(value_a=10, value_b=40), 0)

    def test_negative_when_sample_a_is_higher(self):
        self.assertLess(log2_fold_change(value_a=40, value_b=10), 0)

    def test_zero_when_samples_are_equal(self):
        self.assertEqual(log2_fold_change(value_a=15, value_b=15), 0.0)


class AddTpmTests(unittest.TestCase):
    def test_each_samples_tpm_values_sum_to_one_million(self):
        quantified = [
            {"transcript": "x", "length": 10, "counts": {}, "rpk": {"a": 10, "b": 100}},
            {"transcript": "y", "length": 10, "counts": {}, "rpk": {"a": 30, "b": 100}},
        ]
        normalized = add_tpm(quantified)
        for sample in ("a", "b"):
            total = sum(entry["tpm"][sample] for entry in normalized)
            self.assertAlmostEqual(total, 1_000_000)

    def test_corrects_for_a_sample_having_double_the_total_depth(self):
        # Same relative split (1:3) in both samples, but sample b's RPK
        # totals are 10x sample a's - a pure depth difference. TPM should
        # cancel that out and leave both samples reporting the same ratio.
        quantified = [
            {"transcript": "low", "length": 10, "counts": {}, "rpk": {"a": 10, "b": 100}},
            {"transcript": "high", "length": 10, "counts": {}, "rpk": {"a": 30, "b": 300}},
        ]
        normalized = add_tpm(quantified)
        by_transcript = {entry["transcript"]: entry["tpm"] for entry in normalized}
        self.assertAlmostEqual(by_transcript["low"]["a"], by_transcript["low"]["b"])
        self.assertAlmostEqual(by_transcript["high"]["a"], by_transcript["high"]["b"])

    def test_returns_empty_list_for_no_transcripts(self):
        self.assertEqual(add_tpm([]), [])


class CompareSamplesTests(unittest.TestCase):
    def test_sorts_by_absolute_fold_change_descending(self):
        # A dominant, stable "housekeeping" transcript anchors each sample's
        # total depth - without it, "induced" swinging from 2 to 80 would
        # itself dominate the totals and distort the comparison, the same
        # composition effect a real dataset avoids by having many transcripts.
        quantified = [
            {"transcript": "housekeeping", "length": 10, "counts": {}, "rpk": {"a": 1000, "b": 1000}},
            {"transcript": "flat", "length": 10, "counts": {}, "rpk": {"a": 10, "b": 11}},
            {"transcript": "induced", "length": 10, "counts": {}, "rpk": {"a": 2, "b": 80}},
        ]
        compared = compare_samples(quantified, "a", "b")
        self.assertEqual(compared[0]["transcript"], "induced")

    def test_adds_log2fc_to_every_entry(self):
        quantified = [{"transcript": "t", "length": 10, "counts": {}, "rpk": {"a": 5, "b": 5}}]
        compared = compare_samples(quantified, "a", "b")
        self.assertEqual(compared[0]["log2fc"], 0.0)

    def test_adds_tpm_alongside_log2fc(self):
        quantified = [{"transcript": "t", "length": 10, "counts": {}, "rpk": {"a": 5, "b": 15}}]
        compared = compare_samples(quantified, "a", "b")
        self.assertIn("tpm", compared[0])


if __name__ == "__main__":
    unittest.main()
