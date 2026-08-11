import csv
import tempfile
import unittest
from pathlib import Path

from export_mvp import build_report_rows, write_csv


class BuildReportRowsTests(unittest.TestCase):
    def test_merges_annotation_and_expression_for_a_transcript(self):
        annotations = [{"transcript": "ATGC", "match": "gene_a", "score": 0.9}]
        compared = [
            {
                "transcript": "ATGC",
                "length": 4,
                "counts": {"a": 3, "b": 6},
                "rpk": {"a": 750.0, "b": 1500.0},
                "tpm": {"a": 500000.0, "b": 500000.0},
                "log2fc": 0.0,
            }
        ]
        rows = build_report_rows(annotations, compared)
        self.assertEqual(
            rows,
            [
                {
                    "transcript": "ATGC",
                    "length": 4,
                    "match": "gene_a",
                    "annotation_score": 0.9,
                    "count_a": 3,
                    "count_b": 6,
                    "tpm_a": 500000.0,
                    "tpm_b": 500000.0,
                    "log2fc": 0.0,
                }
            ],
        )

    def test_missing_annotation_reports_no_match(self):
        compared = [
            {
                "transcript": "TTTT",
                "length": 4,
                "counts": {"a": 1},
                "rpk": {"a": 250.0},
                "tpm": {"a": 1000000.0},
                "log2fc": 0.0,
            }
        ]
        rows = build_report_rows([], compared)
        self.assertIsNone(rows[0]["match"])
        self.assertIsNone(rows[0]["annotation_score"])

    def test_returns_empty_list_for_no_transcripts(self):
        self.assertEqual(build_report_rows([], []), [])


class WriteCsvTests(unittest.TestCase):
    def test_writes_a_header_and_one_row_per_transcript(self):
        rows = [
            {"transcript": "ATGC", "length": 4, "match": "gene_a", "annotation_score": 1.0, "log2fc": 0.5},
            {"transcript": "TTTT", "length": 4, "match": None, "annotation_score": None, "log2fc": -0.5},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.csv"
            write_csv(rows, str(path))
            with open(path, newline="") as f:
                reader = list(csv.DictReader(f))

        self.assertEqual(len(reader), 2)
        self.assertEqual(reader[0]["transcript"], "ATGC")
        self.assertEqual(reader[1]["match"], "")

    def test_writes_an_empty_file_for_no_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.csv"
            write_csv([], str(path))
            self.assertEqual(path.read_text(), "")


if __name__ == "__main__":
    unittest.main()
