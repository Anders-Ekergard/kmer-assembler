"""
This module combines annotation and expression results into a single CSV
report - the pipeline's "other end": one row per assembled transcript, with
its best reference match and its per-sample expression alongside it.
"""
import csv


def build_report_rows(annotations: list[dict], compared: list[dict]) -> list[dict]:
    """
    Combine per-transcript annotation and expression results into one flat
    row per transcript, keyed by transcript sequence (unique per assembled
    transcript within a single pipeline run).
    args:
        annotations: list[dict] - output of annotate_transcripts
        compared: list[dict] - output of compare_samples
    returns:
        list[dict] - one row per transcript: transcript, length, match,
        annotation_score, count_<sample>/tpm_<sample> per sample, log2fc
    """
    annotation_by_transcript = {entry["transcript"]: entry for entry in annotations}

    rows = []
    for entry in compared:
        annotation = annotation_by_transcript.get(entry["transcript"], {})
        row = {
            "transcript": entry["transcript"],
            "length": entry["length"],
            "match": annotation.get("match"),
            "annotation_score": annotation.get("score"),
        }
        for sample, count in entry["counts"].items():
            row[f"count_{sample}"] = count
        for sample, tpm in entry["tpm"].items():
            row[f"tpm_{sample}"] = round(tpm, 2)
        row["log2fc"] = entry["log2fc"]
        rows.append(row)
    return rows


def write_csv(rows: list[dict], filepath: str) -> None:
    """
    Write report rows to a CSV file, one row per transcript.
    args:
        rows: list[dict] - output of build_report_rows
        filepath: str - path to write the CSV to
    """
    with open(filepath, "w", newline="") as f:
        if not rows:
            return
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    from annotate_mvp import annotate_transcripts, load_reference_genes
    from inchworm_mvp import load_contigs, load_reads
    from express_mvp import compare_samples, quantify_transcripts

    k = 5
    transcripts = load_contigs("contigs.fasta")
    reference_genes = load_reference_genes("reference_genes.fasta")
    annotations = annotate_transcripts(transcripts, reference_genes, k, min_score=0.5)

    samples = {
        "stimulus_a": load_reads("reads_stimulus_a.fastq"),
        "stimulus_b": load_reads("reads_stimulus_b.fastq"),
    }
    quantified = quantify_transcripts(transcripts, samples)
    compared = compare_samples(quantified, "stimulus_a", "stimulus_b")

    rows = build_report_rows(annotations, compared)
    write_csv(rows, "report.csv")
    print(f"Wrote {len(rows)} row(s) to report.csv")
