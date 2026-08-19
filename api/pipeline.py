"""
Vercel serverless function exposing Inchworm + Chrysalis + Butterfly as a
web API for the interactive teaching demo (see index.html).
"""
import json
import sys
from http.server import BaseHTTPRequestHandler
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from inchworm_mvp import (
    add_reverse_complements,
    expand,
    filter_kmers_by_count,
    find_kmers,
    most_common_kmers,
)
from chrysalis_mvp import cluster_kmers
from butterfly_mvp import run_butterfly
from annotate_mvp import annotate_transcripts, load_reference_genes
from express_mvp import compare_samples, quantify_transcripts
from export_mvp import build_report_rows

_REFERENCE_GENES_PATH = Path(__file__).resolve().parent.parent / "reference_genes_hive.fasta"
_REFERENCE_GENES = load_reference_genes(str(_REFERENCE_GENES_PATH))


def run_pipeline(
    samples: dict[str, list[str]],
    k: int,
    min_kmer_count: int,
    use_reverse_complement: bool,
    min_edge_support: int = 1,
    min_annotation_score: float = 0.5,
    min_shared_kmers: int = 3,
) -> dict:
    if len(samples) != 2:
        raise ValueError("Exactly 2 samples are required to compare expression")

    # Co-assembly: pool every sample's reads into one set of transcripts, so
    # there's a single consistent reference to compare expression against
    # across samples (same approach validated in test_pipeline_end_to_end.py).
    original_reads = [read for reads in samples.values() for read in reads]
    reads = add_reverse_complements(original_reads) if use_reverse_complement else original_reads

    kmers_count = filter_kmers_by_count(most_common_kmers(find_kmers(reads, k)), min_kmer_count)
    contigs = expand(kmers_count, k)
    # min_shared_kmers > 1 guards against merging two unrelated contigs (e.g.
    # different organisms in a metatranscriptomic sample) that happen to
    # share a single k-mer by chance rather than a genuine overlap.
    graphs = cluster_kmers(contigs, k, min_shared_kmers)
    # Butterfly does its own forward/reverse-complement check per read, so it
    # gets the original reads, not the already-doubled list above (doubling
    # first would double-count every matching read's support).
    transcripts = run_butterfly(graphs, original_reads, k, min_edge_support)
    annotations = annotate_transcripts(transcripts, _REFERENCE_GENES, k, min_annotation_score)

    sample_names = list(samples.keys())
    quantified = quantify_transcripts(transcripts, samples)
    expression = compare_samples(quantified, sample_names[0], sample_names[1])
    report_rows = build_report_rows(annotations, expression)

    return {
        "kmer_counts": kmers_count,
        "contigs": contigs,
        "components": [
            {node: sorted(neighbors) for node, neighbors in graph.items()}
            for graph in graphs
        ],
        "butterfly": transcripts,
        "annotations": annotations,
        "expression": expression,
        "report_rows": report_rows,
    }


class handler(BaseHTTPRequestHandler):
    def do_GET(self):
        # This handler is Vercel's single entrypoint for the whole project,
        # so it also has to serve the static demo page itself.
        index_path = Path(__file__).resolve().parent.parent / "index.html"
        body = index_path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length)

        try:
            payload = json.loads(body or b"{}")
            samples = {
                str(name): [str(read).strip().upper() for read in reads if str(read).strip()]
                for name, reads in payload["samples"].items()
            }
            k = int(payload.get("k", 5))
            min_kmer_count = int(payload.get("min_kmer_count", 1))
            use_reverse_complement = bool(payload.get("use_reverse_complement", True))
            min_edge_support = int(payload.get("min_edge_support", 1))
            min_annotation_score = float(payload.get("min_annotation_score", 0.5))
            min_shared_kmers = int(payload.get("min_shared_kmers", 3))

            if len(samples) != 2:
                raise ValueError("Exactly 2 samples are required (e.g. control and test)")
            if any(not reads for reads in samples.values()):
                raise ValueError("Each sample needs at least one read")
            if k < 1:
                raise ValueError("k must be a positive integer")

            result = run_pipeline(
                samples,
                k,
                min_kmer_count,
                use_reverse_complement,
                min_edge_support,
                min_annotation_score,
                min_shared_kmers,
            )
            status = 200
        except Exception as exc:
            result = {"error": str(exc)}
            status = 400

        response = json.dumps(result).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(response)))
        self.end_headers()
        self.wfile.write(response)
