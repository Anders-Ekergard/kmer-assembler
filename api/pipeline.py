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

_REFERENCE_GENES_PATH = Path(__file__).resolve().parent.parent / "reference_genes.fasta"
_REFERENCE_GENES = load_reference_genes(str(_REFERENCE_GENES_PATH))


def run_pipeline(
    reads: list[str],
    k: int,
    min_kmer_count: int,
    use_reverse_complement: bool,
    min_edge_support: int = 1,
    min_annotation_score: float = 0.5,
    min_shared_kmers: int = 3,
) -> dict:
    original_reads = reads
    if use_reverse_complement:
        reads = add_reverse_complements(reads)

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

    return {
        "kmer_counts": kmers_count,
        "contigs": contigs,
        "components": [
            {node: sorted(neighbors) for node, neighbors in graph.items()}
            for graph in graphs
        ],
        "butterfly": transcripts,
        "annotations": annotations,
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
            reads = [str(read).strip().upper() for read in payload["reads"] if str(read).strip()]
            k = int(payload.get("k", 5))
            min_kmer_count = int(payload.get("min_kmer_count", 1))
            use_reverse_complement = bool(payload.get("use_reverse_complement", True))
            min_edge_support = int(payload.get("min_edge_support", 1))
            min_annotation_score = float(payload.get("min_annotation_score", 0.5))
            min_shared_kmers = int(payload.get("min_shared_kmers", 3))

            if not reads:
                raise ValueError("At least one read is required")
            if k < 1:
                raise ValueError("k must be a positive integer")

            result = run_pipeline(
                reads,
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
