"""
Vercel serverless function exposing pathway_mvp's KEGG/STRING lookups as a
web API - deliberately a separate endpoint from api/pipeline.py's
synchronous assembly pipeline. Gene->pathway/PPI mappings are static per
gene (they don't depend on the assembled reads at all), and third-party API
latency shouldn't risk timing out the core DEG-finding request, which
already has no maxDuration headroom to spare (see project plan).
"""
import json
import sys
from http.server import BaseHTTPRequestHandler
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from pathway_mvp import annotate_with_pathways


class handler(BaseHTTPRequestHandler):
    def do_POST(self):
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length)

        try:
            payload = json.loads(body or b"{}")
            genes = sorted({str(name).strip() for name in payload.get("genes", []) if str(name).strip()})
            if not genes:
                raise ValueError("At least one gene name is required")

            result = {"results": annotate_with_pathways(genes)}
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
