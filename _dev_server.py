import sys
from http.server import ThreadingHTTPServer
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "api"))
from pipeline import handler

# ThreadingHTTPServer (not plain HTTPServer): a real browser opens more than
# one connection per page (keep-alive / speculative preconnects) - see
# hybrid_lab/_dev_server.py for the same fix and why it's needed.
#
# This exists because `vercel dev` doesn't work on Windows checkouts of this
# repo: the Vercel CLI embeds this folder's absolute path unescaped into a
# generated .py file, and a `\Users` segment gets parsed as a truncated
# `\UXXXXXXXX` unicode escape (a Vercel CLI bug, unrelated to this repo's
# code). Running api/pipeline.py's real `handler` class directly sidesteps
# that entirely - same code Vercel would run, just served locally.
server = ThreadingHTTPServer(("127.0.0.1", 3111), handler)
print("Serving kmer-assembler on http://127.0.0.1:3111")
server.serve_forever()
