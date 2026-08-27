"""
Pathway (KEGG) and protein-protein interaction (STRING) lookups for
annotated genes, mirroring the manual GEO2R/Enrichr/STRING/KEGG workflow in
the course project (_course_pdf.txt) - e.g. confirming STAT1 sits in the
JAK-STAT signaling pathway and interacts with JAK1/IRF1/PIAS1. Works for any
human gene symbol annotate_mvp.annotate_transcripts might report, not a
curated handful - gene name -> KEGG gene ID resolution is a live NCBI Entrez
lookup (see resolve_kegg_gene_id), not a static table, so this stays as
generic as the rest of the app (see project plan's session-3 addendum on why
a hardcoded gene list would have quietly broken that).

Deliberately separate from api/pipeline.py's synchronous assembly pipeline:
gene->pathway/PPI mappings are static per gene (they don't depend on the
assembled reads at all), and third-party API latency shouldn't risk timing
out the core DEG-finding request, which already has no maxDuration headroom
to spare (see project plan). Stdlib-only (urllib.request), matching the
app's zero-dependency pitch - the cost is KEGG's flat tab-delimited text
format needing manual parsing instead of a JSON client.
"""
import functools
import json
import urllib.error
import urllib.parse
import urllib.request

KEGG_BASE = "https://rest.kegg.jp"
STRING_BASE = "https://string-db.org/api"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@functools.lru_cache(maxsize=64)
def resolve_kegg_gene_id(symbol: str, timeout: float = 10.0) -> str | None:
    """
    Resolve a human gene symbol to a KEGG gene ID via NCBI Entrez `esearch`
    (the same eutils service already used elsewhere in this repo, see
    in_silico_avel/candidate_gene_dossier.py), so pathway lookups work for
    any gene annotate_mvp recognizes rather than a small curated set.
    Cached (per warm process, e.g. a reused serverless container) since the
    same symbol is often looked up repeatedly.
    args:
        symbol: str - a gene symbol, e.g. "STAT1"
        timeout: float - request timeout in seconds
    returns:
        str | None - a KEGG gene ID, e.g. "hsa:6772", or None if Entrez
        found no matching human gene
    """
    query = urllib.parse.urlencode({
        "db": "gene",
        "term": f"{symbol}[sym] AND Homo sapiens[orgn]",
        "retmode": "json",
    })
    url = f"{EUTILS_BASE}/esearch.fcgi?{query}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist", [])
    return f"hsa:{ids[0]}" if ids else None


def parse_kegg_tsv(text: str) -> list[tuple[str, str]]:
    """
    Parse KEGG's flat, tab-delimited response format (one "key\\tvalue" pair
    per line) shared by both the `link` and `list` operations.
    args:
        text: str - raw response body from a KEGG REST API call
    returns:
        list[tuple[str, str]] - (key, value) pairs in response order
    """
    pairs = []
    for line in text.strip().splitlines():
        if not line:
            continue
        key, value = line.split("\t", 1)
        pairs.append((key, value))
    return pairs


def fetch_kegg_pathway_ids(gene_id: str, timeout: float = 10.0) -> list[str]:
    """
    Look up which KEGG pathways a gene belongs to.
    args:
        gene_id: str - a KEGG gene ID, e.g. "hsa:6772" (see resolve_kegg_gene_id)
    returns:
        list[str] - KEGG pathway IDs, e.g. ["path:hsa04630", ...]
    """
    url = f"{KEGG_BASE}/link/pathway/{gene_id}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        text = response.read().decode("utf-8")
    return [pathway_id for _gene, pathway_id in parse_kegg_tsv(text)]


@functools.lru_cache(maxsize=8)
def _fetch_kegg_pathway_list(organism: str, timeout: float) -> tuple[tuple[str, str], ...]:
    """
    Fetch and cache an organism's full pathway ID -> name list. Cached
    (per warm process, e.g. a reused serverless container) since this is
    the same few-hundred-line response regardless of which gene triggered
    the lookup - refetching it per gene, per request would be wasteful once
    this is a live endpoint rather than a one-off local script.
    """
    url = f"{KEGG_BASE}/list/pathway/{organism}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        text = response.read().decode("utf-8")
    return tuple(parse_kegg_tsv(text))


def fetch_kegg_pathway_names(pathway_ids: list[str], organism: str = "hsa", timeout: float = 10.0) -> dict[str, str]:
    """
    Look up human-readable names for KEGG pathway IDs. KEGG's `list`
    operation only resolves whole databases/organisms, not individual
    pathway entries directly, so this fetches the organism's full pathway
    list once (cached, see _fetch_kegg_pathway_list) and filters down to
    the requested IDs (the `list` response keys pathways without the
    "path:" prefix that `link` uses, so that prefix is stripped for
    matching).
    args:
        pathway_ids: list[str] - pathway IDs as returned by fetch_kegg_pathway_ids
        organism: str - KEGG organism code, "hsa" (Homo sapiens) by default
    returns:
        dict[str, str] - pathway ID (with "path:" prefix, as given) to name
    """
    if not pathway_ids:
        return {}
    names_by_bare_id = dict(_fetch_kegg_pathway_list(organism, timeout))
    return {
        pathway_id: names_by_bare_id[pathway_id.removeprefix("path:")]
        for pathway_id in pathway_ids
        if pathway_id.removeprefix("path:") in names_by_bare_id
    }


def fetch_string_ppi_partners(gene_name: str, species: int = 9606, limit: int = 10, timeout: float = 10.0) -> list[dict]:
    """
    Look up STRING's protein-protein interaction partners for a gene.
    args:
        gene_name: str - a gene symbol, e.g. "STAT1"
        species: int - NCBI taxonomy ID, 9606 (Homo sapiens) by default
        limit: int - maximum number of partners to request
    returns:
        list[dict] - one entry per interaction partner: {"partner": str, "score": float}
    """
    url = f"{STRING_BASE}/json/interaction_partners?identifiers={gene_name}&species={species}&limit={limit}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))
    return [{"partner": entry["preferredName_B"], "score": entry["score"]} for entry in data]


def fetch_string_network_image(gene_names: list[str], species: int = 9606, timeout: float = 10.0) -> bytes:
    """
    Fetch a rendered PPI network diagram (PNG) for a set of genes from
    STRING - the same colored node-and-edge figure STRING's own website
    shows (e.g. the STAT1/OAS1/IFI44L/ISG15 hub cluster in the course
    project). Multiple identifiers are separated with "%0d" (a literal
    carriage return), STRING's documented convention for this endpoint.
    args:
        gene_names: list[str] - gene symbols to include in the network
        species: int - NCBI taxonomy ID, 9606 (Homo sapiens) by default
        timeout: float - request timeout in seconds
    returns:
        bytes - PNG image data
    """
    identifiers = "%0d".join(urllib.parse.quote(name) for name in gene_names)
    url = f"{STRING_BASE}/image/network?identifiers={identifiers}&species={species}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read()


def annotate_with_pathways(gene_names: list[str]) -> list[dict]:
    """
    High-level entry point: for each gene name, resolve it to a KEGG gene ID
    (see resolve_kegg_gene_id) and fetch its KEGG pathways and STRING PPI
    partners. Works for any human gene symbol, not a curated set. A gene
    Entrez can't resolve still gets an entry (empty "kegg_pathways") rather
    than being dropped silently - STRING is queried directly by symbol
    regardless, so it can still return partners even when KEGG resolution
    fails, and the caller (index.html) can render "(none found)" either way
    without needing to know which lookup came up empty.
    args:
        gene_names: list[str] - gene names, e.g. from annotate_mvp's "match" field
    returns:
        list[dict] - one entry per gene: {"gene": str,
        "kegg_pathways": [{"id": str, "name": str}, ...],
        "string_partners": [{"partner": str, "score": float}, ...]}
    """
    results = []
    for name in gene_names:
        gene_id = resolve_kegg_gene_id(name)
        if gene_id is not None:
            pathway_ids = fetch_kegg_pathway_ids(gene_id)
            pathway_names = fetch_kegg_pathway_names(pathway_ids)
            kegg_pathways = [{"id": pid, "name": pathway_names.get(pid, pid)} for pid in pathway_ids]
        else:
            kegg_pathways = []
        results.append({
            "gene": name,
            "kegg_pathways": kegg_pathways,
            "string_partners": fetch_string_ppi_partners(name),
        })
    return results


if __name__ == "__main__":
    import sys

    genes = sys.argv[1:] or ["STAT1"]
    for entry in annotate_with_pathways(genes):
        print(f"\n{entry['gene']}:")
        print(f"  KEGG pathways ({len(entry['kegg_pathways'])}):")
        for p in entry["kegg_pathways"]:
            print(f"    {p['id']}: {p['name']}")
        print("  STRING PPI partners:")
        for p in entry["string_partners"]:
            print(f"    {p['partner']} (score {p['score']})")
