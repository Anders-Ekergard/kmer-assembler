"""
This module contains functions for reconstructing transcripts from a
Chrysalis component graph using a greedy, read-supported walk.
"""
from inchworm_mvp import load_contigs, load_reads, reverse_complement
from chrysalis_mvp import cluster_kmers


def edge_support_from_reads(graph: dict[str, set[str]], reads: list[str], k: int) -> dict[tuple[str, str], int]:
    """
    Count how many reads confirm each edge in the graph.
    A read confirms an edge kmer_a -> kmer_b if it contains both k-mers
    back to back (i.e. overlapping by k-1). Each read is checked in its
    given orientation first; only if that finds nothing is its reverse
    complement tried, since a read comes from one physical strand and
    checking both would double-count it.
    args:
        graph: dict[str, set[str]] - the component graph from cluster_kmers
        reads: list[str] - the original reads
        k: int - the k-mer length used throughout
    returns:
        dict[tuple[str, str], int] - number of confirming reads per edge
    """
    def count_matches(read: str) -> dict[tuple[str, str], int]:
        matches: dict[tuple[str, str], int] = {}
        for i in range(len(read) - k):
            kmer_a = read[i:i + k]
            kmer_b = read[i + 1:i + 1 + k]
            if kmer_b in graph.get(kmer_a, ()):
                edge = (kmer_a, kmer_b)
                matches[edge] = matches.get(edge, 0) + 1
        return matches

    support: dict[tuple[str, str], int] = {}
    for read in reads:
        matches = count_matches(read)
        if not matches:
            matches = count_matches(reverse_complement(read))
        for edge, count in matches.items():
            support[edge] = support.get(edge, 0) + count

    return support


def find_start_nodes(graph: dict[str, set[str]]) -> list[str]:
    """
    Find nodes with no incoming edges within the graph.
    args:
        graph: dict[str, set[str]] - the component graph
    returns:
        list[str] - nodes that are never the target of an edge
    """
    has_incoming: set[str] = set()
    for neighbors in graph.values():
        has_incoming |= neighbors
    return [node for node in graph if node not in has_incoming]


def find_end_nodes(graph: dict[str, set[str]]) -> list[str]:
    """
    Find nodes with no outgoing edges within the graph.
    args:
        graph: dict[str, set[str]] - the component graph
    returns:
        list[str] - nodes with no neighbors
    """
    return [node for node, neighbors in graph.items() if not neighbors]


def greedy_walk(
    graph: dict[str, set[str]],
    edge_support: dict[tuple[str, str], int],
    start: str,
    min_edge_support: int = 1,
) -> list[str]:
    """
    Greedily walk the graph from a start node, always following the
    unvisited neighbor with the highest read support. Stops at a dead end,
    when the best remaining candidate is below min_edge_support, or when
    every neighbor has already been visited (repeat/cycle protection).
    args:
        graph: dict[str, set[str]] - the component graph
        edge_support: dict[tuple[str, str], int] - read support per edge
        start: str - the node to walk from
        min_edge_support: int - minimum support required to extend the walk
    returns:
        list[str] - the sequence of nodes visited, starting with start
    """
    path = [start]
    visited = {start}

    current = start
    while True:
        candidates = [n for n in graph.get(current, ()) if n not in visited]
        if not candidates:
            break

        best = max(candidates, key=lambda n: edge_support.get((current, n), 0))
        if edge_support.get((current, best), 0) < min_edge_support:
            break

        path.append(best)
        visited.add(best)
        current = best

    return path


def path_to_transcript(path: list[str]) -> str:
    """
    Reconstruct a transcript sequence from a path of overlapping k-mers.
    args:
        path: list[str] - nodes visited, each overlapping the next by k-1
    returns:
        str - the reconstructed sequence
    """
    if not path:
        return ""
    return path[0] + "".join(node[-1] for node in path[1:])


def assemble_component_transcripts(
    graph: dict[str, set[str]],
    reads: list[str],
    k: int,
    min_edge_support: int = 1,
) -> list[str]:
    """
    Assemble transcripts for a single component graph.
    args:
        graph: dict[str, set[str]] - the component graph from cluster_kmers
        reads: list[str] - the original reads
        k: int - the k-mer length used throughout
        min_edge_support: int - minimum read support required to extend a walk
    returns:
        list[str] - reconstructed transcript sequences, deduplicated
    """
    if not graph:
        return []

    edge_support = edge_support_from_reads(graph, reads, k)
    starts = find_start_nodes(graph) or [min(graph)]

    transcripts: list[str] = []
    seen: set[str] = set()
    for start in starts:
        path = greedy_walk(graph, edge_support, start, min_edge_support)

        is_dead_end = not graph.get(start)
        if len(path) == 1 and not is_dead_end:
            continue

        transcript = path_to_transcript(path)
        if transcript not in seen:
            seen.add(transcript)
            transcripts.append(transcript)

    return transcripts


def run_butterfly(
    graphs: list[dict[str, set[str]]],
    reads: list[str],
    k: int,
    min_edge_support: int = 1,
) -> list[str]:
    """
    Assemble transcripts across all component graphs.
    args:
        graphs: list[dict[str, set[str]]] - component graphs from cluster_kmers
        reads: list[str] - the original reads
        k: int - the k-mer length used throughout
        min_edge_support: int - minimum read support required to extend a walk
    returns:
        list[str] - reconstructed transcript sequences, deduplicated across components
    """
    transcripts: list[str] = []
    seen: set[str] = set()
    for graph in graphs:
        for transcript in assemble_component_transcripts(graph, reads, k, min_edge_support):
            if transcript not in seen:
                seen.add(transcript)
                transcripts.append(transcript)

    return transcripts


if __name__ == "__main__":
    k = 5  # should match inchworm_mvp.py / chrysalis_mvp.py
    min_edge_support = 1

    reads = load_reads("reads.fastq")
    contigs = load_contigs("contigs.fasta")
    graphs = cluster_kmers(contigs, k)

    transcripts = run_butterfly(graphs, reads, k, min_edge_support)
    print(f"Found {len(transcripts)} transcript(s)")
    for transcript in transcripts:
        print(transcript)
