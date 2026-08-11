from inchworm_mvp import kmers_in_contig, load_contigs, reverse_complement, shared_kmer_count


def find(parent: dict[int, int], x: int) -> int:
    while parent[x] != x:
        x = parent[x]
    return x


def union(parent: dict[int, int], x: int, y: int) -> None:
    root_x = find(parent, x)
    root_y = find(parent, y)
    if root_x != root_y:
        parent[root_x] = root_y


def group_contigs_by_kmer_overlap(contigs: list[str], k: int, min_shared_kmers: int = 1) -> list[list[str]]:
    """
    Group contigs into connected components that share at least
    min_shared_kmers k-mers, transitively (union-find over contig indices).
    A contig also overlaps another if it meets that threshold against that
    contig's reverse complement, since the same transcript can be assembled
    from either strand.

    min_shared_kmers defaults to 1 (any shared k-mer merges), matching the
    original behavior. Raising it guards against merging two otherwise
    unrelated contigs - e.g. different organisms' transcripts in a
    metatranscriptomic sample - that happen to share a single k-mer by
    chance: a coincidence touches one k-mer, a genuine overlap spans many.
    args:
        contigs: list[str] - contig sequences to group
        k: int - the k-mer length
        min_shared_kmers: int - minimum number of shared k-mers required to
            merge two contigs into the same component
    """
    parent = {i: i for i in range(len(contigs))}

    for i, contig_a in enumerate(contigs):
        for j, contig_b in enumerate(contigs):
            if i < j and (
                shared_kmer_count(contig_a, contig_b, k) >= min_shared_kmers
                or shared_kmer_count(contig_a, reverse_complement(contig_b), k) >= min_shared_kmers
            ):
                union(parent, i, j)

    groups: dict[int, list[str]] = {}
    for i, contig in enumerate(contigs):
        root = find(parent, i)
        groups.setdefault(root, []).append(contig)

    return list(groups.values())


def orient_component(contigs: list[str], k: int) -> list[str]:
    """
    Normalize a connected component's contigs to one consistent strand
    orientation, so a de Bruijn graph built from them has coherent edges.
    args:
        contigs: list[str] - contigs belonging to a single component
        k: int - the k-mer length
    returns:
        list[str] - the same contigs, each oriented to overlap (by k-mer)
        with the k-mers already collected from earlier contigs in the list
    """
    if not contigs:
        return []

    oriented = [contigs[0]]
    known_kmers = kmers_in_contig(contigs[0], k)
    remaining = contigs[1:]

    while remaining:
        still_remaining: list[str] = []
        progressed = False
        for contig in remaining:
            contig_kmers = kmers_in_contig(contig, k)
            if contig_kmers & known_kmers:
                chosen, chosen_kmers = contig, contig_kmers
            else:
                flipped = reverse_complement(contig)
                flipped_kmers = kmers_in_contig(flipped, k)
                if flipped_kmers & known_kmers:
                    chosen, chosen_kmers = flipped, flipped_kmers
                else:
                    still_remaining.append(contig)
                    continue

            oriented.append(chosen)
            known_kmers |= chosen_kmers
            progressed = True

        if not progressed:
            # Shouldn't happen for a genuinely connected component, but
            # avoid looping forever - keep any leftovers as-is.
            oriented.extend(still_remaining)
            break
        remaining = still_remaining

    return oriented


def build_de_bruijn_graph(contigs: list[str], k: int) -> dict[str, set[str]]:
    """
    Build a de Bruijn-style graph for a group of contigs.
    Nodes are k-mers; a directed edge k-mer A -> k-mer B exists whenever B
    immediately follows A inside one of the contigs (k-1 overlap).
    """
    graph: dict[str, set[str]] = {}
    for contig in contigs:
        kmers = [contig[i:i + k] for i in range(len(contig) - k + 1)]
        for kmer in kmers:
            graph.setdefault(kmer, set())
        
        for kmer_a, kmer_b in zip(kmers, kmers[1:]):
            graph[kmer_a].add(kmer_b)
    return graph


def cluster_kmers(contigs: list[str], k: int, min_shared_kmers: int = 1) -> list[dict[str, set[str]]]:
    """
    Cluster contigs by shared k-mers and build a de Bruijn graph for each
    resulting component.
    args:
        contigs: list[str] - contig sequences (e.g. from Inchworm)
        k: int - the k-mer length
        min_shared_kmers: int - minimum shared k-mers required to merge two
            contigs into the same component (see group_contigs_by_kmer_overlap)
    returns:
        list[dict[str, set[str]]] - one de Bruijn graph per component,
        each graph mapping a k-mer node to the set of k-mers that follow it
    """
    components = group_contigs_by_kmer_overlap(contigs, k, min_shared_kmers)
    oriented_components = [orient_component(component, k) for component in components]
    return [build_de_bruijn_graph(component, k) for component in oriented_components]

if __name__ == "__main__":
    k = 5  # Example k-mer length, should match the value used in inchworm_mvp.py
    
    contigs = load_contigs("contigs.fasta")

    graphs = cluster_kmers(contigs, k)
    print(f"Found {len(graphs)} component(s)")
    for graph in graphs:
        print(graph)