"""
Helpers for building demo read data from real reference gene sequences (e.g.
STAT1), for the HIV-encephalitis (HIVE) DEG demo. Unlike
test_pipeline_end_to_end.py's synthetic genes (engineered to share zero
k-mers - a trustworthy but artificial fixture), these sequences are real
human mRNA, useful for validating the pipeline's k-mer parameters and
runtime at real-gene scale, and for seeding the web demo with a
recognizable example.
"""
from annotate_mvp import load_reference_genes


def load_gene_sequence(filepath: str) -> str:
    """
    Load a single-record FASTA file's sequence (e.g. a raw NCBI RefSeq mRNA
    download).
    args:
        filepath: str - path to a FASTA file with exactly one record
    returns:
        str - the record's sequence
    """
    genes = load_reference_genes(filepath)
    return next(iter(genes.values()))


def tile_reads(sequence: str, read_len: int, depth: int, step: int = 1) -> list[str]:
    """
    Slide a window over sequence and repeat the tiling `depth` times, to
    simulate a given sequencing depth (same approach as
    test_pipeline_end_to_end.py's synthetic fixture, applied here to real
    gene sequences instead).
    args:
        sequence: str - the source sequence to tile reads from
        read_len: int - length of each simulated read
        depth: int - how many times to repeat the full tiling
        step: int - distance to slide the window between reads; step=1
            (the default) covers every position, matching the
            test_pipeline_end_to_end.py fixture's density. A larger step
            gives a smaller, less redundant read set - useful for a
            human-readable UI example - as long as it stays small enough
            that consecutive reads still overlap (step < read_len - k + 1
            for whatever k the pipeline will use), so Inchworm/Chrysalis
            still see unbroken k-mer overlap between them.
    returns:
        list[str] - simulated reads
    """
    tiles = [sequence[i:i + read_len] for i in range(0, len(sequence) - read_len + 1, step)]
    return tiles * depth


if __name__ == "__main__":
    # Quick local check for "does this fragment actually look like STAT1?" -
    # no deployment, no reads/pipeline needed, just the same containment
    # scoring Annotate uses. Useful when manually picking a fragment out of
    # NCBI to paste into the web demo, e.g.:
    #   python hive_demo_data.py ACTGACTG...(your fragment)
    # With no argument, checks a slice of STAT1 against itself as a sanity check.
    import sys

    from annotate_mvp import best_match

    K = 21
    reference_genes = {"STAT1": load_gene_sequence("_stat1_raw.fasta")}
    fragment = sys.argv[1].strip().upper() if len(sys.argv) > 1 else reference_genes["STAT1"][:200]

    name, score = best_match(fragment, reference_genes, K)
    print(f"Best match: {name or '(no match)'} (containment score {score:.3f} at k={K})")
