"""
This module compares assembled transcripts against a reference set of known
genes to find which known gene each transcript most closely matches, using
shared k-mers as a lightweight (BLAST-lite) similarity score.
"""
from inchworm_mvp import kmers_in_contig, load_contigs, reverse_complement


def load_reference_genes(filepath: str) -> dict[str, str]:
    """
    Load a reference gene set from a FASTA file.
    args:
        filepath: str - path to a FASTA file, one header + sequence per gene
    returns:
        dict[str, str] - gene name (header, without '>') to sequence
    """
    genes: dict[str, str] = {}
    name: str | None = None
    sequence_parts: list[str] = []

    def flush() -> None:
        if name is not None:
            genes[name] = "".join(sequence_parts)

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)
        flush()

    return genes


def kmer_containment(query: str, reference: str, k: int) -> float:
    """
    Score how much of a query sequence's k-mer content is explained by a
    reference sequence. Uses containment rather than Jaccard similarity
    because a transcript can be much shorter than the full gene it matches,
    so the reference having "extra" k-mers shouldn't lower the score.
    args:
        query: str - the sequence being annotated (a transcript)
        reference: str - the reference gene sequence
        k: int - the k-mer length
    returns:
        float - fraction of query k-mers also found in reference, 0.0-1.0
    """
    query_kmers = kmers_in_contig(query, k)
    if not query_kmers:
        return 0.0
    reference_kmers = kmers_in_contig(reference, k)
    return len(query_kmers & reference_kmers) / len(query_kmers)


def best_match(transcript: str, reference_genes: dict[str, str], k: int) -> tuple[str | None, float]:
    """
    Find the reference gene whose sequence best explains a transcript's
    k-mer content, checking both the transcript's given orientation and its
    reverse complement since a transcript can be assembled from either
    strand.
    args:
        transcript: str - the assembled transcript sequence
        reference_genes: dict[str, str] - gene name to reference sequence
        k: int - the k-mer length
    returns:
        tuple[str | None, float] - (best-matching gene name, score), or
        (None, 0.0) if no reference gene shares any k-mer with the transcript
    """
    best_name: str | None = None
    best_score = 0.0
    flipped = reverse_complement(transcript)

    for name, reference in reference_genes.items():
        score = max(
            kmer_containment(transcript, reference, k),
            kmer_containment(flipped, reference, k),
        )
        if score > best_score:
            best_name = name
            best_score = score

    return best_name, best_score


def annotate_transcripts(
    transcripts: list[str],
    reference_genes: dict[str, str],
    k: int,
    min_score: float = 0.5,
) -> list[dict]:
    """
    Annotate a list of transcripts against a reference gene set.
    args:
        transcripts: list[str] - assembled transcript sequences
        reference_genes: dict[str, str] - gene name to reference sequence
        k: int - the k-mer length
        min_score: float - minimum containment score required to report a match
    returns:
        list[dict] - one entry per transcript: {"transcript": str,
        "match": str | None, "score": float}; match is None when the best
        score is below min_score
    """
    annotations = []
    for transcript in transcripts:
        name, score = best_match(transcript, reference_genes, k)
        if score < min_score:
            name = None
        annotations.append({"transcript": transcript, "match": name, "score": round(score, 3)})
    return annotations


if __name__ == "__main__":
    k = 5  # should match inchworm_mvp.py / chrysalis_mvp.py / butterfly_mvp.py
    min_score = 0.5

    transcripts = load_contigs("contigs.fasta")
    reference_genes = load_reference_genes("reference_genes.fasta")

    for annotation in annotate_transcripts(transcripts, reference_genes, k, min_score):
        print(annotation)
