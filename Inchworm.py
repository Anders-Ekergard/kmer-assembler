"""Minimal greedy de novo transcript assembler inspired by Inchworm.

The implementation builds contigs by repeatedly extending the current
sequence with the remaining read that shares the largest suffix/prefix
overlap.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple


def _compute_overlap(left: str, right: str, min_overlap: int) -> int:
    """Return the largest overlap length where left suffix matches right prefix."""
    max_len = min(len(left), len(right))
    for overlap in range(min_overlap, max_len + 1):
        if left.endswith(right[:overlap]):
            return overlap
    return 0


def _merge_reads(left: str, right: str, overlap: int) -> str:
    """Merge two reads using a shared overlap."""
    return left + right[overlap:]


def assemble_transcripts(reads: Sequence[str], min_overlap: int = 2) -> List[str]:
    """Assemble a list of reads into one or more transcripts.

    The algorithm is greedy and deterministic:
    1. Start with the first remaining read as a seed contig.
    2. Search the other reads for the largest valid suffix/prefix overlap.
    3. Merge the best match and repeat until no read can extend the contig.
    4. Start a new contig from the remaining reads.
    """
    cleaned_reads = [str(read).strip() for read in reads if str(read).strip()]
    if not cleaned_reads:
        return []

    remaining = list(cleaned_reads)
    transcripts: List[str] = []

    while remaining:
        current = remaining.pop(0)
        extended = True

        while extended:
            extended = False
            best_overlap = 0
            best_index = -1

            for index, candidate in enumerate(remaining):
                overlap = _compute_overlap(current, candidate, min_overlap)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_index = index

            if best_index >= 0 and best_overlap >= min_overlap:
                candidate = remaining.pop(best_index)
                current = _merge_reads(current, candidate, best_overlap)
                extended = True

        transcripts.append(current)

    return transcripts


if __name__ == "__main__":
    example_reads = ["ATG", "TGC", "GCA", "CAT"]
    print(assemble_transcripts(example_reads, min_overlap=2))
