
"""
This module contains functions for assembling contigs from reads using k-mers.
By Anders Ekergård, 2026
For more information, see the README.md file.
"""
import datetime

def find_kmers(reads: list[str], k: int) -> list[str]:
    """
    Find all k-mers in a list of reads.
    args:
        reads: list[str] - list of reads
        k: int - the length of each k-mer
    returns:
        list[str] - a list of all k-mers in the reads
    """
    kmers: list[str] = []
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            kmers.append(kmer)
    return kmers

def kmers_in_contig(contig: str, k: int) -> set[str]:
    """
    Find all k-mers in a contig.
    args:
        contig: str - the contig sequence
        k: int - the length of each k-mer
    returns:
        set[str] - a set of all k-mers in the contig
    """
    result: set[str] = set()
    for i in range(len(contig) - k + 1):
        kmer = contig[i:i+k]
        result.add(kmer)
    return result


def shared_kmer_count(contig_a: str, contig_b: str, k: int) -> int:
    """
    Count how many k-mers two contigs have in common. A single shared k-mer
    can happen by chance between otherwise unrelated sequences (e.g. two
    different organisms' transcripts in a metatranscriptomic sample), but a
    genuine overlap spans many k-mers in a row - so this count is the basis
    for telling a real overlap apart from a coincidence.
    args:
        contig_a: str - the first contig sequence
        contig_b: str - the second contig sequence
        k: int - the k-mer length
    returns:
        int - number of k-mers present in both contigs
    """
    kmers_a = kmers_in_contig(contig_a, k)
    kmers_b = kmers_in_contig(contig_b, k)
    return len(kmers_a & kmers_b)


def share_kmer(contig_a: str, contig_b: str, k: int) -> bool:
    """
    Determine whether two contigs share at least one k-mer.
    args:
        contig_a: str - the first contig sequence
        contig_b: str - the second contig sequence
        k: int - the k-mer length
    returns:
        bool - True when the contigs share a k-mer, otherwise False
    """
    return shared_kmer_count(contig_a, contig_b, k) > 0


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """
    Build the reverse complement of a sequence, since a transcript can be
    sequenced from either strand of the cDNA.
    args:
        seq: str - the sequence
    returns:
        str - the reverse complement of the sequence
    """
    return seq.translate(_COMPLEMENT)[::-1]


def add_reverse_complements(reads: list[str]) -> list[str]:
    """
    Extend a list of reads with each read's reverse complement, so k-mer
    counting also picks up reads that were sequenced from the complementary
    strand.
    args:
        reads: list[str] - list of reads
    returns:
        list[str] - the original reads followed by their reverse complements
    """
    return reads + [reverse_complement(read) for read in reads]


def filter_kmers_by_count(kmers_count: dict[str, int], min_count: int = 1) -> dict[str, int]:
    """
    Drop k-mers seen fewer than min_count times, to filter out k-mers that
    were only created by a sequencing error rather than real coverage.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        min_count: int - the minimum number of occurrences required to keep a k-mer
    returns:
        dict[str, int] - kmers_count with low-support k-mers removed
    """
    return {kmer: count for kmer, count in kmers_count.items() if count >= min_count}


def most_common_kmers(kmers: list[str])-> dict[str, int]:
    """    
    Find the most common kemers from a list of kmers
    args:
        kmers: list[str] - list of kmers
    
    returns:
        dict[str, int] - dictionary of kmers and their counts
    """

    most_common_kmers ={}
    for kmer in kmers:
        if kmer in most_common_kmers:
            most_common_kmers[kmer] += 1
        else:
            most_common_kmers[kmer] = 1
    return most_common_kmers


def expand_forward(kmers_count: dict [str, int], contig: str, k: int)-> str:
    """

    Expand in the contig by finding the next kmer that overlaps with the current kmer.
    Each k-mer is used at most once so the assembly cannot get stuck in a cycle.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        contig: str - the current contig sequence
        k: int - the length of each k-mer
    returns:
        str - the expanded contig sequence
    """
    if not kmers_count:
        return ""

    used_kmers: set[str] = set()

    while True:
        candidates = [
            kmer for kmer in kmers_count
            if kmer not in used_kmers and kmer[:-1] == contig[-(k-1):]
        ]
        if not candidates:
            break

        best = max(candidates, key=kmers_count.get)
        used_kmers.add(best)
        contig += best[-1]
    return contig
def expand_backward(kmers_count: dict [str, int], contig: str, k: int)->str:
    """
    Expand in the contig by finding the next kmer that overlaps with the current kmer.
    Each k-mer is used at most once so the assembly cannot get stuck in a cycle.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        contig: str - the current contig sequence
        k: int - the length of each k-mer
    returns:
        str - the expanded contig sequence
    """
    used_kmers: set[str] = set()
    while True:
        candidates = [
            kmer for kmer in kmers_count
            if kmer not in used_kmers and kmer[1:] == contig[:k-1]
        ]
        if not candidates:
            break

        best = max(candidates, key=kmers_count.get)
        used_kmers.add(best)
        contig = best[0] + contig
    return contig
def expand(kmers_count: dict[str, int], k: int) -> list[str]:
    """
    Expand the contig by finding the next kmer that overlaps with the current kmer.
    args:
        kmers_count: dict[str, int] - dictionary of kmers and their counts
        k: int - the length of each k-mer
    returns:
        list[str] - a list of expanded contig sequences
    """
    used_kmers: set[str] = set()
    contigs: list[str] = []

    available = {kmer: c for kmer, c in kmers_count.items() if kmer not in used_kmers}
    while available:
        start = max(available, key=available.get)
        contig = expand_backward(available, start, k) # expand backward first to find the start of the contig
        
        
        contig = expand_forward(available, contig, k) 
        contigs.append(contig)
        used_kmers |= kmers_in_contig(contig, k)
        available = {kmer: c for kmer, c in kmers_count.items() if kmer not in used_kmers}

    return contigs
def load_reads(filepath: str) -> list[str]:
    """
    Load reads from a FASTQ file.
    """
    reads: list[str] = []
    with open(filepath) as f:
        for index, seq in enumerate(f):
            if index % 4 == 1:
                reads.append(seq.strip())
    return reads

def save_contigs(contigs: list[str], filepath: str) -> None:
    """
    Save contigs to a FASTA file.
    """
    with open(filepath, "w") as f:
        for index, seq in enumerate(contigs):
            f.write(f">contig_{index}\n{seq}\n")
    return None

def load_contigs(filepath: str) -> list[str]:
    """
    Load contigs from a FASTA file (as written by save_contigs).
    """
    contigs: list[str] = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                contigs.append(line)
    return contigs

if __name__ == "__main__":
    k = 5  # Example k-mer length
    min_kmer_count = 2  # drop k-mers only seen once -> likely sequencing errors

    reads = add_reverse_complements(load_reads("reads.fastq"))
    kmers_list = find_kmers(reads, k)
    kmers_count = filter_kmers_by_count(most_common_kmers(kmers_list), min_kmer_count)

    results = expand(kmers_count, k)
    print(f"Results: {results}")

    # Fixed filename so Chrysalis/Butterfly can always find the latest run.
    save_contigs(results, "contigs.fasta")
    # Timestamped copy kept alongside it for history/debugging.
    date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    save_contigs(results, f"contigs_{date}.fasta")