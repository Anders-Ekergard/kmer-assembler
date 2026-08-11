"""
This module quantifies transcript expression per sample by counting how
many reads from each sample support each assembled transcript, and compares
expression between two samples (e.g. two stimuli) as a log2 fold change.
"""
import math

from inchworm_mvp import reverse_complement


def count_read_support(transcript: str, reads: list[str]) -> int:
    """
    Count how many reads from a sample originate from a transcript, checking
    each read in its given orientation first and only trying its reverse
    complement if that finds no match, since a read comes from one physical
    strand and checking both would double-count it.
    args:
        transcript: str - an assembled transcript sequence
        reads: list[str] - reads from a single sample
    returns:
        int - number of reads that are a substring of the transcript
    """
    count = 0
    for read in reads:
        if read in transcript or reverse_complement(read) in transcript:
            count += 1
    return count


def reads_per_kilobase(count: int, transcript_length: int) -> float:
    """
    Normalize a raw read count by transcript length, so a longer transcript
    doesn't look more highly expressed purely because it can contain more
    reads. This does NOT correct for sequencing depth (how many reads the
    sample got in total) - two samples sequenced to different depths still
    aren't comparable on RPK alone. See add_tpm for that correction.
    args:
        count: int - raw read count for a transcript in one sample
        transcript_length: int - length of the transcript in bases
    returns:
        float - reads per kilobase (RPK); 0.0 for a zero-length transcript
    """
    if transcript_length == 0:
        return 0.0
    return count / (transcript_length / 1000)


def quantify_transcripts(transcripts: list[str], samples: dict[str, list[str]]) -> list[dict]:
    """
    Quantify each transcript's expression across samples.
    args:
        transcripts: list[str] - assembled transcript sequences
        samples: dict[str, list[str]] - sample name to that sample's reads
    returns:
        list[dict] - one entry per transcript: {"transcript": str,
        "length": int, "counts": {sample: int}, "rpk": {sample: float}}
    """
    results = []
    for transcript in transcripts:
        counts = {sample: count_read_support(transcript, reads) for sample, reads in samples.items()}
        rpk = {sample: reads_per_kilobase(count, len(transcript)) for sample, count in counts.items()}
        results.append({"transcript": transcript, "length": len(transcript), "counts": counts, "rpk": rpk})
    return results


def add_tpm(quantified: list[dict]) -> list[dict]:
    """
    Add a "tpm" (transcripts per million) column to each quantified entry.
    TPM rescales each transcript's RPK by that sample's total RPK, so it
    corrects for sequencing depth as well as transcript length: a sample
    with twice as many reads doesn't make every transcript look twice as
    expressed, since every transcript in that sample is divided by the same
    inflated total. This is what makes a transcript's expression comparable
    across samples - RPK alone is not, since it only removes the length bias.
    args:
        quantified: list[dict] - output of quantify_transcripts
    returns:
        list[dict] - quantified entries plus a "tpm" dict per entry, keyed
        the same way as "rpk"; each sample's TPM values sum to ~1,000,000
    """
    if not quantified:
        return []

    samples = quantified[0]["rpk"].keys()
    totals = {sample: sum(entry["rpk"][sample] for entry in quantified) for sample in samples}

    return [
        {
            **entry,
            "tpm": {
                sample: (entry["rpk"][sample] / totals[sample] * 1_000_000 if totals[sample] else 0.0)
                for sample in samples
            },
        }
        for entry in quantified
    ]


def log2_fold_change(value_a: float, value_b: float, pseudocount: float = 1.0) -> float:
    """
    Compare a transcript's normalized expression between two samples on a
    log2 scale, since fold-change (not raw difference) is the standard way
    to compare expression magnitudes that can span orders of magnitude. A
    pseudocount avoids taking log2(0) when a transcript is unseen in a sample.
    args:
        value_a: float - normalized expression (e.g. TPM) in sample A
        value_b: float - normalized expression (e.g. TPM) in sample B
        pseudocount: float - added to both values before taking their ratio
    returns:
        float - log2(value_b + pseudocount) - log2(value_a + pseudocount);
        positive means higher in sample B, negative means higher in sample A
    """
    return math.log2(value_b + pseudocount) - math.log2(value_a + pseudocount)


def compare_samples(quantified: list[dict], sample_a: str, sample_b: str) -> list[dict]:
    """
    Add a log2 fold-change column comparing two samples to each quantified
    transcript, sorted by absolute fold-change so the biggest differences
    surface first - a starting point for spotting differentially expressed
    transcripts, not a statistically validated call (that needs replicates).
    Compares TPM rather than raw RPK, so the two samples are comparable even
    if they were sequenced to different depths.
    args:
        quantified: list[dict] - output of quantify_transcripts
        sample_a: str - baseline sample name
        sample_b: str - comparison sample name
    returns:
        list[dict] - quantified entries plus "tpm" and "log2fc", sorted by
        |log2fc| descending
    """
    normalized = add_tpm(quantified)
    compared = [
        {**entry, "log2fc": round(log2_fold_change(entry["tpm"][sample_a], entry["tpm"][sample_b]), 3)}
        for entry in normalized
    ]
    return sorted(compared, key=lambda entry: abs(entry["log2fc"]), reverse=True)


if __name__ == "__main__":
    from inchworm_mvp import load_contigs, load_reads

    transcripts = load_contigs("contigs.fasta")
    samples = {
        "stimulus_a": load_reads("reads_stimulus_a.fastq"),
        "stimulus_b": load_reads("reads_stimulus_b.fastq"),
    }

    quantified = quantify_transcripts(transcripts, samples)
    compared = compare_samples(quantified, "stimulus_a", "stimulus_b")
    for entry in compared:
        print(entry)
