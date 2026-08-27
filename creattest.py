#In

def read_fasta(path: str) -> dict[str, str]:
    """Läser in en fasta-fil och returnerar en dict med sekvenserna."""
    fasta_dict = {}
    with open(path, "r") as f:
        seq_id = ""
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:]
                seq = ""
            else:
                seq += line
        fasta_dict[seq_id] = seq
    return fasta_dict


def create_test_seq(seq: dict[str, str])->list[str]:
    """Simulera restulat av shotgun sekvensering. Returnerar en lista med 100 slumpmässiga substrängar av sekvensen."""
    import random
    test_seq = []
    for _ in range(100):
        start = random.randint(0, len(seq) - 35)
        test_seq.append(seq[start:start + 35])
        
    return test_seq
def create_txt_file(test_seq: list[str], output_path: str) -> None:
    """Skriver testsekvenserna till en textfil."""
    with open(output_path, "w") as f:
        for seq in test_seq:
            f.write(seq + "\n")
    
if __name__ == "__main__":
    fasta_path = "data/test.fasta"
    fasta_dict = read_fasta(fasta_path)
    for seq_id, seq in fasta_dict.items():
        test_seq = create_test_seq(seq)
        print(f"Test sequences for {seq_id}:")
        for s in test_seq:
            print(s)
    create_txt_file(test_seq, "data/test_sequences.txt")