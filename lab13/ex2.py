import json
from pathlib import Path
from typing import Dict, Tuple


NUCLEOTIDES = ("A", "C", "G", "T")


def compute_transition_matrix(sequence: str) -> Tuple[Dict[str, Dict[str, int]], Dict[str, Dict[str, float]]]:

    clean_sequence = sequence.strip().upper()
    if len(clean_sequence) < 2:
        raise ValueError("Sequence must have at least 2 nucleotides.")
    if any(base not in NUCLEOTIDES for base in clean_sequence):
        raise ValueError("Sequence must contain only A, C, G, T characters.")

    counts: Dict[str, Dict[str, int]] = {src: {dst: 0 for dst in NUCLEOTIDES} for src in NUCLEOTIDES}

    for current_base, next_base in zip(clean_sequence, clean_sequence[1:]):
        counts[current_base][next_base] += 1

    probabilities: Dict[str, Dict[str, float]] = {}
    for src in NUCLEOTIDES:
        row_total = sum(counts[src].values())
        probabilities[src] = {
            dst: (counts[src][dst] / row_total if row_total else 0.0)
            for dst in NUCLEOTIDES
        }

    return counts, probabilities


def save_transition_data(sequence: str, output_path: Path) -> None:
    counts, probabilities = compute_transition_matrix(sequence)

    output = {
        "sequence": sequence,
        "length": len(sequence),
        "counts": counts,
        "probabilities": probabilities,
    }

    output_path.write_text(json.dumps(output, indent=2))


if __name__ == "__main__":
    dna_sequence = (
        "ACGTACGTACGTAGCTAGCTTGCATGACGTAGCTAGTCGATACGATCGTA"
    ) 

    destination = Path(__file__).with_name("transition_matrix.json")
    save_transition_data(dna_sequence, destination)
    print(f"Transition matrix saved to {destination}")
