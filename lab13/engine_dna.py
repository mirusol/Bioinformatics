"""DNA sequence generator using transition probabilities from transition_matrix.json."""

import json
import random
from pathlib import Path
from typing import Dict, List


def load_transition_matrix(json_path: Path) -> tuple:
    """Load DNA transition probabilities from JSON file."""
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    probabilities = data['probabilities']
    nucleotides = list(probabilities.keys())
    
    return probabilities, nucleotides


def generate_dna_sequence(probabilities: Dict[str, Dict[str, float]], 
                          nucleotides: List[str], 
                          length: int,
                          start_base: str = None) -> str:
    """Generate DNA sequence using transition probabilities.
    
    Args:
        probabilities: Transition probability matrix
        nucleotides: List of valid nucleotides (A, C, G, T)
        length: Desired length of generated sequence
        start_base: Initial nucleotide (random if None)
    
    Returns:
        Generated DNA sequence as string
    """
    
    if length < 1:
        raise ValueError("Length must be at least 1")
    
    # Start with random base if not specified
    if start_base is None:
        current_base = random.choice(nucleotides)
    else:
        if start_base not in nucleotides:
            raise ValueError(f"Invalid start base: {start_base}")
        current_base = start_base
    
    sequence = [current_base]
    
    # Generate remaining bases using transition probabilities
    for _ in range(length - 1):
        # Get probabilities for next base given current base
        next_base_probs = probabilities[current_base]
        
        # Create weighted choice based on probabilities
        bases = list(next_base_probs.keys())
        weights = [next_base_probs[base] for base in bases]
        
        # Select next base using weighted random choice
        next_base = random.choices(bases, weights=weights, k=1)[0]
        
        sequence.append(next_base)
        current_base = next_base
    
    return ''.join(sequence)


def display_statistics(original_seq: str, generated_seq: str):
    """Display comparison statistics between original and generated sequences."""
    
    def count_bases(seq: str) -> Dict[str, int]:
        return {base: seq.count(base) for base in 'ACGT'}
    
    def count_transitions(seq: str) -> Dict[str, int]:
        transitions = {}
        for i in range(len(seq) - 1):
            trans = f"{seq[i]}→{seq[i+1]}"
            transitions[trans] = transitions.get(trans, 0) + 1
        return transitions
    
    print("\n" + "=" * 80)
    print("SEQUENCE COMPARISON")
    print("=" * 80)
    
    print(f"\nOriginal sequence length: {len(original_seq)}")
    print(f"Generated sequence length: {len(generated_seq)}")
    
    print(f"\nOriginal sequence: {original_seq}")
    print(f"Generated sequence: {generated_seq}")
    
    print("\nBase composition:")
    print(f"{'Base':<10} {'Original':<15} {'Generated':<15}")
    print("-" * 40)
    
    orig_counts = count_bases(original_seq)
    gen_counts = count_bases(generated_seq)
    
    for base in 'ACGT':
        orig_pct = (orig_counts[base] / len(original_seq)) * 100
        gen_pct = (gen_counts[base] / len(generated_seq)) * 100
        print(f"{base:<10} {orig_counts[base]:>5} ({orig_pct:>5.1f}%)  {gen_counts[base]:>5} ({gen_pct:>5.1f}%)")
    
    print("\nTop 10 transitions:")
    print(f"{'Transition':<15} {'Original':<15} {'Generated':<15}")
    print("-" * 45)
    
    orig_trans = count_transitions(original_seq)
    gen_trans = count_transitions(generated_seq)
    
    all_trans = set(list(orig_trans.keys()) + list(gen_trans.keys()))
    trans_list = [(t, orig_trans.get(t, 0), gen_trans.get(t, 0)) for t in all_trans]
    trans_list.sort(key=lambda x: x[1], reverse=True)
    
    for trans, orig_count, gen_count in trans_list[:10]:
        print(f"{trans:<15} {orig_count:<15} {gen_count:<15}")


def main():
    # Load transition matrix
    json_path = Path(__file__).with_name("transition_matrix.json")
    
    if not json_path.exists():
        print(f"Error: {json_path} not found!")
        print("Please run ex2.py first to generate the transition matrix.")
        return
    
    print("Loading DNA transition matrix...")
    probabilities, nucleotides = load_transition_matrix(json_path)
    
    print(f"Loaded transition matrix for nucleotides: {', '.join(nucleotides)}")
    
    # Read original sequence
    with open(json_path, 'r') as f:
        data = json.load(f)
    original_sequence = data['sequence']
    
    # Generate new sequences
    print("\n" + "=" * 80)
    print("GENERATING DNA SEQUENCES")
    print("=" * 80)
    
    # Generate sequence of same length as original
    print(f"\nGenerating sequence of length {len(original_sequence)}...")
    generated = generate_dna_sequence(probabilities, nucleotides, len(original_sequence))
    
    display_statistics(original_sequence, generated)
    
    # Generate additional sequences
    print("\n" + "=" * 80)
    print("ADDITIONAL GENERATED SEQUENCES")
    print("=" * 80)
    
    lengths = [20, 50, 100]
    for length in lengths:
        seq = generate_dna_sequence(probabilities, nucleotides, length)
        print(f"\nLength {length}: {seq}")
    
    # Generate sequences starting with each base
    print("\n" + "=" * 80)
    print("SEQUENCES STARTING WITH EACH BASE (length 30)")
    print("=" * 80)
    
    for base in nucleotides:
        seq = generate_dna_sequence(probabilities, nucleotides, 30, start_base=base)
        print(f"\nStarting with {base}: {seq}")


if __name__ == "__main__":
    main()
