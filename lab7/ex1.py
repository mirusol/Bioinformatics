import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


def read_fasta(filename):
    """Read a FASTA file and extract the DNA sequence."""
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):  # Skip header lines
                sequence += line.upper()
    return sequence


def find_repetitions(sequence, min_length=6, max_length=10):

    all_repetitions = {}
    
    for length in range(min_length, max_length + 1):
        # Extract all subsequences of this length
        repetitions = []
        for i in range(len(sequence) - length + 1):
            subseq = sequence[i:i + length]
            # Only include valid DNA sequences
            if all(base in 'ATGC' for base in subseq):
                repetitions.append(subseq)
        
        # Count occurrences
        counts = Counter(repetitions)
        
        # Only keep sequences that appear more than once
        for subseq, count in counts.items():
            if count > 1:
                if subseq not in all_repetitions:
                    all_repetitions[subseq] = count
    
    return all_repetitions


def plot_repetition_histogram(repetitions):
    """Create a single histogram showing all repetition frequencies."""
    if not repetitions:
        print("No repetitions found!")
        return
    
    # Sort repetitions by frequency (descending) and then by pattern
    sorted_reps = sorted(repetitions.items(), key=lambda x: (-x[1], x[0]))
    
    patterns = [p for p, _ in sorted_reps]
    counts = [c for _, c in sorted_reps]
    
    # Create color map based on pattern length
    colors = []
    for pattern in patterns:
        length = len(pattern)
        if length == 6:
            colors.append('#1f77b4')  # blue
        elif length == 7:
            colors.append('#ff7f0e')  # orange
        elif length == 8:
            colors.append('#2ca02c')  # green
        elif length == 9:
            colors.append('#d62728')  # red
        else:  # 10
            colors.append('#9467bd')  # purple
    
    # Create the figure
    plt.figure(figsize=(16, 10))
    
    bars = plt.bar(range(len(patterns)), counts, color=colors, edgecolor='black', linewidth=0.5)
    
    plt.xlabel('DNA Repetitive Patterns', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency (Number of Occurrences)', fontsize=14, fontweight='bold')
    plt.title('DNA Repetition Frequency Analysis (6-10 bp)', fontsize=16, fontweight='bold', pad=20)
    
    # Set x-axis labels (patterns)
    plt.xticks(range(len(patterns)), patterns, rotation=90, fontfamily='monospace', fontsize=8)
    
    # Add grid for better readability
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label='6 bp'),
        Patch(facecolor='#ff7f0e', label='7 bp'),
        Patch(facecolor='#2ca02c', label='8 bp'),
        Patch(facecolor='#d62728', label='9 bp'),
        Patch(facecolor='#9467bd', label='10 bp')
    ]
    plt.legend(handles=legend_elements, loc='upper right', title='Pattern Length', fontsize=10)
    
    # Add value labels on top of bars (only for taller bars to avoid clutter)
    max_count = max(counts)
    for i, (bar, count) in enumerate(zip(bars, counts)):
        if count > max_count * 0.3:  # Only show labels for bars > 30% of max
            plt.text(i, count, str(count), ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('dna_repetitions_histogram.png', dpi=300, bbox_inches='tight')
    print("✓ Plot saved as: dna_repetitions_histogram.png")
    plt.show()


def print_statistics(sequence, repetitions):
    """Print summary statistics."""
    print("\n" + "="*70)
    print("DNA SEQUENCE ANALYSIS - SUMMARY")
    print("="*70)
    
    print(f"\nSequence Information:")
    print(f"  Length: {len(sequence)} nucleotides")
    print(f"  A: {sequence.count('A')} ({sequence.count('A')/len(sequence)*100:.1f}%)")
    print(f"  T: {sequence.count('T')} ({sequence.count('T')/len(sequence)*100:.1f}%)")
    print(f"  G: {sequence.count('G')} ({sequence.count('G')/len(sequence)*100:.1f}%)")
    print(f"  C: {sequence.count('C')} ({sequence.count('C')/len(sequence)*100:.1f}%)")
    print(f"  GC content: {(sequence.count('G') + sequence.count('C'))/len(sequence)*100:.1f}%")
    
    print(f"\nRepetition Analysis (6-10 bp):")
    print(f"  Total unique repetitive patterns: {len(repetitions)}")
    print(f"  Total repetition occurrences: {sum(repetitions.values())}")
    
    # Group by length
    by_length = {}
    for pattern, count in repetitions.items():
        length = len(pattern)
        if length not in by_length:
            by_length[length] = []
        by_length[length].append((pattern, count))
    
    print(f"\nBreakdown by length:")
    for length in sorted(by_length.keys()):
        patterns = sorted(by_length[length], key=lambda x: x[1], reverse=True)
        total = sum(count for _, count in patterns)
        print(f"  {length}bp: {len(patterns)} unique patterns, {total} total occurrences")
        # Show top 5
        top_5 = patterns[:5]
        print(f"    Top: {', '.join([f'{p}({c}x)' for p, c in top_5])}")
    
    print("\n" + "="*70)


# Main execution
if __name__ == "__main__":
    print("DNA REPETITION ANALYZER")
    print("="*70)
    
    # Read the FASTA file
    filename = "dna.fasta"
    
    try:
        sequence = read_fasta(filename)
        print(f"✓ Loaded sequence from {filename}")
        print(f"  Sequence length: {len(sequence)} bp")
        
        # Find repetitions (6-10 base pairs)
        print("\nSearching for repetitions (6-10 bp)...")
        repetitions = find_repetitions(sequence, min_length=6, max_length=10)
        
        # Print statistics
        print_statistics(sequence, repetitions)
        
        # Create histogram plot
        if repetitions:
            print("\nGenerating histogram...")
            plot_repetition_histogram(repetitions)
            print("\n✓ Analysis complete!")
        else:
            print("\nNo repetitions found in the sequence.")
            
    except FileNotFoundError:
        print(f"Error: Could not find '{filename}'")
        print("Please make sure the file is in the same directory as this script.")
    except Exception as e:
        print(f"Error: {e}")