
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for reliable saving
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import os


def read_fasta(filename: str) -> str:
    """Read DNA sequence from FASTA file."""
    seq_lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq_lines.append(line)
    return ''.join(seq_lines).upper()


def find_repeats(seq: str, min_k: int = 6, max_k: int = 10, min_count: int = 2):
    """Find repeated motifs in DNA sequence."""
    n = len(seq)
    repeats_by_k = {k: defaultdict(int) for k in range(min_k, max_k + 1) if k <= n}
    
    for k in range(min_k, max_k + 1):
        if k > n:
            continue
        for i in range(n - k + 1):
            motif = seq[i:i + k]
            repeats_by_k[k][motif] += 1
    
    # Filter by min_count
    repeated = {k: {m: cnt for m, cnt in d.items() if cnt >= min_count} 
                for k, d in repeats_by_k.items()}
    repeated = {k: d for k, d in repeated.items() if d}
    return repeated


def plot_frequency_distribution(repeated_by_k: dict, filename: str, output_path: str):
    """Create frequency distribution plot for a single file."""
    # Collect all frequencies across all motif lengths
    all_frequencies = []
    for d in repeated_by_k.values():
        all_frequencies.extend(d.values())
    
    if not all_frequencies:
        print(f"  No repeats found in {filename}")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Calculate bins
    min_freq = min(all_frequencies)
    max_freq = max(all_frequencies)
    bins = range(min_freq, max_freq + 2)
    
    # Create histogram
    n, bins_edges, patches = ax.hist(all_frequencies, bins=bins,
                                      color='coral', edgecolor='black', 
                                      alpha=0.75, linewidth=1.2)
    
    # Add value labels on bars
    for i, (count, patch) in enumerate(zip(n, patches)):
        if count > 0:
            height = patch.get_height()
            ax.text(patch.get_x() + patch.get_width()/2., height,
                   f'{int(count)}',
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Styling
    ax.set_xlabel('Number of Occurrences', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of Patterns', fontsize=13, fontweight='bold')
    ax.set_title(f'Frequency Distribution of Repetitive Patterns\n{filename}', 
                 fontsize=15, fontweight='bold', pad=20)
    
    # Grid
    ax.grid(axis='y', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.set_axisbelow(True)
    
    # Tick styling
    ax.tick_params(axis='both', which='major', labelsize=11)
    
    # Statistics box
    total_patterns = len(all_frequencies)
    avg_freq = sum(all_frequencies) / len(all_frequencies)
    max_occ = max(all_frequencies)
    
    stats_text = f'Total patterns: {total_patterns}\nAvg occurrences: {avg_freq:.1f}\nMax occurrences: {max_occ}'
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.8', facecolor='lightyellow', 
                     edgecolor='black', alpha=0.8, linewidth=1.5))
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Plot saved: {output_path}")
    plt.close(fig)


def print_summary(repeated: dict, filename: str):
    """Print summary statistics for a file."""
    if not repeated:
        print(f"  No repeats found in {filename}")
        return
    
    total_motifs = sum(len(d) for d in repeated.values())
    print(f"  Found {total_motifs} repeated motif(s) across lengths {sorted(repeated.keys())}")
    
    for k in sorted(repeated.keys()):
        d = repeated[k]
        print(f"    Length {k}: {len(d)} motif(s)")


def main():
    """Process all 10 influenza FASTA files."""
    # List of influenza files
    influenza_files = [f'influenza_{i}.fasta' for i in range(1, 11)]
    
    print("="*80)
    print("ANALYZING INFLUENZA DNA SEQUENCES FOR REPETITIONS")
    print("="*80)
    
    # Create output directory for plots
    output_dir = 'influenza_plots'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\nCreated output directory: {output_dir}\n")
    
    # Process each file
    for i, filename in enumerate(influenza_files, 1):
        print(f"\n[{i}/10] Processing {filename}...")
        
        if not os.path.exists(filename):
            print(f"  WARNING: File not found - {filename}")
            continue
        
        try:
            # Read sequence
            seq = read_fasta(filename)
            print(f"  Sequence length: {len(seq)} bp")
            
            # Find repeats
            repeated = find_repeats(seq, min_k=6, max_k=10, min_count=2)
            
            # Print summary
            print_summary(repeated, filename)
            
            # Create plot
            output_path = os.path.join(output_dir, f'influenza_{i}_frequency_plot.png')
            plot_frequency_distribution(repeated, filename, output_path)
            
        except Exception as e:
            print(f"  ERROR processing {filename}: {e}")
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print(f"All plots saved in '{output_dir}/' directory")
    print("="*80)


if __name__ == "__main__":
    main()
