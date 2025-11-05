import random
import matplotlib.pyplot as plt
import numpy as np
import re

def read_fasta(file_path):
    """Read DNA sequence from a FASTA file."""
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):  
                sequence += line.strip()
    return sequence

# Common restriction enzymes with their recognition sites
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',      # G^AATTC
    'BamHI': 'GGATCC',      # G^GATCC
    'HindIII': 'AAGCTT',    # A^AGCTT
    'PstI': 'CTGCAG',       # CTGCA^G
    'SmaI': 'CCCGGG',       # CCC^GGG
}

def digest_dna(sequence, enzyme_name, recognition_site):
    """
    Digest DNA sequence with a restriction enzyme.
    Returns list of fragment lengths and positions where cuts occur.
    """
    # Find all positions where the enzyme cuts
    cut_positions = [0]  # Start position
    
    # Find all occurrences of the recognition site
    for match in re.finditer(recognition_site, sequence):
        cut_positions.append(match.start())
    
    cut_positions.append(len(sequence))  # End position
    cut_positions = sorted(set(cut_positions))  # Remove duplicates and sort
    
    # Calculate fragment lengths
    fragments = []
    for i in range(len(cut_positions) - 1):
        fragment_length = cut_positions[i + 1] - cut_positions[i]
        if fragment_length > 0:
            fragments.append(fragment_length)
    
    return fragments, len(cut_positions) - 1

def digest_with_multiple_enzymes(sequence, enzymes_dict):
    """
    Digest DNA with multiple restriction enzymes.
    Returns a dictionary with enzyme names and their resulting fragments.
    """
    results = {}
    
    for enzyme_name, recognition_site in enzymes_dict.items():
        fragments, num_cuts = digest_dna(sequence, enzyme_name, recognition_site)
        results[enzyme_name] = {
            'fragments': fragments,
            'num_cuts': num_cuts,
            'num_fragments': len(fragments)
        }
    
    return results

def draw_gel_electrophoresis_with_enzymes(digestion_results, original_length):
    """
    Create a visual representation of gel electrophoresis for enzyme digestion.
    Each lane represents digestion with a different enzyme.
    """
    fig, ax = plt.subplots(figsize=(12, 10), facecolor='black')
    ax.set_facecolor('black')
    
    enzyme_names = list(digestion_results.keys())
    num_lanes = len(enzyme_names)
    lane_width = 0.6
    
    # Draw the wells (loading points)
    for i in range(num_lanes):
        ax.fill_between([i-lane_width/2, i+lane_width/2], [0, 0], [0.1, 0.1], 
                       color='gray', alpha=0.5)
    
    # Find the maximum fragment length for scaling
    max_length = original_length
    
    # Draw the bands for each enzyme
    for lane_idx, enzyme_name in enumerate(enzyme_names):
        fragments = digestion_results[enzyme_name]['fragments']
        
        # Sort fragments by size (largest first) for better visualization
        sorted_fragments = sorted(fragments, reverse=True)
        
        # Draw each fragment as a band
        for fragment_length in sorted_fragments:
            # Calculate position based on length (longer fragments move slower)
            position = 0.2 + (fragment_length / max_length) * 7
            
            # Draw the band
            ax.fill_between([lane_idx-lane_width/2, lane_idx+lane_width/2], 
                           [position-0.08, position-0.08],
                           [position+0.08, position+0.08],
                           color='white', alpha=0.7, edgecolor='lightgray', linewidth=0.5)
    
    # Add size markers
    markers = [3000, 1500, 1000, 500, 100]
    for size in markers:
        if size <= max_length:
            position = 0.2 + (size / max_length) * 7
            ax.text(-0.8, position, f'{size} bp', color='white', 
                    verticalalignment='center', horizontalalignment='right', fontsize=9)
            ax.axhline(y=position, color='white', alpha=0.2, linestyle='--', linewidth=0.5)
    
    # Customize the plot
    ax.set_ylim(8, -0.5)  # Reverse Y-axis as fragments migrate downward
    ax.set_xlim(-1.2, num_lanes-0.5)
    ax.set_xticks(range(num_lanes))
    ax.set_xticklabels(enzyme_names, rotation=45, ha='right')
    
    # Remove frame
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.title('DNA Gel Electrophoresis - Restriction Enzyme Digestion', 
              color='white', pad=20, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()

def print_digestion_summary(sequence, digestion_results):
    """Print a summary of the digestion results."""
    print("=" * 70)
    print("RESTRICTION ENZYME DIGESTION ANALYSIS")
    print("=" * 70)
    print(f"\nOriginal DNA sequence length: {len(sequence)} bp\n")
    
    for enzyme_name, result in digestion_results.items():
        recognition_site = RESTRICTION_ENZYMES[enzyme_name]
        fragments = result['fragments']
        num_cuts = result['num_cuts']
        
        print(f"\n{enzyme_name} (Recognition site: {recognition_site})")
        print("-" * 70)
        print(f"  Number of cuts: {num_cuts}")
        print(f"  Number of fragments: {len(fragments)}")
        
        if len(fragments) > 0:
            print(f"  Fragment sizes (bp): {sorted(fragments, reverse=True)}")
            print(f"  Largest fragment: {max(fragments)} bp")
            print(f"  Smallest fragment: {min(fragments)} bp")
            print(f"  Average fragment size: {sum(fragments)/len(fragments):.1f} bp")
        else:
            print("  No fragments (enzyme did not cut)")
    
    print("\n" + "=" * 70)

def main():
    # Read the DNA sequence from FASTA file
    sequence = read_fasta(r'dna.fasta')
    
    print(f"\nLoaded DNA sequence: {len(sequence)} bp")
    print(f"First 60 nucleotides: {sequence[:60]}...\n")
    
    # Perform digestion with 5 restriction enzymes
    digestion_results = digest_with_multiple_enzymes(sequence, RESTRICTION_ENZYMES)
    
    # Print detailed summary
    print_digestion_summary(sequence, digestion_results)
    
    # Visualize on gel electrophoresis
    draw_gel_electrophoresis_with_enzymes(digestion_results, len(sequence))

if __name__ == "__main__":
    main()
