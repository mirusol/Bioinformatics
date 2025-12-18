import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np

# Log-likelihoods matrix from previous exercise
log_matrix = {
    'A': [-0.080,  0.613, -0.486, -1.179, -1.179,  0.613,  0.767, -0.080, -0.486],
    'C': [-0.080, -0.080, -0.486, -1.179, -1.179, -0.080, -1.179, -0.486, -0.080],
    'G': [-0.486, -0.486,  0.767,  1.124, -1.179, -0.080, -0.080,  0.431, -0.486],
    'T': [ 0.431, -0.486, -0.486, -1.179,  1.124, -1.179, -0.486, -0.080,  0.613]
}

def parse_fasta(filename):
    """Parse FASTA file and extract sequences"""
    sequences = []
    current_seq = ""
    current_id = ""
    current_desc = ""

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences.append({
                        'id': current_id,
                        'description': current_desc,
                        'sequence': current_seq.upper()
                    })
                # Start new sequence
                parts = line[1:].split(' ', 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = ""
            else:
                current_seq += line
        if current_id:
            sequences.append({
                'id': current_id,
                'description': current_desc,
                'sequence': current_seq.upper()
            })

    return sequences

def scan_sequence(sequence, log_matrix, motif_length=9):
    scores = []

    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i:i + motif_length]
        score = 0
        valid = True

        for pos, nuc in enumerate(window):
            if nuc in log_matrix:
                score += log_matrix[nuc][pos]
            else:
                valid = False
                break

        if valid:
            scores.append({
                'position': i,
                'window': window,
                'score': score
            })
        else:
            scores.append({
                'position': i,
                'window': window,
                'score': -10
            })

    return scores

def identify_significant_motifs(scores, threshold=0):
    significant = []

    for score_data in scores:
        if score_data['score'] > threshold:
            significant.append(score_data)

    return significant

def create_genome_chart(genome_data, scores, output_filename):
    sequence = genome_data['sequence']
    genome_id = genome_data['id']
    description = genome_data['description']

    positions = [s['position'] for s in scores]
    score_values = [s['score'] for s in scores]

    significant = identify_significant_motifs(scores, threshold=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'Motif Analysis: {genome_id}\n{description[:80]}', fontsize=12, fontweight='bold')

    ax1.plot(positions, score_values, linewidth=0.5, color='steelblue', alpha=0.7)
    ax1.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Threshold (0)')
    ax1.fill_between(positions, 0, score_values, where=[s > 0 for s in score_values],
                      color='green', alpha=0.3, label='Positive scores')
    ax1.fill_between(positions, 0, score_values, where=[s <= 0 for s in score_values],
                      color='gray', alpha=0.2, label='Negative scores')

    if significant:
        top_motifs = sorted(significant, key=lambda x: x['score'], reverse=True)[:10]
        top_positions = [m['position'] for m in top_motifs]
        top_scores = [m['score'] for m in top_motifs]
        ax1.scatter(top_positions, top_scores, color='red', s=50, zorder=5,
                   label=f'Top {len(top_motifs)} motifs', marker='*')

    ax1.set_xlabel('Position in Genome (bp)', fontsize=10)
    ax1.set_ylabel('Log-likelihood Score', fontsize=10)
    ax1.set_title('Motif Score Distribution Across Genome', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=8)

    ax2.set_xlim(0, len(sequence))
    ax2.set_ylim(0, 1)

    ax2.plot([0, len(sequence)], [0.5, 0.5], 'k-', linewidth=3, label='Genome')

    
    if significant:
        for motif in significant:
            pos = motif['position']
            score = motif['score']
            if score > 2:
                color = 'darkred'
                height = 0.3
            elif score > 1:
                color = 'red'
                height = 0.25
            else:
                color = 'orange'
                height = 0.2

            ax2.add_patch(Rectangle((pos, 0.5 - height/2), 9, height,
                                    facecolor=color, edgecolor='none', alpha=0.6))

    ax2.set_xlabel('Genome Position (bp)', fontsize=10)
    ax2.set_yticks([])
    ax2.set_title(f'Functional Motif Locations (n={len(significant)} motifs with score > 0)', fontsize=11)

    red_patch = mpatches.Patch(color='darkred', label='High (score > 2)')
    orange_patch = mpatches.Patch(color='red', label='Medium (1-2)')
    yellow_patch = mpatches.Patch(color='orange', label='Low (0-1)')
    ax2.legend(handles=[red_patch, orange_patch, yellow_patch], loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    plt.close()

    return significant

def generate_summary_report(all_results):

    print("\n" + "="*80)
    print("SUMMARY REPORT - ALL GENOMES")
    print("="*80)

    for result in all_results:
        genome_id = result['genome_id']
        total_windows = result['total_windows']
        significant_count = result['significant_count']
        max_score = result['max_score']
        top_motif = result['top_motif']

        print(f"\n{genome_id}:")
        print(f"  Genome length: {result['genome_length']:,} bp")
        print(f"  Windows analyzed: {total_windows:,}")
        print(f"  Significant motifs (score > 0): {significant_count}")
        print(f"  Max score: {max_score:.3f}")
        if top_motif:
            print(f"  Best motif: '{top_motif['window']}' at position {top_motif['position']:,}")
            print(f"  Chart saved: {result['chart_file']}")

def main():
    print("INFLUENZA GENOME MOTIF SCANNER")
    genomes = parse_fasta("influenza_genomes.fasta")
    print(f"    [OK] Loaded {len(genomes)} genomes")
    all_results = []

    for i, genome in enumerate(genomes, 1):
        genome_id = genome['id']
        sequence = genome['sequence']

        print(f"\n    [{i}/{len(genomes)}] Analyzing {genome_id} ({len(sequence):,} bp)...")

        scores = scan_sequence(sequence, log_matrix)
        significant = identify_significant_motifs(scores, threshold=0)
        if scores:
            max_score_data = max(scores, key=lambda x: x['score'])
            max_score = max_score_data['score']
            top_motif = max_score_data if max_score > 0 else None
        else:
            max_score = 0
            top_motif = None

        print(f"        - Total windows: {len(scores):,}")
        print(f"        - Significant motifs: {len(significant)}")
        print(f"        - Max score: {max_score:.3f}")

        chart_filename = f"genome_{i}_{genome_id}_motif_chart.png"
        print(f"        - Creating chart: {chart_filename}")
        create_genome_chart(genome, scores, chart_filename)

        all_results.append({
            'genome_id': genome_id,
            'genome_length': len(sequence),
            'total_windows': len(scores),
            'significant_count': len(significant),
            'max_score': max_score,
            'top_motif': top_motif,
            'chart_file': chart_filename
        })

    generate_summary_report(all_results)
    
if __name__ == "__main__":
    main()
