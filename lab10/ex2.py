import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List, Dict
import os
import zipfile
from pathlib import Path
from multiprocessing import Pool, cpu_count
import time



def read_fasta(filename: str) -> Dict[str, str]:
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def calculate_cg_content(sequence: str) -> float:
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    if len(sequence) == 0: return 0.0
    return ((c_count + g_count) / len(sequence)) * 100

def calculate_index_of_coincidence(sequence: str) -> float:
    n = len(sequence)
    if n == 0:
        return 0.0
    
    counts = {'A': sequence.count('A'), 'C': sequence.count('C'), 
              'G': sequence.count('G'), 'T': sequence.count('T')}
    
    expected = n / 4
    chi_sq = sum((count - expected) ** 2 / expected for count in counts.values())
    
    return chi_sq

def sliding_window_analysis(sequence: str, window_size: int) -> Tuple[List[float], List[float]]:
    cg_values = []
    ic_values = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        cg_values.append(calculate_cg_content(window))
        ic_values.append(calculate_index_of_coincidence(window))
    
    return cg_values, ic_values

def calculate_center_of_weight(cg_values: List[float], ic_values: List[float]) -> Tuple[float, float]:
    weights = np.arange(1, len(cg_values) + 1)
    if len(weights) == 0: return 0.0, 0.0
    cg_center = np.average(cg_values, weights=weights)
    ic_center = np.mean(ic_values)
    return cg_center, ic_center

def plot_pattern(cg_values: List[float], ic_values: List[float], 
                 title: str, filename: str):
    plt.figure(figsize=(8, 5))
    
    colors = np.arange(len(cg_values))
    scatter = plt.scatter(cg_values, ic_values, c=colors, cmap='viridis', 
                          alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    plt.colorbar(scatter, label='Window Position')
    plt.xlabel('C+G Content (%)', fontsize=10)
    plt.ylabel('Index of Coincidence (χ²)', fontsize=10)
    plt.title(title, fontsize=11, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
    plt.scatter([center_cg], [center_ic], color='red', s=200, 
               marker='X', linewidths=2, label=f'Center ({center_cg:.2f}, {center_ic:.2f})', 
               edgecolors='darkred', zorder=5)
    plt.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

def analyze_promoter(sequence: str, window_size: int = 30) -> dict:
    cg_values, ic_values = sliding_window_analysis(sequence, window_size)
    center = calculate_center_of_weight(cg_values, ic_values)
    
    return {
        'cg_values': cg_values,
        'ic_values': ic_values,
        'center': center,
        'avg_cg': center[0],
        'avg_ic': center[1],
        'sequence_length': len(sequence),
        'window_count': len(cg_values)
    }

def process_single_sequence(args):
    idx, seq_id, sequence, window_size, output_dir = args
    
    if len(sequence) < window_size:
        return None
    
    results = analyze_promoter(sequence, window_size)
    
    safe_id = "".join(c if c.isalnum() or c in (' ', '', '-') else '' for c in seq_id)
    safe_id = safe_id[:50] # Shorten filename
    
    output_file = "Skipped (Limit Reached)"
    
    if idx <= 50:
        output_file = os.path.join(output_dir, f"pattern_{idx:05d}_{safe_id}.png")
        plot_pattern(results['cg_values'], results['ic_values'],
                    title=f"DNA Pattern - Sequence {idx}\n{seq_id[:40]}...",
                    filename=output_file)
    
    return {
        'idx': idx,
        'seq_id': seq_id,
        'center': results['center'],
        'results': results,
        'output_file': output_file
    }

def plot_all_centers(centers: List[Tuple[float, float]], filename: str):
    plt.figure(figsize=(12, 8))
    
    centers_cg = [c[0] for c in centers]
    centers_ic = [c[1] for c in centers]
    
    plt.scatter(centers_cg, centers_ic, s=30, alpha=0.5, c='blue', edgecolors='black', linewidth=0.3)
    
    plt.xlabel('C+G Content (%)', fontsize=12)
    plt.ylabel('Index of Coincidence (χ²)', fontsize=12)
    plt.title(f'Centers of Weight - All Sequences (n={len(centers)})', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_zip_archive(source_dir, output_filename):
    print(f"Creating ZIP archive: {output_filename}...")
    with zipfile.ZipFile(output_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, os.path.dirname(source_dir))
                zipf.write(file_path, arcname)
    
    size_mb = os.path.getsize(output_filename) / (1024 * 1024)
    print(f"✓ ZIP created successfully. Size: {size_mb:.2f} MB")

def process_fasta_file_parallel(fasta_file: str, window_size: int = 30, output_dir: str = "dna_patterns_output", n_cores: int = None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Loading sequences from {fasta_file}...")
    try:
        sequences = read_fasta(fasta_file)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    print(f"Loaded {len(sequences)} sequences\n")
    
    if n_cores is None:
        n_cores = cpu_count()
    
    print(f"Using {n_cores} CPU cores for parallel processing")
    print("NOTE: Individual plot generation limited to first 50 sequences to save space.\n")
    
    args_list = [(idx, seq_id, sequence, window_size, output_dir) 
                 for idx, (seq_id, sequence) in enumerate(sequences.items(), 1)]
    
    start_time = time.time()
    
    with Pool(processes=n_cores) as pool:
        results_list = []
        total = len(args_list)
        
        for i, result in enumerate(pool.imap_unordered(process_single_sequence, args_list), 1):
            if result:
                results_list.append(result)
            
            if i % 100 == 0:
                print(f"Processed {i}/{total} sequences ({i/total*100:.1f}%)")
    
    elapsed_time = time.time() - start_time
    print(f"\n✓ Processed {len(results_list)} sequences in {elapsed_time:.2f} seconds")
    
    
    all_centers = [r['center'] for r in results_list]
    centers_file = os.path.join(output_dir, "all_centers_summary.png")
    print(f"Generating aggregate centers plot...")
    plot_all_centers(all_centers, centers_file)
    
    summary_file = os.path.join(output_dir, "summary_data.txt")
    print(f"Writing summary file...")
    with open(summary_file, 'w') as f:
        f.write("DNA PATTERN ANALYSIS SUMMARY\n")
        f.write("="*70 + "\n")
        f.write(f"Total sequences: {len(results_list)}\n")
        f.write("Note: Only first 50 sequences have individual plot images.\n")
        f.write("="*70 + "\n\n")
        
        for result in sorted(results_list, key=lambda x: x['idx']):
            f.write(f"Seq {result['idx']}: {result['seq_id'][:50]}\n")
            f.write(f"  Center: ({result['center'][0]:.2f}, {result['center'][1]:.2f})\n")
            if result['idx'] <= 50:
                 f.write(f"  Image: {os.path.basename(result['output_file'])}\n")
            f.write("\n")


if __name__ == "__main__":
    fasta_file = "geme_promoters_complete.fasta"
    
    
    if not os.path.exists(fasta_file):
        fasta_file = input("Enter FASTA file path: ").strip()
    
    if os.path.exists(fasta_file):
        process_fasta_file_parallel(fasta_file)
    else:
        print("File not found.")