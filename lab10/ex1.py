import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List


def calculate_cg_percentage(sequence: str) -> float:
    if len(sequence) == 0:
        return 0.0
    
    c_count = sequence.upper().count('C')
    g_count = sequence.upper().count('G')
    total = len(sequence)
    
    cg_percentage = ((c_count + g_count) / total) * 100
    return round(cg_percentage, 2)


def calculate_kappa_ic(sequence: str) -> float:
    if len(sequence) <= 1:
        return 0.0
    
    sequence = sequence.upper()
    N = len(sequence)
    
    counts = {
        'A': sequence.count('A'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'T': sequence.count('T')
    }
    
    sum_ni_ni_minus_1 = sum(n * (n - 1) for n in counts.values())
    ic_observed = (sum_ni_ni_minus_1 / (N * (N - 1))) * 100
   
    kappa_ic = ic_observed - 25.0
    
    return round(kappa_ic, 2)


def sliding_window_analysis(sequence: str, window_size: int) -> Tuple[List[float], List[float], List[int]]:
    cg_values = []
    ic_values = []
    positions = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        cg = calculate_cg_percentage(window)
        ic = calculate_kappa_ic(window)
        
        cg_values.append(cg)
        ic_values.append(ic)
        positions.append(i)
    
    return cg_values, ic_values, positions


def calculate_center_of_weight(cg_values: List[float], ic_values: List[float]) -> Tuple[float, float]:
    if len(cg_values) == 0 or len(ic_values) == 0:
        return (0.0, 0.0)
    
    center_cg = np.mean(cg_values)
    center_ic = np.mean(ic_values)
    
    return (round(center_cg, 2), round(center_ic, 2))


def plot_pattern(cg_values: List[float], ic_values: List[float], 
                 title: str = "DNA Pattern (CG% vs IC)", 
                 show_center: bool = True):
    plt.figure(figsize=(10, 8))
    plt.scatter(cg_values, ic_values, alpha=0.6, s=30, c='blue', label='Windows')
    
    if show_center:
        center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
        plt.scatter([center_cg], [center_ic], c='red', s=200, marker='.', 
                   label=f'Center ({center_cg}, {center_ic})', zorder=5)
    
    plt.xlabel('C+G %', fontsize=12)
    plt.ylabel('Kappa IC', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_pattern_centers(centers: List[Tuple[float, float]], labels: List[str] = None):
    plt.figure(figsize=(10, 8))
    
    if labels is None:
        labels = [f"Pattern {i+1}" for i in range(len(centers))]
    
    for i, (cg, ic) in enumerate(centers):
        plt.scatter([cg], [ic], s=200, marker='.', label=labels[i])
        plt.annotate(labels[i], (cg, ic), xytext=(5, 5), 
                    textcoords='offset points', fontsize=10)
    
    plt.xlabel('C+G %', fontsize=12)
    plt.ylabel('Kappa IC', fontsize=12)
    plt.title('Pattern Centers', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    test_sequence = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    window_size = 30

    print("DNA PATTERN ANALYSIS")
    print(f"\nSequence: {test_sequence}")
    print(f"Length: {len(test_sequence)} bp")
    print(f"Window size: {window_size} bp")
    
    print("\n" + "-" * 60)
    print("Step 3: Testing CG% calculation")
    print("-" * 60)
    test_cg = calculate_cg_percentage(test_sequence)
    print(f"CG% for entire sequence: {test_cg}%")
    print(f"Expected: 29.27%")
    
    print("\n" + "-" * 60)
    print("Step 4: Testing Kappa IC calculation")
    print("-" * 60)
    test_ic = calculate_kappa_ic(test_sequence)
    print(f"IC for entire sequence: {test_ic}")
    print(f"Expected: 27.53")
    
    print("\n" + "-" * 60)
    print("Performing sliding window analysis...")
    print("-" * 60)
    cg_values, ic_values, positions = sliding_window_analysis(test_sequence, window_size)
    
    print(f"Number of windows: {len(cg_values)}")
    print(f"CG% range: {min(cg_values):.2f}% - {max(cg_values):.2f}%")
    print(f"IC range: {min(ic_values):.2f} - {max(ic_values):.2f}")
    
    print("\n" + "-" * 60)
    print("Step 6: Calculating center of weight")
    print("-" * 60)
    center_cg, center_ic = calculate_center_of_weight(cg_values, ic_values)
    print(f"Center of weight: CG% = {center_cg}%, IC = {center_ic}")
    
    print("\n" + "-" * 60)
    print("Step 5: Plotting pattern...")
    print("-" * 60)
    plot_pattern(cg_values, ic_values, 
                title=f"DNA Pattern - Test Sequence (Window={window_size}bp)",
                show_center=True)
    
    print("\n" + "-" * 60)
    print("Step 7: Plotting pattern centers...")
    print("-" * 60)
    centers = [(center_cg, center_ic)]
    labels = ["Test Sequence"]
    plot_pattern_centers(centers, labels)
    
    
    promoter_example = "ATATATGCGCGCGCATATATATAGCGCGCGCATATATATATGCGCGCGCATATATAT"
    print(f"\nExample promoter: {promoter_example}")
    cg_p, ic_p, pos_p = sliding_window_analysis(promoter_example, window_size)
    if len(cg_p) > 0:
        center_p = calculate_center_of_weight(cg_p, ic_p)
        print(f"Promoter center: CG% = {center_p[0]}%, IC = {center_p[1]}")
        plot_pattern(cg_p, ic_p, 
                    title=f"DNA Pattern - Example Promoter (Window={window_size}bp)",
                    show_center=True)
        
      
        centers_all = [(center_cg, center_ic), center_p]
        labels_all = ["Test Sequence", "Example Promoter"]
        plot_pattern_centers(centers_all, labels_all)
    
if __name__ == "__main__":
    main()
