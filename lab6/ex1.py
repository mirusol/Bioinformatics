import random
import matplotlib.pyplot as plt
import numpy as np

def read_fasta(file_path):
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):  
                sequence += line.strip()
    return sequence

def generate_random_samples(sequence, num_samples=10, min_length=100, max_length=3000):
    samples = []
    seq_length = len(sequence)
    
    for _ in range(num_samples):
        length = random.randint(min_length, min(max_length, seq_length))
        start = random.randint(0, seq_length - length)
        sample = sequence[start:start + length]
        samples.append((sample, length))
    
    return sorted(samples, key=lambda x: x[1], reverse=True)  

def draw_gel_electrophoresis(samples):
    fig, ax = plt.subplots(figsize=(8, 10), facecolor='black')
    
    ax.set_facecolor('black')
    
    num_samples = len(samples)
    lane_width = 0.6
    
    for i in range(num_samples):
        ax.fill_between([i-lane_width/2, i+lane_width/2], [0, 0], [0.1, 0.1], 
                       color='gray', alpha=0.5)
    
    max_length = max(length for _, length in samples)
    for i, (_, length) in enumerate(samples):
        position = 0.2 + (length / max_length) * 7
        
        ax.fill_between([i-lane_width/2, i+lane_width/2], 
                       [position-0.1, position-0.1],
                       [position+0.1, position+0.1],
                       color='white', alpha=0.8)
    
    markers = [3000, 1500, 500]
    for i, size in enumerate(markers):
        position = 0.2 + (size / max_length) * 7
        ax.text(-1, position, f'{size} bp', color='white', 
                verticalalignment='center', horizontalalignment='right')
        ax.axhline(y=position, color='white', alpha=0.3, linestyle='--', xmin=-0.1, xmax=1.1)
    

    ax.set_ylim(8, -0.5) 
    ax.set_xlim(-1.5, num_samples-0.5)
    ax.set_xticks(range(num_samples))
    ax.set_xticklabels([f'Sample {i+1}' for i in range(num_samples)], rotation=45)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.title('DNA Gel Electrophoresis Simulation', color='white', pad=20)
    plt.tight_layout()
    plt.show()

def main():
    
    sequence = read_fasta(r'D:\facultate\anul4\sem1\bioinformatics\lab\lab6\dna.fasta')
    
    print(f"Loaded sequence length: {len(sequence)} bp")
    
    samples = generate_random_samples(sequence)
    
    print("\nSample lengths (bp):")
    for i, (_, length) in enumerate(samples, 1):
        print(f"Sample {i}: {length} bp")
    
    draw_gel_electrophoresis(samples)

if __name__ == "__main__":
    main()