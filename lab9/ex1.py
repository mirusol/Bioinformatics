from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


class RestrictionEnzyme:
    
    def __init__(self, name: str, recognition_seq: str, cut_position: int):
        self.name = name
        self.recognition_seq = recognition_seq.upper()
        self.cut_position = cut_position
        self.length = len(recognition_seq)
    
    def find_sites(self, dna_sequence: str) -> List[int]:
        dna_sequence = dna_sequence.upper()
        sites = []
        
        for i in range(len(dna_sequence) - self.length + 1):
            if dna_sequence[i:i + self.length] == self.recognition_seq:
                sites.append(i)
        
        return sites
    
    def digest(self, dna_sequence: str) -> Tuple[int, List[int], List[int]]:
        sites = self.find_sites(dna_sequence)
        num_cuts = len(sites)
        
        if num_cuts == 0:
            return 0, [], [len(dna_sequence)]
        cut_positions = [site + self.cut_position for site in sites]
        
        fragment_lengths = []
        prev_pos = 0
        
        for cut_pos in cut_positions:
            fragment_lengths.append(cut_pos - prev_pos)
            prev_pos = cut_pos
        
        fragment_lengths.append(len(dna_sequence) - prev_pos)
        
        return num_cuts, cut_positions, fragment_lengths


def read_fasta(filename: str) -> Tuple[str, str]:
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    header = ""
    sequence = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:]
        else:
            sequence += line.upper()
    
    return header, sequence


def create_gel_simulation(enzyme_results: Dict[str, Tuple[int, List[int], List[int]]], 
                          dna_length: int,
                          output_file: str = "gel.png"):
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    num_lanes = len(enzyme_results) + 1
    lane_width = 0.8
    lane_spacing = 1.0
    gel_height = 10
    
    gel_rect = patches.Rectangle((0, 0), num_lanes * lane_spacing, gel_height, 
                                 linewidth=2, edgecolor='black', 
                                 facecolor='lightgray', alpha=0.3)
    ax.add_patch(gel_rect)

    well_y = gel_height - 0.5
    for i in range(num_lanes):
        well_x = i * lane_spacing + (lane_spacing - lane_width) / 2
        well = patches.Rectangle((well_x, well_y), lane_width, 0.3,
                                linewidth=1, edgecolor='black', 
                                facecolor='white')
        ax.add_patch(well)
    
    marker_lane = num_lanes - 1
    marker_x = marker_lane * lane_spacing + lane_spacing / 2
    
    markers = [3000, 2000, 1500, 1000, 750, 500, 250, 100]
    markers = [m for m in markers if m <= dna_length * 1.2]
    
    for marker_size in markers:
        y_pos = gel_height - 1.5 - (np.log10(marker_size) / np.log10(max(markers))) * 7
        
        band = patches.Rectangle((marker_x - lane_width/2, y_pos - 0.05), 
                                lane_width, 0.1,
                                linewidth=0, facecolor='blue', alpha=0.7)
        ax.add_patch(band)
        
        ax.text(marker_x + lane_width/2 + 0.1, y_pos, f"{marker_size} bp", 
               fontsize=8, va='center')
    for idx, (enzyme_name, (num_cuts, cut_positions, fragment_lengths)) in enumerate(enzyme_results.items()):
        lane_x = idx * lane_spacing + lane_spacing / 2
        
        ax.text(lane_x, gel_height - 0.3, enzyme_name, 
               ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        for fragment_length in fragment_lengths:
            if fragment_length > 0:
                max_size = max(markers) if markers else dna_length
                y_pos = gel_height - 1.5 - (np.log10(max(fragment_length, 10)) / np.log10(max_size)) * 7
                intensity = min(0.9, 0.3 + (fragment_length / dna_length) * 0.6)
                
                band = patches.Rectangle((lane_x - lane_width/2, y_pos - 0.05), 
                                        lane_width, 0.1,
                                        linewidth=0, facecolor='red', alpha=intensity)
                ax.add_patch(band)
    
    ax.text(marker_x, gel_height - 0.3, "Size\nMarker", 
           ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_xlim(-0.5, num_lanes * lane_spacing + 0.5)
    ax.set_ylim(0, gel_height + 0.5)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.annotate('', xy=(0.2, 0.5), xytext=(0.2, gel_height - 1),
               arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.text(0.1, gel_height/2, 'Migration', rotation=90, 
           va='center', ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved gel image: {output_file}")
    plt.show()


def analyze_dna_digestion(dna_sequence: str, enzymes: List[RestrictionEnzyme]):
   
    print("DNA Analysis")
    print(f"\nDNA Sequence Length: {len(dna_sequence)} bp")
    print(f"Number of Enzymes: {len(enzymes)}")
    print()
    
    enzyme_results = {}
    
    for enzyme in enzymes:
        print(f"\n{'-'*80}")
        print(f"Enzyme: {enzyme.name}")
        print(f"Recognition Sequence: 5'-{enzyme.recognition_seq}-3'")
        print(f"Cut Position: After base {enzyme.cut_position}")
        
        num_cuts, cut_positions, fragment_lengths = enzyme.digest(dna_sequence)
        enzyme_results[enzyme.name] = (num_cuts, cut_positions, fragment_lengths)
        
        print(f"\nNumber of Cleavages: {num_cuts}")
        
        if num_cuts > 0:
            print(f"Cleavage Positions: {', '.join(map(str, cut_positions))}")
            print(f"\nFragment Lengths ({len(fragment_lengths)} fragments):")
            for i, length in enumerate(sorted(fragment_lengths, reverse=True), 1):
                print(f"  Fragment {i}: {length} bp")
        else:
            print("No cleavage sites found.")
            print(f"Fragment Length: {fragment_lengths[0]} bp (uncut)")
    
    print(f"\n{'='*80}\n")
    
    return enzyme_results


def main():
    enzymes = [
        RestrictionEnzyme("EcoRI", "GAATTC", 1),      
        RestrictionEnzyme("BamHI", "GGATCC", 1),     
        RestrictionEnzyme("HindIII", "AAGCTT", 1), 
        RestrictionEnzyme("TaqI", "TCGA", 1),         
        RestrictionEnzyme("HaeIII", "GGCC", 2),       
    ]
    
    fasta_file = "dna.fasta"
    
    try:
        header, dna_sequence = read_fasta(fasta_file)
        print(f"Loaded sequence: {header}")
        print(f"Sequence length: {len(dna_sequence)} bp\n")
        
        if len(dna_sequence) < 1000 or len(dna_sequence) > 3000:
            print(f"Warning: Sequence length ({len(dna_sequence)} bp) is outside recommended range (1000-3000 bp)")
        
        enzyme_results = analyze_dna_digestion(dna_sequence, enzymes)
        
        create_gel_simulation(enzyme_results, len(dna_sequence))
        
        total_cuts = sum(result[0] for result in enzyme_results.values())
        print(f"Total cleavages across all enzymes: {total_cuts}")
        
        for enzyme_name, (num_cuts, _, fragment_lengths) in enzyme_results.items():
            if num_cuts > 0:
                avg_fragment = sum(fragment_lengths) / len(fragment_lengths)
                print(f"{enzyme_name}: {num_cuts} cuts, {len(fragment_lengths)} fragments, "
                      f"avg fragment size: {avg_fragment:.1f} bp")
            else:
                print(f"{enzyme_name}: No cuts (DNA remains intact)")
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{fasta_file}'")
        return


if __name__ == "__main__":
    main()
