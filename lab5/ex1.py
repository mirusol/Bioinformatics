import random
from collections import defaultdict

def load_sequence():
    sequence = """
AGCGAAAGCAGGTCAATTATATTCAATATGGAAAGAATAAAAGAACTAAGAAATCTAATGTCGCAGTCTC
GCACCCGCGAGATACTCACAAAAACCACCGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAG
ACAGGAGAAGAACCCAGCACTTAGGATGAAATGGATGATGGCAATGAAATATCCAATTACAGCAGACAAG
AGGATAACGGAAATGATTCCTGAGAGAAATGAGCAAGGACAAACTTTATGGAGTAAAATGAATGATGCCG
GATCAGACCGAGTGATGGTATCACCTCTGGCTGTGACATGGTGGAATAGGAATGGACCAATGACAAATAC
AGTTCATTATCCAAAAATCTACAAAACTTATTTTGAAAGAGTCGAAAGGCTAAAGCATGGAACCTTTGGC
CCTGTCCATTTTAGAAACCAAGTCAAAATACGTCGGAGAGTTGACATAAATCCTGGTCATGCAGATCTCA
GTGCCAAGGAGGCACAGGATGTAATCATGGAAGTTGTTTTCCCTAACGAAGTGGGAGCCAGGATACTAAC
ATCGGAATCGCAACTAACGATAACCAAAGAGAAGAAAGAAGAACTCCAGGATTGCAAAATTTCTCCTTTG
ATGGTTGCATACATGTTGGAGAGAGAACTGGTCCGCAAAACGAGATTCCTCCCAGTGGCTGGTGGAACAA
GCAGTGTGTACATTGAAGTGTTGCATTTGACTCAAGGAACATGCTGGGAACAGATGTATACTCCAGGAGG
GGAAGTGAAGAATGATGATGTTGATCAAAGCTTGATTATTGCTGCTAGGAACATAGTGAGAAGAGCTGCA
GTATCAGCAGACCCACTAGCATCTTTATTGGAGATGTGCCACAGCACACAGATTGGTGGAATTAGGATGG
TAGACATCCTTAAGCAGAACCCAACAGAAGAGCAAGCCGTGGGTATATGCAAGGCTGCAATGGGACTGAG
AATTAGCTCATCCTTCAGTTTTGGTGGATTCACATTTAAGAGAACAAGCGGATCATCAGTCAAGAGAGAG
GAAGAGGTGCTTACGGGCAATCTTCAAACATTGAAGATAAGAGTGCATGAGGGATATGAAGAGTTCACAA

    """
    return sequence.replace("\n", "").replace(" ", "")

def generate_samples(sequence, num_samples=2000, sample_size=100):
    samples = []
    seq_len = len(sequence)
    for _ in range(num_samples):
        start = random.randint(0, seq_len - sample_size)
        sample = sequence[start:start + sample_size]
        samples.append((start, sample))
    return samples

def reconstruct_sequence(samples, original_length):
    consensus = ['N'] * original_length
    coverage = defaultdict(int)
    
    for start, sample in samples:
        for i, base in enumerate(sample):
            pos = start + i
            if pos < original_length:
                if consensus[pos] == 'N':
                    consensus[pos] = base
                elif consensus[pos] != base:
                    coverage[pos] += 1
    
    return ''.join(consensus)

def main():
    original_sequence = load_sequence()
    print(f"Original sequence length: {len(original_sequence)}")
    
    samples = generate_samples(original_sequence)
    print(f"Generated {len(samples)} samples")
    
    reconstructed_sequence = reconstruct_sequence(samples, len(original_sequence))
    print(f"Reconstructed sequence length: {len(reconstructed_sequence)}")
    
    matches = sum(1 for a, b in zip(original_sequence, reconstructed_sequence) 
                 if a == b)
    accuracy = (matches / len(original_sequence)) * 100
    print(f"Reconstruction accuracy: {accuracy:.2f}%")
    

if __name__ == "__main__":
    main()