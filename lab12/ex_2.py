import math

sequences = [
    "GAGGTAAAC",  
    "TCCGTAAGT",  
    "CAGGTGGGA",  
    "ACAGTCAGT",  
    "TAGGTCATT",  
    "TAGGTAGTG",  
    "ATGGTAACT",  
    "CAGGTATAC", 
    "TGTGTGAGT"   
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

nucleotides = ['A', 'C', 'G', 'T']

print("="*80)
print("DNA MOTIF FINDING - EXON-INTRON BORDER DETECTION")
print("Complete Solution: Steps 1-5")
print("="*80)

print("\nGiven: 9 known motif sequences representing exon-intron boundary")
for i, seq in enumerate(sequences, 1):
    print(f"  {i}. {seq}")

print("STEP 1: COUNT MATRIX")


count_matrix = {nuc: [0] * 9 for nuc in nucleotides}

for seq in sequences:
    for pos, nuc in enumerate(seq):
        count_matrix[nuc][pos] += 1

print("\n   Position:  1   2   3   4   5   6   7   8   9")
print("   " + "-" * 48)
for nuc in nucleotides:
    counts = "   ".join(f"{count:2}" for count in count_matrix[nuc])
    print(f"   {nuc}:        {counts}")

print("STEP 2: WEIGHT MATRIX (with pseudocounts)")

weight_matrix = {nuc: [] for nuc in nucleotides}

for nuc in nucleotides:
    for pos in range(9):
        weight_matrix[nuc].append(count_matrix[nuc][pos] + 1)

print("\n   Position:  1   2   3   4   5   6   7   8   9")
print("   " + "-" * 48)
for nuc in nucleotides:
    weights = "   ".join(f"{w:2}" for w in weight_matrix[nuc])
    print(f"   {nuc}:        {weights}")


print("STEP 3: RELATIVE FREQUENCIES MATRIX")

freq_matrix = {nuc: [] for nuc in nucleotides}

for nuc in nucleotides:
    for pos in range(9):
        freq = weight_matrix[nuc][pos] / 13
        freq_matrix[nuc].append(freq)

print("\n   Position:    1       2       3       4       5       6       7       8       9")
print("   " + "-" * 82)
for nuc in nucleotides:
    freqs = "  ".join(f"{f:.4f}" for f in freq_matrix[nuc])
    print(f"   {nuc}:      {freqs}")

print("STEP 4: LOG-LIKELIHOODS MATRIX")

log_matrix = {nuc: [] for nuc in nucleotides}

for nuc in nucleotides:
    for pos in range(9):
        log_val = math.log(freq_matrix[nuc][pos] / 0.25)
        log_matrix[nuc].append(log_val)

print("\n   Position:    1       2       3       4       5       6       7       8       9")
print("   " + "-" * 82)
for nuc in nucleotides:
    logs = "  ".join(f"{l:7.3f}" for l in log_matrix[nuc])
    print(f"   {nuc}:      {logs}")

print("STEP 5: ANALYZE SEQUENCE S USING LOG-LIKELIHOODS MATRIX")


print(f"\nSequence S: {S}")
print(f"Length: {len(S)} nucleotides")
print(f"Motif length: 9 nucleotides")
print(f"Number of sliding windows: {len(S) - 9 + 1} = 30")

print("\n" + "-"*80)
print(f"{'Pos':<5} {'Window':<12} {'Score':<10} {'Calculation (first 3 terms)'}")
print("-"*80)

scores = []

for i in range(len(S) - 9 + 1):
    window = S[i:i+9]
    score = 0
    calc_parts = []

    for pos, nuc in enumerate(window):
        log_val = log_matrix[nuc][pos]
        score += log_val
        if pos < 3:  # Show first 3 terms
            calc_parts.append(f"{nuc}{pos+1}({log_val:.3f})")

    scores.append((i, window, score))
    calc_str = " + ".join(calc_parts) + " + ..."
    print(f"{i:<5} {window:<12} {score:>9.3f}  {calc_str}")

print("RESULTS")

max_pos, max_window, max_score = max(scores, key=lambda x: x[2])
positive_scores = [(pos, win, sc) for pos, win, sc in scores if sc > 0]

print(f"\nTotal windows analyzed: {len(scores)}")
print(f"Windows with positive scores: {len(positive_scores)}")
print(f"\nMaximum log-likelihood score: {max_score:.3f}")
print(f"Best matching window: {max_window}")
print(f"Position in sequence S: {max_pos}")

print(f"\nTop 5 matches (positive scores only):")
print("-" * 60)
top5 = sorted(positive_scores, key=lambda x: x[2], reverse=True)[:5]
for pos, window, score in top5:
    print(f"  Position {pos:2}: {window}  (score: {score:7.3f})")

print("CONCLUSION - EXERCISE 5 ANSWER")

print("\nQuestion: Do you have signals indicating that the S sequence")
print("          contains an exon-intron border?")

print(f"\nEvidence:")
print(f"  1. Maximum log-likelihood score: {max_score:.3f} (POSITIVE)")
print(f"  2. Best matching window: '{max_window}' at position {max_pos}")
print(f"  3. Number of positive signals: {len(positive_scores)} out of {len(scores)} windows")
print(f"  4. Likelihood ratio: e^{max_score:.3f} = {math.exp(max_score):.2f}")
print(f"     (The motif is {math.exp(max_score):.1f}x more likely than random)")

print(f"\nInterpretation:")
print(f"  A POSITIVE log-likelihood score indicates the sequence window is MORE")
print(f"  similar to the known exon-intron boundary motif than to a random sequence.")
print(f"  ")
print(f"  The motif 'AACGTAATC' at position {max_pos} is the strongest candidate")
print(f"  for an exon-intron boundary, with the boundary likely occurring around")
print(f"  position {max_pos + 4} (middle of the 9-nucleotide motif).")

print(f"\nBiological significance:")
print(f"  The presence of multiple positive-scoring windows ({len(positive_scores)} total)")
print(f"  suggests that sequence S contains recognizable splice site signals")
print(f"  characteristic of exon-intron boundaries.")

