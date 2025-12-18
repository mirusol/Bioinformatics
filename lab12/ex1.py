import math

bases = ['A', 'C', 'G', 'T']
base_to_row = {b: i for i, b in enumerate(bases)}

motif_len = 9
motif_sequences = None
pseudocount = 1.0

background = {b: 0.25 for b in bases}

P_given = [
    [0.3, 0.6, 0.1, 0.0, 0.0, 0.6, 0.7, 0.2, 0.1],
    [0.2, 0.2, 0.1, 0.0, 0.0, 0.0, 0.2, 0.1, 0.2],
    [0.1, 0.0, 0.7, 0.0, 0.0, 0.1, 0.1, 0.5, 0.1],
    [0.4, 0.1, 0.1, 1.0, 1.0, 0.1, 0.0, 0.2, 0.6],
]

def counts_from_sequences(seqs, L):
    counts = [[0] * L for _ in bases]
    for s in seqs:
        assert len(s) == L, "All sequences must have equal length L"
        for i, ch in enumerate(s):
            counts[base_to_row[ch]][i] += 1
    return counts

def probs_from_counts(counts, alpha=1.0):
    L = len(counts[0])
    probs = [[0.0] * L for _ in bases]
    for j in range(L):
        col_total = sum(counts[i][j] + alpha for i in range(4))
        for i in range(4):
            probs[i][j] = (counts[i][j] + alpha) / col_total
    return probs

def log_likelihood_matrix_from_pwm(pwm, bg):
    eps = 1e-12
    ll = []
    for i in range(4):
        row = []
        for j in range(len(pwm[0])):
            p = max(pwm[i][j], eps)
            row.append(math.log(p / bg[bases[i]]))
        ll.append(row)
    return ll

def print_matrix(title, matrix, rounding=3):
    print(title)
    for b, row in zip(bases, matrix):
        if isinstance(row[0], float):
            print(b, [round(x, rounding) for x in row])
        else:
            print(b, row)

if motif_sequences:
    num_seqs = len(motif_sequences)
    motif_len = len(motif_sequences[0])
    count_matrix = counts_from_sequences(motif_sequences, motif_len)
    print("1) COUNT MATRIX (from provided sequences)")
    print_matrix("", count_matrix)

    weight_matrix = probs_from_counts(count_matrix, alpha=pseudocount)
    print("\n2) WEIGHT MATRIX (probabilities with pseudocounts)")
    print_matrix("", weight_matrix)
else:
    num_seqs = 10
    count_matrix = [[int(round(p * num_seqs)) for p in row] for row in P_given]
    print("1) COUNT MATRIX (derived by rounding from given PWM)")
    print_matrix("", count_matrix)

    weight_matrix = P_given
    print("\n2) WEIGHT MATRIX (given PWM)")
    print_matrix("", weight_matrix)

rel_freq_matrix = weight_matrix
print("\n3) RELATIVE FREQUENCIES (pwm)")
print_matrix("", rel_freq_matrix)

ll_matrix = log_likelihood_matrix_from_pwm(rel_freq_matrix, background)
print("\n4) LOG-LIKELIHOOD MATRIX ln(P / bg)")
print_matrix("", ll_matrix)

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

def score_window(window: str) -> float:
    score = 0.0
    for i, base in enumerate(window):
        score += ll_matrix[base_to_row[base]][i]
    return score

def scan_sequence(seq: str, k: int):
    scores = []
    for start in range(0, len(seq) - k + 1):
        w = seq[start:start + k]
        scores.append((start, w, score_window(w)))
    return scores

print("\n5) SLIDING-WINDOW SCORES ON S")
window_scores = scan_sequence(S, motif_len)
for start, w, s in window_scores:
    print(f"pos {start:2d}-{start+motif_len-1:2d}  {w}  score = {s:.3f}")

score_threshold = 0.0
top_n = 10

sorted_hits = sorted([x for x in window_scores if x[2] >= score_threshold],
                     key=lambda x: x[2], reverse=True)

if sorted_hits:
    print("\nHITS (score >= threshold, sorted):")
    for i, (start, w, s) in enumerate(sorted_hits[:top_n], 1):
        print(f"{i:2d}. pos {start}-{start+motif_len-1}  {w}  score = {s:.3f}")
else:
    print("\nNo windows passed the threshold.")

best_start, best_window, best_score = max(window_scores, key=lambda x: x[2])
print("\nBEST WINDOW:")
print(f"pos {best_start}-{best_start+motif_len-1}, {best_window}, score = {best_score:.3f}")

if best_score > 0:
    print("\nConclusion: YES, there is a strong signal that S contains an exon–intron border.")
else:
    print("\nConclusion: NO clear exon–intron border signal above random background.")