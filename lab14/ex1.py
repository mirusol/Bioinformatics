import math

# Known sequences
S1 = "ATCGATTCGATATCATACACGTAT"  # CpG+ island
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"  # CpG- non-island
S_test = "CAGGTTGGAAACGTAA"  # Test sequence

# Nucleotides
nucleotides = ['A', 'C', 'G', 'T']

def count_transitions(sequence):
    """Count transition frequencies from a sequence"""
    # Initialize transition count matrix
    transitions = {n1: {n2: 0 for n2 in nucleotides} for n1 in nucleotides}
    
    # Count transitions
    for i in range(len(sequence) - 1):
        current = sequence[i]
        next_base = sequence[i + 1]
        transitions[current][next_base] += 1
    
    return transitions

def calculate_probabilities(transitions):
    """Calculate transition probabilities from counts"""
    probabilities = {n1: {n2: 0.0 for n2 in nucleotides} for n1 in nucleotides}
    
    for n1 in nucleotides:
        total = sum(transitions[n1].values())
        if total > 0:
            for n2 in nucleotides:
                probabilities[n1][n2] = transitions[n1][n2] / total
        else:
            # If no transitions from this nucleotide, uniform distribution
            for n2 in nucleotides:
                probabilities[n1][n2] = 0.25
    
    return probabilities

def print_matrix(matrix, title):
    """Pretty print a transition matrix"""
    print(f"\n{title}")
    print(f"{'':>8}", end='')
    for n in nucleotides:
        print(f"{n:>10}", end='')
    print()
    
    for n1 in nucleotides:
        print(f"{n1:>8}", end='')
        for n2 in nucleotides:
            print(f"{matrix[n1][n2]:>10.3f}", end='')
        print()

def calculate_log_likelihood_matrix(prob_plus, prob_minus):
    """Calculate log-likelihood matrix: β = log2(P+/P-)"""
    beta = {n1: {n2: 0.0 for n2 in nucleotides} for n1 in nucleotides}
    
    for n1 in nucleotides:
        for n2 in nucleotides:
            p_plus = prob_plus[n1][n2]
            p_minus = prob_minus[n1][n2]
            
            # Avoid division by zero and log(0)
            if p_plus > 0 and p_minus > 0:
                ratio = p_plus / p_minus
                # log2(x) = ln(x) / ln(2)
                beta[n1][n2] = math.log(ratio) / math.log(2)
            elif p_plus > 0 and p_minus == 0:
                beta[n1][n2] = float('inf')  # Strongly favors CpG+
            elif p_plus == 0 and p_minus > 0:
                beta[n1][n2] = float('-inf')  # Strongly favors CpG-
            else:
                beta[n1][n2] = 0.0  # Both zero, neutral
    
    return beta

def calculate_log_likelihood_score(sequence, beta_matrix):
    """Calculate the log-likelihood score for a sequence"""
    score = 0.0
    
    for i in range(len(sequence) - 1):
        current = sequence[i]
        next_base = sequence[i + 1]
        
        if current in nucleotides and next_base in nucleotides:
            beta_value = beta_matrix[current][next_base]
            if not math.isinf(beta_value):
                score += beta_value
                print(f"  {current} -> {next_base}: β = {beta_value:.4f}, cumulative score = {score:.4f}")
            else:
                print(f"  {current} -> {next_base}: β = {beta_value}, cumulative score = {score:.4f}")
    
    return score

# Step 1: Count transitions for CpG+ model (S1)
print("="*70)
print("STEP 1: Count transitions for CpG+ model (S1)")
print("="*70)
print(f"S1 (CpG+ island): {S1}")
transitions_plus = count_transitions(S1)
print("\nTransition counts for S1:")
print_matrix(transitions_plus, "CpG+ Transition Counts")

# Calculate probabilities
prob_plus = calculate_probabilities(transitions_plus)
print_matrix(prob_plus, "CpG+ Transition Probabilities")

# Step 2: Count transitions for CpG- model (S2)
print("\n" + "="*70)
print("STEP 2: Count transitions for CpG- model (S2)")
print("="*70)
print(f"S2 (CpG- non-island): {S2}")
transitions_minus = count_transitions(S2)
print("\nTransition counts for S2:")
print_matrix(transitions_minus, "CpG- Transition Counts")

# Calculate probabilities
prob_minus = calculate_probabilities(transitions_minus)
print_matrix(prob_minus, "CpG- Transition Probabilities")

# Step 3: Calculate log-likelihood matrix
print("\n" + "="*70)
print("STEP 3: Calculate log-likelihood matrix β = log2(P+/P-)")
print("="*70)
beta_matrix = calculate_log_likelihood_matrix(prob_plus, prob_minus)
print_matrix(beta_matrix, "Log-Likelihood Matrix (β)")

# Step 4: Test the new sequence
print("\n" + "="*70)
print("STEP 4: Test new sequence")
print("="*70)
print(f"Test sequence S: {S_test}")
print("\nCalculating log-likelihood score:")
score = calculate_log_likelihood_score(S_test, beta_matrix)

print("\n" + "="*70)
print("RESULTS")
print("="*70)
print(f"Test sequence: {S_test}")
print(f"Total log-likelihood score: {score:.4f}")
print("\nInterpretation:")
if score > 0:
    print(f"  Score > 0: The sequence LIKELY belongs to a CpG island (+)")
    print(f"  The positive score of {score:.4f} indicates the sequence is more")
    print(f"  consistent with the CpG+ model than the CpG- model.")
elif score < 0:
    print(f"  Score < 0: The sequence LIKELY DOES NOT belong to a CpG island (-)")
    print(f"  The negative score of {score:.4f} indicates the sequence is more")
    print(f"  consistent with the CpG- model than the CpG+ model.")
else:
    print(f"  Score = 0: The sequence is NEUTRAL")
    print(f"  The sequence is equally consistent with both models.")

print("\n" + "="*70)
