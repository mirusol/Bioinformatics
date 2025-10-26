genetic_code = {
    # U row
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    
    # C row
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    
    # A row
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    
    # G row
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def dna_to_rna(sequence):
    return sequence.replace('T', 'U')

def translate_sequence(sequence):
    sequence = sequence.upper().strip()
    
    if 'T' in sequence:
        sequence = dna_to_rna(sequence)
    
    amino_acids = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
    
            if amino_acid == 'Stop':
                break
            
            amino_acids.append(amino_acid)
        else:
            print(f"Warning: Unknown codon '{codon}' at position {i}")
    
    return amino_acids
def main():    
    print("Enter your own sequence (DNA or RNA):")
    user_sequence = input("Sequence: ")
    
    if user_sequence:
        result = translate_sequence(user_sequence)
        print(f"\nAmino Acid Sequence: {'-'.join(result)}")
        print(f"Full Result: {result}")

if __name__ == "__main__":
    main()