def find_dinucleotides_trinucleotides(sequence):

    dinucleotides = set()
    trinucleotides = set()
    
    for i in range(len(sequence) - 1):
        dinucleotides.add(sequence[i:i+2])
    
    
    for i in range(len(sequence) - 2):
        trinucleotides.add(sequence[i:i+3])
    
    return dinucleotides, trinucleotides


S = "ABAA"
di, tri = find_dinucleotides_trinucleotides(S)

print(f"Dinucleotides: {sorted(di)}")
print(f"Trinucleotides: {sorted(tri)}")
