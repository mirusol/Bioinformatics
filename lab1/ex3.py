def parse_fasta(filename):
    """
    Parse a FASTA file and return the sequence(s).
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequences = []
    current_sequence = ""
    current_header = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_sequence:
                sequences.append({
                    'header': current_header,
                    'sequence': current_sequence
                })
            current_header = line
            current_sequence = ""
        else:
            current_sequence += line.replace(" ", "").upper()
    
   
    if current_sequence:
        sequences.append({
            'header': current_header,
            'sequence': current_sequence
        })
    
    return sequences


def analyze_sequence(s):
    """
    Analyze a sequence and return alphabet and percentages.
    Based on the original algorithm provided.
    """
   
    s = s.replace(" ", "")
    
    
    alphabet = []
    for c in s:
        if c not in alphabet:
            alphabet.append(c)
    
    print("Alphabet:", ''.join(alphabet))
    print(f"Sequence length: {len(s)}")
    print("-" * 40)
    
    
    results = []
    for letter in alphabet:
        count = s.count(letter)
        percent = count * 100 / len(s)
        results.append((letter, count, percent))
        print(f"{letter}: {percent:.2f}% (count: {count})")
    
    return alphabet, results


def main():
    """
    Main function to run the FASTA analyzer.
    """
    
    filename = input("Enter FASTA filename (or press Enter for 'sample.fasta'): ").strip()
    if not filename:
        filename = "sample.fasta"
    
    try:
        # Parse FASTA file
        sequences = parse_fasta(filename)
        
        print(f"\nFound {len(sequences)} sequence(s) in the file\n")
        print("=" * 60)
        
        # Analyze each sequence
        for idx, seq_data in enumerate(sequences, 1):
            print(f"\nSequence {idx}:")
            print(f"Header: {seq_data['header']}")
            print(f"\nFirst 60 characters: {seq_data['sequence'][:60]}...")
            print()
            
            analyze_sequence(seq_data['sequence'])
            print("=" * 60)
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()