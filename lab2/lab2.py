import tkinter as tk

def generate_combinations(length):
    bases = ['A', 'C', 'G', 'T']
    combinations = []
    
    def recursive_combine(current, depth):
        if depth == length:
            combinations.append(current)
            return
        for base in bases:
            recursive_combine(current + base, depth + 1)
    
    recursive_combine('', 0)
    return combinations

def calculate_percentage(sequence, pattern):
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            count += 1
    
    total_possible_positions = len(sequence) - len(pattern) + 1
    percentage = (count / total_possible_positions) * 100
    return percentage

def main():
    
    S = "ATTGTCCCAATCTGTTG"
    
    print("Dinucleotide Percentages:")
    dinucleotides = generate_combinations(2)
    for di in dinucleotides:
        percentage = calculate_percentage(S, di)
        print(f"{di}: {percentage:.2f}%")
    
    print("\nTrinucleotide Percentages:")
    trinucleotides = generate_combinations(3)
    for tri in trinucleotides:
        percentage = calculate_percentage(S, tri)
        print(f"{tri}: {percentage:.2f}%")

if __name__ == "__main__":
    main()