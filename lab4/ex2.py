from collections import Counter
import matplotlib.pyplot as plt

genetic_code = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def read_fasta(filename):
    with open(filename, 'r') as file:
        sequence = ''.join(file.readlines()[1:]).replace('\n', '')
    return sequence

def get_codons(sequence):
    sequence = sequence.upper()
    codons = []
    for i in range(0, len(sequence)-2):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    return codons

def plot_frequencies(counter, title):
    top_10 = counter.most_common(10)
    codons, counts = zip(*top_10)
    
    plt.figure(figsize=(10, 6))
    plt.bar(codons, counts)
    plt.title(f'Top 10 Most Frequent Codons in {title}')
    plt.xlabel('Codons')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def main():
    covid_seq = read_fasta('covid.fasta')
    influenza_seq = read_fasta('influenza.fasta')

    covid_codons = get_codons(covid_seq)
    influenza_codons = get_codons(influenza_seq)
    
    covid_counter = Counter(covid_codons)
    influenza_counter = Counter(influenza_codons)

    #a
    plot_frequencies(covid_counter, 'COVID-19')

    #b
    plot_frequencies(influenza_counter, 'Influenza')

    #c
    covid_top = set([codon for codon, _ in covid_counter.most_common(10)])
    influenza_top = set([codon for codon, _ in influenza_counter.most_common(10)])

    print("\nCommon frequent codons between viruses:")
    common = covid_top.intersection(influenza_top)
    print(f"Common codons: {common}")

    # d
    print("\nTop 3 codons in COVID-19 and their amino acids:")
    for codon, count in covid_counter.most_common(3):
        amino_acid = genetic_code.get(codon, 'Unknown')
        print(f"Codon {codon} ({amino_acid}): {count} times")

    print("\nTop 3 codons in Influenza and their amino acids:")
    for codon, count in influenza_counter.most_common(3):
        amino_acid = genetic_code.get(codon, 'Unknown')
        print(f"Codon {codon} ({amino_acid}): {count} times")

if __name__ == "__main__":
    main()