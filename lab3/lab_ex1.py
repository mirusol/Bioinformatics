import math

def calculate_tm_basic(dna_sequence):
    
    
    sequence = dna_sequence.upper()
    
    
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    
    
    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    
    return tm
def calculate_tm_advanced(dna_sequence, na_concentration=0.05):
    
    
    sequence = dna_sequence.upper()
    length = len(sequence)
    
    gc_count = sequence.count('G') + sequence.count('C')
    gc_percentage = (gc_count / length) * 100
    

    tm = 81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc_percentage - (600 / length)
    
    return tm

def main():
    
    dna_sequence = input("DNA seq: ").strip()
    
    valid_nucleotides = set('ATCG')
        
    basic_tm = calculate_tm_basic(dna_sequence)
    advanced_tm = calculate_tm_advanced(dna_sequence)
    

    print(f"Basic Tm: {basic_tm:.2f}°C")
    print(f"Advanced Tm: {advanced_tm:.2f}°C")

if __name__ == "__main__":
    main()