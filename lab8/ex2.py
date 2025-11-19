from Bio import SeqIO
import re
import os

BACTERIAL_GENOMES = [
    ("Escherichia_coli.fasta", "Escherichia coli"),
    ("Bacillus_subtilis.fasta", "Bacillus subtilis"),
    ("Mycoplasma_genitalium.fasta", "Mycoplasma genitalium")
]

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def load_genome(filename):
    print(f"\nLoading genome from {filename}...")
    try:
        if not os.path.exists(filename):
            print(f"Error: File {filename} not found")
            return None
        
        record = SeqIO.read(filename, "fasta")
        print(f"Loaded: {record.description}")
        print(f"Sequence length: {len(record.seq)} bp")
        return str(record.seq).upper()
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def find_inverted_repeats(sequence, min_repeat_len=4, max_repeat_len=6, 
                         max_spacer_len=50, max_results=100):
    inverted_repeats = []
    seq_len = len(sequence)
    
    print(f"\nSearching for inverted repeats ({min_repeat_len}-{max_repeat_len} bp)...")
    
    positions_checked = 0
    total_positions = min(seq_len, 10000)
    
    for i in range(min(seq_len, 10000)):
        if i % 1000 == 0:
            print(f"Progress: {i}/{total_positions} positions checked...")
        
        for repeat_len in range(min_repeat_len, max_repeat_len + 1):
            if i + repeat_len >= seq_len:
                break
                
            left_repeat = sequence[i:i + repeat_len]
            
            if 'N' in left_repeat:
                continue
            
            right_repeat_target = reverse_complement(left_repeat)
            
            search_start = i + repeat_len
            search_end = min(i + repeat_len + max_spacer_len + repeat_len, seq_len)
            
            search_region = sequence[search_start:search_end]
            pos = search_region.find(right_repeat_target)
            
            if pos != -1:
                actual_pos = search_start + pos
                spacer_length = actual_pos - (i + repeat_len)
                transposon_start = i
                transposon_end = actual_pos + repeat_len - 1
                transposon_length = transposon_end - transposon_start + 1
                
                inverted_repeats.append({
                    'start': transposon_start,
                    'end': transposon_end,
                    'left_repeat': left_repeat,
                    'right_repeat': right_repeat_target,
                    'repeat_length': repeat_len,
                    'spacer_length': spacer_length,
                    'total_length': transposon_length
                })
                
                if len(inverted_repeats) >= max_results:
                    print(f"\nReached maximum of {max_results} results.")
                    return inverted_repeats
    
    return inverted_repeats

def filter_overlapping_repeats(repeats):
    if not repeats:
        return []
    
    repeats.sort(key=lambda x: (x['start'], -x['total_length']))
    
    filtered = []
    for repeat in repeats:
        is_duplicate = False
        for existing in filtered:
            if (repeat['start'] == existing['start'] and 
                repeat['end'] == existing['end']):
                is_duplicate = True
                break
        
        if not is_duplicate:
            filtered.append(repeat)
    
    return filtered

def print_inverted_repeats(repeats, genome_name, top_n=20):
    print(f"INVERTED REPEATS FOUND IN: {genome_name}")
    print(f"Total inverted repeats detected: {len(repeats)}")
    
    if not repeats:
        print("No inverted repeats found.")
        return
    
    repeats_by_length = {}
    for repeat in repeats:
        rl = repeat['repeat_length']
        if rl not in repeats_by_length:
            repeats_by_length[rl] = 0
        repeats_by_length[rl] += 1
    
    print(f"\nDistribution by repeat length:")
    for length in sorted(repeats_by_length.keys()):
        print(f"  {length} bp repeats: {repeats_by_length[length]}")
    
    print(f"\nShowing top {min(top_n, len(repeats))} inverted repeats:")
    print(f"{'-'*80}")
    
    for idx, repeat in enumerate(repeats[:top_n], 1):
        print(f"\n#{idx}")
        print(f"  Position: {repeat['start']} - {repeat['end']} ({repeat['total_length']} bp total)")
        print(f"  Left repeat:  {repeat['left_repeat']} ({repeat['repeat_length']} bp)")
        print(f"  Right repeat: {repeat['right_repeat']} ({repeat['repeat_length']} bp)")
        print(f"  Spacer length: {repeat['spacer_length']} bp")

def analyze_genome(filename, name):
    
    sequence = load_genome(filename)
    
    if sequence is None:
        print(f"Failed to load genome {filename}")
        return None
    
    repeats = find_inverted_repeats(
        sequence, 
        min_repeat_len=4, 
        max_repeat_len=6,
        max_spacer_len=50,
        max_results=100
    )
    
    filtered_repeats = filter_overlapping_repeats(repeats)
    
    print_inverted_repeats(filtered_repeats, name, top_n=20)
    
    return filtered_repeats

def main():

    print("BACTERIAL GENOME TRANSPOSON DETECTION")
    print("Searching for inverted repeats (4-6 bp) in bacterial genomes")
    all_results = {}
    
    for filename, name in BACTERIAL_GENOMES:
        try:
            results = analyze_genome(filename, name)
            all_results[name] = results
        except Exception as e:
            print(f"\nError analyzing {name}: {e}")
            continue
    
    for name, results in all_results.items():
        if results:
            print(f"{name}: {len(results)} potential transposons detected")
        else:
            print(f"{name}: Analysis failed or no results")

if __name__ == "__main__":
    main()
