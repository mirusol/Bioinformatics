import random
import re

def generate_random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def create_transposon(internal_length=20, repeat_length=8):
    left_repeat = generate_random_dna(repeat_length)
    
    internal_seq = generate_random_dna(internal_length)
    
    right_repeat = reverse_complement(left_repeat)
    
    has_mutation = False
    if random.random() < 0.3:
        has_mutation = True
        right_repeat_list = list(right_repeat)
        num_mutations = random.randint(1, 3)
        for _ in range(num_mutations):
            pos = random.randint(0, len(right_repeat_list) - 1)
            right_repeat_list[pos] = random.choice('ACGT')
        right_repeat = ''.join(right_repeat_list)
    
    transposon = left_repeat + internal_seq + right_repeat
    
    return transposon, left_repeat, right_repeat, has_mutation

def create_dna_with_transposons(total_length=300, num_transposons=3, allow_overlap=True):
    dna_sequence = generate_random_dna(total_length)
    transposons_info = []
    
    for i in range(num_transposons):
        transposon, left_repeat, right_repeat, has_mutation = create_transposon(
            internal_length=random.randint(15, 25),
            repeat_length=random.randint(7, 10)
        )
        
        if allow_overlap and i > 0 and random.random() < 0.4:
            existing = random.choice(transposons_info)
            
            if random.random() < 0.5:
                available_space = (existing['end'] - existing['start'] + 1) - len(transposon)
                if available_space > 10:
                    max_start = existing['end'] - len(transposon)
                    start_pos = random.randint(existing['start'] + 5, max_start)
                    overlap_type = 'nested'
                else:
                    start_pos = random.randint(0, max(0, total_length - len(transposon)))
                    overlap_type = 'none'
            else:
                overlap_amount = random.randint(5, min(15, (existing['end'] - existing['start']) // 2))
                start_pos = max(0, existing['end'] - overlap_amount)
                overlap_type = 'overlapping'
        else:
            start_pos = random.randint(0, max(0, total_length - len(transposon)))
            overlap_type = 'none'
        
        end_pos = start_pos + len(transposon) - 1
        
        if end_pos < total_length:
            dna_sequence = dna_sequence[:start_pos] + transposon + dna_sequence[end_pos + 1:]
            
            transposons_info.append({
                'id': i + 1,
                'start': start_pos,
                'end': end_pos,
                'sequence': transposon,
                'left_repeat': left_repeat,
                'right_repeat': right_repeat,
                'length': len(transposon),
                'has_mutation': has_mutation,
                'overlap_type': overlap_type
            })
    
    transposons_info.sort(key=lambda x: x['start'])
    for idx, trans in enumerate(transposons_info):
        trans['id'] = idx + 1
    
    return dna_sequence, transposons_info

def detect_inverted_repeats(dna_sequence, min_repeat_length=6, max_internal_gap=50):
    detected_transposons = []
    seq_length = len(dna_sequence)
    
    for i in range(seq_length):
        for repeat_len in range(min_repeat_length, min(15, seq_length - i)):
            left_repeat = dna_sequence[i:i + repeat_len]
            right_repeat_target = reverse_complement(left_repeat)
            
            search_start = i + repeat_len
            search_end = min(i + repeat_len + max_internal_gap + repeat_len, seq_length)
            
            pos = dna_sequence.find(right_repeat_target, search_start, search_end)
            
            if pos != -1:
                start = i
                end = pos + repeat_len - 1
                transposon_seq = dna_sequence[start:end + 1]
                internal_length = pos - (i + repeat_len)
                
                detected_transposons.append({
                    'start': start,
                    'end': end,
                    'length': end - start + 1,
                    'left_repeat': left_repeat,
                    'right_repeat': right_repeat_target,
                    'repeat_length': repeat_len,
                    'internal_length': internal_length,
                    'sequence': transposon_seq
                })
    
    detected_transposons.sort(key=lambda x: (x['start'], -x['length']))
    
    filtered = []
    for trans in detected_transposons:
        is_duplicate = False
        for existing in filtered:
            if (trans['start'] == existing['start'] and 
                trans['end'] == existing['end'] and
                trans['repeat_length'] == existing['repeat_length']):
                is_duplicate = True
                break
        
        if not is_duplicate:
            filtered.append(trans)
    
    return filtered

def print_transposon_info(transposons, title="Transposons"):
    
    
    for trans in transposons:
        print(f"\nTransposon {trans.get('id', '?')}:")
        print(f"Position: {trans['start']} - {trans['end']}")
        print(f"Length: {trans['length']} bp")
        print(f"Left Repeat:  {trans['left_repeat']}")
        print(f"Right Repeat: {trans['right_repeat']}")
        if trans.get('has_mutation'):
            print(f"Contains mutations")
        if trans.get('overlap_type') and trans['overlap_type'] != 'none':
            print(f"Relationship: {trans['overlap_type'].upper()}")
        if 'internal_length' in trans:
            print(f"Internal Length: {trans['internal_length']} bp")
        print(f"Sequence: {trans['sequence'][:50]}{'...' if len(trans['sequence']) > 50 else ''}")

def main():
    
    dna_length = random.randint(200, 400)
    num_transposons = random.randint(3, 4)

    print(f"Total length: {dna_length} bp")
    print(f"Number of transposons: {num_transposons}")
    
    dna_sequence, original_transposons = create_dna_with_transposons(dna_length, num_transposons)
    
    print(f"\nGenerated DNA sequence ({len(dna_sequence)} bp):")
    print(f"{dna_sequence[:100]}")
    print(f"{dna_sequence[-100:]}")
    
    print_transposon_info(original_transposons, "ORIGINAL TRANSPOSONS")
    
    detected = detect_inverted_repeats(dna_sequence, min_repeat_length=8, max_internal_gap=50)
    
    if detected:
        print_transposon_info(detected, "DETECTED TRANSPOSONS")
        
        print(f"Original transposons: {len(original_transposons)}")
        print(f"Detected transposons: {len(detected)}")
        
        has_overlap = False
        for i, trans1 in enumerate(original_transposons):
            for trans2 in original_transposons[i+1:]:
                if trans1['start'] <= trans2['start'] <= trans1['end']:
                    has_overlap = True
                    if trans2['end'] <= trans1['end']:
                        print(f"  🔹 Transposon {trans2['id']} is NESTED within Transposon {trans1['id']}")
                    else:
                        overlap = trans1['end'] - trans2['start'] + 1
                        print(f"Transposon {trans1['id']} and {trans2['id']} OVERLAP by {overlap} bp")
        if not has_overlap:
            print("No overlapping or nested transposons detected.")
       
        for orig in original_transposons:
            matched = False
            for det in detected:
                if det['start'] == orig['start'] and det['end'] == orig['end']:
                    print(f"Transposon {orig['id']} at position {orig['start']}-{orig['end']} - EXACTLY DETECTED")
                    matched = True
                    break
            
            if not matched:
                for det in detected:
                    if not (det['end'] < orig['start'] or det['start'] > orig['end']):
                        overlap_start = max(det['start'], orig['start'])
                        overlap_end = min(det['end'], orig['end'])
                        overlap_percent = ((overlap_end - overlap_start + 1) / orig['length']) * 100
                        print(f"Transposon {orig['id']} at position {orig['start']}-{orig['end']} - PARTIALLY DETECTED ({det['start']}-{det['end']}, {overlap_percent:.1f}% overlap)")
                        matched = True
                        break
                
                if not matched:
                    print(f"Transposon {orig['id']} at position {orig['start']}-{orig['end']} - NOT DETECTED")
    else:
        print("\nNo transposons detected")
    

if __name__ == "__main__":
    random.seed() 
    main()
