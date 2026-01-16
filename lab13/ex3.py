import numpy as np
import json
import re
from collections import defaultdict

def get_sample_text():
    text = """
		Autumn sunlight filters through the quiet library windows while students whisper about 
		upcoming exams. In the corner, a curious fox plush sits beside a stack of notebooks,
		reminding everyone to pause, breathe, and trust the steady pace of learning each day.Calm review turns work into steady recall.
    """
    return text.strip()


def preprocess_text(text):
    words = re.findall(r'\b[a-z]+\b', text.lower())
    return words


def create_word_to_symbol_mapping(words):
    unique_words = sorted(set(words)) 
    
    word_to_symbol = {}
    symbol_to_word = {}
    
    for idx, word in enumerate(unique_words):
        symbol = f"S{idx}" 
        word_to_symbol[word] = symbol
        symbol_to_word[symbol] = word
    
    return word_to_symbol, symbol_to_word


def calculate_word_transition_matrix(words, word_to_symbol):
    symbols = sorted(word_to_symbol.values())
    symbol_to_idx = {symbol: idx for idx, symbol in enumerate(symbols)}
    n = len(symbols)
    
    transition_counts = defaultdict(lambda: defaultdict(int))
    
    for i in range(len(words) - 1):
        current_word = words[i]
        next_word = words[i + 1]
        
        current_symbol = word_to_symbol[current_word]
        next_symbol = word_to_symbol[next_word]
        
        transition_counts[current_symbol][next_symbol] += 1
    
    transition_matrix = defaultdict(lambda: defaultdict(float))
    
    for from_symbol in symbols:
        total = sum(transition_counts[from_symbol].values())
        
        if total > 0:
            for to_symbol in symbols:
                count = transition_counts[from_symbol].get(to_symbol, 0)
                transition_matrix[from_symbol][to_symbol] = count / total
        else:
            for to_symbol in symbols:
                transition_matrix[from_symbol][to_symbol] = 1.0 / n
    
    return transition_matrix, transition_counts, symbols


def save_to_json(text, words, word_to_symbol, symbol_to_word, transition_matrix, 
                 transition_counts, filename='word_transition_matrix.json'):
    data = {
        'original_text': text,
        'text_length': len(text),
        'word_count': len(words),
        'unique_words': len(word_to_symbol),
        'word_to_symbol_mapping': word_to_symbol,
        'symbol_to_word_mapping': symbol_to_word,
        'transition_matrix': {k: dict(v) for k, v in transition_matrix.items()},
        'transition_counts': {k: dict(v) for k, v in transition_counts.items()},
        'word_frequencies': {}
    }
    
    for word in word_to_symbol.keys():
        count = words.count(word)
        data['word_frequencies'][word] = {
            'symbol': word_to_symbol[word],
            'count': count,
            'frequency': count / len(words)
        }
    
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    
    print(f"Word transition matrix saved to: {filename}")


def display_results(text, words, word_to_symbol, symbol_to_word, 
                   transition_matrix, transition_counts):
    
    print("WORD TRANSITION MATRIX ANALYSIS")
    print(f"\nOriginal Text ({len(text)} characters):")
    print(text[:200] + "..." if len(text) > 200 else text)
   
    print(f"TEXT STATISTICS")
    print(f"Total characters: {len(text)}")
    print(f"Total words: {len(words)}")
    print(f"Unique words: {len(word_to_symbol)}")
    
    print("WORD TO SYMBOL MAPPING")
    print(f"\n{'Word':<15} {'Symbol':<10} {'Frequency':<10}")
    

    word_freq = [(word, words.count(word)) for word in word_to_symbol.keys()]
    word_freq.sort(key=lambda x: x[1], reverse=True)
    
    for word, freq in word_freq[:20]:  # Show top 20
        symbol = word_to_symbol[word]
        print(f"{word:<15} {symbol:<10} {freq:<10}")
    
    if len(word_freq) > 20:
        print(f"... and {len(word_freq) - 20} more words")
    
    print(f"\n" + "=" * 80)
    print("SAMPLE TRANSITION PROBABILITIES")
    print("=" * 80)
    print("\nShowing transitions with probability > 0.1:")
    print(f"\n{'From Word':<15} {'To Word':<15} {'Symbol':<15} {'Prob':<10} {'Count':<10}")
    print("-" * 80)
    
    transitions = []
    for from_symbol, to_dict in transition_matrix.items():
        from_word = symbol_to_word[from_symbol]
        for to_symbol, prob in to_dict.items():
            if prob > 0.1:  # Only show significant transitions
                to_word = symbol_to_word[to_symbol]
                count = transition_counts[from_symbol][to_symbol]
                transitions.append((from_word, to_word, f"{from_symbol}→{to_symbol}", prob, count))
    
    # Sort by probability
    transitions.sort(key=lambda x: x[3], reverse=True)
    
    for from_word, to_word, symbol_trans, prob, count in transitions[:30]:
        print(f"{from_word:<15} {to_word:<15} {symbol_trans:<15} {prob:<10.4f} {count:<10}")
    
    if len(transitions) > 30:
        print(f"... and {len(transitions) - 30} more transitions")
    
    print(f"\n" + "=" * 80)
    print("MATRIX VALIDATION")
    print("=" * 80)
    
    # Check if each row sums to 1
    all_valid = True
    for symbol in sorted(word_to_symbol.values())[:5]:  # Check first 5
        row_sum = sum(transition_matrix[symbol].values())
        is_valid = abs(row_sum - 1.0) < 0.0001
        print(f"Symbol {symbol} ({symbol_to_word[symbol]}): sum = {row_sum:.6f} - {'✓' if is_valid else '✗'}")
        all_valid = all_valid and is_valid
    
    print(f"\nAll rows sum to 1.0: {'Yes' if all_valid else 'No'}")


def main():
    # Get sample text (approximately 300 letters)
    text = get_sample_text()
    
    # Option: Use custom text (uncomment to use)
    # text = "Your custom text here..."
    
    # Preprocess text to extract words
    words = preprocess_text(text)
    
    print(f"Processing text with {len(text)} characters and {len(words)} words...")
    
    # Create word to symbol mapping
    word_to_symbol, symbol_to_word = create_word_to_symbol_mapping(words)
    
    # Calculate transition matrix
    transition_matrix, transition_counts, symbols = calculate_word_transition_matrix(
        words, word_to_symbol
    )
    
    display_results(text, words, word_to_symbol, symbol_to_word, 
                   transition_matrix, transition_counts)
    
    json_filename = 'word_transition_matrix.json'
    save_to_json(text, words, word_to_symbol, symbol_to_word, 
                transition_matrix, transition_counts, json_filename)


if __name__ == "__main__":
    main()