"""English text generator using transition probabilities from word_transition_matrix.json."""

import json
import random
from pathlib import Path
from typing import Dict, List, Tuple
from collections import Counter


def load_word_transition_matrix(json_path: Path) -> Tuple[Dict, Dict, Dict]:
    """Load word transition probabilities and mappings from JSON file."""
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    transition_matrix = data['transition_matrix']
    symbol_to_word = data['symbol_to_word_mapping']
    word_to_symbol = data['word_to_symbol_mapping']
    
    return transition_matrix, symbol_to_word, word_to_symbol


def generate_text(transition_matrix: Dict[str, Dict[str, float]], 
                  symbol_to_word: Dict[str, str],
                  word_to_symbol: Dict[str, str],
                  num_words: int,
                  start_word: str = None) -> str:
    """Generate text using word transition probabilities.
    
    Args:
        transition_matrix: Symbol-to-symbol transition probabilities
        symbol_to_word: Mapping from symbols to words
        word_to_symbol: Mapping from words to symbols
        num_words: Number of words to generate
        start_word: Starting word (random if None)
    
    Returns:
        Generated text as string
    """
    
    if num_words < 1:
        raise ValueError("Number of words must be at least 1")
    
    # Select starting word
    if start_word is None:
        current_word = random.choice(list(word_to_symbol.keys()))
    else:
        if start_word.lower() not in word_to_symbol:
            print(f"Warning: '{start_word}' not in vocabulary. Using random start.")
            current_word = random.choice(list(word_to_symbol.keys()))
        else:
            current_word = start_word.lower()
    
    words = [current_word]
    current_symbol = word_to_symbol[current_word]
    
    # Generate remaining words using transition probabilities
    for _ in range(num_words - 1):
        # Get probabilities for next word given current symbol
        next_symbol_probs = transition_matrix[current_symbol]
        
        # Filter out zero probabilities for efficiency
        valid_symbols = [(sym, prob) for sym, prob in next_symbol_probs.items() if prob > 0]
        
        if not valid_symbols:
            # If no valid transitions, choose random word
            current_symbol = random.choice(list(symbol_to_word.keys()))
        else:
            # Create weighted choice based on probabilities
            symbols = [sym for sym, _ in valid_symbols]
            weights = [prob for _, prob in valid_symbols]
            
            # Select next symbol using weighted random choice
            next_symbol = random.choices(symbols, weights=weights, k=1)[0]
            current_symbol = next_symbol
        
        next_word = symbol_to_word[current_symbol]
        words.append(next_word)
    
    # Capitalize first word and join
    text = ' '.join(words)
    return text.capitalize()


def display_statistics(original_text: str, generated_text: str, 
                       word_to_symbol: Dict[str, str]):
    """Display comparison statistics between original and generated text."""
    
    import re
    
    def tokenize(text: str) -> List[str]:
        return re.findall(r'\b[a-z]+\b', text.lower())
    
    orig_words = tokenize(original_text)
    gen_words = tokenize(generated_text)
    
    print("\n" + "=" * 80)
    print("TEXT COMPARISON")
    print("=" * 80)
    
    print(f"\nOriginal text length: {len(original_text)} characters, {len(orig_words)} words")
    print(f"Generated text length: {len(generated_text)} characters, {len(gen_words)} words")
    
    print(f"\nOriginal text:\n{original_text[:200]}...")
    print(f"\nGenerated text:\n{generated_text}")
    
    print("\nWord frequency comparison (top 10):")
    print(f"{'Word':<15} {'Original':<15} {'Generated':<15}")
    print("-" * 45)
    
    orig_counter = Counter(orig_words)
    gen_counter = Counter(gen_words)
    
    all_words = set(orig_counter.keys()) | set(gen_counter.keys())
    word_freqs = [(w, orig_counter.get(w, 0), gen_counter.get(w, 0)) for w in all_words]
    word_freqs.sort(key=lambda x: x[1], reverse=True)
    
    for word, orig_count, gen_count in word_freqs[:10]:
        print(f"{word:<15} {orig_count:<15} {gen_count:<15}")
    
    print(f"\nUnique words - Original: {len(orig_counter)}, Generated: {len(gen_counter)}")
    
    # Calculate vocabulary overlap
    orig_vocab = set(orig_words)
    gen_vocab = set(gen_words)
    overlap = orig_vocab & gen_vocab
    
    print(f"Vocabulary overlap: {len(overlap)}/{len(orig_vocab)} " +
          f"({(len(overlap)/len(orig_vocab)*100):.1f}% of original)")


def main():
    # Load transition matrix
    json_path = Path(__file__).with_name("word_transition_matrix.json")
    
    if not json_path.exists():
        print(f"Error: {json_path} not found!")
        print("Please run ex3.py first to generate the word transition matrix.")
        return
    
    print("Loading word transition matrix...")
    transition_matrix, symbol_to_word, word_to_symbol = load_word_transition_matrix(json_path)
    
    print(f"Loaded transition matrix with {len(word_to_symbol)} unique words")
    
    # Read original text
    with open(json_path, 'r') as f:
        data = json.load(f)
    original_text = data['original_text']
    original_word_count = data['word_count']
    
    # Generate new text
    print("\n" + "=" * 80)
    print("GENERATING TEXT")
    print("=" * 80)
    
    # Generate text of same length as original
    print(f"\nGenerating text with {original_word_count} words...")
    generated = generate_text(transition_matrix, symbol_to_word, word_to_symbol, 
                             original_word_count)
    
    display_statistics(original_text, generated, word_to_symbol)
    
    # Generate additional texts
    print("\n" + "=" * 80)
    print("ADDITIONAL GENERATED TEXTS")
    print("=" * 80)
    
    word_counts = [10, 20, 30]
    for count in word_counts:
        text = generate_text(transition_matrix, symbol_to_word, word_to_symbol, count)
        print(f"\n{count} words: {text}")
    
    # Generate texts starting with specific words
    print("\n" + "=" * 80)
    print("TEXTS STARTING WITH SPECIFIC WORDS (15 words each)")
    print("=" * 80)
    
    start_words = ['the', 'autumn', 'students', 'learning']
    for start_word in start_words:
        if start_word in word_to_symbol:
            text = generate_text(transition_matrix, symbol_to_word, word_to_symbol, 
                               15, start_word=start_word)
            print(f"\nStarting with '{start_word}':\n{text}")
        else:
            print(f"\nWord '{start_word}' not in vocabulary, skipping...")


if __name__ == "__main__":
    main()
