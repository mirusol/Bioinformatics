import math
import re
from collections import defaultdict


# Mihai Eminescu 
EMINESCU_TEXT = """
In zori sub cerul plin de stele argintii
Când luna-și varsă razele pe ape-ntinse lin,
O zână-i-acolo iar și iar mi-apare mie
Și-al ei priviri purtă-o-nțeles divin.
Ființa ei plutind în văluri de lumini
Cu glas de flori și-al ochilor o-nțeles pus
Așa mi-apare-n visuri mele sfinte
Atunci când noaptea veghează-mi-a sus.
"""

# Nichita Stanescu
STANESCU_TEXT = """
Cuvintele mele sunt din piatra și din apă
Și-ntre ele-n golul vostru curg ploi negre și bărci
Vine tăcerea după fiecare vers
Ca o umbră care-și caută pe tine-ntru întuneric.
Să vorbesc cu dumnezeu prin tine-ntru codru
Și să-nțeleg că viața mea e doar un semn
Că toată dragostea mea-i doar o durere
Și toată lumina mea e doar cuvânt.
"""

MIHAI_TEXT = """
În zori sub cerul plin de stele argintii
Cuvintele mele sunt din piatra și din apă
Când luna-și varsă razele pe ape-ntinse lin
O zână-i-acolo iar și iar mi-apare mie
Vine tăcerea după fiecare vers
Și-al ei priviri purtă-o-nțeles divin
Ca o umbră care-și caută pe tine-ntru întuneric
Ființa ei plutind în văluri de lumini
Să vorbesc cu dumnezeu prin tine-ntru codru
Cu glas de flori și-al ochilor o-nțeles pus
Și să-nțeleg că viața mea e doar un semn
Că toată dragostea mea-i doar o durere
"""

def preprocess_text(text):
    """Convert text to lowercase and extract words"""
    text = text.lower()
    text = re.sub(r'[^\w\s]', '', text)
    words = [w for w in text.split() if w]
    return words

def build_transition_matrix(words):
    """Build a transition matrix from word sequences"""
    transitions = defaultdict(lambda: defaultdict(int))
    

    for i in range(len(words) - 1):
        current_word = words[i]
        next_word = words[i + 1]
        transitions[current_word][next_word] += 1
    
    prob_matrix = {}
    for word in transitions:
        total = sum(transitions[word].values())
        prob_matrix[word] = {}
        for next_word in transitions[word]:
            prob_matrix[word][next_word] = transitions[word][next_word] / total
    
    return prob_matrix

def get_transition_probability(matrix, current_word, next_word, smoothing=0.001):
    """Get probability of transition, with smoothing for unseen transitions"""
    if current_word in matrix and next_word in matrix[current_word]:
        return matrix[current_word][next_word]
    else:
        return smoothing

def calculate_log_likelihood_ratio(prob_eminescu, prob_stanescu, current_word, next_word):
    """Calculate log-likelihood ratio: log2(P_E / P_S)"""
    p_e = get_transition_probability(prob_eminescu, current_word, next_word, smoothing=0.001)
    p_s = get_transition_probability(prob_stanescu, current_word, next_word, smoothing=0.001)
    
    if p_e > 0 and p_s > 0:
        ratio = p_e / p_s
        return math.log2(ratio)
    elif p_e > 0:
        return float('inf')
    elif p_s > 0:
        return float('-inf')
    else:
        return 0.0

def analyze_with_sliding_window(words, llr_matrix, window_size=3):
    """Analyze text using sliding window"""
    scores = []
    positions = []
    
    for i in range(len(words) - window_size):
        window = words[i:i+window_size]
        window_score = 0.0
        
        # Sum LLR for transitions in the window
        for j in range(len(window) - 1):
            current = window[j]
            next_w = window[j + 1]
            llr = calculate_log_likelihood_ratio(llr_matrix['eminescu'], 
                                                 llr_matrix['stanescu'],
                                                 current, next_w)
            if not math.isinf(llr):
                window_score += llr
        
        scores.append(window_score)
        positions.append(i)
    
    return scores, positions

def classify_segment(score, threshold=0.3):
    """Classify a segment based on its score"""
    if score > threshold:
        return "EMINESCU", score
    elif score < -threshold:
        return "STANESCU", score
    else:
        return "NEITHER", score

def print_colored_analysis(words, scores, positions, threshold=0.3):
    """Print the analysis with classification"""
    print("\n" + "="*80)
    print("PLAGIARISM ANALYSIS - SLIDING WINDOW RESULTS")
    print("="*80)
    
    classification_map = {}
    for i, pos in enumerate(positions):
        window_size = 3
        for word_idx in range(pos, min(pos + window_size, len(words))):
            if word_idx not in classification_map:
                classification_map[word_idx] = []
            classification_map[word_idx].append(scores[i])
    
    # Print the text with classifications
    print("\nText Analysis (word by word):\n")
    current_author = None
    segment_scores = []
    segment_words = []
    
    for idx, word in enumerate(words):
        if idx in classification_map:
            scores_list = classification_map[idx]
            # Use average score
            avg_score = sum(scores_list) / len(scores_list)
            author = classify_segment(avg_score, threshold)[0]
        else:
            author = "UNKNOWN"
            avg_score = 0
        
        if author != current_author:
            if segment_words:
                print_segment(segment_words, current_author, segment_scores)
            segment_words = [word]
            segment_scores = []
            current_author = author
        else:
            segment_words.append(word)
            segment_scores.append(avg_score)
    
    if segment_words:
        print_segment(segment_words, current_author, segment_scores)
    
    print("\n" + "="*80)
    print("SLIDING WINDOW SCORES")
    print("="*80)
    for i, (pos, score) in enumerate(zip(positions, scores)):
        author, score_val = classify_segment(score, threshold)
        window_text = " ".join(words[pos:pos+3])
        print(f"Position {pos:3d}: [{window_text:40s}] | Score: {score:7.3f} | Author: {author}")

def print_segment(words, author, scores):
    """Print a text segment with color coding"""
    text = " ".join(words)
    if author == "EMINESCU":
        color_code = "\033[92m"  # Green
        display = "[EMINESCU]"
    elif author == "STANESCU":
        color_code = "\033[94m"  # Blue
        display = "[STANESCU]"
    else:
        color_code = "\033[93m"  # Yellow
        display = "[NEITHER]"
    
    reset_code = "\033[0m"
    
    if scores:
        avg_score = sum(scores) / len(scores)
        print(f"{color_code}{display} {text} (avg score: {avg_score:.3f}){reset_code}")
    else:
        print(f"{color_code}{display} {text}{reset_code}")

def print_matrix_info(matrix, title, top_n=10):
    """Print information about a transition matrix"""
    print(f"\n{title}")
    print("-" * 80)
    
    print(f"Total unique words: {len(matrix)}")
    print(f"\nTop {top_n} most common starting words:")
    
    # Sort by number of outgoing transitions
    sorted_words = sorted(matrix.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (word, transitions) in enumerate(sorted_words[:top_n], 1):
        print(f"  {i:2d}. '{word}' -> {len(transitions)} different next words")

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("PLAGIARISM DETECTION SYSTEM - LITERARY ANALYSIS")
print("Case: The People vs. Mihai - Plagiarism Investigation")
print("="*80)

# Step 1: Preprocess texts
print("\n" + "="*80)
print("STEP 1: TEXT PREPROCESSING")
print("="*80)

eminescu_words = preprocess_text(EMINESCU_TEXT)
stanescu_words = preprocess_text(STANESCU_TEXT)
mihai_words = preprocess_text(MIHAI_TEXT)

print(f"\nEminescu text: {len(eminescu_words)} words")
print(f"First 10 words: {' '.join(eminescu_words[:10])}")

print(f"\nStanescu text: {len(stanescu_words)} words")
print(f"First 10 words: {' '.join(stanescu_words[:10])}")

print(f"\nMihai text: {len(mihai_words)} words")
print(f"First 10 words: {' '.join(mihai_words[:10])}")

# Step 2: Build transition matrices
print("\n" + "="*80)
print("STEP 2: BUILD TRANSITION MATRICES")
print("="*80)

prob_eminescu = build_transition_matrix(eminescu_words)
prob_stanescu = build_transition_matrix(stanescu_words)

print_matrix_info(prob_eminescu, "Eminescu Transition Matrix")
print_matrix_info(prob_stanescu, "Stanescu Transition Matrix")

# Step 3: Create LLR matrix
print("\n" + "="*80)
print("STEP 3: LOG-LIKELIHOOD RATIO MATRIX")
print("="*80)

llr_matrix = {
    'eminescu': prob_eminescu,
    'stanescu': prob_stanescu
}

print("\nLog-Likelihood Ratio = log2(P_Eminescu / P_Stanescu)")
print("Positive scores favor Eminescu")
print("Negative scores favor Stanescu")
print("Scores near 0 favor neither")

# Step 4: Analyze Mihai's text
print("\n" + "="*80)
print("STEP 4: ANALYZE MIHAI'S TEXT")
print("="*80)

scores, positions = analyze_with_sliding_window(mihai_words, llr_matrix, window_size=3)

print(f"\nAnalyzed {len(positions)} windows of 3 consecutive words")
print(f"Score statistics:")
print(f"  Mean score: {sum(scores) / len(scores):.4f}")
print(f"  Min score: {min(scores):.4f}")
print(f"  Max score: {max(scores):.4f}")

# Step 5: Print detailed analysis
print_colored_analysis(mihai_words, scores, positions, threshold=0.3)


eminescu_count = sum(1 for s in scores if s > 0.3)
stanescu_count = sum(1 for s in scores if s < -0.3)
neither_count = len(scores) - eminescu_count - stanescu_count

total = len(scores)

print(f"\nTotal windows analyzed: {total}")
print(f"Windows attributed to Eminescu: {eminescu_count} ({100*eminescu_count/total:.1f}%)")
print(f"Windows attributed to Stanescu: {stanescu_count} ({100*stanescu_count/total:.1f}%)")
print(f"Windows attributed to neither: {neither_count} ({100*neither_count/total:.1f}%)")

print(f"\nCourt Decision:")
if eminescu_count > 0 and stanescu_count > 0:
    print(f"✓ The text contains material from BOTH poets:")
    print(f"  - Approximately {100*eminescu_count/total:.1f}% matches Eminescu's style")
    print(f"  - Approximately {100*stanescu_count/total:.1f}% matches Stanescu's style")
    print(f"  - Mihai appears to have plagiarized from both sources")
elif eminescu_count > stanescu_count:
    print(f"✓ The text is primarily from Eminescu ({100*eminescu_count/total:.1f}%)")
    print(f"  - Mihai is likely guilty of plagiarizing Eminescu's work")
elif stanescu_count > eminescu_count:
    print(f"✓ The text is primarily from Stanescu ({100*stanescu_count/total:.1f}%)")
    print(f"  - Mihai is likely guilty of plagiarizing Stanescu's work")
else:
    print(f"✓ The text does not match either poet's typical style patterns")
    print(f"  - Mihai appears to have written original work or mixed sources unclear")

print("\n" + "="*80)
