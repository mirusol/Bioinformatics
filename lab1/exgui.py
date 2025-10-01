import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox


def parse_fasta(filename):
    """Parse FASTA file and return sequences."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequences = []
    current_sequence = ""
    current_header = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_sequence:
                sequences.append({'header': current_header, 'sequence': current_sequence})
            current_header = line
            current_sequence = ""
        else:
            current_sequence += line.replace(" ", "").upper()
    
    if current_sequence:
        sequences.append({'header': current_header, 'sequence': current_sequence})
    
    return sequences


def analyze_sequence(s):
    """Analyze sequence - original algorithm."""
    s = s.replace(" ", "")
    
    # Build alphabet
    alphabet = []
    for c in s:
        if c not in alphabet:
            alphabet.append(c)
    
    # Calculate percentages
    results = []
    for letter in alphabet:
        count = s.count(letter)
        percent = count * 100 / len(s)
        results.append((letter, count, percent))
    
    return alphabet, results


def load_file():
    """Load and analyze FASTA file."""
    filename = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    
    if not filename:
        return
    
    try:
        sequences = parse_fasta(filename)
        
        if not sequences:
            messagebox.showerror("Error", "No sequences found!")
            return
        
        # Clear previous results
        text_area.delete(1.0, tk.END)
        
        # Display results
        text_area.insert(tk.END, f"File: {filename}\n")
        text_area.insert(tk.END, f"Found {len(sequences)} sequence(s)\n")
        text_area.insert(tk.END, "="*60 + "\n\n")
        
        for idx, seq_data in enumerate(sequences, 1):
            sequence = seq_data['sequence']
            header = seq_data['header']
            
            text_area.insert(tk.END, f"Sequence {idx}: {header}\n")
            text_area.insert(tk.END, f"Length: {len(sequence)}\n\n")
            
            # Analyze
            alphabet, results = analyze_sequence(sequence)
            
            text_area.insert(tk.END, f"Alphabet: {''.join(alphabet)}\n\n")
            
            # Display percentages
            for letter, count, percent in results:
                text_area.insert(tk.END, f"{letter}: {percent:.2f}% (count: {count})\n")
            
            text_area.insert(tk.END, "\n" + "="*60 + "\n\n")
        
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")


def clear_results():
    """Clear the text area."""
    text_area.delete(1.0, tk.END)


# Create main window
root = tk.Tk()
root.title("FASTA Sequence Analyzer")
root.geometry("700x500")

# Buttons frame
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

load_btn = tk.Button(button_frame, text="Load FASTA File", command=load_file, width=15)
load_btn.pack(side=tk.LEFT, padx=5)

clear_btn = tk.Button(button_frame, text="Clear", command=clear_results, width=15)
clear_btn.pack(side=tk.LEFT, padx=5)

# Text area for results
text_area = scrolledtext.ScrolledText(root, font=("Courier", 10), wrap=tk.WORD)
text_area.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Initial message
text_area.insert(tk.END, "Click 'Load FASTA File' to begin...\n")

root.mainloop()