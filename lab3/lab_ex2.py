import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def read_fasta(file_path):
    sequence = ""
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if not line.startswith('>'): 
                    sequence += line.strip()
        return sequence
    except Exception as e:
        messagebox.showerror("Error", f"Error reading FASTA file: {str(e)}")
        return None

def calculate_tm_basic(sequence):
    seq = sequence.upper()
    return 4 * (seq.count('G') + seq.count('C')) + 2 * (seq.count('A') + seq.count('T'))

def calculate_tm_advanced(sequence, na_conc=0.05):
    seq = sequence.upper()
    length = len(seq)
    gc_pct = (seq.count('G') + seq.count('C')) / length * 100
    return 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_pct - (600 / length)

def sliding_window_analysis(sequence, window_size=8):
    results = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        results.append({
            'position': i+1,
            'window': window,
            'basic_tm': calculate_tm_basic(window),
            'advanced_tm': calculate_tm_advanced(window)
        })
    return results

class app:
    def __init__(self, root):
        self.root = root
        self.root.geometry("1000x700")
        
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(2, weight=1)
        
        self.create_file_frame(main_frame)
        
        self.create_controls_frame(main_frame)
        
        self.create_results_notebook(main_frame)

    def create_file_frame(self, parent):
        frame = ttk.LabelFrame(parent, text="FASTA File Selection", padding="5")
        frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5)
        
        self.file_path = tk.StringVar()
        ttk.Entry(frame, textvariable=self.file_path, width=70).grid(row=0, column=0, padx=5)
        ttk.Button(frame, text="Browse", command=self.browse_file).grid(row=0, column=1, padx=5)

    def create_controls_frame(self, parent):
        frame = ttk.LabelFrame(parent, text="Analysis Controls", padding="5")
        frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(frame, text="Window Size:").grid(row=0, column=0, padx=5)
        self.window_size = ttk.Entry(frame, width=10)
        self.window_size.insert(0, "8")
        self.window_size.grid(row=0, column=1, padx=5)
        
        ttk.Button(frame, text="Analyze", command=self.analyze).grid(row=0, column=2, padx=5)
        ttk.Button(frame, text="Clear", command=self.clear).grid(row=0, column=3, padx=5)

    def create_results_notebook(self, parent):
        notebook = ttk.Notebook(parent)
        notebook.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        results_frame = ttk.Frame(notebook, padding="5")
        notebook.add(results_frame, text="Data Table")
        
        self.results_text = scrolledtext.ScrolledText(results_frame, width=80, height=25)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        
        plot_frame = ttk.Frame(notebook, padding="5")
        notebook.add(plot_frame, text="Plot")
        
        self.figure = plt.Figure(figsize=(8, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )
        if filename:
            self.file_path.set(filename)

    def analyze(self):
        if not self.file_path.get():
            messagebox.showerror("Error", "Please select a FASTA file")
            return
        
        try:
            window_size = int(self.window_size.get())
            if window_size < 1:
                raise ValueError("Window size must be positive")
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid window size: {e}")
            return
        
        sequence = read_fasta(self.file_path.get())
        if not sequence:
            return
        
        if len(sequence) < window_size:
            messagebox.showerror("Error", f"Sequence length ({len(sequence)}) is shorter than window size ({window_size})")
            return
        
        results = sliding_window_analysis(sequence, window_size)
        
        self.display_results(results)
        
        self.plot_results(results)

    def display_results(self, results):
        self.results_text.delete(1.0, tk.END)
        header = f"{'Position':<10}{'Window':<12}{'Basic Tm':<15}{'Advanced Tm':<15}\n"
        self.results_text.insert(tk.END, header)
        self.results_text.insert(tk.END, "-" * 50 + "\n")
        
        for r in results:
            line = f"{r['position']:<10}{r['window']:<12}{r['basic_tm']:.2f}°C{'':<7}{r['advanced_tm']:.2f}°C\n"
            self.results_text.insert(tk.END, line)

    def plot_results(self, results):
        self.ax.clear()
        
        positions = [r['position'] for r in results]
        basic_tms = [r['basic_tm'] for r in results]
        advanced_tms = [r['advanced_tm'] for r in results]
        
        self.ax.plot(positions, basic_tms, 'b-', label='Basic Tm', linewidth=2)
        self.ax.plot(positions, advanced_tms, 'r-', label='Advanced Tm', linewidth=2)
        
        self.ax.set_xlabel('Position', fontsize=11)
        self.ax.set_ylabel('Melting Temperature (°C)', fontsize=11)
        self.ax.set_title('DNA Melting Temperature - Sliding Window Analysis', fontsize=12, fontweight='bold')
        self.ax.legend(loc='best')
        self.ax.grid(True, alpha=0.3)
        
        self.figure.tight_layout()
        self.canvas.draw()

    def clear(self):
        self.file_path.set("")
        self.window_size.delete(0, tk.END)
        self.window_size.insert(0, "8")
        self.results_text.delete(1.0, tk.END)
        self.ax.clear()
        self.canvas.draw()

def main():
    root = tk.Tk()
    result = app(root)
    root.mainloop()

if __name__ == "__main__":
    main()