import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class FastaAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sliding Window Analyzer")
        self.root.geometry("1000x700")
        
        self.sequence = ""
        self.window_size = 30
        
        self.setup_gui()
    
    def setup_gui(self):
        # Top frame for controls
        control_frame = tk.Frame(self.root, padx=10, pady=10)
        control_frame.pack(side=tk.TOP, fill=tk.X)
        
        # File selection button
        self.select_btn = tk.Button(
            control_frame, 
            text="Select FASTA File", 
            command=self.select_file,
            font=("Arial", 12),
            bg="#4CAF50",
            fg="white",
            padx=20,
            pady=10
        )
        self.select_btn.pack(side=tk.LEFT, padx=5)
        
        # Window size input
        tk.Label(control_frame, text="Window Size:", font=("Arial", 10)).pack(side=tk.LEFT, padx=5)
        self.window_entry = tk.Entry(control_frame, width=10)
        self.window_entry.insert(0, "30")
        self.window_entry.pack(side=tk.LEFT, padx=5)
        
        # Analyze button
        self.analyze_btn = tk.Button(
            control_frame,
            text="Analyze",
            command=self.analyze_sequence,
            font=("Arial", 12),
            bg="#2196F3",
            fg="white",
            padx=20,
            pady=10,
            state=tk.DISABLED
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=5)
        
        # Info label
        self.info_label = tk.Label(
            control_frame,
            text="No file selected",
            font=("Arial", 10),
            fg="gray"
        )
        self.info_label.pack(side=tk.LEFT, padx=10)
        
        # Frame for matplotlib figure
        self.plot_frame = tk.Frame(self.root)
        self.plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)
    
    def select_file(self):
        filepath = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
        )
        
        if filepath:
            try:
                self.sequence = self.read_fasta(filepath)
                if self.sequence:
                    self.info_label.config(
                        text=f"Loaded: {len(self.sequence)} nucleotides",
                        fg="green"
                    )
                    self.analyze_btn.config(state=tk.NORMAL)
                else:
                    messagebox.showerror("Error", "No sequence found in FASTA file")
            except Exception as e:
                messagebox.showerror("Error", f"Error reading file: {str(e)}")
    
    def read_fasta(self, filepath):
        """Read FASTA file and extract the sequence"""
        sequence = ""
        with open(filepath, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):  # Skip header lines
                    sequence += line.upper()
        return sequence
    
    def compute_frequencies(self):
        """Compute relative frequencies using sliding window"""
        try:
            self.window_size = int(self.window_entry.get())
            if self.window_size < 1:
                raise ValueError("Window size must be positive")
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid window size: {str(e)}")
            return None
        
        if len(self.sequence) < self.window_size:
            messagebox.showerror(
                "Error", 
                f"Sequence length ({len(self.sequence)}) is smaller than window size ({self.window_size})"
            )
            return None
        
        # Initialize frequency vectors for each nucleotide
        freq_A = []
        freq_C = []
        freq_G = []
        freq_T = []
        
        # Sliding window analysis
        num_windows = len(self.sequence) - self.window_size + 1
        
        for i in range(num_windows):
            window = self.sequence[i:i + self.window_size]
            
            # Count nucleotides in current window
            count_A = window.count('A')
            count_C = window.count('C')
            count_G = window.count('G')
            count_T = window.count('T')
            
            # Compute relative frequencies
            freq_A.append(count_A / self.window_size)
            freq_C.append(count_C / self.window_size)
            freq_G.append(count_G / self.window_size)
            freq_T.append(count_T / self.window_size)
        
        return {
            'A': freq_A,
            'C': freq_C,
            'G': freq_G,
            'T': freq_T
        }
    
    def analyze_sequence(self):
        """Analyze sequence and plot results"""
        frequencies = self.compute_frequencies()
        
        if frequencies is None:
            return
        
        self.plot_frequencies(frequencies)
    
    def plot_frequencies(self, frequencies):
        """Plot the frequency signals"""
        # Clear previous plot
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Position vector (x-axis)
        positions = range(len(frequencies['A']))
        
        # Plot each nucleotide signal
        ax.plot(positions, frequencies['A'], label='A (Adenine)', color='red', linewidth=1.5, alpha=0.8)
        ax.plot(positions, frequencies['C'], label='C (Cytosine)', color='blue', linewidth=1.5, alpha=0.8)
        ax.plot(positions, frequencies['G'], label='G (Guanine)', color='green', linewidth=1.5, alpha=0.8)
        ax.plot(positions, frequencies['T'], label='T (Thymine)', color='orange', linewidth=1.5, alpha=0.8)
        
        # Customize plot
        ax.set_xlabel('Window Position', fontsize=12)
        ax.set_ylabel('Relative Frequency', fontsize=12)
        ax.set_title(f'Nucleotide Frequency Analysis (Window Size: {self.window_size})', fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)
        
        # Embed plot in tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


def main():
    root = tk.Tk()
    app = FastaAnalyzerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()