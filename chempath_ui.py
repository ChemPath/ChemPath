import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rdkit import Chem
from rdkit.Chem import Draw
from chempath_api import ChemPathAPI

class ChemPathUI:
    def __init__(self, master):
        self.master = master
        self.master.title("ChemPath Advanced UI")
        self.api = ChemPathAPI("chempath_database.db")
        
        self.setup_ui()

    def setup_ui(self):
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.search_frame = ttk.Frame(self.notebook)
        self.optimize_frame = ttk.Frame(self.notebook)
        self.retrosynthesis_frame = ttk.Frame(self.notebook)

        self.notebook.add(self.search_frame, text="Search")
        self.notebook.add(self.optimize_frame, text="Optimize")
        self.notebook.add(self.retrosynthesis_frame, text="Retrosynthesis")

        self.setup_search_tab()
        self.setup_optimize_tab()
        self.setup_retrosynthesis_tab()

    def setup_search_tab(self):
        ttk.Label(self.search_frame, text="Search Compounds:").pack(pady=10)
        self.search_entry = ttk.Entry(self.search_frame, width=50)
        self.search_entry.pack(pady=5)
        ttk.Button(self.search_frame, text="Search", command=self.perform_search).pack(pady=5)

        self.search_results = ttk.Treeview(self.search_frame, columns=("Name", "SMILES", "Molecular Weight", "LogP"), show="headings")
        self.search_results.heading("Name", text="Name")
        self.search_results.heading("SMILES", text="SMILES")
        self.search_results.heading("Molecular Weight", text="Molecular Weight")
        self.search_results.heading("LogP", text="LogP")
        self.search_results.pack(pady=10, fill=tk.BOTH, expand=True)



    def setup_optimize_tab(self):
        ttk.Label(self.optimize_frame, text="Enter SMILES to Optimize:").pack(pady=10)
        self.optimize_entry = ttk.Entry(self.optimize_frame, width=50)
        self.optimize_entry.pack(pady=5)
        ttk.Button(self.optimize_frame, text="Optimize", command=self.perform_optimization).pack(pady=5)

        self.optimization_result = tk.Text(self.optimize_frame, height=10, width=50)
        self.optimization_result.pack(pady=10)

        self.structure_canvas = tk.Canvas(self.optimize_frame, width=400, height=400)
        self.structure_canvas.pack(pady=10)

    def setup_retrosynthesis_tab(self):
        ttk.Label(self.retrosynthesis_frame, text="Enter SMILES for Retrosynthesis:").pack(pady=10)
        self.retro_entry = ttk.Entry(self.retrosynthesis_frame, width=50)
        self.retro_entry.pack(pady=5)
        ttk.Button(self.retrosynthesis_frame, text="Analyze", command=self.perform_retrosynthesis).pack(pady=5)

        self.retro_result = tk.Text(self.retrosynthesis_frame, height=10, width=50)
        self.retro_result.pack(pady=10)

        self.retro_canvas = tk.Canvas(self.retrosynthesis_frame, width=600, height=400)
        self.retro_canvas.pack(pady=10)

    
    def perform_search(self):
        print("Performing search...")
        query = self.search_entry.get()
        results = self.api.search_compounds(query)
        self.search_results.delete(*self.search_results.get_children())
        for result in results:
            display_values = result[:4]  # Ensure we're using the first 4 elements
            print(f"Display values: {display_values}")
            self.search_results.insert("", "end", values=display_values)
            print(f"Inserted values: {display_values}")
        print(f"Total results inserted: {len(results)}")



    def update_status(self, message):
        if not hasattr(self, 'status_label'):
            self.status_label = ttk.Label(self.search_frame, text="")
            self.status_label.pack(pady=5)
        self.status_label.config(text=message)

    def validate_query(self, query):
        return query.replace("'", "''")  # Basic SQL injection prevention


    def perform_optimization(self):
        smiles = self.optimize_entry.get()
        optimized_smiles = self.api.optimize_structure(smiles)
        self.optimization_result.delete(1.0, tk.END)
        self.optimization_result.insert(tk.END, f"Optimized SMILES: {optimized_smiles}")
        self.draw_structure(smiles, optimized_smiles)

    def perform_retrosynthesis(self):
        smiles = self.retro_entry.get()
        retro_steps = self.api.perform_retrosynthesis(smiles)
        self.retro_result.delete(1.0, tk.END)
        self.retro_result.insert(tk.END, "Retrosynthesis Steps:\n" + "\n".join(retro_steps))
        self.draw_retrosynthesis(smiles, retro_steps)

    def draw_structure(self, original_smiles, optimized_smiles):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        
        mol1 = Chem.MolFromSmiles(original_smiles)
        mol2 = Chem.MolFromSmiles(optimized_smiles)
        
        img1 = Draw.MolToImage(mol1)
        img2 = Draw.MolToImage(mol2)
        
        ax1.imshow(img1)
        ax1.axis('off')
        ax1.set_title("Original")
        
        ax2.imshow(img2)
        ax2.axis('off')
        ax2.set_title("Optimized")
        
        canvas = FigureCanvasTkAgg(fig, master=self.structure_canvas)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def draw_retrosynthesis(self, target_smiles, retro_steps):
        fig, ax = plt.subplots(figsize=(10, 6))
        
        mol = Chem.MolFromSmiles(target_smiles)
        img = Draw.MolToImage(mol)
        
        ax.imshow(img)
        ax.axis('off')
        ax.set_title("Target Molecule")
        
        for i, step in enumerate(retro_steps):
            ax.text(0.1, 0.9 - i*0.1, f"Step {i+1}: {step}", transform=ax.transAxes)
        
        canvas = FigureCanvasTkAgg(fig, master=self.retro_canvas)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def main():
    root = tk.Tk()
    app = ChemPathUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
