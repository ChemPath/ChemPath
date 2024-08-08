import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rdkit import Chem
from rdkit.Chem import Draw
from chempath_api import ChemPathAPI
import logging
from src.core.chempath_core import optimize_structure, predict_properties, retrosynthesize
from chempath_api import ChemPathAPI
from src.core.chempath_core import create_connection, create_tables, create_indexes, insert_compound
from models.advanced_retrosynthesis import advanced_retrosynthetic_analysis, get_advanced_retrosynthesis_result

class ChemPathUI:
    def __init__(self, master, db_path):
        self.master = master
        self.db_path = db_path
        self.api = ChemPathAPI(db_path)
        
        self.setup_ui()

    def setup_ui(self):
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.search_frame = ttk.Frame(self.notebook)
        self.optimize_frame = ttk.Frame(self.notebook)
        self.retrosynthesis_frame = ttk.Frame(self.notebook)
        self.advanced_retro_frame = ttk.Frame(self.notebook)

        self.notebook.add(self.search_frame, text="Search")
        self.notebook.add(self.optimize_frame, text="Optimize")
        self.notebook.add(self.retrosynthesis_frame, text="Retrosynthesis")
        self.notebook.add(self.advanced_retro_frame, text="Advanced Retrosynthesis")

        self.setup_search_tab()
        self.setup_optimize_tab()
        self.setup_retrosynthesis_tab()
        self.setup_advanced_retro_tab()

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
        ttk.Button(self.retrosynthesis_frame, text="Analyze", command=lambda: self.perform_retrosynthesis(self.retro_entry.get())).pack(pady=5)

        self.retro_result = tk.Text(self.retrosynthesis_frame, height=10, width=50)
        self.retro_result.pack(pady=10)

        self.retro_canvas = tk.Canvas(self.retrosynthesis_frame, width=600, height=400)
        self.retro_canvas.pack(pady=10)

    def setup_advanced_retro_tab(self):
        ttk.Label(self.advanced_retro_frame, text="Enter SMILES for Advanced Retrosynthesis:").pack(pady=10)
        self.adv_retro_entry = ttk.Entry(self.advanced_retro_frame, width=50)
        self.adv_retro_entry.pack(pady=5)

        ttk.Label(self.advanced_retro_frame, text="Depth:").pack()
        self.depth_var = tk.StringVar(value="3")
        self.depth_entry = ttk.Entry(self.advanced_retro_frame, textvariable=self.depth_var, width=5)
        self.depth_entry.pack()

        ttk.Label(self.advanced_retro_frame, text="Strategy:").pack()
        self.strategy_var = tk.StringVar(value="dfs")
        self.strategy_combo = ttk.Combobox(self.advanced_retro_frame, textvariable=self.strategy_var, values=["dfs", "bfs", "dijkstra", "astar"])
        self.strategy_combo.pack()

        ttk.Label(self.advanced_retro_frame, text="Weights:").pack()
        self.weight_frame = ttk.Frame(self.advanced_retro_frame)
        self.weight_frame.pack()

        self.feasibility_weight = self.create_weight_input(self.weight_frame, "Feasibility:", 0, 0, 0.4)
        self.availability_weight = self.create_weight_input(self.weight_frame, "Availability:", 0, 1, 0.3)
        self.cost_weight = self.create_weight_input(self.weight_frame, "Cost:", 0, 2, 0.3)

        ttk.Label(self.advanced_retro_frame, text="Criteria:").pack()
        self.criteria_frame = ttk.Frame(self.advanced_retro_frame)
        self.criteria_frame.pack()

        self.min_feasibility = self.create_criteria_input(self.criteria_frame, "Min Feasibility:", 0, 0, 0.6)
        self.min_availability = self.create_criteria_input(self.criteria_frame, "Min Availability:", 0, 1, 0.5)
        self.max_cost = self.create_criteria_input(self.criteria_frame, "Max Cost:", 0, 2, 100)

        ttk.Button(self.advanced_retro_frame, text="Analyze", command=self.perform_advanced_retrosynthesis).pack(pady=5)

        self.adv_retro_result = tk.Text(self.advanced_retro_frame, height=10, width=50)
        self.adv_retro_result.pack(pady=10)

        self.adv_retro_canvas = tk.Canvas(self.advanced_retro_frame, width=600, height=400)
        self.adv_retro_canvas.pack(pady=10)

    def create_weight_input(self, parent, label, row, col, default):
        ttk.Label(parent, text=label).grid(row=row, column=col, padx=5, pady=2)
        var = tk.DoubleVar(value=default)
        entry = ttk.Entry(parent, textvariable=var, width=5)
        entry.grid(row=row+1, column=col, padx=5, pady=2)
        return var

    def create_criteria_input(self, parent, label, row, col, default):
        ttk.Label(parent, text=label).grid(row=row, column=col, padx=5, pady=2)
        var = tk.DoubleVar(value=default)
        entry = ttk.Entry(parent, textvariable=var, width=5)
        entry.grid(row=row+1, column=col, padx=5, pady=2)
        return var

    def perform_search(self):
        print("Performing search...")
        query = self.search_entry.get()
        results, total_count = self.api.search_compounds(query)
        self.search_results.delete(*self.search_results.get_children())
        for result in results:
            display_values = result[:4]  # Ensure we're using the first 4 elements
            print(f"Display values: {display_values}")
            self.search_results.insert("", "end", values=display_values)
            print(f"Inserted values: {display_values}")
        print(f"Total results inserted: {len(results)}")
        self.update_status(f"Found {total_count} results")

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

    def perform_retrosynthesis(self, smiles):
        retro_steps = self.api.perform_retrosynthesis(smiles)
        self.retro_result.delete(1.0, tk.END)
        if retro_steps and 'steps' in retro_steps:
            steps_text = "\n".join([f"{step['reaction_name']}: {' + '.join(step['products'])}" for step in retro_steps['steps']])
            self.retro_result.insert(tk.END, f"Retrosynthesis Steps:\n{steps_text}")
        else:
            self.retro_result.insert(tk.END, "No retrosynthesis steps found.")
        self.draw_retrosynthesis(smiles, retro_steps)

    def perform_advanced_retrosynthesis(self):
        smiles = self.adv_retro_entry.get()
        depth = int(self.depth_var.get())
        strategy = self.strategy_var.get()

        weights = {
            'feasibility': self.feasibility_weight.get(),
            'availability': self.availability_weight.get(),
            'cost': self.cost_weight.get()
        }

        criteria = {
            'min_feasibility': self.min_feasibility.get(),
            'min_availability': self.min_availability.get(),
            'max_cost': self.max_cost.get()
        }

        result = advanced_retrosynthetic_analysis(1, smiles, depth, strategy, weights, criteria)
        self.display_advanced_retro_result(result)

    def display_advanced_retro_result(self, result):
        self.adv_retro_result.delete(1.0, tk.END)
        self.adv_retro_result.insert(tk.END, f"Result ID: {result['result_id']}\n")
        self.adv_retro_result.insert(tk.END, f"Number of Nodes: {result['num_nodes']}\n")
        self.adv_retro_result.insert(tk.END, f"Number of Edges: {result['num_edges']}\n")
        self.adv_retro_result.insert(tk.END, f"Number of Pathways: {result['num_pathways']}\n")
        self.adv_retro_result.insert(tk.END, "\nTop Pathways:\n")
        
        for i, pathway in enumerate(result['top_pathways'], 1):
            self.adv_retro_result.insert(tk.END, f"Pathway {i}:\n")
            for step in pathway:
                self.adv_retro_result.insert(tk.END, f"  {' + '.join(step['reactants'])} -> {' + '.join(step['products'])}\n")

        # Display the retrosynthesis tree image
        img = plt.imread(result['tree_image'])
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(img)
        ax.axis('off')
        
        canvas = FigureCanvasTkAgg(fig, master=self.adv_retro_canvas)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

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
    
        if retro_steps and 'steps' in retro_steps:
            for i, step in enumerate(retro_steps['steps']):
                ax.text(0.1, 0.9 - i*0.1, f"Step {i+1}: {step['reaction_name']}", transform=ax.transAxes)
    
        canvas = FigureCanvasTkAgg(fig, master=self.retro_canvas)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def main():
    root = tk.Tk()
    db_path = "chempath_database.db"  # or use a configuration file to set this
    app = ChemPathUI(root, db_path)
    root.mainloop()

if __name__ == '__main__':
    main()

