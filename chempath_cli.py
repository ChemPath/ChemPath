from chempath_api import ChemPathAPI
from database_operations import get_compound_by_smiles, get_retrosynthesis_data, store_retrosynthesis_informed_optimization
from chempath_core import retrosynthesis_informed_optimization

def add_compound(api):
    name = input("Enter compound name: ")
    smiles = input("Enter SMILES string: ")
    molecular_weight = float(input("Enter molecular weight: "))
    logp = float(input("Enter LogP value: "))
    h_bond_donors = int(input("Enter number of H-bond donors: "))
    h_bond_acceptors = int(input("Enter number of H-bond acceptors: "))

    compound = {
        'name': name,
        'smiles': smiles,
        'molecular_weight': molecular_weight,
        'logp': logp,
        'h_bond_donors': h_bond_donors,
        'h_bond_acceptors': h_bond_acceptors
    }

    api.add_compound(compound)
    print(f"Compound {name} added successfully.")

def main():
    api = ChemPathAPI("C:/Users/Dr. Contessa Petrini/ChemPath/chempath_database.db")

    while True:
        print("\nChemPath Menu:")
        print("1. Search compounds")
        print("2. Add compound")
        print("3. Get therapeutic areas")
        print("4. Predict therapeutic areas")
        print("5. Optimize structure")
        print("6. Explore chemical space")
        print("7. Generate analogs")
        print("8. Perform retrosynthesis")
        print("9. Analyze drug-likeness")
        print("10. Predict ADMET properties")
        print("11. Perform molecular docking")
        print("12. Analyze protein-ligand interactions")
        print("13. Generate pharmacophore model")
        print("14. Perform virtual screening")
        print("15. Analyze structure-activity relationships")
        print("16. Retrosynthesis-Informed Optimization")
        print("17. Train ML models")
        print("18. Predict retrosynthesis feasibility")
        print("19. Predict reaction class")
        print("20. AI-driven Optimization Prioritization")
        print("21. Exit")

        choice = input("Enter your choice (1-21): ")

        if choice == '1':
            # Implement search compounds
            pass
        elif choice == '2':
            add_compound(api)
        elif choice == '3':
            # Implement get therapeutic areas
            pass
        elif choice == '4':
            # Implement predict therapeutic areas
            pass
        elif choice == '5':
            # Implement optimize structure
            pass
        elif choice == '6':
            # Implement explore chemical space
            pass
        elif choice == '7':
            # Implement generate analogs
            pass
        elif choice == '8':
            # Implement perform retrosynthesis
            pass
        elif choice == '9':
            # Implement analyze drug-likeness
            pass
        elif choice == '10':
            # Implement predict ADMET properties
            pass
        elif choice == '11':
            # Implement perform molecular docking
            pass
        elif choice == '12':
            # Implement analyze protein-ligand interactions
            pass
        elif choice == '13':
            # Implement generate pharmacophore model
            pass
        elif choice == '14':
            # Implement perform virtual screening
            pass
        elif choice == '15':
            # Implement analyze structure-activity relationships
            pass
        elif choice == '16':
            smiles = input("Enter the SMILES string of the compound: ")
            compound = get_compound_by_smiles(api.conn, smiles)
            if compound:
                retrosynthesis_data = get_retrosynthesis_data(api.conn, compound[0])
                optimized_smiles = retrosynthesis_informed_optimization(smiles, retrosynthesis_data)
                print(f"Retrosynthesis-informed optimized structure: {optimized_smiles}")
                store_retrosynthesis_informed_optimization(api.conn, compound[0], optimized_smiles)
            else:
                print("Compound not found in the database.")
        elif choice == '17':
            api.train_ml_models()
            print("ML models trained successfully.")
        elif choice == '18':
            smiles = input("Enter SMILES string: ")
            feasibility = api.predict_retrosynthesis_feasibility(smiles)
            print(f"Retrosynthesis feasibility: {feasibility}")
        elif choice == '19':
            smiles = input("Enter SMILES string: ")
            reaction_class = api.predict_reaction_class(smiles)
            print(f"Predicted reaction class: {reaction_class}")
        elif choice == '20':
            compounds = api.get_all_compounds()
            optimized_compounds = api.optimize_compounds(compounds)
            print("Top 5 compounds prioritized for optimization:")
            for i, compound in enumerate(optimized_compounds[:5], 1):
                print(f"\n{i}. {compound['name']}:")
                print(f"   Priority Score: {compound['optimization_priority']:.2f}")
                print(f"   Molecular Weight: {compound['molecular_weight']:.2f}")
                print(f"   LogP: {compound['logp']:.2f}")
                print(f"   H-Bond Donors: {compound['h_bond_donors']}")
                print(f"   H-Bond Acceptors: {compound['h_bond_acceptors']}")
        elif choice == '21':
            print("Thank you for using ChemPath. Goodbye!")
            api.close_connection()
            break
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()
