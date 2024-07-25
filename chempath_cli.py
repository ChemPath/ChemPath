from chempath_api import ChemPathAPI

def main():
    api = ChemPathAPI('chempath_database.db')

    while True:
        print("\nChemPath Menu:")
        print("1. Search compounds")
        print("2. Add compound")
        print("3. Predict therapeutic areas")
        print("4. Exit")

        choice = input("Enter your choice (1-4): ")

        if choice == '1':
            query = input("Enter search query: ")
            compounds, total_count = api.search(query)
            # Display compounds here

        elif choice == '2':
            # Collect compound data from user
            compound_data = {}  # Fill this dictionary with user input
            api.add_compound(compound_data)

        elif choice == '3':
            smiles = input("Enter SMILES string: ")
            areas = api.predict_therapeutic_areas(smiles)
            print("Predicted areas:", areas)

        elif choice == '4':
            break

    api.close_connection()

if __name__ == '__main__':
    main()
