import csv
import pandas as pd
from rdkit import Chem

def read_aurora_kinase_b_interactions(filepath):
    """
    Reads and parses the aurora_kinase_B_interactions.csv file,
    and returns a list of RDKit Mol objects from the 'smiles' column.

    Args:
        filepath (str): Path to the CSV file.

    Returns:
        list[rdkit.Chem.rdchem.Mol]: List of RDKit Mol objects.
    """
    data = []
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            cleaned_row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            data.append(cleaned_row)
    df = pd.DataFrame(data)
    mols = [Chem.MolFromSmiles(smiles) for smiles in df["smiles"]]
    return mols

# Load mol objects
file_path = "/vol/data/airs/AIRS/OpenMI/GraphBP/GraphBP/post_hoc_filtering/data/aurora_kinase_B_interactions.csv"
mols = read_aurora_kinase_b_interactions(file_path)

# Print first 5 mol objects
print(mols[:5])