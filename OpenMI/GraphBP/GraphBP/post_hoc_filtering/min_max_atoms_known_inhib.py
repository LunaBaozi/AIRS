import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def read_aurora_kinase_b_interactions(filepath, smiles_col='smiles'):
    data = []
    # fps = []
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            cleaned_row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            data.append(cleaned_row)
    df = pd.DataFrame(data)
    # mols = [Chem.MolFromSmiles(sm) for sm in df[smiles_col] if Chem.MolFromSmiles(sm) is not None]
    atom_counts = []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            atom_counts.append(mol.GetNumAtoms())

    if atom_counts:
        print(f"Minimum number of atoms: {min(atom_counts)}")
        print(f"Maximum number of atoms: {max(atom_counts)}")
    else:
        print("No valid molecules found.")
    return atom_counts

read_aurora_kinase_b_interactions('data/aurora_kinase_B_interactions.csv', smiles_col='smiles')
