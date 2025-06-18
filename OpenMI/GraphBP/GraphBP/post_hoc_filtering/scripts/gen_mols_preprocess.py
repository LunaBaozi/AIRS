import os
from rdkit import Chem

def load_mols_from_sdf_folder(folder_path):
    mols = []
    smiles = []
    filenames = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            mol = Chem.SDMolSupplier(os.path.join(folder_path, filename))[0]
            if mol is not None:
                mols.append(mol)
                smiles.append(Chem.MolToSmiles(mol))
                filenames.append(filename)
    return mols, smiles, filenames