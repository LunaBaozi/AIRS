import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


df = pd.read_csv('SMILES.csv')

mols = [Chem.MolFromSmiles(smi) for smi in df.SMILES]
hmols = [Chem.AddHs(m) for m in mols]

for mol  in hmols:
    AllChem.EmbedMolecule(mol,AllChem.ETKDG())
    print(AllChem.UFFOptimizeMolecule(mol,1000))

smiles = list(df.SMILES)
sid = list(df.SOURCE_ID)
libs = df[df.columns[0]]
writer = Chem.SDWriter('TEST.sdf')

for n in range(len(libs)):
    hmols[n].SetProp("_Library","%s"%libs[n])
    hmols[n].SetProp("_Name","%s"%sid[n])
    hmols[n].SetProp("_SourceID","%s"%sid[n])
    hmols[n].SetProp("_SMILES","%s"%smiles[n])
    writer.write(hmols[n])
writer.close()