# NOTE: there is a problem. OpenBabel cannot be installed by pip, 
# but it only works from conda. Have to try similarity calculation procedure
# suggested by Todeschini on colab maybe, or we need to create a venv

import os, copy, csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from cinfony import pybel, rdk, cdk


known_aurk_inhibitors_csv = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_inhibitors_copy.csv', header=0)
known_aurk_inhibitors = dict(zip(known_aurk_inhibitors_csv.inhibitor, known_aurk_inhibitors_csv.smiles))


def tanimoto_calc(smi1, smi2):
    
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)

    return s




# ------------------------------------------------------------------------------------------------------------------
# LET'S TRY TO GENERATE FINGERPRINTS WITH CINFONY

def similarity_cinfony():

    tanimoto = [['knowninhib_name', 'knowninhib_smiles', 'genmol_name', 'genmol_smiles', 'tanimoto_simil', 'tanimoto_simil_cinfony']]

    gen_mols_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/gen_mols_epoch_33_1k/content/GraphBP/GraphBP/trained_model/gen_mols_epoch_33'

    # for gen_mol_dir in os.listdir(gen_mols_dir):
    #     rel_path = os.path.join(gen_mols_dir, gen_mol_dir)
        
    for filename in os.listdir(gen_mols_dir):
        if 'sdf' in filename:
            sppl = Chem.SDMolSupplier(os.path.join(gen_mols_dir, filename)) 
                
            for mol in sppl:
                if mol is not None:  # some compounds cannot be loaded.
                    smiles_pattern = Chem.MolToSmiles(mol)
                        
                    for name, smile in known_aurk_inhibitors.items():
                        # print(compound)
                        smiles_test = str(smile)

                        smiles_list = [smiles_pattern, smiles_test]
                        mols = [pybel.readstring('smi', x) for x in smiles_list]
                        fps = [x.calcfp() for x in mols]
                        tanimoto_cinfony = fps[0] | fps[1]

                        # using basic smiles
                        simil = tanimoto_calc(smiles_pattern, smiles_test)
                        tanimoto.append([name, smiles_test, smiles_pattern, simil, tanimoto_cinfony])

                            
    return tanimoto #matches


hello = similarity_cinfony()

with open('./hello_output.csv', 'w', newline='') as csvfile1:
    writer = csv.writer(csvfile1)
    writer.writerows(hello)