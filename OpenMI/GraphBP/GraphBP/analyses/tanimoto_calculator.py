import os, copy, csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from openbabel import openbabel
from cinfony import pybel   #, cdk  #rdk, 

# AURORA KINASE B KNOWN INHIBITORS following https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=1937
# interessante la presenza del fluoro in molti inibitori
# interessante anche come scompaia la @ nella forma canonica delle smiles

def addHs2molfromsmiles(smiles, filename='Hmol.sdf', output_type=None):  

    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter(filename)  
    writer.write(hmol)
    writer.close()

    if output_type == 'smiles':
        return Chem.MolToSmiles(hmol)
    elif output_type == 'mol':
        return hmol
    #elif output_type == 'sdf'
    return  



# print(known_aurkb_inhibitors.iloc[:, 0])    # target id
# print(known_aurkb_inhibitors.iloc[:, 1])    # target uniprot
# print(known_aurkb_inhibitors.iloc[:, 2])    # target species
# print(known_aurkb_inhibitors.iloc[:, 3])    # ligand
# print(known_aurkb_inhibitors.iloc[:, 4])    # ligand id
# print(known_aurkb_inhibitors.iloc[:, 5])    # ligand species
# print(known_aurkb_inhibitors.iloc[:, 6])    # ligand pubchem cid
# print(known_aurkb_inhibitors.iloc[:, 7])    # smiles

aurkA_inhibitors = pd.read_csv('./aurorakinaseAinteractions.csv')
aurkB_inhibitors = pd.read_csv('./aurorakinaseBinteractions.csv')
aurkA_inhibitors_dict = dict(zip(aurkA_inhibitors.iloc[:, 3], aurkA_inhibitors.iloc[:, 7]))
aurkB_inhibitors_dict = dict(zip(aurkB_inhibitors.iloc[:, 3], aurkB_inhibitors.iloc[:, 7]))


def tanimoto_calc(smi1, smi2):
    
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)

    return s



def compare_mol_smiles(gen_mols_dir, inhibitors_dict):

    tanimoto_similarity = [['knowninhib_name', 'knowninhib_smiles', 'genmol_name', 'genmol_smiles', 'simil', 'cinfony-simil']]

    for gen_complex_dir in os.listdir(gen_mols_dir):
        rel_path = os.path.join(gen_mols_dir, gen_complex_dir)
        
        for filename in os.listdir(rel_path):
            if 'ligand' in filename:
                sppl = Chem.SDMolSupplier(os.path.join(rel_path, filename)) 
                
                for mol in sppl:
                    if mol is not None:  
                        basic_smiles_pattern = Chem.MolToSmiles(mol)
                        canon_smiles_pattern = Chem.CanonSmiles(basic_smiles_pattern)
                        
                        for name, compound in inhibitors_dict.items():
                            basic_smiles_test = str(compound)
                            canon_smiles_test = Chem.CanonSmiles(compound)

                            # Using canon smiles
                            simil = tanimoto_calc(canon_smiles_pattern, canon_smiles_test)

                            # Using cinfony fingerprints
                            smiles_list = [basic_smiles_pattern, basic_smiles_test]
                            mols = [pybel.readstring('smi', x) for x in smiles_list]
                            fps = [x.calcfp() for x in mols]
                            tanimoto_cin = fps[0] | fps[1]

                            tanimoto_similarity.append([name, canon_smiles_test, filename, canon_smiles_pattern, simil, tanimoto_cin])

    return tanimoto_similarity   



aurkA_gen_mols = './results/AURKA_GEN_COMPLEX'
aurkB_gen_mols = './results/AURKB_GEN_COMPLEX'
aurkA_NOBS_gen_mols = './results/AURKA_NOBS_GEN_COMPLEX'
aurkB_NOBS_gen_mols = './results/AURKB_NOBS_GEN_COMPLEX'

tanimoto_canon_A = compare_mol_smiles(aurkA_gen_mols, aurkA_inhibitors_dict)
tanimoto_canon_B = compare_mol_smiles(aurkB_gen_mols, aurkB_inhibitors_dict)
tanimoto_canon_A_nobs = compare_mol_smiles(aurkA_NOBS_gen_mols, aurkA_inhibitors_dict)
tanimoto_canon_B_nobs = compare_mol_smiles(aurkB_NOBS_gen_mols, aurkB_inhibitors_dict)

with open('./results/tanimoto_simil_AURKA_4byi.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(tanimoto_canon_A)

with open('./results/tanimoto_simil_AURKB.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(tanimoto_canon_B)

with open('./results/tanimoto_simil_AURKA_NOBS.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(tanimoto_canon_A_nobs)

with open('./results/tanimoto_simil_AURKB_NOBS.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(tanimoto_canon_B_nobs)



# # INTERESSANTE! A e B hanno tanimoto coeff di solo 0.153
# A = 'CN1CCN(CC1)c1nc(Sc2ccc(cc2)NC(=O)C2CC2)nc(c1)Nc1[nH]nc(c1)C'
# B = 'CCS(=O)(=O)Nc1ccc2c(c1)C(=C(c1ccccc1)Nc1ccc(cc1)CN1CCCCC1)C(=O)N2'
# C = 'CCS(=O)(=O)Nc1ccc2c(c1)C(=C(c1ccccc1)Nc1ccc(cc1)CN1CCCCC1)C(=O)N2'

# print(tanimoto_calc(B, C))