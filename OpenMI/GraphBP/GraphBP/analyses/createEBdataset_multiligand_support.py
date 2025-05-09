# BEFORE DOING THIS WE SHOULD PULL FROM EQUIBIND (since multiligand support is a new feature)
# so for the time being we still use the old method for generating the EB dataset

import os 
import shutil 
from rdkit import Chem


aur_protein_pdb = '4af3A'
rec_name = '4af3_A_rec'
aurkb_file = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
gen_ligands_path = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/gen_mols_epoch_33_1k/content/GraphBP/GraphBP/trained_model/gen_mols_epoch_33/'
aurkb_gen_complex = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex_1k'
write_dir = '/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_gen_complex_1k/molsHs_added'



def converter(file_name):
    sppl = Chem.SDMolSupplier(file_name)
    outname = file_name.replace(".sdf", ".txt")
    out_file = open(outname, "w")
    for mol in sppl:
        if mol is not None:  # some compounds cannot be loaded.
            smi = Chem.MolToSmiles(mol)
            name = mol.GetProp("_Name")
            out_file.write(f"{smi}\t{name}\n")

            out_file.close()
            return smi
    return None
        # else: break
    # out_file.close()
    # return smi


def addHs2molfromsmiles(smiles):  

    # Generate a 3D structure from smiles
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter('Hmol.sdf')  #(os.path.join(f'{write_dir}', 'Hmol.sdf'))
    writer.write(hmol)
    writer.close()

    return hmol 


def createEquiBind_dataset():

    os.makedirs(os.path.dirname(write_dir), exist_ok=True)

    for i in range(1, 1001):   # NEED TO MAKE IT UNIVERSAL
        os.makedirs(os.path.join(aurkb_gen_complex, f'{aur_protein_pdb}{i}'))

        subfolder = os.path.join(aurkb_gen_complex, f'{aur_protein_pdb}{i}')
        noH_mol = os.path.join(gen_ligands_path, f'{i}.sdf')

        shutil.copy(aurkb_file, subfolder)
        shutil.copy(f'{gen_ligands_path}/{i}.sdf', subfolder)

        # Returns
        sppl = Chem.SDMolSupplier(noH_mol)
        outname = noH_mol.replace(".sdf", ".txt")
        out_file = open(outname, "w")
        for mol in sppl:
            if mol is not None:  # some compounds cannot be loaded.
                smi = Chem.MolToSmiles(mol)
                name = mol.GetProp("_Name")
                out_file.write(f"{smi}\t{name}\n")
            else: break

        out_file.close()
        # return smi
        # smi = converter(noH_mol)
        if smi is not None:
            mol_addedHs = addHs2molfromsmiles(smi)
            shutil.move('Hmol.sdf', subfolder)

            os.rename(f'{subfolder}/{rec_name}.pdb', f'{subfolder}/{aur_protein_pdb}{i}_protein.pdb')
            os.rename(f'{subfolder}/Hmol.sdf', f'{subfolder}/{aur_protein_pdb}{i}_ligand.sdf')

        
createEquiBind_dataset()