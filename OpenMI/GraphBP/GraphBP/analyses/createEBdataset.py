import os 
import shutil 
from rdkit import Chem
from config_aurora import conf



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


def addHs2molfromsmiles(smiles):  

    # Generate a 3D structure from smiles
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    writer = Chem.SDWriter('Hmol.sdf')  
    writer.write(hmol)
    writer.close()

    return hmol 



def createEquiBindDataset(aurora_type, known_binding_site=True):

    if aurora_type == 'A' and known_binding_site:
        ligands_path = conf[1]['gen_ligands_path']
        os.makedirs(os.path.dirname(conf[1]['gen_complexes']), exist_ok=True)

        for i in range(1, 200):  #(len(os.listdir(ligands_path))-2)//2
            os.makedirs(os.path.join(conf[1]['gen_complexes'], f"{conf[1]['pdb']}_{i}"))
            subfolder = os.path.join(conf[1]['gen_complexes'], f"{conf[1]['pdb']}_{i}")
            noH_mol = os.path.join(ligands_path, f"{i}.sdf")
            shutil.copy(conf[1]['pdb_file_path'], subfolder)
            shutil.copy(f"{ligands_path}/{i}.sdf", subfolder)
            
            sppl = Chem.SDMolSupplier(noH_mol)
            outname = noH_mol.replace(".sdf", ".txt")
            out_file = open(outname, "w")
            for mol in sppl:
                if mol is not None:  
                    smi = Chem.MolToSmiles(mol)
                    name = mol.GetProp("_Name")
                    out_file.write(f"{smi}\t{name}\n")
                else: break

            out_file.close()
            if smi is not None:
                mol_addedHs = addHs2molfromsmiles(smi)
                shutil.move('Hmol.sdf', subfolder)

                os.rename(f"{subfolder}/{conf[1]['rec_name']}.pdb", f"{subfolder}/{conf[1]['pdb']}_{i}_protein.pdb")
                os.rename(f"{subfolder}/Hmol.sdf", f"{subfolder}/{conf[1]['pdb']}_{i}_ligand.sdf")

    elif aurora_type == 'B' and known_binding_site:
        ligands_path = conf[2]['gen_ligands_path']
        os.makedirs(os.path.dirname(conf[2]['gen_complexes']), exist_ok=True)

        for i in range(1, 200):
            os.makedirs(os.path.join(conf[2]['gen_complexes'], f"{conf[2]['pdb']}_{i}"))
            subfolder = os.path.join(conf[2]['gen_complexes'], f"{conf[2]['pdb']}_{i}")
            noH_mol = os.path.join(ligands_path, f"{i}.sdf")
            shutil.copy(conf[2]['pdb_file_path'], subfolder)
            shutil.copy(f"{ligands_path}/{i}.sdf", subfolder)
            
            sppl = Chem.SDMolSupplier(noH_mol)
            outname = noH_mol.replace(".sdf", ".txt")
            out_file = open(outname, "w")
            for mol in sppl:
                if mol is not None:  
                    smi = Chem.MolToSmiles(mol)
                    name = mol.GetProp("_Name")
                    out_file.write(f"{smi}\t{name}\n")
                else: break

            out_file.close()
            if smi is not None:
                mol_addedHs = addHs2molfromsmiles(smi)
                shutil.move('Hmol.sdf', subfolder)

                os.rename(f"{subfolder}/{conf[2]['rec_name']}.pdb", f"{subfolder}/{conf[2]['pdb']}_{i}_protein.pdb")
                os.rename(f"{subfolder}/Hmol.sdf", f"{subfolder}/{conf[2]['pdb']}_{i}_ligand.sdf")

    elif aurora_type == 'A' and not known_binding_site:
        ligands_path = conf[3]['gen_ligands_path']
        os.makedirs(os.path.dirname(conf[3]['gen_complexes']), exist_ok=True)

        for i in range(1, 200):
            os.makedirs(os.path.join(conf[3]['gen_complexes'], f"{conf[3]['pdb']}_{i}"))
            subfolder = os.path.join(conf[3]['gen_complexes'], f"{conf[3]['pdb']}_{i}")
            noH_mol = os.path.join(ligands_path, f"{i}.sdf")
            shutil.copy(conf[3]['pdb_file_path'], subfolder)
            shutil.copy(f"{ligands_path}/{i}.sdf", subfolder)
            
            sppl = Chem.SDMolSupplier(noH_mol)
            outname = noH_mol.replace(".sdf", ".txt")
            out_file = open(outname, "w")
            for mol in sppl:
                if mol is not None:  
                    smi = Chem.MolToSmiles(mol)
                    name = mol.GetProp("_Name")
                    out_file.write(f"{smi}\t{name}\n")
                else: break

            out_file.close()
            if smi is not None:
                mol_addedHs = addHs2molfromsmiles(smi)
                shutil.move('Hmol.sdf', subfolder)

                os.rename(f"{subfolder}/{conf[3]['rec_name']}.pdb", f"{subfolder}/{conf[3]['pdb']}_{i}_protein.pdb")
                os.rename(f"{subfolder}/Hmol.sdf", f"{subfolder}/{conf[3]['pdb']}_{i}_ligand.sdf")

    elif aurora_type == 'B' and not known_binding_site:
        ligands_path = conf[4]['gen_ligands_path']
        os.makedirs(os.path.dirname(conf[4]['gen_complexes']), exist_ok=True)

        for i in range(1, 200):
            os.makedirs(os.path.join(conf[4]['gen_complexes'], f"{conf[4]['pdb']}_{i}"))
            subfolder = os.path.join(conf[4]['gen_complexes'], f"{conf[4]['pdb']}_{i}")
            noH_mol = os.path.join(ligands_path, f"{i}.sdf")
            shutil.copy(conf[4]['pdb_file_path'], subfolder)
            shutil.copy(f"{ligands_path}/{i}.sdf", subfolder)
            
            sppl = Chem.SDMolSupplier(noH_mol)
            outname = noH_mol.replace(".sdf", ".txt")
            out_file = open(outname, "w")
            for mol in sppl:
                if mol is not None:  
                    smi = Chem.MolToSmiles(mol)
                    name = mol.GetProp("_Name")
                    out_file.write(f"{smi}\t{name}\n")
                else: break

            out_file.close()
            if smi is not None:
                mol_addedHs = addHs2molfromsmiles(smi)
                shutil.move('Hmol.sdf', subfolder)

                os.rename(f"{subfolder}/{conf[4]['rec_name']}.pdb", f"{subfolder}/{conf[4]['pdb']}_{i}_protein.pdb")
                os.rename(f"{subfolder}/Hmol.sdf", f"{subfolder}/{conf[4]['pdb']}_{i}_ligand.sdf")

    else:
        raise Exception('aurora kinase type or binding site flag not recognized.')

            
createEquiBindDataset(aurora_type='A', known_binding_site=True)
createEquiBindDataset(aurora_type='B', known_binding_site=True)
createEquiBindDataset(aurora_type='A', known_binding_site=False)
createEquiBindDataset(aurora_type='B', known_binding_site=False)
