from rdkit import Chem
from rdkit.Chem import Draw
import os

def sdf_folder_to_smiles_and_images(sdf_folder, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    smiles_list = []
    sdf_files = [f for f in os.listdir(sdf_folder) if f.lower().endswith('.sdf')]
    for sdf_file in sdf_files:
        sdf_path = os.path.join(sdf_folder, sdf_file)
        suppl = Chem.SDMolSupplier(sdf_path)
        for idx, mol in enumerate(suppl):
            if mol is None:
                continue
            # Compute 2D coordinates for better depiction
            Chem.rdDepictor.Compute2DCoords(mol)
            smiles = Chem.MolToSmiles(mol)
            smiles_list.append(smiles)
            print(sdf_file)
            print(smiles)
            # Use MolDraw2DCairo for higher quality images
            drawer = Draw.MolDraw2DCairo(500, 500)
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.addStereoAnnotation = True
            opts.legendFontSize = 20
            opts.bondLineWidth = 2.0
            drawer.DrawMolecule(mol, legend=smiles)
            drawer.FinishDrawing()
            img_bytes = drawer.GetDrawingText()
            img_path = os.path.join(output_dir, f"{os.path.splitext(sdf_file)[0]}_mol_{idx+1}.png")
            with open(img_path, "wb") as img_file:
                img_file.write(img_bytes)
            # break
        # break
    # Save SMILES to a text file
    smiles_file = os.path.join(output_dir, "smiles.txt")
    with open(smiles_file, "w") as f:
        for smi, sdf in zip(smiles_list, sdf_files):
            f.write(f"{sdf}\t{smi}\n")
    print(f"Processed {len(smiles_list)} molecules from {len(sdf_files)} SDF files. Images and SMILES saved to {output_dir}")

# Example usage:
sdf_folder_to_smiles_and_images("4af3/experiment_epoch_0_mols_0_bs_0_pdbid_4af3/ligands/", 
                                "4af3/experiment_epoch_0_mols_0_bs_0_pdbid_4af3/output_images")