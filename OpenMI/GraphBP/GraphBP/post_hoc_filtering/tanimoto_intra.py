import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import matplotlib.pyplot as plt

epoch = 99
num_gen = 1000
known_binding_site = True

def get_mol_fps_and_smiles(sdf_folder):
    fps = []
    smiles_list = []
    for fname in os.listdir(sdf_folder):
        if fname.endswith('.sdf'):
            sdf_path = os.path.join(sdf_folder, fname)
            suppl = Chem.SDMolSupplier(sdf_path)
            for mol in suppl:
                if mol is not None:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    smiles = Chem.MolToSmiles(mol)
                    fps.append(fp)
                    smiles_list.append(smiles)
    return smiles_list, fps

def compute_tanimoto_scores(fps, smiles_list):
    n = len(fps)
    scores = []
    for i in range(n):
        for j in range(i+1, n):
            score = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            scores.append((smiles_list[i], smiles_list[j], score))
    return scores

def compute_tanimoto_matrix(fps):
    n = len(fps)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1):
            mat[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    return mat

if __name__ == "__main__":
    sdf_folder = f"../trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{str(known_binding_site)}/sdf/"
    smiles_list, fps = get_mol_fps_and_smiles(sdf_folder)
    scores = compute_tanimoto_scores(fps, smiles_list)
    # csv_path = "results/tanimoto_results_intra.csv"
    results_folder = f"results_bs_{str(known_binding_site)}"
    os.makedirs(results_folder, exist_ok=True)
    output_csv = os.path.join(results_folder, "tanimoto_results_intra.csv")

    with open(output_csv, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["mol_1", "mol_2", "tanimoto"])
        for smi1, smi2, score in scores:
            writer.writerow([smi1, smi2, f"{score:.3f}"])
    print(f"Tanimoto scores saved to {output_csv}")

    

    # Generate lower triangular heat map
    mat = compute_tanimoto_matrix(fps)
    print(mat.max())
    mask = np.tril(np.ones_like(mat, dtype=bool))
    plt.figure(figsize=(8, 8))
    plt.imshow(np.where(mask, mat, np.nan), cmap='viridis', vmin=0, vmax=1)
    plt.colorbar(label='Tanimoto Similarity')
    plt.title('Lower Triangular Tanimoto Similarity Heatmap')
    plt.xlabel('Ligand Index')
    plt.ylabel('Ligand Index')
    plt.tight_layout()
    plt.savefig(os.path.join(results_folder, "tanimoto_heatmap_intra.png"))
    plt.close()
    print("Lower triangular heatmap saved to tanimoto_heatmap_intra.png")