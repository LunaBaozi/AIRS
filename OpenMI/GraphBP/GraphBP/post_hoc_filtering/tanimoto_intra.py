import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import matplotlib.pyplot as plt

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

# def get_fps_from_smiles(smiles_list):
#     fps = []
#     for smi in smiles_list:
#         mol = Chem.MolFromSmiles(smi)
#         if mol is not None:
#             fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
#             fps.append(fp)
#     return fps


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
    sdf_folder = "/vol/data/airs/AIRS/OpenMI/GraphBP/GraphBP/trained_model_reduced_dataset_100_epochs/gen_mols_epoch_99_10000/sdf"
    smiles_list, fps = get_mol_fps_and_smiles(sdf_folder)
    scores = compute_tanimoto_scores(fps, smiles_list)
    csv_path = "tanimoto_results_intra.csv"
    
    # smiles_list = ["C=C1[C@@H](C)O[C@H]2C[C@H]2N1C",
    #                "C=C1[C@@H](C)O[C@H]2C[C@H]2N1N"]
    # fps = get_fps_from_smiles(smiles_list)
    # scores = compute_tanimoto_scores(fps, smiles_list)
    # csv_path = "tanimoto_scores_selected.csv"
    with open(csv_path, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["smiles ligand 1", "smiles ligand 2", "tanimoto similarity"])
        for smi1, smi2, score in scores:
            writer.writerow([smi1, smi2, f"{score:.3f}"])
    print(f"Tanimoto scores saved to {csv_path}")

    

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
    plt.savefig("tanimoto_heatmap_intra.png", dpi=300)
    plt.close()
    print("Lower triangular heatmap saved to tanimoto_heatmap_intra.png")