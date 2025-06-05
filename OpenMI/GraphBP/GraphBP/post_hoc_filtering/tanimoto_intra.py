import os
import argparse
# import csv
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



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
    parser = argparse.ArgumentParser(description="Compute intra-set Tanimoto similarity for generated molecules.")
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--num_gen', type=int, required=True, help='Number of molecules generated')
    parser.add_argument('--known_binding_site', type=str2bool, required=True, help='Known binding site (True/False)')
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site

    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")    
    results_folder = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}")

    smiles_list, fps = get_mol_fps_and_smiles(sdf_folder)
    scores = compute_tanimoto_scores(fps, smiles_list)

    df = pd.DataFrame(scores, columns=["mol_1", "mol_2", "tanimoto"])
    # df["tanimoto"] = df["tanimoto"].round(3)
    os.makedirs(results_folder, exist_ok=True)
    output_csv = os.path.join(results_folder, "tanimoto_results_intra.csv")
    df.to_csv(output_csv, index=False)
    print(f"Tanimoto scores saved to {output_csv}")

    # Generate lower triangular heat map
    # mat = compute_tanimoto_matrix(fps)
    # print(mat.max())
    # mask = np.tril(np.ones_like(mat, dtype=bool))
    # plt.figure(figsize=(8, 8))
    # plt.imshow(np.where(mask, mat, np.nan), cmap='viridis', vmin=0, vmax=1)
    # plt.colorbar(label='Tanimoto Similarity')
    # plt.title('Lower Triangular Tanimoto Similarity Heatmap')
    # plt.xlabel('Ligand Index')
    # plt.ylabel('Ligand Index')
    # plt.tight_layout()
    # plt.savefig(os.path.join(results_folder, "tanimoto_heatmap_intra.png"))
    # plt.close()
    # print("Lower triangular heatmap saved to tanimoto_heatmap_intra.png")