import os
import csv
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator
import numpy as np
import seaborn as sns

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def load_all_mols_from_sdf_folder(sdf_folder):
    fps = []
    filenames = []
    smiles_list = []
    for fname in os.listdir(sdf_folder):
        if fname.endswith('.sdf'):
            sdf_path = os.path.join(sdf_folder, fname)
            suppl = Chem.SDMolSupplier(sdf_path)
            for mol in suppl:
                if mol is not None:
                    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                    fp = generator.GetFingerprint(mol)
                    smiles = Chem.MolToSmiles(mol)
                    fps.append(fp)
                    smiles_list.append(smiles)
                    filenames.append(fname)
                    # break
    return smiles_list, fps, filenames

def read_aurora_kinase_b_interactions(filepath, smiles_col='smiles'):
    data = []
    # fps = []
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            cleaned_row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            data.append(cleaned_row)
    df = pd.DataFrame(data)
    smiles = df['smiles']
    mols = [Chem.MolFromSmiles(sm) for sm in df[smiles_col] if Chem.MolFromSmiles(sm) is not None]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = [generator.GetFingerprint(mol) for mol in mols]
    # fps.append(fp)
    # smiles = [sm for sm in df[smiles_col] if Chem.MolFromSmiles(sm) is not None]
    return mols, fp, smiles

def tanimoto_similarity_matrix(fps1, fps2):
    sim_matrix = []
    for fp1 in fps1:
        sims = DataStructs.BulkTanimotoSimilarity(fp1, fps2)
        sim_matrix.append(sims)
    return sim_matrix

if __name__ == "__main__":
    # import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Calculate Tanimoto similarity between generated and known molecules.")
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--num_gen', type=int, required=True, help='Number of molecules generated')
    parser.add_argument('--known_binding_site', type=str2bool, required=True, help='Known binding site (True/False)')
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site

    # Paths to input files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")
    results_folder = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}")
    csv_path = os.path.join(script_dir, f"data/aurora_kinase_B_interactions.csv")
    output_csv = "tanimoto_results_inter.csv"

    # Load all molecules from all SDF files in the directory
    smiles_sdf, fps_sdf, filenames = load_all_mols_from_sdf_folder(sdf_folder)
    if not fps_sdf:
        raise FileNotFoundError("No molecules found in any SDF file in the specified directory.")

    mols_csv, fps_csv, smiles = read_aurora_kinase_b_interactions(csv_path)

    sim_matrix = tanimoto_similarity_matrix(fps_sdf, fps_csv)
    # Prepare data for CSV
    rows = []
    for i, smiles_sdf_i in enumerate(smiles_sdf):
        for j, smile in enumerate(smiles):
            rows.append({
                "filename": filenames[i],  # add filename column
                "smiles": smiles_sdf_i,  # generated molecules
                "mol_2": smile, #Chem.MolToSmiles(mol_csv),  # known inhibitors
                "tanimoto": sim_matrix[i][j]
            })

    # Save to CSV
    df_out = pd.DataFrame(rows)
    df_out.to_csv(os.path.join(results_folder, output_csv), index=False)
    print(f"Results saved to {output_csv}")

    # Plot full similarity heatmap
    # sim_arr = np.array(sim_matrix)

    # plt.figure(figsize=(12, 10))
    # sns.heatmap(sim_arr, cmap="viridis", cbar_kws={"label": "Tanimoto similarity"})
    # plt.title("Tanimoto Similarity Heatmap")
    # plt.xlabel("CSV molecules")
    # plt.ylabel("SDF molecules")
    # plt.tight_layout()
    # plt.savefig(os.path.join(results_folder, "tanimoto_heatmap_inter.png"))
    # plt.close()
    # print("Heatmap saved to tanimoto_heatmap_inter.png")
