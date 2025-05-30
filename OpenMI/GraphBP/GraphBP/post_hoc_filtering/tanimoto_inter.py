import os
import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
import seaborn as sns

# Paths to input files
sdf_path = "/vol/data/airs/AIRS/OpenMI/GraphBP/GraphBP/trained_model_reduced_dataset_100_epochs/gen_mols_epoch_99_10000/sdf"
csv_path = "data/aurora_kinase_B_interactions.csv"
output_csv = "tanimoto_results_inter.csv"

# def load_mols_from_sdf(sdf_file):
#     suppl = Chem.SDMolSupplier(sdf_file)
#     return [mol for mol in suppl if mol is not None]

# def load_all_mols_from_sdf_folder(folder):
#     mols = []
#     for f in os.listdir(folder):
#         if f.endswith('.sdf'):
#             mols.extend(load_mols_from_sdf(os.path.join(folder, f)))
#     return mols

def load_all_mols_from_sdf_folder(sdf_folder):
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

def read_aurora_kinase_b_interactions(filepath, smiles_col='smiles'):
    data = []
    # fps = []
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            cleaned_row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            data.append(cleaned_row)
    df = pd.DataFrame(data)
    mols = [Chem.MolFromSmiles(sm) for sm in df[smiles_col] if Chem.MolFromSmiles(sm) is not None]
    fp = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols]
    # fps.append(fp)
    # smiles = [sm for sm in df[smiles_col] if Chem.MolFromSmiles(sm) is not None]
    return mols, fp  #, smiles

def tanimoto_similarity_matrix(fps1, fps2):
    sim_matrix = []
    for fp1 in fps1:
        sims = DataStructs.BulkTanimotoSimilarity(fp1, fps2)
        sim_matrix.append(sims)
    return sim_matrix

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Load all molecules from all SDF files in the directory
    smiles_sdf, fps_sdf = load_all_mols_from_sdf_folder(sdf_path)
    if not fps_sdf:
        raise FileNotFoundError("No molecules found in any SDF file in the specified directory.")

    mols_csv, fps_csv = read_aurora_kinase_b_interactions(csv_path)

    sim_matrix = tanimoto_similarity_matrix(fps_sdf, fps_csv)

    # Prepare data for CSV
    rows = []
    for i, smiles_sdf_i in enumerate(smiles_sdf):
        for j, mol_csv in enumerate(mols_csv):
            rows.append({
                "sdf smiles": smiles_sdf_i,
                "csv smiles": Chem.MolToSmiles(mol_csv),
                "tanimoto": sim_matrix[i][j]
            })

    # Save to CSV
    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    # Plot full similarity heatmap
    sim_arr = np.array(sim_matrix)

    plt.figure(figsize=(12, 10))
    sns.heatmap(sim_arr, cmap="viridis", cbar_kws={"label": "Tanimoto similarity"})
    plt.title("Tanimoto Similarity Heatmap")
    plt.xlabel("CSV molecules")
    plt.ylabel("SDF molecules")
    plt.tight_layout()
    plt.savefig("tanimoto_heatmap.png")
    plt.close()
    print("Heatmap saved to tanimoto_heatmap_inter.png")
