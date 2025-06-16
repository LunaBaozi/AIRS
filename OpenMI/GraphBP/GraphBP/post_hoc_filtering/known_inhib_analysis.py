import os, sys
import csv
import pandas as pd
import numpy as np
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFingerprintGenerator
from rdkit.Chem import Lipinski, Descriptors, Crippen

from SA_Score import sascorer 
from NP_Score import npscorer

from scripts import scscorer_standalone 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Initialize SCScorer
scscorer = scscorer_standalone.SCScorer()
scscorer.restore()


def read_aurora_kinase_interactions(filepath, smiles_col='smiles'):
    data = []
    with open(filepath, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', skipinitialspace=True)
        for row in reader:
            cleaned_row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            data.append(cleaned_row)
    df = pd.DataFrame(data)
    atom_counts = []
    for smi in df[smiles_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mol = Chem.AddHs(mol)
            atom_counts.append(mol.GetNumAtoms())
    if atom_counts:
        print(f"Minimum number of atoms: {min(atom_counts)}")
        print(f"Maximum number of atoms: {max(atom_counts)}")
    else:
        print("No valid molecules found.")
    return df, atom_counts

def log_partition_coefficient(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{smiles} returns a None molecule")
    return Crippen.MolLogP(mol)

def lipinski_trial_on_smiles(smiles_list):
    passed = []
    failed = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failed.append({'smiles': smiles, 'reason': 'Invalid SMILES'})
            continue
        num_hdonors = Lipinski.NumHDonors(mol)
        num_hacceptors = Lipinski.NumHAcceptors(mol)
        mol_weight = Descriptors.MolWt(mol)
        mol_logp = Crippen.MolLogP(mol)
        reasons = []
        if num_hdonors > 5:
            reasons.append(f'Over 5 H-bond donors, found {num_hdonors}')
        if num_hacceptors > 10:
            reasons.append(f'Over 10 H-bond acceptors, found {num_hacceptors}')
        if mol_weight >= 500:
            reasons.append(f'Molecular weight over 500, calculated {mol_weight}')
        if mol_logp >= 5:
            reasons.append(f'Log partition coefficient over 5, calculated {mol_logp}')
        if reasons:
            failed.append({'smiles': smiles, 'reason': '; '.join(reasons)})
        else:
            passed.append({'smiles': smiles})
    return passed, failed

def calculate_sa_score(mol):
    return sascorer.calculateScore(mol) 

def calculate_np_score(mol):
    fscore = npscorer.readNPModel()
    score = npscorer.scoreMol(mol, fscore)
    confidence = npscorer.scoreMolWConfidence(mol, fscore)
    return score, confidence

def calculate_all_scores_and_save(df, smiles_col, csv_path):
    """
    Calculates all scores and Lipinski results for all SMILES in a dataframe and writes them to a single CSV file.
    Columns: index, smiles, SA_score, SCScore, NP_score, passed, failed
    """
    # Initialize SCScorer and Syba if needed
    scscorer = scscorer_standalone.SCScorer()
    scscorer.restore()
    # syba = SybaClassifier()
    # syba.fitDefaultScore()

    # Prepare Lipinski results
    smiles_list = df[smiles_col].tolist()
    passed, failed = lipinski_trial_on_smiles(smiles_list)
    lipinski_dict = {}
    for entry in passed:
        lipinski_dict[entry['smiles']] = {'passed': 'Yes', 'failed': ''}
    for entry in failed:
        lipinski_dict[entry['smiles']] = {'passed': '', 'failed': entry['reason']}

    # Calculate scores
    data = []
    for idx, row in df.iterrows():
        ligand = row['ligand']
        smi = row[smiles_col]
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            sa_score = calculate_sa_score(mol)
            np_score, _ = calculate_np_score(mol)
            (_, sc_score) = scscorer.get_score_from_smi(smi)
            # syba_score = syba.predict(smi)  # Not included in output as per request
            lipinski = lipinski_dict.get(smi, {'passed': '', 'failed': 'Not checked'})
            data.append({
                'index': idx,
                'ligand': ligand,
                'smiles': smi,
                'SA_score': sa_score,
                'SCScore': sc_score,
                'NP_score': np_score,
                'passed': lipinski['passed'],
                'failed': lipinski['failed']
            })
    result_df = pd.DataFrame(data)
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    result_df.to_csv(csv_path, index=False)




def plot_sa_scscore(csv_path):
    df = pd.read_csv(csv_path)
    plot_name = str(csv_path.replace('_scores.csv', '_scatterplot.png'))
    print(plot_name)
    df['SA_score'] = pd.to_numeric(df['SA_score'], errors='coerce')
    df['SCScore'] = pd.to_numeric(df['SCScore'], errors='coerce')
    df = df.dropna(subset=['SA_score', 'SCScore'])

    colors = df['passed'].apply(lambda x: 'blue' if x == 'Yes' else 'red')
    plt.figure(figsize=(8,6))
    plt.scatter(df['SA_score'], df['SCScore'], c=colors, alpha=0.7, edgecolor='k')
    plt.xlim(df['SA_score'].min()-0.5, df['SA_score'].max()+0.5)
    plt.ylim(df['SCScore'].min()-0.5, df['SCScore'].max()+0.5)
    plt.xlabel('SA_score')
    plt.ylabel('SCScore')
    plt.title('SA_score vs SCScore (Lipinski color-coded)')
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Lipinski Passed', markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Lipinski Failed', markerfacecolor='red', markersize=8)
    ]
    plt.legend(handles=legend_elements)
    plt.tight_layout()

    # Find the point with both minimum SA_score and minimum SCScore
    min_sa = df['SA_score'].min()
    min_sc = df['SCScore'].min()
    mask = (df['SA_score'] == min_sa) & (df['SCScore'] == min_sc)
    if mask.any():
        row = df[mask].iloc[0]
        plt.annotate(
            f"{row['ligand']} ({row['SA_score']:.2f}, {row['SCScore']:.2f})",
            (row['SA_score'], row['SCScore']),
            textcoords="offset points",
            xytext=(10,10),
            ha='left',
            fontsize=9,
            color='black',
            bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.3)
        )

    plt.show()
    plt.savefig(plot_name, dpi=300)
    plt.close()



def get_mol_fps_and_smiles(aurora_file):
    # df = pd.read_csv(aurora_file)
    fps = []
    smiles_list = df['smiles']
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            fp = generator.GetFingerprint(mol)
            # smiles = Chem.MolToSmiles(mol)
            fps.append(fp)
            # smiles_list.append(smiles)
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

# Generate lower triangular heat map
def generate_lower_triangular_heatmap(fps, filename):
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
    plt.savefig(os.path.join('data', filename))
    plt.close()
    print("Lower triangular heatmap saved to tanimoto_heatmap_intra.png")







# Example usage:
df, _ = read_aurora_kinase_interactions('data/aurora_kinase_B_interactions.csv', smiles_col='smiles')
# calculate_all_scores_and_save(df=df, smiles_col='smiles', csv_path='data/aurkb_scores.csv')
# plot_sa_scscore(csv_path='data/aurkb_scores.csv')
smiles_list, fps = get_mol_fps_and_smiles(df)
scores = compute_tanimoto_scores(fps, smiles_list)

df = pd.DataFrame(scores, columns=["mol_1", "mol_2", "tanimoto"])
# df["tanimoto"] = df["tanimoto"].round(3)
output_csv = os.path.join('data', "aurkb_tanimoto_results_intra.csv")
df = df.sort_values(by="tanimoto", ascending=False)
df.to_csv(output_csv, index=False)

print(f"Tanimoto scores saved to {output_csv}")
generate_lower_triangular_heatmap(fps, filename='tanimoto_heatmap_aurkb.png')



df, _ = read_aurora_kinase_interactions('data/aurora_kinase_A_interactions.csv', smiles_col='smiles')
# calculate_all_scores_and_save(df=df, smiles_col='smiles', csv_path='data/aurka_scores.csv')
# plot_sa_scscore(csv_path='data/aurka_scores.csv')
smiles_list, fps = get_mol_fps_and_smiles(df)
scores = compute_tanimoto_scores(fps, smiles_list)

df = pd.DataFrame(scores, columns=["mol_1", "mol_2", "tanimoto"])
# df["tanimoto"] = df["tanimoto"].round(3)
output_csv = os.path.join('data', "aurka_tanimoto_results_intra.csv")
df = df.sort_values(by="tanimoto", ascending=False)
df.to_csv(output_csv, index=False)
print(f"Tanimoto scores saved to {output_csv}")
generate_lower_triangular_heatmap(fps, filename='tanimoto_heatmap_aurka.png')