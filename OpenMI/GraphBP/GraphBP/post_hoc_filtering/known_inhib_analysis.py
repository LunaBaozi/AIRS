import os, sys
import csv
import pandas as pd
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import Lipinski, Descriptors, Crippen

from SA_Score import sascorer 
from NP_Score import npscorer
# from syba.syba import SybaClassifier

from scripts import scscorer_standalone 
from matplotlib.lines import Line2D

# Initialize SCScorer
scscorer = scscorer_standalone.SCScorer()
scscorer.restore()

# Initialize SybaClassifier
# syba = SybaClassifier()
# syba.fitDefaultScore()

def read_aurora_kinase_b_interactions(filepath, smiles_col='smiles'):
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




import matplotlib.pyplot as plt

def plot_sa_scscore(csv_path):
    df = pd.read_csv('data/aurkb_scores.csv')
    df['SA_score'] = pd.to_numeric(df['SA_score'], errors='coerce')
    df['SCScore'] = pd.to_numeric(df['SCScore'], errors='coerce')
    df = df.dropna(subset=['SA_score', 'SCScore'])

    colors = df['passed'].apply(lambda x: 'blue' if x == 'Yes' else 'red')
    # # print(colors)
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
    plt.show()
    plt.savefig('data/aurkb_synth_scores_scatterplot.png', dpi=300)
    plt.close()


# Example usage:
# df, _ = read_aurora_kinase_b_interactions('data/aurora_kinase_B_interactions.csv', smiles_col='smiles')
# calculate_all_scores_and_save(df=df, smiles_col='smiles', csv_path='data/aurkb_scores.csv')
plot_sa_scscore(csv_path='data/aurkb_scores.csv')