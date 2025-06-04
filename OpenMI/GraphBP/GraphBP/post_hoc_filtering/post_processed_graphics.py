import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

results_folder = "results_bs_True"
graphics_folder = os.path.join(results_folder, 'graphics')
os.makedirs(graphics_folder, exist_ok=True)

# 1. Histogram: SMILES length
df_smiles = pd.read_csv(os.path.join(results_folder, "top100_synthesizable_len_smiles.csv"))
plt.figure(figsize=(6,4))
plt.hist(df_smiles['len_smiles'], bins=20, color='skyblue', edgecolor='black')
plt.xlabel('SMILES Length')
plt.ylabel('Count')
plt.title('Histogram of SMILES Length')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_len_smiles.png"))
plt.close()

# 2. Histogram: NP Score
df_np = pd.read_csv(os.path.join(results_folder, "top100_synthesizable_NP_score.csv"))
plt.figure(figsize=(6,4))
plt.hist(df_np['NP_score'], bins=20, color='lightgreen', edgecolor='black')
plt.xlabel('NP Score')
plt.ylabel('Count')
plt.title('Histogram of NP Score')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_NP_score.png"))
plt.close()

# 3. Histogram: SA Score
df_sa = pd.read_csv(os.path.join(results_folder, "top100_synthesizable_SA_score.csv"))
plt.figure(figsize=(6,4))
plt.hist(df_sa['SA_score'], bins=20, color='salmon', edgecolor='black')
plt.xlabel('SA Score')
plt.ylabel('Count')
plt.title('Histogram of SA Score')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_SA_score.png"))
plt.close()

# 4. Histogram: SCScore
df_sc = pd.read_csv(os.path.join(results_folder, "top100_synthesizable_SCScore.csv"))
plt.figure(figsize=(6,4))
plt.hist(df_sc['SCScore'], bins=20, color='plum', edgecolor='black')
plt.xlabel('SCScore')
plt.ylabel('Count')
plt.title('Histogram of SCScore')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_SCScore.png"))
plt.close()

# 5. Histogram: Sybe Score
df_syba = pd.read_csv(os.path.join(results_folder, "top100_synthesizable_Syba_score.csv"))
plt.figure(figsize=(6,4))
plt.hist(df_syba['Syba_score'], bins=20, color='gold', edgecolor='black')
plt.xlabel('Syba Score')
plt.ylabel('Count')
plt.title('Histogram of Syba Score')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_Syba_score.png"))
plt.close()

# 6. Heatmap: Tanimoto pairs inter
df_inter = pd.read_csv(os.path.join(results_folder, "top100_tanimoto_pairs_inter.csv"))
df_inter_agg = df_inter.groupby(['mol_1', 'mol_2'], as_index=False)['tanimoto'].mean()
matrix_inter = df_inter_agg.pivot(index="mol_1", columns="mol_2", values="tanimoto")
plt.figure(figsize=(8,6))
sns.heatmap(matrix_inter, cmap='viridis', annot=False)
plt.title('Tanimoto Similarity (Inter)')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "heatmap_tanimoto_inter.png"))
plt.close()

# Alternative visualization: Histogram of Tanimoto values (Inter)
plt.figure(figsize=(6,4))
plt.hist(df_inter_agg['tanimoto'], bins=20, color='deepskyblue', edgecolor='black')
plt.xlabel('Tanimoto Similarity')
plt.ylabel('Count')
plt.title('Histogram of Tanimoto Similarity (Inter)')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_tanimoto_inter.png"))
plt.close()

# 7. Heatmap: Tanimoto pairs intra
df_intra = pd.read_csv(os.path.join(results_folder, "top100_tanimoto_pairs_intra.csv"))
df_intra_agg = df_intra.groupby(['mol_1', 'mol_2'], as_index=False)['tanimoto'].mean()
matrix_intra = df_intra_agg.pivot(index="mol_1", columns="mol_2", values="tanimoto")
plt.figure(figsize=(8,6))
sns.heatmap(matrix_intra, cmap='viridis', annot=False)
plt.title('Tanimoto Similarity (Intra)')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "heatmap_tanimoto_intra.png"))
plt.close()

# Alternative visualization: Histogram of Tanimoto values (Intra)
plt.figure(figsize=(6,4))
plt.hist(df_intra_agg['tanimoto'], bins=20, color='mediumorchid', edgecolor='black')
plt.xlabel('Tanimoto Similarity')
plt.ylabel('Count')
plt.title('Histogram of Tanimoto Similarity (Intra)')
plt.tight_layout()
plt.savefig(os.path.join(graphics_folder, "hist_tanimoto_intra.png"))
plt.close()
