import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

epoch = 99
num_gen = 10000

# Load the data
df = pd.read_csv(f"results/molecule_scores_{epoch}_{num_gen}.csv")

# List of score types and their labels
score_types = [
    ("len_smiles", "SMILES Length"),
    ("SA_score", "SA Score"),
    ("NP_score", "NP Score"),
    ("SCScore", "SC Score"),
    ("Syba_score", "Syba Score")
]



# SMILES
fig, ax = plt.subplots(figsize=(8, 5))

bin_centers = np.arange(df["len_smiles"].min(), df["len_smiles"].max()+1, 10)
bins = np.arange(df["len_smiles"].min()+0.75, df["len_smiles"].max()+1, 5.5)
ax.hist(df["len_smiles"].dropna(), bins=bins, color='green', edgecolor='black', rwidth=0.8)

ax.set_xlabel("SMILES length")
ax.set_ylabel("Count of Molecules")
ax.set_xticks(bin_centers)
ax.set_xlim([df["len_smiles"].min()+0.5, df["len_smiles"].max()+0.5])
ax.set_title("Histogram of SMILES Length")
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
fig.savefig(f'results/histogram_len_smiles_{epoch}_{num_gen}.png')


# SA_SCORE
fig, ax = plt.subplots(figsize=(8, 5))

bin_centers = np.arange(1, 11)
bins = np.arange(0.75, 10.25, 0.5)
ax.hist(df["SA_score"].dropna(), bins=bins, color='orange', edgecolor='black', rwidth=0.8)

ax.set_xlabel("SA_Score (centered bins 1-10)")
ax.set_ylabel("Count of Molecules")
ax.set_xticks(bin_centers)
ax.set_xlim([0.5, 10.5])
ax.set_title("Histogram of SA_Score (centered bins 1-10)")
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
fig.savefig(f'results/histogram_SA_score_{epoch}_{num_gen}.png')



# SC_SCORE
fig, ax = plt.subplots(figsize=(8, 5))

bin_centers = np.arange(1, 6)
bins = np.arange(0.75, 5.25, 0.5)
ax.hist(df["SCScore"].dropna(), bins=bins, color='blue', edgecolor='black', rwidth=0.8)

ax.set_xlabel("SCScore")
ax.set_ylabel("Count of Molecules")
ax.set_xticks(bin_centers)
ax.set_xlim([0.5, 5.5])
ax.set_title("Histogram of SCScore")
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
fig.savefig(f'results/histogram_SCScore_{epoch}_{num_gen}.png')



# NP_SCORE
fig, ax = plt.subplots(figsize=(8, 5))

bin_centers = np.arange(-5, 6)
bins = np.arange(-5.75, 5.25, 0.5)
ax.hist(df["NP_score"].dropna(), bins=bins, color='pink', edgecolor='black', rwidth=0.8)

ax.set_xlabel("NP_score")
ax.set_ylabel("Count of Molecules")
ax.set_xticks(bin_centers)
ax.set_xlim([-5.5, 5.5])
ax.set_title("Histogram of NP_score")
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
fig.savefig(f'results/histogram_NPscore_{epoch}_{num_gen}.png')
