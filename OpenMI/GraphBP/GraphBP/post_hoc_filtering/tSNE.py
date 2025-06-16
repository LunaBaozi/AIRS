import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


df = pd.read_csv('results_epoch_99_mols_1000_bs_False/molecule_scores_99_1000.csv')

# Generate Morgan fingerprints (as numpy arrays)
def smiles_to_fp(smiles, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = generator.GetFingerprint(mol)
    arr = np.zeros((n_bits,), dtype=int)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fps = df['smiles'].apply(smiles_to_fp)
fps = fps.dropna()
fps_matrix = np.stack(fps.values)

tsne = TSNE(n_components=2, random_state=42, perplexity=30)
tsne_results = tsne.fit_transform(fps_matrix)

# plt.figure(figsize=(8,6))
# plt.scatter(tsne_results[:,0], tsne_results[:,1], alpha=0.7, s=30)
# plt.xlabel('t-SNE 1')
# plt.ylabel('t-SNE 2')
# plt.title('t-SNE of Chemical Space')
# plt.tight_layout()
# plt.show()


plt.figure(figsize=(8,6))
plt.scatter(tsne_results[:,0], tsne_results[:,1], c=df.loc[fps.index, 'SA_score'], cmap='viridis', alpha=0.7, s=30)
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.title('t-SNE of Chemical Space (colored by SA_score)')
plt.colorbar(label='SA_score')
plt.tight_layout()
plt.savefig("tsne.png")
# plt.show()