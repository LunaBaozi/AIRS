import os
import argparse
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

from scripts.aurk_int_preprocess import read_aurora_kinase_interactions
from scripts.gen_mols_preprocess import load_mols_from_sdf_folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper for CADD pipeline targeted to Aurora protein kinases.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=True, help='Aurora kinase type (str, A, B)')
    args = parser.parse_args()

    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    epoch = args.epoch
    aurora = args.aurora

    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")
    known_inhib_file = os.path.join(script_dir, f"data/aurora_kinase_{aurora}_interactions.csv")
    results_dir = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}_aurora_{aurora}")
    output_png = os.path.join(results_dir, f"tSNE_{epoch}_{num_gen}_{known_binding_site}_{aurora}.png")
    
    os.makedirs(os.path.dirname(output_png), exist_ok=True)

    if epoch != 0:
        # Calculating scores for generated molecules
        mols, smiles, filenames, fps = load_mols_from_sdf_folder(sdf_folder)

    else:
        # Calculating scores for Aurora inhibitors
        mols, smiles, filenames, fps = read_aurora_kinase_interactions(known_inhib_file)
    
    synth_df = pd.read_csv(os.path.join(results_dir, f"synthesizability_scores_{epoch}_{num_gen}_{known_binding_site}_{aurora}.csv"))
    n_samples = len(fps)
    if n_samples <= 1:
        raise ValueError("At least 2 samples are required for t-SNE.")
    if n_samples <= 30:
        perplexity = max(1, n_samples - 1)
    else:
        perplexity = 30
    # Ensure perplexity is strictly less than n_samples
    if perplexity >= n_samples:
        perplexity = n_samples - 1
    print(perplexity, n_samples)
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    tsne_results = tsne.fit_transform(np.stack(fps))

    plt.figure(figsize=(8,6))
    plt.scatter(tsne_results[:,0], tsne_results[:,1], c=synth_df['SA_score'].values, cmap='viridis', alpha=0.7, s=30)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('t-SNE of Chemical Space (colored by SA_score)')
    plt.colorbar(label='SA_score')
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, output_png))

    print(f"tSNE plot saved to {results_dir}")
    