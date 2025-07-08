# import os
import numpy as np
import scikit_posthocs as sp
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit.Chem.QED import qed

from aurk_int_preprocess import read_aurora_kinase_interactions
from gen_mols_preprocess import load_mols_from_sdf_folder
import pandas as pd

def calculate_qed(mol):
    """
    Calculates the Quantitative Estimate of Drug-likeness (QED) score for a given molecule.
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    Returns:
        tuple: A tuple containing the QED score.
    """
    return qed(mol)

def merge_data(aurora_file, mols_folder):
    """
    Merges Aurora kinase interactions with molecular data from SDF files into a single DataFrame.

    Args:
        aurora_file (str): Path to the Aurora kinase interactions file.
        mols_folder (str): Path to the folder containing SDF files.

    Returns:
        pd.DataFrame: DataFrame with columns ['filename', 'source', 'mols', 'smiles', 'fps'].
    """
    mols_aur, smiles_aur, filenames_aur, fps_aur = read_aurora_kinase_interactions(aurora_file)
    mols_gen, smiles_gen, filenames_gen, fps_gen = load_mols_from_sdf_folder(mols_folder)

    df_aur = pd.DataFrame({
        'filename': filenames_aur,
        'source': ['aurora'] * len(filenames_aur),
        'mols': mols_aur,
        'smiles': smiles_aur,
        'fps': fps_aur
    })

    df_gen = pd.DataFrame({
        'filename': filenames_gen,
        'source': ['gen'] * len(filenames_gen),
        'mols': mols_gen,
        'smiles': smiles_gen,
        'fps': fps_gen
    })

    df = pd.concat([df_aur, df_gen], ignore_index=True)
    df['QED'] = df['mols'].apply(calculate_qed)
    df.to_csv("merged_data.csv", index=False)
    return df

def kdeplotting(df):
    sns.set_theme(rc={'figure.figsize':(6,6)})
    sns.set_style('white')
    ax = sns.kdeplot(x="QED", hue="source", data=df, common_norm=False)
    ax.set_ylabel("QED Score")
    # plt.xlabel("QED Score")
    # plt.ylabel("Density")
    plt.savefig("kdeplotting.png")
    plt.clf()
    # plt.show()
    return

def boxplotting(df):
    sns.set_theme(rc={'figure.figsize':(6,6)})
    sns.set_style('white')
    # sns.set_context('notebook')
    ax = sns.boxenplot(x="source",y="QED",data=df)
    ax.set_ylabel("QED Score")
    plt.savefig("boxplotting.png")
    plt.clf()
    # plt.show()
    return

def sig_plot_maybe(df):
    sns.set_theme(rc={'figure.figsize':(8,6)})
    sns.set_context("notebook", font_scale=1.5)
    sns.set_style('white')
    pc = sp.posthoc_mannwhitney(df, val_col="QED", group_col="source", p_adjust='holm')
    heatmap_args = {'linewidths':0.25, 'linecolor':'0.5', 'clip_on':False, 'square':True, 'cbar_ax_bbox':[0.80, 0.35, 0.04, 0.3]}
    _ = sp.sign_plot(pc, **heatmap_args)
    plt.savefig("sig_plot.png")
    plt.clf()
    return

merged_data = merge_data("/vol/data/airs_unbroken/AIRS/OpenMI/GraphBP/GraphBP/post_hoc_filtering/data/aurora_kinase_B_interactions.csv", 
           "/vol/data/airs_unbroken/AIRS/OpenMI/GraphBP/GraphBP/trained_model_reduced_dataset_100_epochs/gen_mols_epoch_99_mols_1000_bs_False_aurora_B/sdf")

boxplotting(merged_data)
    
sig_plot_maybe(merged_data)

kdeplotting(merged_data)

# =====================================
# 2. MOLECULAR PROPERTY DISTRIBUTIONS
# =====================================
def plot_molecular_properties(df):
    """Plot distributions of molecular properties"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # properties = ['synthesizability', 'qed', 'tanimoto_similarity', 
    #              'molecular_weight', 'logp', 'num_rotatable_bonds']
    properties = ["SA_score", "SCScore", "NP_score", "QED", 
                  "tanimoto", "len_smiles"]
    
    for i, prop in enumerate(properties):
        row, col = i // 3, i % 3
        ax = axes[row, col]
        
        if prop in df.columns:
            # Histogram with KDE
            ax.hist(df[prop], bins=30, alpha=0.7, density=True, color='skyblue', edgecolor='black')
            
            # Add KDE curve
            from scipy import stats
            x = np.linspace(df[prop].min(), df[prop].max(), 100)
            kde = stats.gaussian_kde(df[prop])
            ax.plot(x, kde(x), 'r-', linewidth=2, label='KDE')
            
            # Add vertical line for mean
            mean_val = df[prop].mean()
            ax.axvline(mean_val, color='red', linestyle='--', alpha=0.8, 
                      label=f'Mean: {mean_val:.3f}')
            
            ax.set_title(f'{prop.replace("_", " ").title()} Distribution', fontsize=12, weight='bold')
            ax.set_xlabel(prop.replace("_", " ").title())
            ax.set_ylabel('Density')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig