import pandas as pd
import os

def merge_on_smiles(scores_path="results_epoch_99_mols_1000_bs_False/molecule_scores_99_1000.csv", 
                    lipinski_path="results_epoch_99_mols_1000_bs_False/lipinski_pass_epoch_99_mols_1000_bs_False.csv",
                    tanimoto_path="results_epoch_99_mols_1000_bs_False/tanimoto_results_inter.csv",
                    output_dir="results_epoch_99_mols_1000_bs_False",
                    output_filename="merged_molecule_data.csv"):
    """
    Merges three CSV files on the 'smiles' column, sorts by 'tanimoto' in decreasing order, and saves the result.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_filename)
    scores_df = pd.read_csv(scores_path)
    lipinski_df = pd.read_csv(lipinski_path)
    tanimoto_df = pd.read_csv(tanimoto_path)
    # merged_df = pd.merge(scores_df, lipinski_df, on='smiles', how='outer')
    merged_df = pd.merge(scores_df, tanimoto_df, on='filename', how='outer')
    merged_df = pd.merge(merged_df, lipinski_df, on='filename', how='left')
    merged_df = merged_df.sort_values(by='tanimoto', ascending=False)

    print(set(scores_df['filename']) == set(tanimoto_df['filename']))
    # print(set(lipinski_df['filename']) == set(tanimoto_df['filename']))

    # Specify the desired column order
    desired_columns = [
        'filename', 'smiles', 'len_smiles', 'SA_score', 'NP_score', 'SCScore', 'Syba_score',
        'mol_2', 'tanimoto', 'passed', 'failed'
    ] # 'passed', 'failed',
    # Only keep columns that exist in the merged DataFrame
    columns_to_use = [col for col in desired_columns if col in merged_df.columns]
    merged_df = merged_df[columns_to_use]

    merged_df.to_csv(output_path, index=False)
    return merged_df

merge_on_smiles()


