import os
import argparse
import pandas as pd

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_top100_per_metric(results_dir, input_filename):
    input_path = os.path.join(results_dir, input_filename)
    df = pd.read_csv(input_path)

    # "filename" is now the first column, "smiles" is the second
    filename_col = "filename"
    mol_id_col = "smiles"

    # Load Lipinski violations
    lipinski_path = os.path.join(results_dir, f"lipinski_pass_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}.csv")
    lipinski_df = pd.read_csv(lipinski_path)
    # Assume molecule id is the "smiles" column in both files
    violations_col = 'failed'
    if violations_col not in lipinski_df.columns:
        raise ValueError(f"Column '{violations_col}' not found in Lipinski file.")

    # Compute rankings for all metrics (skip first two columns: filename, smiles)
    metric_columns = df.columns[3:]
    rankings = {}
    for metric in metric_columns:
        rankings[metric] = df[metric].rank(method='min', ascending=True).astype(int)

    for metric in metric_columns:
        top100 = df.nsmallest(100, metric).copy()
        # Add ranking columns for all metrics
        for other_metric in metric_columns:
            rank_col = f"rank_in_{other_metric}"
            top100[rank_col] = rankings[other_metric].loc[top100.index].values

        # Merge Lipinski violations
        top100 = top100.merge(
            lipinski_df[[mol_id_col, violations_col]],
            on=mol_id_col,
            how='left'
        )

        output_filename = f"top100_synthesizable_{metric}.csv"
        output_path = os.path.join(results_dir, output_filename)
        # Output columns: filename, smiles, metric, all ranks, violations
        cols_to_save = [filename_col, mol_id_col, metric] + [f"rank_in_{m}" for m in metric_columns] + [violations_col]
        top100.to_csv(output_path, columns=cols_to_save, index=False)
        print(f"Saved top 100 for {metric} to {output_filename}")

    # --- New code: Save top 100 Tanimoto pairs with scores for mol_1 ---
    tanimoto_path = os.path.join(results_dir, "tanimoto_results_inter.csv")
    if os.path.exists(tanimoto_path):
        tanimoto_df = pd.read_csv(tanimoto_path)
        # Get top 100 pairs by tanimoto
        top100_pairs = tanimoto_df.nlargest(100, "tanimoto").copy()
        # Merge scores for mol_1 from df (input_filename)
        # Assume mol_1 column contains the smiles string
        scores_cols = [col for col in df.columns if col not in [filename_col, mol_id_col]]
        top100_pairs = top100_pairs.merge(
            df[[mol_id_col] + list(scores_cols)],
            left_on="mol_1",
            right_on=mol_id_col,
            how="left",
            suffixes=('', '_mol1')
        )
        # Save to CSV
        output_tanimoto_filename = "top100_tanimoto_pairs_with_scores.csv"
        output_tanimoto_path = os.path.join(results_dir, output_tanimoto_filename)
        top100_pairs.to_csv(output_tanimoto_path, index=False)
        print(f"Saved top 100 Tanimoto pairs with mol_1 scores to {output_tanimoto_filename}")
    else:
        print(f"Tanimoto file {tanimoto_path} not found, skipping top100 Tanimoto pairs output.")

def get_top_tanimoto_pairs(results_dir, tanimoto_filename, top_n=100, group="inter"):
    """
    Save the top N pairs of molecules with the highest Tanimoto similarity.
    """
    tanimoto_path = os.path.join(results_dir, tanimoto_filename)
    df = pd.read_csv(tanimoto_path)
    # Columns: 'filename', 'mol_1', 'mol_2', 'tanimoto'
    required_cols = {'filename', 'mol_1', 'mol_2', 'tanimoto'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Tanimoto file must contain columns: {required_cols}. Found columns: {list(df.columns)}")
    top_pairs = df.nlargest(top_n, 'tanimoto').copy()
    output_filename = f"top{top_n}_tanimoto_pairs_{group}.csv"
    output_path = os.path.join(results_dir, output_filename)
    # Save all columns, including 'filename'
    top_pairs.to_csv(output_path, index=False)
    print(f"Saved top {top_n} Tanimoto pairs to {output_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and rank molecules based on metrics and Tanimoto similarity.")
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--num_gen', type=int, required=True, help='Number of molecules generated')
    parser.add_argument('--known_binding_site', type=str2bool, required=True, help='Known binding site (True/False)')
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site

    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_folder = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}")
    
    # Uncomment if you want to run top100 per metric as well:
    get_top100_per_metric(results_folder, f"molecule_scores_{epoch}_{num_gen}.csv")
    get_top_tanimoto_pairs(results_folder, "tanimoto_results_inter.csv", top_n=100, group="inter")
    get_top_tanimoto_pairs(results_folder, "tanimoto_results_intra.csv", top_n=100, group="intra")