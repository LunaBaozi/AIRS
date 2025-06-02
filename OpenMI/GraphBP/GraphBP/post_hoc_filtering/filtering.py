import os
import pandas as pd

def get_top100_per_metric(results_dir, input_filename):
    input_path = os.path.join(results_dir, input_filename)
    df = pd.read_csv(input_path)

    # Assume the first column is SMILES, second is molecule name/id, rest are metrics
    smiles_col = df.columns[1]
    metric_columns = df.columns[2:]

    # Compute rankings for all metrics
    rankings = {}
    for metric in metric_columns:
        rankings[metric] = df[metric].rank(method='min', ascending=True).astype(int)

    for metric in metric_columns:
        top100 = df.nsmallest(100, metric).copy()
        # Add ranking columns for all metrics
        for other_metric in metric_columns:
            rank_col = f"rank_in_{other_metric}"
            # Map the index of top100 to the original df's rankings
            top100[rank_col] = rankings[other_metric].loc[top100.index].values

        output_filename = f"top100_synthesizable_{metric}.csv"
        output_path = os.path.join(results_dir, output_filename)
        # Save SMILES, metric, and all rank columns
        cols_to_save = [smiles_col, metric] + [f"rank_in_{m}" for m in metric_columns]
        top100.to_csv(output_path, columns=cols_to_save, index=False)
        print(f"Saved top 100 for {metric} to {output_filename}")

if __name__ == "__main__":
    results_folder = os.path.join(os.path.dirname(__file__), "results")
    get_top100_per_metric(results_folder, "molecule_scores_99_10000.csv")