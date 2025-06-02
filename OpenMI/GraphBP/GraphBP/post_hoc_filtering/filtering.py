import os
import pandas as pd

def get_top100_per_metric(results_dir, input_filename):
    input_path = os.path.join(results_dir, input_filename)
    df = pd.read_csv(input_path)

    # Load Lipinski violations
    lipinski_path = os.path.join(results_dir, "lipinski_pass_99_10000.csv")
    lipinski_df = pd.read_csv(lipinski_path)
    # Assume molecule id is the second column in both files
    mol_id_col = df.columns[1]
    # Assume violations column is named 'violations'
    violations_col = 'failed'
    if violations_col not in lipinski_df.columns:
        raise ValueError(f"Column '{violations_col}' not found in Lipinski file.")

    # Compute rankings for all metrics
    metric_columns = df.columns[2:]
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
        cols_to_save = [mol_id_col, metric] + [f"rank_in_{m}" for m in metric_columns] + [violations_col]
        top100.to_csv(output_path, columns=cols_to_save, index=False)
        print(f"Saved top 100 for {metric} to {output_filename}")

if __name__ == "__main__":
    results_folder = os.path.join(os.path.dirname(__file__), "results")
    get_top100_per_metric(results_folder, "molecule_scores_99_10000.csv")