import os, argparse
import csv
import pandas as pd

import matplotlib.pyplot as plt

def postprocess_vina_results(input_path, output_path):
    # Read and filter lines with at least one comma
    with open(input_path, 'r') as infile:
        lines = [line for line in infile if ',' in line]

    # Parse CSV and sort by affinity_kcal/mol
    reader = csv.DictReader(lines)
    rows = list(reader)
    rows.sort(key=lambda x: float(x['affinity_kcal/mol']))

    # Write sorted rows to output
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def plot_sa_vs_affinity(sa_score_csv, vina_csv):
    # Read CSVs
    sa_df = pd.read_csv(sa_score_csv)
    vina_df = pd.read_csv(vina_csv)

    # Merge on a common column, assuming 'name' is present in both
    merged = pd.merge(sa_df, vina_df, left_on='filename', right_on='ligand', suffixes=('_sa', '_vina'))

    # Prepare data
    x = merged['SA_score']
    y = merged['affinity_kcal/mol']
    colors = merged['tanimoto']
    edge_colors = merged['failed'].apply(lambda v: 'red' if pd.notna(v) and str(v).strip() else 'blue')

    # Plot
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        x, y, c=colors, cmap='viridis', edgecolors=edge_colors, linewidths=1, s=60
    )
    plt.xlabel('SA_score')
    plt.ylabel('affinity_kcal/mol')
    plt.title('SA_score vs. affinity_kcal/mol')
    cbar = plt.colorbar(scatter)
    cbar.set_label('Tanimoto')
    plt.tight_layout()
    plt.savefig("4af3/experiment_epoch_0_mols_0_bs_0_aurora_B/sa_vs_affinity_plot.png")
    plt.show()





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper for CADD pipeline targeted to Aurora protein kinases.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=True, help='Aurora kinase type (str, A, B)')
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    aurora = args.aurora

    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    # sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")
    # known_inhib_file = os.path.join(script_dir, f"data/aurora_kinase_{aurora}_interactions.csv")
    results_dir = os.path.join(parent_dir, f"post_hoc_filtering/results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}_aurora_{aurora}")
    synth_csv = os.path.join(results_dir, f"top_50_sascore_{epoch}_{num_gen}_{known_binding_site}_{aurora}.csv")
    # lipinski_csv = os.path.join(results_dir, f"lipinski_pass_{epoch}_{num_gen}_{known_binding_site}_{aurora}.csv")
    # tanimoto_inter_csv = os.path.join(results_dir, f"tanimoto_inter_{epoch}_{num_gen}_{known_binding_site}_{aurora}.csv")
    
    if aurora == "A":
        aur_type = '4cfg'
    elif aurora == "B":
        aur_type = '4af3'
    else:
        raise ValueError("Aurora type must be 'A' or 'B'.")
   
    experiment_dir = os.path.join(script_dir, f"{aur_type}/experiment_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}_aurora_{aurora}")
    vina_csv = os.path.join(experiment_dir, f"vina_results.csv")
    vina_postprocess = os.path.join(experiment_dir, "vina_results_postprocessed.csv")

    os.makedirs(os.path.dirname(vina_postprocess), exist_ok=True)

    postprocess_vina_results(
        vina_csv,
        vina_postprocess
    )
    
    plot_sa_vs_affinity(
        synth_csv,
        vina_postprocess
    )

    print("Results saved.")

    