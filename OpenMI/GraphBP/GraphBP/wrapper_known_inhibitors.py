import argparse
import subprocess
import os, sys, csv, time


def run_script(script_path, args=None):
    cmd = [sys.executable, script_path]
    if args:
        cmd += args
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Wrapper for CADD pipeline targeted to Aurora protein kinases.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=True, help='Aurora kinase type (str, A, B)')
    args = parser.parse_args()

    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    post_hoc_dir = os.path.join(base_dir, 'post_hoc_filtering')
    scripts = [
        # os.path.join(post_hoc_dir, 'synthesizability_scores.py'),
        os.path.join(post_hoc_dir, 'lipinski.py'),
        # os.path.join(post_hoc_dir, 'tanimoto_inter.py'),
        # os.path.join(post_hoc_dir, 'tanimoto_intra.py'),
        # os.path.join(post_hoc_dir, 'lipinski.py'),
        # os.path.join(post_hoc_dir, 'filtering.py'),
    #     os.path.join(post_hoc_dir, 'post_processed_graphics.py'),
    ]

    # Arguments to pass to main_gen.py and main_eval.py
    param_args = [
        '--num_gen', str(args.num_gen),
        '--epoch', str(args.epoch),
        '--known_binding_site', str(args.known_binding_site),
        '--aurora', str(args.aurora).upper()
    ]

    # Run analyses and measure execution time
    start_time = time.time()

    results_dir = os.path.join(
        base_dir,
        "post_hoc_filtering",
        f"results_epoch_{args.epoch}_mols_{args.num_gen}_bs_{args.known_binding_site}_aurora_{args.aurora}"
    )
    os.makedirs(results_dir, exist_ok=True)
    
    # Run post-hoc filtering scripts (no extra args)
    for script in scripts:
        run_script(script, param_args)

    end_time = time.time() - start_time
    # print(f"main_gen.py executed in {elapsed_time:.2f} seconds")
    print(f"Whole pipeline executed in {end_time:.2f} seconds")

if __name__ == '__main__':
    main()