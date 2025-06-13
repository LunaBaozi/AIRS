import argparse
import subprocess
import sys
import os
import time

def run_script(script_path, args=None):
    cmd = [sys.executable, script_path]
    if args:
        cmd += args
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Wrapper for AIRS pipeline.")
    parser.add_argument('--num_gen', type=int, required=True, help='Number of generations')
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--known_binding_site', type=str, required=True, help='Known binding site')
    args = parser.parse_args()

    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    main_gen_path = os.path.join(base_dir, 'main_gen.py')
    main_eval_path = os.path.join(base_dir, 'main_eval.py')
    post_hoc_dir = os.path.join(base_dir, 'post_hoc_filtering')
    scripts = [
        os.path.join(post_hoc_dir, 'synthesizability_scores.py'),
        os.path.join(post_hoc_dir, 'tanimoto_inter.py'),
        os.path.join(post_hoc_dir, 'tanimoto_intra.py'),
        os.path.join(post_hoc_dir, 'lipinski.py'),
        os.path.join(post_hoc_dir, 'filtering.py'),
    #     os.path.join(post_hoc_dir, 'post_processed_graphics.py'),
    ]

    # Arguments to pass to main_gen.py and main_eval.py
    param_args = [
        '--num_gen', str(args.num_gen),
        '--epoch', str(args.epoch),
        '--known_binding_site', str(args.known_binding_site)
    ]

    # Run main_gen.py and measure execution time
    start_time = time.time()
    run_script(main_gen_path, param_args)
    elapsed_time = time.time() - start_time

    # Run main_eval.py
    run_script(main_eval_path, param_args)

    results_dir = os.path.join(
        base_dir,
        "post_hoc_filtering",
        f"results_epoch_{args.epoch}_mols_{args.num_gen}_bs_{args.known_binding_site}"
    )
    os.makedirs(results_dir, exist_ok=True)
    
    # Run post-hoc filtering scripts (no extra args)
    for script in scripts:
        run_script(script, param_args)

    end_time = time.time() - start_time
    print(f"main_gen.py executed in {elapsed_time:.2f} seconds")
    print(f"Whole pipeline executed in {end_time:.2f} seconds")

if __name__ == '__main__':
    main()