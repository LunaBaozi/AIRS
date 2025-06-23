import argparse
import subprocess
import os, sys, csv, time


def run_script(script_path, args=None):
    # Determine the correct executable
    if script_path.endswith(".sh"):
        cmd = ["bash", script_path]
    elif script_path.endswith(".py"):
        cmd = [sys.executable, script_path]
    else:
        raise ValueError("Unsupported script type. Only .sh and .py are supported.")
    if args:
        cmd += args
    subprocess.run(cmd, check=True)

def run_bash_script_in_conda(script_path, args, conda_env):
    arg_str = ' '.join(args) if args else ''
    command = (
        f"source ~/../../vol/data/miniconda3/etc/profile.d/conda.sh && "
        f"conda activate {conda_env} && "
        f"bash {script_path} {arg_str}"
    )
    subprocess.run(["bash", "-c", command], check=True)



def main():
    parser = argparse.ArgumentParser(description="Wrapper for CADD pipeline targeted to Aurora protein kinases.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=False, default='A', help='Aurora kinase type (str, A, B)')
    args = parser.parse_args()

    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    post_hoc_dir = os.path.join(base_dir, 'post_hoc_filtering')
    docking_dir = os.path.join(base_dir, 'docking')
    
    bash_scripts = [
    os.path.join(docking_dir, 'run_pipeline.sh'),
    ]

    python_scripts = [
        # os.path.join(post_hoc_dir, 'synthesizability_scores.py'),
        # os.path.join(post_hoc_dir, 'lipinski.py'),
        # os.path.join(post_hoc_dir, 'tanimoto_intra.py'),
        # os.path.join(post_hoc_dir, 'tanimoto_inter.py'),
        # os.path.join(post_hoc_dir, 'graphics.py'),
        os.path.join(post_hoc_dir, 'post_processing.py'),
    ]

    bash_param_args = [
        str(args.num_gen),
        str(args.epoch),
        str(args.known_binding_site),
        str(args.aurora).upper()
    ]

    python_param_args = [
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
    
    for script in python_scripts:
        run_script(script, python_param_args)

    # for script in bash_scripts:
    #     run_bash_script_in_conda(script, bash_param_args, conda_env='vina')

    run_script(os.path.join(docking_dir, 'top_scoring_docking.py'), python_param_args)
    
    end_time = time.time() - start_time
    print(f"Whole pipeline executed in {end_time:.2f} seconds")

if __name__ == '__main__':
    main()