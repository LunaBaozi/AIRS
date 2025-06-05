import argparse
import sys
import os
import pandas as pd
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))

from rdkit import Chem
from SA_Score import sascorer 
from NP_Score import npscorer
from syba.syba import SybaClassifier


from scripts import scscorer_standalone 


# Initialize SCScorer
scscorer = scscorer_standalone.SCScorer()
scscorer.restore()

# Initialize SybaClassifier
syba = SybaClassifier()
syba.fitDefaultScore()


def str2bool(v):
    """
    Converts a string representation of a boolean to a boolean.
    Args:
        v (str): The string to convert.
    Returns:
        bool: The converted boolean value.
    """
    if isinstance(v, bool):
        return v
    if v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def sdf_to_mol(sdf_path, mol_path):
    """
    Converts an SDF file to a MOL file.
    Args:
        sdf_path (str): Path to the input SDF file.
        mol_path (str): Path to the output MOL file.
    """
    suppl = Chem.SDMolSupplier(sdf_path)
    for mol in suppl:
        if mol is not None:
            Chem.MolToMolFile(mol, mol_path)
            break  # Only write the first molecule


def calculate_sa_score(mol):
    """
    Calculates the Synthetic Accessibility Score (SAS) for a given molecule.
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    Returns:
        float: The SAS score of the molecule.
    """
    return sascorer.calculateScore(mol) 


def calculate_np_score(mol):
    """
    Calculates the Novelty and Predictability (NP) score for a given molecule.
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    Returns:
        tuple: A tuple containing the NP score and its confidence.
    """
    fscore = npscorer.readNPModel()
    score = npscorer.scoreMol(mol, fscore)
    confidence = npscorer.scoreMolWConfidence(mol, fscore)
    return score, confidence
    
    

def calculate_scores_to_csv(folder_path, csv_path):
    """
    Calculates scores for all SDF files in a folder and writes them to a CSV file using Pandas.

    Args:
        folder_path (str): Path to the folder containing SDF files.
        csv_path (str): Path to the output CSV file.
    """
    data = []
    idx = 0
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            sdf_path = os.path.join(folder_path, filename)
            suppl = Chem.SDMolSupplier(sdf_path)
            for mol in suppl:
                if mol is not None:
                    sa_score = calculate_sa_score(mol)
                    np_score, _ = calculate_np_score(mol)
                    smi = Chem.MolToSmiles(mol)
                    (smi, sc_score) = scscorer.get_score_from_smi(smi)
                    syba_score = syba.predict(smi)
                    data.append({
                        'index': idx,
                        'smiles': smi,
                        'len_smiles': len(smi),
                        'SA_score': sa_score,
                        'NP_score': np_score,
                        'SCScore': sc_score,
                        'Syba_score': syba_score
                    })
                    idx += 1
                break  # Only process the first molecule in each file
    df = pd.DataFrame(data)
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    df.to_csv(csv_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate synthesizability scores for generated molecules.")
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--num_gen', type=int, required=True, help='Number of molecules generated')
    parser.add_argument('--known_binding_site', type=str2bool, required=True, help='Known binding site (True/False)')
    args = parser.parse_args()

    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    epoch = args.epoch

    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")
    results_dir = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}")
    csv_path = os.path.join(results_dir, f"molecule_scores_{epoch}_{num_gen}.csv")
    calculate_scores_to_csv(
        sdf_folder,
        csv_path
    )


# NOTE: The SA_Score ranges from 1 to 10 with 1 being easy to make and 10 being hard to make.
# NOTE: The NP_Score ranges from -5 to 5 with -5 being easy to make and 5 being hard to make.
# NOTE: The SC_Score ranges from 1 to 5 with 1 being easy to make and 5 being hard to make.
# NOTE: While SYBA score can theoretically assume values between plus and minus infinity, 
#       a majority of compounds will have SYBA score between − 100 and +100 in real applications. 
#       It must be stressed here that the absolute value of the SYBA score is the measure of the 
#       confidence of the prediction and not of the degree of the synthetic accessibility.