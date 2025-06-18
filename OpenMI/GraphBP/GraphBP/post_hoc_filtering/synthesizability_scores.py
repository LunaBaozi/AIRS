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
from scripts.aurk_int_preprocess import read_aurora_kinase_interactions


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

def calculate_sc_score(smi):
    """
    Calculates the Synthetic Accessibility Score (SAS) for a given molecule.
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    Returns:
        float: The SAS score of the molecule.
    """
    return scscorer.get_score_from_smi(smi)


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


def calculate_syba_score(smi):
    """
    Calculates the Synthetic Accessibility Score (SAS) for a given molecule.
    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
    Returns:
        float: The SAS score of the molecule.
    """
    return syba.predict(smi)


def load_mols_from_sdf_folder(folder_path):
    mols = []
    smiles = []
    filenames = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            mol = Chem.SDMolSupplier(os.path.join(folder_path, filename))[0]
            if mol is not None:
                mols.append(mol)
                smiles.append(Chem.MolToSmiles(mol))
                filenames.append(filename)
    return mols, smiles, filenames
    

def calculate_scores(mols, smiles, filenames):
    results = []
    for mol, smi, fn in zip(mols, smiles, filenames):
        if mol is None:
            continue
        sa_score = calculate_sa_score(mol)
        np_score, _ = calculate_np_score(mol)
        (smi, sc_score) = calculate_sc_score(smi)
        syba_score = calculate_syba_score(smi)
        results.append({
            'filename': fn,
            'smiles': smi,
            'len_smiles': len(smi),
            'SA_score': sa_score,
            'SCScore': sc_score,
            'NP_score': np_score,
            'Syba_score': syba_score
        })
    return pd.DataFrame(results)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper for AIRS pipeline.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Number of generations')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Known binding site')
    parser.add_argument('--aurora', type=str, required=True, help='Aurora kinase')
    args = parser.parse_args()

    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    epoch = args.epoch
    aurora = args.aurora

    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
    sdf_folder = os.path.join(parent_dir, f"trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}/sdf")
    known_inhib_file = os.path.join(script_dir, f"data/aurora_kinase_{aurora}_interactions.csv")
    results_dir = os.path.join(script_dir, f"results_epoch_{epoch}_mols_{num_gen}_bs_{known_binding_site}_aurora_{aurora}")
    csv_path = os.path.join(results_dir, f"synthesizability_scores_{epoch}_{num_gen}_{known_binding_site}_{aurora}.csv")
    
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    
    if epoch != 0:
        # Calculating scores for generated molecules
        mols, smiles, filenames = load_mols_from_sdf_folder(sdf_folder)
        synth = calculate_scores(mols, smiles, filenames)
        synth.to_csv(csv_path, index=False)

    else:
        # Calculating scores for Aurora inhibitors
        mols, smiles, filenames, fps = read_aurora_kinase_interactions(known_inhib_file)
        synth = calculate_scores(mols, smiles, filenames)
        synth.to_csv(csv_path, index=False)



# NOTE: The SA_Score ranges from 1 to 10 with 1 being easy to make and 10 being hard to make.
# NOTE: The NP_Score ranges from -5 to 5 with -5 being easy to make and 5 being hard to make.
# NOTE: The SC_Score ranges from 1 to 5 with 1 being easy to make and 5 being hard to make.
# NOTE: While SYBA score can theoretically assume values between plus and minus infinity, 
#       a majority of compounds will have SYBA score between − 100 and +100 in real applications. 
#       It must be stressed here that the absolute value of the SYBA score is the measure of the 
#       confidence of the prediction and not of the degree of the synthetic accessibility.