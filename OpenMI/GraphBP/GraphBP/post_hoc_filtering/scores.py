import sys
import os
import csv
sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'share','RDKit','Contrib'))

from rdkit import Chem
from SA_Score import sascorer 
from NP_Score import npscorer
from syba.syba import SybaClassifier


from scscorer_standalone import SCScorer

# Initialize SCScorer
scscorer = SCScorer()
scscorer.restore()

# Initialize SybaClassifier
syba = SybaClassifier()
syba.fitDefaultScore()

epoch = 99
num_gen = 10000

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
    Calculates scores for all SDF files in a folder and writes them to a CSV file.

    Args:
        folder_path (str): Path to the folder containing SDF files.
        csv_path (str): Path to the output CSV file.
    """
    with open(csv_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['index', 'smiles', 'len_smiles', 'SA_score', 'NP_score', 'SCScore', 'Syba_score'])
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
                        writer.writerow([idx, smi, len(smi), sa_score, np_score, sc_score, syba_score])
                        idx += 1
                    break  # Only process the first molecule in each file

# Example usage:
calculate_scores_to_csv(
    f"../trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_{num_gen}/sdf",
    f"molecule_scores_{epoch}_{num_gen}.csv"
)

# NOTE: The SA_Score ranges from 1 to 10 with 1 being easy to make and 10 being hard to make.
# NOTE: The NP_Score ranges from -5 to 5 with -5 being easy to make and 5 being hard to make.
# NOTE: The SC_Score ranges from 1 to 5 with 1 being easy to make and 5 being hard to make.
# NOTE: While SYBA score can theoretically assume values between plus and minus infinity, 
#       a majority of compounds will have SYBA score between − 100 and +100 in real applications. 
#       It must be stressed here that the absolute value of the SYBA score is the measure of the 
#       confidence of the prediction and not of the degree of the synthetic accessibility.