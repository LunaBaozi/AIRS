#!/usr/bin/env python

import os
import csv

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

epoch = 99
num_gen = 10000

class SmilesError(Exception): pass

def log_partition_coefficient(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        raise SmilesError('%s returns a None molecule' % smiles)
    return Crippen.MolLogP(mol)

def lipinski_trial(sdf_path):
    passed = []
    failed = []
    smiles = []
    suppl = Chem.SDMolSupplier(sdf_path)
    for mol in suppl:
        if mol is None:
            raise Exception('Not a valid mol')
        smiles.append(Chem.MolToSmiles(mol))
        num_hdonors = Lipinski.NumHDonors(mol)
        num_hacceptors = Lipinski.NumHAcceptors(mol)
        mol_weight = Descriptors.MolWt(mol)
        mol_logp = Crippen.MolLogP(mol)
        if num_hdonors > 5:
            failed.append('Over 5 H-bond donors, found %s' % num_hdonors)
        else:
            passed.append('Found %s H-bond donors' % num_hdonors)
        if num_hacceptors > 10:
            failed.append('Over 10 H-bond acceptors, found %s' % num_hacceptors)
        else:
            passed.append('Found %s H-bond acceptors' % num_hacceptors)
        if mol_weight >= 500:
            failed.append('Molecular weight over 500, calculated %s' % mol_weight)
        else:
            passed.append('Molecular weight: %s' % mol_weight)
        if mol_logp >= 5:
            failed.append('Log partition coefficient over 5, calculated %s' % mol_logp)
        else:
            passed.append('Log partition coefficient: %s' % mol_logp)
    return passed, failed, smiles

def lipinski_pass(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            sdf_path = os.path.join(folder_path, filename)
    passed, failed, smiles = lipinski_trial(sdf_path)
    if failed:
        return False
    else:
        return True

def save_lipinski_results_to_csv(folder_path, csv_path):
    '''
    Runs Lipinski trial on all SDF files in a folder and saves results to a CSV.
    '''
    results = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.sdf'):
            sdf_path = os.path.join(folder_path, filename)
            try:
                passed, failed, smiles = lipinski_trial(sdf_path)
                results.append({
                    'filename': filename,
                    'smiles': '; '.join(smiles),
                    'passed': '; '.join(passed),
                    'failed': '; '.join(failed)
                })
            except Exception as e:
                results.append({
                    'filename': filename,
                    'smiles': smiles,
                    'passed': '',
                    'failed': f'Error: {str(e)}'
                })
    with open(csv_path, 'w', newline='') as csvfile:
        fieldnames = ['filename', 'smiles', 'passed', 'failed']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

# Example usage:
save_lipinski_results_to_csv(f'../trained_model_reduced_dataset_100_epochs/gen_mols_epoch_{epoch}_{num_gen}/sdf', 
                             f'results/lipinski_pass_{epoch}_{num_gen}.csv')
