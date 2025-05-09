
# ciau2 = pd.read_csv('/home/luna/Documents/Coding/GraphBP/GraphBP/aurdata/aurkb_inhibitors_copy.csv', header=0)
# known_aurk_inhibitors2 = dict(zip(known_aurk_inhibitors_csv.inhibitor, known_aurk_inhibitors_csv.smiles))

# def compare_ciaus(dict1, dict2):

#     tanimoto = list()

#     for name, smile in dict1.items():
#         for name2, smile2 in dict2.items():
            
#             if name != name2:
#                 simil = tanimoto_calc(smile, smile2)
#                 if [name2, name, simil] not in tanimoto:
#                     tanimoto.append([name, name2, simil])

#     return tanimoto



# with open('./aurdata/results/output_tanimoto_knownaurkinhibitors.csv', 'w', newline='') as csvfile1:
#     writer = csv.writer(csvfile1)
#     writer.writerows(compare_ciaus(known_aurk_inhibitors, known_aurk_inhibitors2))


# def compare_mol_smiles():

#     tanimoto_basic = list()
#     tanimoto_canon = list()

#     gen_mols_dir = './aurdata/aurkb_gen_complex_1k'

#     for gen_mol_dir in os.listdir(gen_mols_dir):
#         rel_path = os.path.join(gen_mols_dir, gen_mol_dir)
        
#         for filename in os.listdir(rel_path):
#             if 'ligand' in filename:
#                 sppl = Chem.SDMolSupplier(os.path.join(rel_path, filename)) 
                
#                 for mol in sppl:
#                     if mol is not None:  # some compounds cannot be loaded.
#                         basic_smiles_pattern = Chem.MolToSmiles(mol)
#                         canon_smiles_pattern = Chem.CanonSmiles(basic_smiles_pattern)
                        
#                         for name, compound in inhibitors_dict.items():
#                             # print(compound)
#                             basic_smiles_test = str(compound)
#                             canon_smiles_test = Chem.CanonSmiles(compound)

#                             # using basic smiles
#                             simil = tanimoto_calc(basic_smiles_pattern, basic_smiles_test)
#                             tanimoto_basic.append([name, basic_smiles_test, basic_smiles_pattern, simil])

#                             # using canon smiles
#                             simil = tanimoto_calc(canon_smiles_pattern, canon_smiles_test)
#                             tanimoto_canon.append([name, canon_smiles_pattern, canon_smiles_test, simil])

#     return tanimoto_basic, tanimoto_canon #matches




# tanimoto_basic, tanimoto_canon = compare_mol_smiles()

# with open('./aurdata/results/output_tanimoto_basic_1k.csv', 'w', newline='') as csvfile1:
#     writer = csv.writer(csvfile1)
#     writer.writerows(tanimoto_basic)

# with open('./aurdata/results/output_tanimoto_canon_1k.csv', 'w', newline='') as csvfile2:
#     writer = csv.writer(csvfile2)
#     writer.writerows(tanimoto_canon)
