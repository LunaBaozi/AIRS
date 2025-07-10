conf = {}
conf_aurkA = {}
conf_aurkB = {}
conf_aurkA_NOBS = {}
conf_aurkB_NOBS = {}


conf_aurkA['pdb'] = '4ztq'
conf_aurkA['rec_name'] = '4ztq_A_rec'
conf_aurkA['pdb_file_path'] = '../datav11/crossdock2020/AURKA_HUMAN_122_397_0/4ztq_A_rec.pdb'
conf_aurkA['gen_ligands_path'] = '../trained_model/gen_mols_epoch_33_aurkA/'
conf_aurkA['gen_complexes'] = './results/AURKA_GEN_COMPLEX_4byi'

conf_aurkB['pdb'] = '4af3'
conf_aurkB['rec_name'] = '4af3_A_rec'
conf_aurkB['pdb_file_path'] = '../datav11/crossdock2020/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
conf_aurkB['gen_ligands_path'] = '../trained_model/gen_mols_epoch_33_aurkB/'
conf_aurkB['gen_complexes'] = './results/AURKB_GEN_COMPLEX'

conf_aurkA_NOBS['pdb'] = '4ztq'
conf_aurkA_NOBS['rec_name'] = '4ztq_A_rec'
conf_aurkA_NOBS['pdb_file_path'] = '../datav11/crossdock2020/AURKA_HUMAN_122_397_0/4ztq_A_rec.pdb'
conf_aurkA_NOBS['gen_ligands_path'] = '../trained_model/gen_mols_epoch_33_aurkA_NOKNOWNBS/'
conf_aurkA_NOBS['gen_complexes'] = './results/AURKA_NOBS_GEN_COMPLEX'

conf_aurkB_NOBS['pdb'] = '4af3'
conf_aurkB_NOBS['rec_name'] = '4af3_A_rec'
conf_aurkB_NOBS['pdb_file_path'] = '../datav11/crossdock2020/AURKB_HUMAN_54_344_0/4af3_A_rec.pdb'
conf_aurkB_NOBS['gen_ligands_path'] = '../trained_model/gen_mols_epoch_33_aurkB_NOKNOWNBS/'
conf_aurkB_NOBS['gen_complexes'] = './results/AURKB_NOBS_GEN_COMPLEX'

conf[1] = conf_aurkA
conf[2] = conf_aurkB
conf[3] = conf_aurkA_NOBS
conf[4] = conf_aurkB_NOBS