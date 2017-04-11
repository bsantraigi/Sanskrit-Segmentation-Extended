import os

folders = ['skt_dcs_DS.bz2', 'skt_dcs_DS.bz2_10K', 'skt_dcs_DS.bz2_1L_bigram_10K', \
'skt_dcs_DS.bz2_1L_bigram_heldout', 'skt_dcs_DS.bz2_1L_pmi_10K', 'skt_dcs_DS.bz2_1L_pmi_heldout', \
'skt_dcs_DS.bz2_4K_pmi_10K', 'skt_dcs_DS.bz2_4K_pmi_heldout', 'skt_dcs_DS.bz2_heldout_mifeats', 'skt_dcs_DS.bz2_1L_bigram_10K_d2K',\
          'skt_dcs_DS.bz2_1L_bigram_heldout_d2K']

folders = ['skt_dcs_DS.bz2_1L_bigram_10K_d2K', 'skt_dcs_DS.bz2_1L_bigram_heldout_d2K']
folders = ['skt_dcs_DS.bz2_4K_bigram_rfe_10K', 'skt_dcs_DS.bz2_4K_bigram_rfe_heldout', 'skt_dcs_DS.bz2_1L_pmi_mir_10K_again', 'skt_dcs_DS.bz2_1L_pmi_mir_heldout_again', 'skt_dcs_DS.bz2_1L_bigram_rfe_10K', 'skt_dcs_DS.bz2_1L_bigram_rfe_heldout']
for folder in folders:
	path = os.path.join('../NewData', folder)
	c = len(os.listdir(path))
	print('Folder: {:35s} ------> File_Count: {}\n'.format(folder, c))
