# Copyright (C) 2020 Joseph Young - see GPLv2_note.txt for full notice
# IMPORTS
from CCMI import CCMI
import scipy.io as sio
import numpy as np
import sys

# ESTIMATE CDI
load_file_name = 'ccdi_mat/data_for_ccdi_%s.mat' % sys.argv[1]
CCMI_output = CCMI(sio.loadmat(load_file_name)['X'],
               sio.loadmat(load_file_name)['Y'],
               sio.loadmat(load_file_name)['Z'],
               tester = 'Classifier', metric = 'donsker_varadhan',
               num_boot_iter = int(sys.argv[2]), h_dim = 64, max_ep = 20).get_cmi_est()

est_CDI = CCMI_output[0]
est_CDI_list = CCMI_output[1]

# SAVE ESTIMATE
save_file_name = 'ccdi_mat/ccdi_est_CDI_%s.mat' % sys.argv[1]
sio.savemat(save_file_name,{'est_CDI':est_CDI, 'est_CDI_list':est_CDI_list})
