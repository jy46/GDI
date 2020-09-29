# Copyright (C) 2020 Joseph Young - see GPLv2_note.txt for full notice

################################################################################
# IMPORTS
################################################################################
from CCMI import CCMI
import numpy as np
import numpy.matlib
import copy
import itertools
from multiprocessing import Pool
from functools import partial
################################################################################





################################################################################
# pair_GDI: COMPUTE GDI BETWEEN TWO CHANNELS
#   INPUTS:
#       X_past_win: 2D array where rows are samples and columns are dimensions, 
#                   i.e. variables representing past values of each channel.
#       X_current: Vector containing samples of current value of target channel.
#       M: History length, i.e. number of past samples to use.
#       B: Number of bootstrap iterations to use for training classifiers
#       chan_pair: Tuple specifying channels to compute GDI between.
#                  GDI is computed from first element to second element.
#   OUTPUTS:
#       GDI_estimate: Estimate of the GDI from first channel of chan_pair to 
#                     second channel of chan_pair.
################################################################################
def pair_GDI(X_past_win,X_current,M,B,chan_pair):
    chan1 = chan_pair[0]
    chan2 = chan_pair[1]
    if chan1!=chan2:
        # CREATE PAST WINDOW TO CONDITION ON
        X_past_win_cond = np.delete(copy.deepcopy(X_past_win),
                            np.arange((chan1*M),((chan1+1)*M)).tolist(),
                            axis=1)

        # ESTIMATE GDI
        GDI_estimate = CCMI(X_past_win[:,(chan1*M):((chan1+1)*M)],
            X_current[:,[chan2]],
            X_past_win_cond,
            tester = 'Classifier',
            metric = 'donsker_varadhan',
            num_boot_iter = B,
            h_dim = 64, max_ep = 20).get_cmi_est()

    else:
        GDI_estimate = 0

    return GDI_estimate
################################################################################
# END OF pair_GDI
################################################################################





################################################################################
# GDI: COMPUTE GDI BETWEEN COLUMNS OF X
#   INPUTS:
#       X: Input data with dim (sample)x(channel)
#       M: History length, i.e. number of past samples to use
#       B: Number of bootstrap iterations to use for training classifiers
#   OUTPUTS:
#       GDI_estimate: Estimate of the GDI from rows to columns. Shape is:
#                     (channel)x(channel)
################################################################################
def GDI(X,M,B):

    # MAKE SURE INTEGERS
    M = int(M)
    B = int(B)

    # REMOVE EXCESS SAMPLES AT END
    num_channels = X.shape[1]
    num_samples = X.shape[0]
    num_windows = int(np.floor(num_samples/(M+1)))
    num_samples_to_keep = int(num_windows*(M+1))
    X_trim = X[:num_samples_to_keep,:]

    # RESHAPE TO GET ARRAY OF PAST VALUES FOR EACH CHANNEL AND 
    # TO GET ANOTHER ARRAY OF CURRENT VALUES FOR EACH CHANNEL
    X_past_win = np.zeros((num_windows,num_channels*M))
    X_current  = np.zeros((num_windows,num_channels))
    for chan in range(num_channels):
        current_chan_mat = np.reshape(X_trim[:,chan],(num_windows,(M+1)))
        X_past_win[:,(chan*M):((chan+1)*M)] = current_chan_mat[:,:M]
        X_current[:,chan] = current_chan_mat[:,M]

    # ESTIMATE GDI
    chan_list = np.arange(num_channels).tolist()
    chan_list_pairs = list(itertools.product(chan_list, chan_list))

    pool = Pool()
    func = partial(pair_GDI,X_past_win,X_current,M,B)
    GDI_estimate_list = pool.map(func,chan_list_pairs)

    # UNWRAP
    GDI_estimate = np.zeros((num_channels,num_channels))
    for ii in range(num_channels*num_channels):
        GDI_estimate[chan_list_pairs[ii][0],chan_list_pairs[ii][1]] = GDI_estimate_list[ii]

    # RETURN GDI ESTIMATE
    return GDI_estimate
################################################################################
# END OF GDI
################################################################################





################################################################################
# pair_GDI_mask: COMPUTE GDI BETWEEN TWO CHANNELS WITH ONLY CONDITIONING ON 
#                CHANNELS SPECIFIED IN MASK
#   INPUTS:
#       X_past_win: 2D array where rows are samples and columns are dimensions, 
#                   i.e. variables representing past values of each channel.
#       X_current: Vector containing samples of current value of target channel.
#       M: History length, i.e. number of past samples to use.
#       B: Number of bootstrap iterations to use for training classifiers
#       mask: 2D array with dimensions (channel)x(channel). This function looks
#             at the column mask[:,chan_pair[1]] and only conditions GDI on 
#             the other channels of that column which are 1.
#       chan_pair: Tuple specifying channels to compute GDI between.
#                  GDI is computed from first element to second element.
#   OUTPUTS:
#       GDI_estimate: Estimate of the GDI from first channel of chan_pair to 
#                     second channel of chan_pair.
################################################################################
def pair_GDI_mask(X_past_win,X_current,M,B,mask,chan_pair):
    chan1 = chan_pair[0]
    chan2 = chan_pair[1]
    
    # CHECK TO MAKE SURE NOT DIAGONAL & THAT MASK SPECIFIED TO DO THIS COMPUTATION
    if (chan1!=chan2) and (mask[chan1,chan2]!=0):

        # USE MASK TO KNOW WHICH CHANNELS TO CONDITION ON
        mask_col = mask[:,chan2]
        col_of_X_past_win_to_delete = np.arange((chan1*M),((chan1+1)*M)).tolist()
        for ii in range(mask.shape[0]):
            if ii!=chan1:
                if ii!=chan2:
                    if mask_col[ii]==0:
                        col_of_X_past_win_to_delete = col_of_X_past_win_to_delete+np.arange((ii*M),((ii+1)*M)).tolist()

        # CREATE PAST WINDOW TO CONDITION ON
        X_past_win_cond = np.delete(copy.deepcopy(X_past_win),
                            col_of_X_past_win_to_delete,
                            axis=1)

        # ESTIMATE GDI
        GDI_estimate = CCMI(X_past_win[:,(chan1*M):((chan1+1)*M)],
            X_current[:,[chan2]],
            X_past_win_cond,
            tester = 'Classifier',
            metric = 'donsker_varadhan',
            num_boot_iter = B,
            h_dim = 64, max_ep = 20).get_cmi_est()

    else:
        GDI_estimate = 0

    return GDI_estimate
################################################################################
# END OF pair_GDI_mask
################################################################################





################################################################################
# GDI_mask: COMPUTE GDI BETWEEN COLUMNS OF X & CONDITION AS SPECIFIED BY MASK
#   INPUTS:
#       X: Input data with dim (sample)x(channel)
#       M: History length, i.e. number of past samples to use
#       B: Number of bootstrap iterations to use for training classifiers
#       mask: Square matrix containing zeros and ones which specifies which 
#             time series / columns GDI is to be computed for as well as which 
#             time series are to be conditioned on in such GDI computations. 
#             For example, a mask with all zeros except for ones at elements 
#             (2,3), (4,3), and (5,4) would mean that GDI would only be computed 
#             from columns 2 to 3, 4 to 3, and 5 to 4 of X. Furthermore, the 
#             GDI computation from column 2 to 3 would only be conditioned on 
#             column 4, while the GDI from 4 to 3 would only be conditioned on 2, 
#             and finally the GDI from 5 to 4 would not be conditioned on any 
#             other column, making it equivalent to the DI from 5 to 4. 
#   OUTPUTS:
#       GDI_estimate: Estimate of the GDI from rows to columns. Shape is:
#                     (channel)x(channel)
################################################################################
def GDI_mask(X,M,B,mask):

    # MAKE SURE INTEGERS
    M = int(M)
    B = int(B)

    # REMOVE EXCESS SAMPLES AT END
    num_channels = X.shape[1]
    num_samples = X.shape[0]
    num_windows = int(np.floor(num_samples/(M+1)))
    num_samples_to_keep = int(num_windows*(M+1))
    X_trim = X[:num_samples_to_keep,:]

    # RESHAPE TO GET ARRAY OF PAST VALUES FOR EACH CHANNEL AND 
    # TO GET ANOTHER ARRAY OF CURRENT VALUES FOR EACH CHANNEL
    X_past_win = np.zeros((num_windows,num_channels*M))
    X_current  = np.zeros((num_windows,num_channels))
    for chan in range(num_channels):
        current_chan_mat = np.reshape(X_trim[:,chan],(num_windows,(M+1)))
        X_past_win[:,(chan*M):((chan+1)*M)] = current_chan_mat[:,:M]
        X_current[:,chan] = current_chan_mat[:,M]

    # ESTIMATE GDI
    chan_list = np.arange(num_channels).tolist()
    chan_list_pairs = list(itertools.product(chan_list, chan_list))

    pool = Pool()
    func = partial(pair_GDI_mask,X_past_win,X_current,M,B,mask)
    GDI_estimate_list = pool.map(func,chan_list_pairs)

    # UNWRAP
    GDI_estimate = np.zeros((num_channels,num_channels))
    for ii in range(num_channels*num_channels):
        GDI_estimate[chan_list_pairs[ii][0],chan_list_pairs[ii][1]] = GDI_estimate_list[ii]

    # RETURN GDI ESTIMATE
    return GDI_estimate
################################################################################
# END OF GDI_mask
################################################################################





################################################################################
# pair_DI: COMPUTE (NON-GRAPHICAL) DI BETWEEN TWO CHANNELS
#   INPUTS:
#       X_past_win: 2D array where rows are samples and columns are dimensions, 
#                   i.e. variables representing past values of each channel.
#       X_current: Vector containing samples of current value of target channel.
#       M: History length, i.e. number of past samples to use.
#       B: Number of bootstrap iterations to use for training classifiers
#       chan_pair: Tuple specifying channels to compute DI between.
#                  DI is computed from first element to second element.
#   OUTPUTS:
#       DI_estimate: Estimate of the DI from first channel of chan_pair to 
#                     second channel of chan_pair.
################################################################################
def pair_DI(X_past_win,X_current,M,B,chan_pair):
    chan1 = chan_pair[0]
    chan2 = chan_pair[1]
    if chan1!=chan2:

        # ESTIMATE DI
        DI_estimate = CCMI(X_past_win[:,(chan1*M):((chan1+1)*M)],
            X_current[:,[chan2]],
            X_past_win[:,(chan2*M):((chan2+1)*M)],
            tester = 'Classifier',
            metric = 'donsker_varadhan',
            num_boot_iter = B,
            h_dim = 64, max_ep = 20).get_cmi_est()

    else:
        DI_estimate = 0

    return DI_estimate
################################################################################
# END OF pair_DI
################################################################################





################################################################################
# DI: FUNCTION TO COMPUTE DI BETWEEN COLUMNS OF X
#   INPUTS:
#       X: Input data with dim (sample)x(channel)
#       M: History length, i.e. number of past samples to use
#       B: Number of bootstrap iterations to use for training classifiers
#   OUTPUTS:
#       DI_estimate: Estimate of the DI from rows to columns. Shape is:
#                     (channel)x(channel)
################################################################################
def DI(X,M,B):

    # MAKE SURE INTEGERS
    M = int(M)
    B = int(B)

    # REMOVE EXCESS SAMPLES AT END
    num_channels = X.shape[1]
    num_samples = X.shape[0]
    num_windows = int(np.floor(num_samples/(M+1)))
    num_samples_to_keep = int(num_windows*(M+1))
    X_trim = X[:num_samples_to_keep,:]

    # RESHAPE SO THAT MATRIX HAS DIM:
    #   (samples)x(X1(-M+i),X1(-M+1+i),...,X1(i),X2(-M+1),X2(-M+1+i),...)
    # WHERE EACH SAMPLE CORRESPONDS TO A WINDOW
    X_past_win = np.zeros((num_windows,num_channels*M))
    X_current  = np.zeros((num_windows,num_channels))
    for chan in range(num_channels):
        current_chan_mat = np.reshape(X_trim[:,chan],(num_windows,(M+1)))
        X_past_win[:,(chan*M):((chan+1)*M)] = current_chan_mat[:,:M]
        X_current[:,chan] = current_chan_mat[:,M]

    # ESTIMATE DI
    chan_list = np.arange(num_channels).tolist()
    chan_list_pairs = list(itertools.product(chan_list, chan_list))

    pool = Pool()
    func = partial(pair_DI,X_past_win,X_current,M,B)
    DI_estimate_list = pool.map(func,chan_list_pairs)

    # UNWRAP
    DI_estimate = np.zeros((num_channels,num_channels))
    for ii in range(num_channels*num_channels):
        DI_estimate[chan_list_pairs[ii][0],chan_list_pairs[ii][1]] = DI_estimate_list[ii]

    # RETURN DI ESTIMATE
    return DI_estimate
################################################################################
# END OF DI
################################################################################





################################################################################
# partial_corr: FUNCTION TO PERFORM SIGN INFERENCE OF RELATIONSHIPS BETWEEN
#                 COLUMNS OF X
#   INPUTS:
#       X: Input data with dim (sample)x(channel)
#       M: History length, i.e. number of past samples to use
#   OUTPUTS:
#       DI_estimate: Estimate of the DI from rows to columns. Shape is:
#                     (channel)x(channel)
################################################################################
def partial_corr(X):

    # COMPUTE CORRELATION COEFFICIENT MATRIX
    corr_matrix = np.corrcoef(X.T)

    # COMPUTE PRECISION MATRIX
    precision_matrix = np.linalg.inv(corr_matrix)

    # NORMALIZE PRECISION MATRIX BY DIAGONAL VALUES TO GET PARTIAL CORRELATIONS
    precision_matrix_diagonal = np.diag(precision_matrix)
    precision_col = np.matlib.repmat(precision_matrix_diagonal,corr_matrix.shape[0],1)
    precision_row = precision_col.T
    normalization_matrix = np.sqrt(precision_row*precision_col)
    partial_corr_matrix = -precision_matrix/normalization_matrix

    # RETURN PARTIAL & REGULAR CORRELATION MATRICES
    return partial_corr_matrix, corr_matrix

################################################################################
# END OF partial_corr
################################################################################





################################################################################
# sign_inference: FUNCTION TO PERFORM SIGN INFERENCE OF RELATIONSHIPS BETWEEN
#                 COLUMNS OF X
#   INPUTS:
#       X: Input data with dim (sample)x(channel)
#       M: History length, i.e. number of past samples to use
#   OUTPUTS:
#       DI_estimate: Estimate of the DI from rows to columns. Shape is:
#                     (channel)x(channel)
################################################################################
def sign_inference(X,M):

    # MAKE SURE INTEGERS
    M = int(M)

    # PARAMETER(S)
    num_channels = X.shape[1]

    # LOOP THROUGH TAUS & CHANNELS
    X_pcorr_sign  = np.zeros((num_channels,num_channels))
    X_corr_sign   = np.zeros((num_channels,num_channels))
    current_pcorr = np.zeros((num_channels,num_channels,M))
    current_corr  = np.zeros((num_channels,num_channels,M))
    for chan2 in range(num_channels):
        for tt in range(M):
            X_copy = copy.deepcopy(X)
            X_copy[:(-(tt+1)),chan2] = X_copy[(tt+1):,chan2]
            X_copy = X_copy[:(-(tt+1)),:]
            current_pcorr[:,:,tt], current_corr[:,:,tt] = partial_corr(X_copy)
        # EXTRACT RELEVANT COLUMN
        current_pcorr_col = np.squeeze(current_pcorr[:,chan2,:])
        current_corr_col  = np.squeeze(current_corr[:,chan2,:])

        # DEBUG: PRINT
        print(chan2)
        print(' ')
        print(current_pcorr_col)
        print(' ')
        print(current_corr_col)

        # TAKE MAX ABS VAL ACROSS TAUS
        max_pcorr_val = np.zeros((num_channels,))
        max_corr_val  = np.zeros((num_channels,))
        for chan1 in range(num_channels):
            max_pcorr_val[chan1] = current_pcorr_col[chan1,np.argmax(np.absolute(current_pcorr_col[chan1,:]))]
            max_corr_val[chan1]  = current_corr_col[chan1,np.argmax(np.absolute(current_corr_col[chan1,:]))]

        # DEBUG: PRINT
        print(' ')
        print(max_pcorr_val)
        print(max_corr_val)
        print(' ')
#        if 1:# length of vals more than 1, then check equality
#            1
        X_pcorr_sign[:,chan2] = np.sign(max_pcorr_val)
        X_corr_sign[:,chan2]  = np.sign(max_corr_val)

    # RETURN CORR & PCORR SIGNS
    return X_pcorr_sign, X_corr_sign

################################################################################
# END OF DI
################################################################################
