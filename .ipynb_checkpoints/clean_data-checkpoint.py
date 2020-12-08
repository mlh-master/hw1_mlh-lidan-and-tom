# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 17:14:23 2019

@author: smorandv
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def rm_ext_and_nan(CTG_features, extra_feature):
    """

    :param CTG_features: Pandas series of CTG features
    :param extra_feature: A feature to be removed
    :return: A dictionary of clean CTG called c_ctg
    """
#     ------------------ IMPLEMENT YOUR CODE HERE:------------------------------
#     c_ctg = CTG_features.apply(pd.to_numeric, errors='coerce').to_dict('dict')
#     c_ctg = {k: {a: b for a, b in v.items() if not np.isnan(b)} for k, v in CTG_features.loc[:, CTG_features.columns != extra_feature].apply(pd.to_numeric, errors='coerce')} #.to_dict('series').items()
#     CTG_features.loc[:, CTG_features.columns != extra_feature]
    # c_ctg = {k: {a: b for a, b in v.items() if not np.isnan(b)} for k, v in c_ctg.items()}
    
    #all above are code pices used to try it in one row
    
    CTG_features=CTG_features.drop(extra_feature, axis=1).apply(pd.to_numeric, errors='coerce')
    c_ctg = {}
    keys = CTG_features.columns;
    for i in keys:
        c_ctg[i] = CTG_features[i].dropna()
    return c_ctg
#     --------------------------------------------------------------------------


def nan2num_samp(CTG_features, extra_feature):
    """

    :param CTG_features: Pandas series of CTG features
    :param extra_feature: A feature to be removed
    :return: A pandas dataframe of the dictionary c_cdf containing the "clean" features
    """
    c_cdf = {}
    # ------------------ IMPLEMENT YOUR CODE HERE:------------------------------
    CTG_features=CTG_features.drop(extra_feature, axis=1).apply(pd.to_numeric, errors='coerce')
    keys = CTG_features.columns;
    for i in keys:
        c_cdf[i] = CTG_features[i].fillna(pd.Series(np.random.choice(CTG_features[i], size=len(CTG_features.index))))
    # -------------------------------------------------------------------------
    return pd.DataFrame(c_cdf)


def sum_stat(c_feat):
    """

    :param c_feat: Output of nan2num_cdf
    :return: Summary statistics as a dicionary of dictionaries (called d_summary) as explained in the notebook
    """
    # ------------------ IMPLEMENT YOUR CODE HERE:------------------------------
    d_summary = {}
    for x in c_feat.columns:
        tmp = c_feat.describe()[x].to_dict()
        d_summary[x] =  {'min': tmp['min'], 'Q1': tmp['25%'], 'median': tmp['50%'], 'Q3': tmp['75%'], 'max': tmp['max']  }  
    # -------------------------------------------------------------------------
    return d_summary


def rm_outlier(c_feat, d_summary):
    """

    :param c_feat: Output of nan2num_cdf
    :param d_summary: Output of sum_stat
    :return: Dataframe of the dictionary c_no_outlier containing the feature with the outliers removed
    """
    c_no_outlier = {}
    # ------------------ IMPLEMENT YOUR CODE HERE:------------------------------
    column_names = c_feat.columns;

    print ()
    
    for i in column_names:
        Q1_tmp = d_summary[i]['Q1']
        Q3_tmp = d_summary[i]['Q3']
        IQR_tmp = Q3_tmp - Q1_tmp
        low_thresh = Q1_tmp - 1.5 * IQR_tmp
        up_thresh = Q3_tmp + 1.5 * IQR_tmp
        no_low_out = c_feat[i][c_feat[i] > low_thresh]
        no_up_ou = no_low_out[c_feat[i] < up_thresh]
        c_no_outlier[i] = no_up_ou
    # -------------------------------------------------------------------------
    return pd.DataFrame(c_no_outlier)


def phys_prior(c_cdf, feature, thresh):
    """

    :param c_cdf: Output of nan2num_cdf
    :param feature: A string of your selected feature
    :param thresh: A numeric value of threshold
    :return: An array of the "filtered" feature called filt_feature
    """
    # ------------------ IMPLEMENT YOUR CODE HERE:----------------------------- 
    thresh_low = thresh[0]
    thresh_up = thresh[1]
    no_low_outlier = c_cdf[feature][c_cdf[feature] > thresh_low]
    filt_feature = no_low_outlier[c_cdf[feature] < thresh_up]#filt_feature = no_up_outlier
    # -------------------------------------------------------------------------
    return filt_feature


def norm_standard(CTG_features, selected_feat=('LB', 'ASTV'), mode='none', flag=False):
    """

    :param CTG_features: Pandas series of CTG features
    :param selected_feat: A two elements tuple of strings of the features for comparison
    :param mode: A string determining the mode according to the notebook
    :param flag: A boolean determining whether or not plot a histogram
    :return: Dataframe of the normalized/standardazied features called nsd_res
    """
    x, y = selected_feat
    # ------------------ IMPLEMENT YOUR CODE HERE:------------------------------

    if mode is "standard":
        Mean=np.mean(CTG_features)
        Std= np.std(CTG_features)
        nsd_res=(CTG_features-Mean)/Std

    elif mode is 'MinMax':
        Min=np.amin(CTG_features)
        Max=np.amax(CTG_features)
        nsd_res = (CTG_features - Min) / (Max - Min)

    elif (mode == 'mean'):

        Mean = np.mean(CTG_features)
        Min=np.amin(CTG_features)
        Max=np.amax(CTG_features)
        nsd_res = (CTG_features-Mean)/(Max - Min)

    elif (mode == 'none'):
        nsd_res = pd.DataFrame(CTG_features)

    else: print('Choose correct mode')

    if flag:
        axarr = nsd_res.hist(column=[x,y], bins=100,layout = (1, 2),figsize=(20, 10))
        for i, ax in enumerate(axarr.flatten()):
            ax.set_xlabel('Values')
            ax.set_ylabel("Counts")
            
        plt.show()
    
    # -------------------------------------------------------------------------
    return pd.DataFrame(nsd_res)
