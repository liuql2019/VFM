#!/usr/bin/python

# VFM predicts bins as phages or bacteria

import os
import re
import argparse
import numpy as np
import pandas as pd
import methods as ms


# Print results on the screen
def show_bin_result(id, seq_num, len_ave, label, prob):
    pd.set_option('max_rows', None)
    pd.set_option('max_columns', None)
    res = pd.DataFrame(id, columns=["ID"])
    res["Seq Num"] = seq_num
    res["Len Ave"] = len_ave
    pred = ["Y" if val else "N" for val in label]
    res["Phage"] = pred
    res["Score"] = prob
    print(res)
    print("\n")
    return res


# VFM predicting program
def VFM_predict(bins_dir, out_dir, dims, threads, model_name=None, training_csv=None):
    if threads == None:
        threads = 1
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    bin_infor = ms.test_preprocess(bins_dir, out_dir, False, training_csv, dims, threads)
    # load model
    kmer_model = ms.load_model("LR_0.5k_4mer.m")
    if model_name is None:
        hybrid_model = ms.load_model("model_24d_bins.m")
    else:
        hybrid_model = ms.load_model(model_name, True)
    # result
    kmer_pred, hybrid_pred, kmer_prob, hybrid_prob = ([], [], [], [])
    if bin_infor["kmer"][0] != []:
        kmer_prob = kmer_model.predict_proba(bin_infor["kmer"][1])
        kmer_pred = kmer_prob[:, 1] >= 0.7
    if bin_infor["hybrid"][0] != []:
        hybrid_prob = hybrid_model.predict_proba(bin_infor["hybrid"][1])
        hybrid_pred = hybrid_prob[:, 1] >= 0.7

    # show and save the result
    ids = bin_infor["kmer"][0] + bin_infor["hybrid"][0]
    seq_nums = bin_infor["kmer"][2] + bin_infor["hybrid"][2]
    len_aves = bin_infor["kmer"][3] + bin_infor["hybrid"][3]
    pred = np.append(kmer_pred, hybrid_pred)
    if kmer_prob == []:
        prob = hybrid_prob
    elif hybrid_prob == []:
        prob = kmer_prob
    else:
        prob = np.append(kmer_prob, hybrid_prob, axis=0)
    print("The result is as follows:")
    res = show_bin_result(ids, seq_nums, len_aves, pred, prob[:, 1])
    print("The columns mean bin ID, sequnece number in the bin, length average for contigs in the bin, phage or not and probability to be a phage bin respectively.")
    res.to_csv(os.path.join(out_dir, "bin-VFM_result.csv"))
    print("The result has been stored in the folder "+os.path.join(bins_dir, "pred_results")+".")


#####
# Predict bins
#####

# Parse the param
parser = argparse.ArgumentParser(description='Predict phage draft genomes in a binned metagenome.')
parser.add_argument('-d', action="store", required=True, dest="bins_dir", help='The folder of the bins to be predicted.')
parser.add_argument('-t', action="store", required=False, dest="threads", help='The CPU number to run program.')
parser.add_argument('-m', action="store", required=False, dest="model", help='The model to be used for prediction.')

args = parser.parse_args()
bins_dir = args.bins_dir
threads = args.threads
if args.model != None:
    if not re.search('\.m$', args.model):
        model = args.model+".m"
        training_features = args.model+".csv"
    else:
        model = args.model
        (filepath, tempfilename) = os.path.split(model)
        (filename, extension) = os.path.splitext(tempfilename)
        training_features = filename+".csv"
else:
    # use original model
    model = None
    training_features = None

dims = 24
out_dir = os.path.join(bins_dir, 'pred_results/')
VFM_predict(bins_dir, out_dir, dims, threads, model, training_features)
