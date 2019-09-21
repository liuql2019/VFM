#!/usr/bin/python

# VFM predicts unbinned phages or bacteria

import os
import re
import argparse
import numpy as np
import pandas as pd
import methods as ms


# VFM predicting method for unbinned contigs
def VFM_unbin_predict(fa_folder, out_dir, out_file, dims, threads, model_names=None, training_csv=None):
    if threads == None:
        threads = 1
    contig_ids = [[], [], []]
    contig_features = [[], [], []]
    contig_lengths = [[], [], []]

    contig_infor = ms.test_preprocess(fa_folder, out_dir, True, training_csv, dims, threads)
    lengths = contig_infor["hybrid"][3]
    # Distribute contig information by lengths
    for index, seq_len in enumerate(lengths):
        contig_id = contig_infor["hybrid"][0][index]
        features = contig_infor["hybrid"][1][index]
        # if seq_len < ms.lens[0]:
        #     print("Contig "+contig_infor["hybrid"][0][index]+" shorter than "+str(ms.lens[0])+" is not processed.")
        #     continue
        if 0 <= seq_len < ms.lens[1]:
            contig_ids[0] += [contig_id]
            contig_features[0] += [features]
            contig_lengths[0] += [seq_len]
        elif ms.lens[1] <= seq_len < ms.lens[2]:
            contig_ids[1] += [contig_id]
            contig_features[1] += [features]
            contig_lengths[1] += [seq_len]
        elif ms.lens[2] <= seq_len:
            contig_ids[2] += [contig_id]
            contig_features[2] += [features]
            contig_lengths[2] += [seq_len]
    if model_names == None:
        model_names = ["2k_model_unbin.m", "4k_model_unbin.m", "8k_model_unbin.m"]
        user = False
    else:
        user = True
    prob = [[], [], []]
    pred = [[], [], []]
    for index, features in enumerate(contig_features):
        model = ms.load_model(model_names[index], user)
        if features != []:
            features = np.array(features)
            prob[index] = model.predict_proba(features)
            pred[index] = prob[index][:, 1] >= 0.7
    print("The result is as follows:")
    res = show_contig_result(contig_ids, contig_lengths, pred, prob)
    print("The columns mean contig ID, the length of the contig, phage or not and probability to be a phage respectively.")
    res.to_csv(os.path.join(out_dir, out_file))
    print("The result has been saved in "+out_dir)


# Show the prediction result
def show_contig_result(contig_ids, contig_lengths, contig_pred, contig_prob):
    ids = []
    lens = []
    pred = []
    prob = []
    for id_group in contig_ids:
        ids += id_group
    for len_group in contig_lengths:
        lens += len_group
    for pred_group in contig_pred:
        pred += list(pred_group)
    for prob_group in contig_prob:
        prob += list(prob_group)

    res = pd.DataFrame(ids, columns=["ID"])
    res["Length"] = lens
    pred = ["Y" if val else "N" for val in pred]
    res["Phage"] = pred
    res["Score"] = np.array(prob)[:, 1]
    print(res)
    print("\n")
    return res


# Parse the param
parser = argparse.ArgumentParser(description='Predict phage contigs in a metagenome.')
parser.add_argument('-f', action="store", required=True, dest="fa_file", help='The fasta file for predicting.')
parser.add_argument('-t', action="store", required=False, dest="threads", default=1, help='The CPU number for running this program.')
parser.add_argument('-m', action="store", required=False, dest="model", help='The selected model for prediction.')
args = parser.parse_args()
fa_file = args.fa_file
threads = args.threads
if args.model != None:
    if re.search('\.m$', args.model):
        model_id = re.findall('^(.+)\.m$', args.model)[0]
    else:
        model_id = args.model
    model_names = [model_id+"_2k_model_unbin.m", model_id+"_4k_model_unbin.m", model_id+"_8k_model_unbin.m"]
    training_csv = [model_id+"_training_unbin_2k.csv", model_id+"_training_unbin_4k.csv", model_id+"_training_unbin_8k.csv"]
else:
    model_names = None
    training_csv = None
dims = 24

# For release
fa_folder = ms.split_fasta(fa_file)
(filepath, tempfilename) = os.path.split(fa_file)
out_file = ms.get_prefix(fa_file)+"_result.csv"
VFM_unbin_predict(fa_folder, filepath, out_file, dims, threads, model_names, training_csv)
