#!/usr/bin/python

# Train the random forest model with  unbinned contigs longer than 2kbp by 24D features

import argparse
import os
import re
import methods as ms
import VFM_train as vt

parser = argparse.ArgumentParser(description='Train a model for predicting phage draft genomes in metagenomic contigs.')
parser.add_argument('-vir', action="store", required=True, dest="vir_file", default=None, help='The fasta file of virus contigs for training.')
parser.add_argument('-bac', action="store", required=True, dest="bac_file", default=None, help='The fasta file of bacterium contigs for training.')
parser.add_argument('-cpus', action="store", required=False, dest="threads", default=1, help='The thread number for computing')
parser.add_argument('-model', action="store", required=True, dest="model_name", help='Model name.')

args = parser.parse_args()
bac_file = args.bac_file
vir_file = args.vir_file
threads = args.threads

# Params
if re.search('\.m$', args.model_name):
    model_id = re.findall('^(.+)\.m$', args.model_name)[0]
else:
    model_id = args.model_name
model_name = [model_id+"_2k_model_unbin.m", model_id+"_4k_model_unbin.m", model_id+"_8k_model_unbin.m"]
unbin_training_csv =[model_id+"_training_unbin_2k.csv", model_id+"_training_unbin_4k.csv", model_id+"_training_unbin_8k.csv"]
lengths = [2000, 4000, 8000]

for index, length in enumerate(lengths):
    print("Generating "+str(length)+" bp training set...")
    bac_temp_folder, vir_temp_folder = ms.split_data(bac_file, vir_file, length)

    vir_dir, vir_file_name = os.path.split(vir_file)
    vir_folder = os.path.join(vir_dir, ms.get_prefix(vir_file)+"_"+str(length))
    bac_dir, bac_file_name = os.path.split(bac_file)
    bac_folder = os.path.join(bac_dir, ms.get_prefix(bac_file)+"_"+str(length))

    ms.sample_fasta(vir_temp_folder, 1, vir_folder)
    ms.balance(bac_temp_folder, vir_folder, bac_folder)

    data_dir = [bac_folder, vir_folder]
    vt.VFM_train(model_name[index], unbin_training_csv[index], threads, data_dir)
