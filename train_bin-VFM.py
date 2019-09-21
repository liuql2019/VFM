#!/usr/bin/python

import argparse
import re
import os
import VFM_train as vt

# Parameters
parser = argparse.ArgumentParser(description='Train a model for predicting phage draft genomes in metagenomic bins.')
parser.add_argument('-vir', action="store", required=True, dest="vir_dir", default=None, help='The directory of the virus bins for training.')
parser.add_argument('-bac', action="store", required=True, dest="bac_dir", default=None, help='The directory of the bacterium bins for training.')
parser.add_argument('-cpus', action="store", required=False, dest="threads", default=1, help='The thread number.')
parser.add_argument('-model', action="store", required=True, dest="model_name", help='The model name.')

# Parse the args
args = parser.parse_args()
data_dir = [args.bac_dir, args.vir_dir]
threads = args.threads
if not re.search('\.m$', args.model_name):
    model_name = args.model_name+".m"
    features_csv = args.model_name+".csv"
else:
    model_name = args.model_name
    (filepath, tempfilename) = os.path.split(model_name)
    (filename, extension) = os.path.splitext(tempfilename)
    features_csv = filename+".csv"

vt.VFM_train(model_name, features_csv, threads, data_dir)
