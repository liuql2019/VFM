#!/usr/bin/python

import os
import subprocess
from collections import Counter
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.externals import joblib
import pickle
import random
import shutil
from ctypes import *  # package for c/c++ calling


# Sequence lengths
lens = [2000, 4000, 8000]

# Training features stored in the form of csv
unbin_training_csv =["training_unbin_2k.csv", "training_unbin_4k.csv", "training_unbin_8k.csv"]
bin_training_csv = "training_data.csv"

# kmer counter
kmer_counter = "count_k-mer.so"
dll_dir, dll_file = os.path.split(os.path.abspath(__file__))
dll = CDLL(dll_dir + "/" + kmer_counter)


# Calculate dsDNA kmer number
def kmer_num(k):
    dll.kmerNum.restype = c_ulong
    res = dll.kmerNum(k)
    return res


# Pointer to kmer freqs for c++ function return
class StructPointer(Structure):
    _fields_ = [("kmer_freq", c_double * 65536)]


# Calculate dsDNA kmer frequencies
def ds_kmer_frequency(ds_dna, k):
    kmer_dim = kmer_num(k)
    dll.countSeqFeatureCpp.restype = POINTER(StructPointer)
    kmer = dll.countSeqFeatureCpp(ds_dna.encode(), k)
    res = []
    counter = 0
    while counter < kmer_dim:
        res = res + [kmer.contents.kmer_freq[counter]]
        counter += 1
    return res


# Calculate kmer frequency dimension
def kmer_freq_dim(k):
    if k % 2 != 0:
        res = int(pow(4, k)/2)
    else:
        res = int(pow(4, k/2) + (pow(4, k) - pow(4, k/2))/2)
    return res


# Calculate ssDNA kmer frequencies
def kmer_frequency(dna, kmer):
    freq = {}
    freq.update(Counter(chunks(dna, 3)))  # update with counter of dinucleotide
    if kmer in freq:
        for each in freq:
            freq[each] = freq[each] / len(dna)
        return (freq[kmer])
    else:
        return 0


# Seq num in a bin
def seq_num(bin_dir, binn):
    records = list(SeqIO.parse(os.path.join(bin_dir, binn), "fasta"))
    return len(records)


# The average length of sequences in the bin
def len_ave(bin_dir, binn):
    lengths = []
    for record in SeqIO.parse(os.path.join(bin_dir, binn), "fasta"):
        lengths += [len(record.seq)]
    return np.mean(lengths)


# Auxiliary function to kmer_frequencies
def chunks(l, n):
    for i in range(0, len(l) - (n - 1)):
        yield l[i:i + n]


# Run Prokka
def run_prokka(binn, bins_dir, threads):
    # Bin name
    prefix = get_prefix(binn)
    bin_file = os.path.join(bins_dir, binn)
    out_dir = os.path.join(bins_dir, 'results/prokka/' + prefix)
    command_line = ('prokka --kingdom Viruses --gcode 11 --cpus ' + str(threads) + ' --force --quiet --prefix prokka_results_' + prefix + ' --fast --centre X --compliant --norrna --notrna --outdir ' + out_dir + ' ' + bin_file).split()

    return_code = subprocess.call(command_line, stderr=subprocess.PIPE)
    # Check whether Prokka run smoothly
    if return_code == 1:
        print("Prokka did not end well.\nExiting...")
        quit()


# Calculate coden usage bias
def count_coden_usage(record, single=True):
    gene_length_sum = 0
    number = 0
    # iterate each gene to sum all gene lengths
    for feature in record.features:
        if feature.type == "CDS":  # coding seq
            if single and feature.location.strand != -1:  # ssDNA
                break
            number += 1
            seq = record.seq
            start = int(feature.location.start)
            end = int(feature.location.end)
            if feature.location.strand == -1:
                cds_seq = seq[start:end].reverse_complement().upper()  # Reverse complementary coding sequence
            else:
                cds_seq = seq[start:end].upper()
            gene_length = len(cds_seq) - cds_seq.count('N')
            gene_length_sum += gene_length

    contig_coden_usage = [0, 0, 0, 0, 0]
    # iterate each gene for coden usage bias
    for feature in record.features:
        if feature.type == "CDS":  # coding seq
            if single and feature.location.strand != -1:
                break
            n_TC2 = 0
            n_C2 = 0
            n_AG2 = 0
            n_G2 = 0
            n_N4 = 0
            n_C4 = 0
            n_G4 = 0
            n_A4 = 0
            pC2 = 0
            pG2 = 0
            pC4 = 0
            pG4 = 0
            pA4 = 0
            start = feature.location.start
            end = feature.location.end
            # gene_length = end - start
            # process CDS seq
            seq = record.seq
            if feature.location.strand == -1:
                cds_seq = seq[start:end].reverse_complement().upper()  # Reverse complementary coding sequence
            else:
                cds_seq = seq[start:end].upper()
            pos = 0
            while pos+3 <= len(cds_seq):
                coden = cds_seq[pos:pos+3]
                if(coden == "TTT" or coden == "TTC" or
                   coden == "TTA" or coden == "CTA" or
                   coden == "TTG" or coden == "CTG" or
                   coden == "TAT" or coden == "TAC" or
                   coden == "TGC" or coden == "TGT" or
                   coden == "CAT" or coden == "CAC" or
                   coden == "ATT" or coden == "ATC" or
                   coden == "AAT" or coden == "AAC" or
                   coden == "AGT" or coden == "AGC" or
                   coden == "GAT" or coden == "GAC"):
                    n_TC2 += 1
                    if(coden == "TTC" or
                       coden == "CTA" or
                       coden == "CTG" or
                       coden == "TAC" or
                       coden == "TGC" or
                       coden == "CAC" or
                       coden == "ATC" or
                       coden == "AAC" or
                       coden == "AGC" or
                       coden == "GAC"):
                        n_C2 += 1
                if(coden == "TAA" or coden == "TAG" or coden == "TGA" or
                   coden == "CAA" or coden == "CAG" or
                   coden == "AAA" or coden == "AAG" or
                   coden == "AGA" or coden == "AGG" or
                   # coden == "AAT" or coden == "GAT" or
                   # coden == "AAC" or coden == "GAC" or
                   coden == "GAA" or coden == "GAG" or
                   coden == "TTA" or coden == "TTG"):
                    n_AG2 += 1
                    if(coden == "TAG" or coden == "TGA" or
                       coden == "CAG" or
                       coden == "AAG" or
                       coden == "AGG" or
                       # coden == "GAT" or
                       # coden == "GAC" or
                       coden == "GAG" or
                       coden == "TTG"):
                        n_G2 += 1
                if(coden == "CTT" or coden == "CTC" or coden == "CTA" or coden == "CTG" or
                   coden == "TCT" or coden == "TCC" or coden == "TCA" or coden == "TCG" or
                   coden == "CCT" or coden == "CCC" or coden == "CCA" or coden == "CCG" or
                   coden == "CGT" or coden == "CGC" or coden == "CGA" or coden == "CGG" or
                   coden == "ACT" or coden == "ACC" or coden == "ACA" or coden == "ACG" or
                   coden == "GTT" or coden == "GTC" or coden == "GTA" or coden == "GTG" or
                   coden == "GCT" or coden == "GCC" or coden == "GCA" or coden == "GCG" or
                   coden == "GGT" or coden == "GGC" or coden == "GGA" or coden == "GGG"):
                    n_N4 += 1
                    if(coden == "CTC" or
                       coden == "TCC" or
                       coden == "CCC" or
                       coden == "CGC" or
                       coden == "ACC" or
                       coden == "GTC" or
                       coden == "GCC" or
                       coden == "GGC"):
                        n_C4 += 1
                    if(coden == "CTG" or
                       coden == "TCG" or
                       coden == "CCG" or
                       coden == "CGG" or
                       coden == "ACG" or
                       coden == "GTG" or
                       coden == "GCG" or
                       coden == "GGG"):
                        n_G4 += 1
                    if(coden == "CTA" or
                       coden == "TCA" or
                       coden == "CCA" or
                       coden == "CGA" or
                       coden == "ACA" or
                       coden == "GTA" or
                       coden == "GCA" or
                       coden == "GGA"):
                        n_A4 += 1
                pos += 3
            if n_TC2 != 0:
                pC2 = n_C2 / n_TC2
            if n_AG2 != 0:
                pG2 = n_G2 / n_AG2
            if n_N4 != 0:
                pC4 = n_C4 / n_N4
                pG4 = n_G4 / n_N4
                pA4 = n_A4 / n_N4
            gene_length = len(cds_seq) - cds_seq.count('N')
            length_ratio = gene_length / gene_length_sum
            contig_coden_usage[0] += pC2 * length_ratio
            contig_coden_usage[1] += pG2 * length_ratio
            contig_coden_usage[2] += pC4 * length_ratio
            contig_coden_usage[3] += pG4 * length_ratio
            contig_coden_usage[4] += pA4 * length_ratio
    return contig_coden_usage, gene_length_sum


# Calculate feature averages
def average(feature_mat, weight_mat):
    weight_list = list(weight_mat.T)
    feature_list = list(feature_mat.T)
    count = 0
    res = []
    for weight in weight_list:
        # process weight 0
        if np.sum(weight) == 0:
            res += [np.nan]
        else:
            res += [np.average(feature_list[count], weights=weight)]
        count += 1
    return res


# HMM scan
def run_hmmscan(binn, bins_dir, db_dir, threads):
    # Make hmmscan the directory of hmmscan results
    hmmscan_dir = os.path.join(bins_dir, 'results/hmmscan/')
    if not os.path.exists(hmmscan_dir):
        os.makedirs(hmmscan_dir)
    prefix = get_prefix(binn)
    contain_gene, only_one = contain_genes(bins_dir, binn)
    if contain_gene:
        query_faa = os.path.join(bins_dir, 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa')
        hmm_db = os.path.join(db_dir, 'all_vogs_hmm_profiles_feb2018.hmm')
        out_hmm = os.path.join(hmmscan_dir, prefix + '_hmmscan.out')
        tbl_hmm = os.path.join(hmmscan_dir, prefix + '_hmmscan.tbl')
        command_line_hmmscan = 'hmmscan -o ' + out_hmm + ' --cpu ' + str(threads) + ' --tblout ' + tbl_hmm + ' --noali ' + hmm_db + ' ' + query_faa
        # In case hmmscan returns an error
        try:
            subprocess.call(command_line_hmmscan, shell=True)
        except:
            print('Error when calling HMMscan:', command_line_hmmscan+'\nExiting...')
            quit()


# Run Blastp
def run_blastp(binn, bins_dir, db_dir, db_name, threads):
    bins_dir = os.path.abspath(bins_dir)
    # Change the working directory to the directory of the COG database
    old_dir = os.getcwd()
    os.chdir(db_dir)
    # Make Blastp directory of Blastp results
    blastp_dir = os.path.join(bins_dir, 'results/blastp')
    if not os.path.exists(blastp_dir):
        os.makedirs(blastp_dir)
    # Construct the blastp command line
    prefix = get_prefix(binn)
    contain_gene, only_one = contain_genes(bins_dir, binn)
    if contain_gene:
        bin_faa = os.path.join(bins_dir, 'results/prokka/'+prefix+'/prokka_results_'+prefix+'.faa')
        output = os.path.join(blastp_dir, prefix+'_out.blastp')
        command_line_blastp = "blastp -query "+bin_faa+" -out "+output+" -db "+db_name+" -outfmt 6 -evalue 1e-3 -max_target_seqs 1 -num_threads "+str(threads)
        # In case blastp returns an error
        try:
            subprocess.call(command_line_blastp, shell=True)
        except:
            print('Error when calling blastp.\nExiting...')
            quit()
    os.chdir(old_dir)


# Calculate the proportion of pVOG hits
def VOG_ratio(input_folder, binn):
    # Prefix for naming results
    prefix = get_prefix(binn)

    num_proteins_bin = 0
    with open(os.path.join(input_folder, 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa'), 'r') as faa:
        for line in faa:
            if re.search('^>', line):
                num_proteins_bin += 1
    # Parse hmmscan output files and store pVOG hits of proteins within the binn
    dic_matches = {}
    hmm_res = os.path.join(input_folder, 'results/hmmscan/' + prefix + '_hmmscan.tbl')
    if os.path.exists(hmm_res):
        with open(hmm_res, 'r') as hmmscan_out:
            for line in hmmscan_out:
                match = re.search('^VOG\d+\s+-\s+(\S+)\s+-\s+(\S+)\s+.+$', line)
                if match:
                    if match.group(1) not in dic_matches:
                        dic_matches[match.group(1)] = float(match.group(2))
        # Get the proportion of proteins matching the pVOGs
        i_sig = 0
        for key in dic_matches:
            if dic_matches[key] <= 1e-10:
                i_sig += 1
        res = i_sig / num_proteins_bin
    else:
        res = np.nan
    return res


# Calculate the proportion of COG hits
def COG_ratio(input_folder, binn):
    # Prefix for naming results
    prefix = get_prefix(binn)
    # Protein number of the bin
    num_proteins_bin = 0
    with open(os.path.join(input_folder, 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa'), 'r') as faa:
        for line in faa:
            if re.search('^>', line):
                num_proteins_bin += 1
    # COG hit number of the bin
    num_COG_hit = 0
    blastp_res = os.path.join(input_folder, 'results/blastp/' + prefix + '_out.blastp')
    if os.path.exists(blastp_res):
        with open(blastp_res, 'r') as out_blastp:
            for line in out_blastp:
                if line != '':
                    num_COG_hit += 1
        res = num_COG_hit/num_proteins_bin
    else:
        res = np.nan
    return res


# Extract features from genbank record and prepare the vector of features
def extract_features(record, single=True):

    count = 0
    sum_cds_length = 0
    strand_shift = 0
    non_coding_spacing = []
    # coden usage
    (coden_usage, gene_length_sum) = count_coden_usage(record, single)
    # 1-mer,2-mer frequency
    if not single:
        one_mer = ds_kmer_frequency(str(record.seq), 1)
        two_mer = ds_kmer_frequency(str(record.seq), 2)
    for feature in record.features:
        if feature.type == "CDS":
            if single and feature.location.strand != -1:
                break
            if 'translation' in feature.qualifiers:
                    count += 1
                    start = feature.location.start
                    end = feature.location.end
                    sum_cds_length += (end - start)
                    if count == 1:
                        strand_prev = feature.location.strand
                        end_prev = end
                    else:
                        non_coding_spacing.append(start - end_prev)
                        if not single and strand_prev != feature.location.strand:
                            strand_shift += 1
                    end_prev = end
                    strand_prev = feature.location.strand
            else:
                print('WARNING: Prokka predicted a CDS, but there is no translation. Record ID: ', record.id)
    # initial value
    mean_gene_size = 0
    mean_spacing_size = 0
    density = count / (len(record.seq) / 1000)
    non_coding_percent = 0
    strand_shift_frq = 0
    ATG_freq = 0
    if count != 0:
        mean_gene_size = sum_cds_length / count
        ATG_freq = kmer_frequency(str(record.seq), 'ATG')
    if len(non_coding_spacing) >= 1:
        sum_spacing = 0
        for i in non_coding_spacing:
            sum_spacing += i
        mean_spacing_size = sum_spacing / len(non_coding_spacing)
        non_coding_percent = sum_spacing / len(record.seq)
        strand_shift_frq = strand_shift / len(non_coding_spacing)
    if single:
        features = [mean_gene_size, mean_spacing_size, density, non_coding_percent, ATG_freq] + coden_usage
        return features
    else:
        features = [mean_gene_size, mean_spacing_size, density, strand_shift_frq, ATG_freq] + coden_usage + one_mer + two_mer
        # feature weights
        weights = [count, len(non_coding_spacing), len(record.seq), len(non_coding_spacing), len(record.seq) - 3 + 1] + [gene_length_sum] * 5 + [len(record.seq) - 1 + 1] * kmer_num(1) + [len(record.seq) - 2 + 1] * kmer_num(2)
        return features, weights


# Get prefix from bins
def get_prefix(fa_file):
    # prefix = re.sub('.fasta', '', binn)
    (filepath, tempfilename) = os.path.split(fa_file)
    (filename, extension) = os.path.splitext(tempfilename)
    return filename


# Load data from csv files
def load_data(csv_file, for_train, dims):
    data = np.loadtxt(csv_file, delimiter=',')
    features = data[:, 0:dims]
    if for_train:
        labels = data[:, dims]
        return features, labels
    else:
        return features
    # classes = np.array([data[:, data.shape[1]-3]]).T
    # gene_seq_sum = np.array([data[:, data.shape[1]-2]]).T
    # seq_num = np.array([data[:, data.shape[1]-1]]).T


# Load k-mer data from pkl files
def load_kmer_data(csv_label, seq_length, sparse):
    file = open(csv_label+'_feature_'+seq_length+'.pkl', 'rb')
    data = pickle.load(file)
    if sparse:
        data = data.tocsc()
    features = data[:, 0:data.shape[1] - 2]
    classes = data[:, data.shape[1] - 2]
    seq_num = data[:, data.shape[1] - 1]
    return features, classes, seq_num


# Save predicted and true labels
def save_label(dir, seq_length, true_label, pred_label, pred_prob=None):
    result = pd.DataFrame(true_label, columns=["true_label"])
    result["pred_label"] = pred_label
    result["pred_prob"] = pred_prob
    result.to_csv(dir+"/test_label_"+seq_length+".csv")


# Fill nan values
def fill_nan(data, ref=None):
    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    if ref is not None:
        f = imp.fit(ref)
        res = imp.transform(data)
    else:
        res = imp.fit_transform(data)
    return res


# Fill nan values of s unbin feature vector
def fill_unbin_nan(feature, seq_len, training_csv, dims):
    res = feature
    if 0 <= seq_len < lens[1]:
        if training_csv is None:
            ref = load_data(os.path.join("models/", unbin_training_csv[0]), False, dims)
        else:
            ref = load_data(os.path.join("user_models/", training_csv[0]), False, dims)
        res = fill_nan([feature], ref)
    elif lens[1] <= seq_len < lens[2]:
        if training_csv is None:
            ref = load_data(os.path.join("models/", unbin_training_csv[1]), False, dims)
        else:
            ref = load_data(os.path.join("user_models/",training_csv[1]), False, dims)
        res = fill_nan([feature], ref)
    elif lens[2] <= seq_len:
        if training_csv is None:
            ref = load_data(os.path.join("models/", unbin_training_csv[2]), False, dims)
        else:
            ref = load_data(os.path.join("user_models/", training_csv[2]), False, dims)
        res = fill_nan([feature], ref)
    return res[0]


# Find fasta files in the input directory
def find_fasta(bins_dir):
    list_bins_temp = os.listdir(bins_dir)
    list_bins = []
    for binn in list_bins_temp:
        if not binn.startswith('.'):
            if re.search('.fasta$', binn, re.IGNORECASE):
                list_bins.append(binn)
            elif re.search('.fa$', binn, re.IGNORECASE):
                list_bins.append(binn)
            elif re.search('.fna$', binn, re.IGNORECASE):
                list_bins.append(binn)
    return list_bins


# Split sequences for a certain length and store each one into a fasta file
def split_data(bac_file, vir_file, length):
    print("Splitting the training data...")
    files = [bac_file, vir_file]
    out_dirs = []
    for file in files:
        (filepath, tempfilename) = os.path.split(file)
        prefix = get_prefix(file)
        out_dir = os.path.join(filepath, prefix)+"_temp_"+str(length)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for record in SeqIO.parse(file, "fasta"):
            substr_num = int(len(record)/length)
            start = 0
            for i in range(substr_num):
                end = start + length
                new_id = record.id+"_"+str(i)
                new_seq = record.seq[start:end]
                new_record = SeqRecord(new_seq, id=new_id)
                SeqIO.write(new_record, os.path.join(out_dir, new_record.id+".fasta"), "fasta")
                start = end
        out_dirs += [out_dir]
    return out_dirs


# Sample fasta files by ratio
def sample_fasta(fa_dir, ratio, out_dir):
    fa_list = find_fasta(fa_dir)
    total = len(fa_list)
    num = int(total * ratio)
    index_list = random.sample(range(total), num)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for i in index_list:
        fa_file = fa_list[i]
        out_file = os.path.join(out_dir, fa_file)
        shutil.copy(os.path.join(fa_dir, fa_file), out_file)


# Balance the file number of the source folder with ref_folder
def balance(source_folder, ref_folder, out_folder):
    print("Balancing the number of phage and bacteria...")
    source_files = find_fasta(source_folder)
    ref_files = find_fasta(ref_folder)
    ratio = len(ref_files)/len(source_files)
    sample_fasta(source_folder, ratio, out_folder)


# Split a fasta file into some files, each containing one sequence
def split_fasta(fa_file):
    filepath, tempfilename = os.path.split(fa_file)
    prefix = get_prefix(fa_file)
    out_dir = os.path.join(filepath, prefix)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for record in SeqIO.parse(fa_file, "fasta"):
        SeqIO.write(record, os.path.join(out_dir, record.id+".fasta"), "fasta")
    return out_dir


# Preprocess training data in bin format
def train_preprocess(bins_dir, label, threads):
    # List all bins contained inside the bins directory
    list_bins = find_fasta(bins_dir)
    if list_bins == []:
        print('There is no bin in fasta format inside the input bins\' directory (%s).\nExiting...' % bins_dir)
        quit()
    # Gene predict
    predict_gene(bins_dir, list_bins, threads)
    # Current file path
    current_dir, file_name = os.path.split(os.path.abspath(__file__))
    db_dir = os.path.join(current_dir, 'models')
    # pVOGs scan
    scan_pVOGs(bins_dir, list_bins, db_dir, threads)
    # COGs scan
    COG_db = 'COG_sample_db'
    scan_COGs(bins_dir, list_bins, db_dir, COG_db, threads)
    # Extract features
    print('Extracting features ...\n')
    feature_array = None
    ids = []
    count = 0
    for binn in list_bins:
        prefix = get_prefix(binn)
        ids += [prefix]
        feature_vec = extract_bin_features(bins_dir, binn, only_kmer=False)
        if feature_array is None:
            feature_array = [feature_vec]
        else:
            feature_array += [feature_vec]
        count += 1
        if count % 10 == 0:
            print('Done with %d bins/contigs features extracting.' % count)
    feature_array = np.array(feature_array)
    labels = [label]*feature_array.shape[0]
    labels = np.array([labels]).T
    return ids, feature_array, labels


# Preprocess testing data in the form of bins including one-contig bin(unbinned contig)
def test_preprocess(bins_dir, out_dir, unbin, training_csv, dims, threads):
    # List all bins contained inside the bins directory
    list_bins = find_fasta(bins_dir)
    if list_bins == []:
        print('There is no bin in fasta format inside the input bin directory (%s).\nExiting...' % bins_dir)
        quit()

    # Gene predict
    predict_gene(bins_dir, list_bins, threads)

    # Current file path
    current_dir, file_name = os.path.split(os.path.abspath(__file__))
    db_dir = os.path.join(current_dir, 'models')

    # pVOGs scan
    scan_pVOGs(bins_dir, list_bins, db_dir, threads)

    # COGs scan
    COG_db = 'COG_sample_db'
    scan_COGs(bins_dir, list_bins, db_dir, COG_db, threads)

    # Extract features
    print('Extracting features ...')
    # k-mer feature information
    kmer_ids = []
    kmer_array = None
    kmer_seq_nums = []
    kmer_len_aves = []
    # hybrid feature information
    hf_ids = []
    hf_array = None
    hf_seq_nums = []
    hf_len_aves = []
    count_features = 0

    for binn in list_bins:
        prefix = get_prefix(binn)
        (contain_gene, only_one_gene) = contain_genes(bins_dir, binn)
        if unbin:
            only_kmer = False
        else:
            only_kmer = not contain_gene

        feature_vec = extract_bin_features(bins_dir, binn, only_kmer)
        # Fill nan for the unbinned contig
        if unbin and np.nan in feature_vec:
            feature_vec = fill_unbin_nan(feature_vec, len_ave(bins_dir, binn), training_csv, dims)
        if only_kmer:
            kmer_ids += [prefix]
            kmer_seq_nums += [seq_num(bins_dir, binn)]
            kmer_len_aves += [len_ave(bins_dir, binn)]
            if kmer_array is None:
                kmer_array = [feature_vec]
            else:
                kmer_array = np.concatenate((kmer_array, [feature_vec]), axis=0)
        else:
            hf_ids += [prefix]
            hf_seq_nums += [seq_num(bins_dir, binn)]
            hf_len_aves += [len_ave(bins_dir, binn)]
            if hf_array is None:
                hf_array = [feature_vec]
            else:
                hf_array = np.concatenate((hf_array, [feature_vec]), axis=0)
        count_features += 1
        if count_features % 10 == 0:
            print('Done with %d bins for feature extraction...' % count_features)

    # Fill nan for bins
    if not unbin and np.isnan(hf_array).sum() > 0:
        if training_csv is None:
            training_features = load_data(os.path.join("models/", bin_training_csv), False, dims)
            hf_array = fill_nan(hf_array, training_features)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # store the bin names and features
    kmer_id_file = open(os.path.join(out_dir, 'kmer_id.pkl'), 'wb')
    pickle.dump(kmer_ids, kmer_id_file)
    kmer_id_file.close()
    hf_id_file = open(os.path.join(out_dir, 'hf_id.pkl'), 'wb')
    pickle.dump(hf_ids, hf_id_file)
    hf_id_file.close()
    kmer_file = open(os.path.join(out_dir, 'kmer_features.pkl'), 'wb')
    pickle.dump(kmer_array, kmer_file)
    kmer_file.close()
    hf_file = open(os.path.join(out_dir, 'hf_features.pkl'), 'wb')
    pickle.dump(hf_array, hf_file)
    hf_file.close()
    # print('Features have been generated and stored in the bin folder.\n')
    res = {"kmer": [kmer_ids, kmer_array, kmer_seq_nums, kmer_len_aves],
           "hybrid": [hf_ids, hf_array, hf_seq_nums, hf_len_aves]}
    return res


# Whether contain genes in the bin
def contain_genes(bins_dir, binn):
    contain = False
    max_gene_num = 0
    prefix = get_prefix(binn)
    file_name = os.path.join(bins_dir, 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.gbk')
    try:
        for record in SeqIO.parse(file_name, "genbank"):
            gene_count = 0
            for feature in record.features:
                if feature.type == "CDS":
                    gene_count += 1
            if gene_count > max_gene_num:
                max_gene_num = gene_count
    except:
        print("Fault within "+file_name)
    if max_gene_num > 0:
        contain = True
    if max_gene_num == 1:
        only_one = True
    else:
        only_one = False
    return contain, only_one


# Predict genes
def predict_gene(bins_dir, list_bins, threads):
    count_prokka = 0
    print('Gene prediction starts, please wait...')
    for binn in list_bins:
        run_prokka(binn, bins_dir, threads)
        count_prokka += 1
        if count_prokka % 10 == 0:
            print('%d bins/contigs done...' % count_prokka)
    print('Gene predition has finished.\n')


# Scan pVOGs for bins
def scan_pVOGs(bins_dir, list_bins, db_dir, threads):
    print('Starting pVOGs scan, this may take awhile.')
    count_hmm = 0
    for binn in list_bins:
        run_hmmscan(binn, bins_dir, db_dir, threads)
        count_hmm += 1
        if count_hmm % 10 == 0:
            print('Done with %d bins/contigs for HMM search...' % count_hmm)
    print('Finished scan to the pVOG database.\n')


# Scan COGs for bins
def scan_COGs(bins_dir, list_bins, db_dir, db_name, threads):
    print('Starting COGs scan, please wait awhile...')
    count_blastp = 0
    # makeblastdb -in prot2003-2014.fa -dbtype prot -out COG_db -parse_seqids
    for binn in list_bins:
        run_blastp(binn, bins_dir, db_dir, db_name, threads)
        count_blastp += 1
        if count_blastp % 10 == 0:
            print('Done with %d bins/contigs for blastp search...' % count_blastp)
    print('Finished scan to the COG database.\n')


# Extract the feature vector of a bin
def extract_bin_features(bin_dir, binn, only_kmer):
    k = 4
    sub_data_bin = []
    weights_list = []
    prefix = get_prefix(binn)
    file_name = os.path.join(bin_dir, 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.gbk')
    for record in SeqIO.parse(file_name, "genbank"):
        if only_kmer:
            weights = [len(record) - k + 1] * kmer_freq_dim(k)
            temp_vec = ds_kmer_frequency(str(record.seq), k)
        else:
            (temp_vec, weights) = extract_features(record, single=False)
        sub_data_bin.append(temp_vec)
        weights_list.append(weights)
    bin_features = np.array(sub_data_bin)
    bin_weights = np.array(weights_list)
    bin_features_ave = average(bin_features, bin_weights)
    if not only_kmer:
        prop_hmms_hits = VOG_ratio(bin_dir, binn)
        bin_features_ave.insert(5, prop_hmms_hits)
        prop_blastp_hits = COG_ratio(bin_dir, binn)
        bin_features_ave.insert(6, prop_blastp_hits)
    return bin_features_ave


# Load model
def load_model(model_name, user=False):
    program_dir = os.path.split(os.path.realpath(__file__))[0]
    if not user:
        model_file = program_dir+"/models/"+model_name
    else:
        model_file = program_dir+"/user_models/"+model_name
    model = joblib.load(model_file)
    return model

