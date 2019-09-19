Mine phages from metagenomes in the form of bins and contigs  
============================================================
VFM is developed for identifing phages from metagenomic bins or contigs, which is designed in two versions： bin-VFM for metagenomic bins and unbin-VFM for metagenomic contigs.  
This project is composed of predicting and training scripts.The details are as follows.
* Predicting scripts - Predict bins or contigs as phages or bacteria    
bin-VFM_predict.py    
unbin-VFM_pred.py
* Training scripts - Train a new model by user's own data  
train_bin-VFM.py  
train_unbin-VFM.py  
## Requirement and Dependency
The system must be Linux, with Python3 installed on. Some python packages and modules should be installed in this way:  
```
pip3 install numpy pandas scipy biopython scikit-learn
``` 
VFM depends on some bioinformatic tools: 
* [Prokka](https://github.com/tseemann/prokka) - For gene prediction
* [HMM tools](http://www.hmmer.org/) - For hidden Markov models associated with gene/protein families
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - For sequence alignments

## Setup
VFM can be set up by:
```
git clone https://github.com/liuql2019/VFM
```
Pressing the download button can also fix it.

## Usage
To predct metagenomic bins as phages or bacteria, bin-VFM is used as follows:
```
python3 bin-VFM_predict.py -d BINS_DIR [-t THREADS] [-m MODEL]
```
where the parameter -d means the directory of the bins, -t means cpu number and -m means user-trained model for predicting. If the parameter -t is omitted, default value 1 will be used. Omitting parameter -m means choosing the default model that has been stored in VFM before release.  
Similarly, unbin-VFM is used to predict metagenomic contigs by:
```
python3 unbin-VFM_pred.py -f FA_FILE [-t THREADS] [-m MODEL]
```
where the parameter -f means the fasta file in which the contigs are stored. The meanings of -t and -m are the same as above. 
## Demo
The data in the folder data/test will be used to demenstrate how to run VFM.  
Run bin-VFM as follows:
```
python3 bin-VFM_predict.py -d ./data/test/bins -t 4
```
Run unbin-VFM as follows:
```
python3 unbin-VFM_pred.py -f ./data/test/contigs.fasta -t 4
```
## Train user's model
Users could train their own model(s) by scripts train_bin-VFM.py and train_unbin-VFM.py.   
Run train_bin-VFM.py as follows:
```
python3 train_bin-VFM.py -vir VIR_DIR -bac BAC_DIR [-cpus THREADS] -model MODEL_NAME
```
where the parameters -vir and -bac are folders of phage and bacterium bins for training, -model is the model name. The paremeter -cpus is cpu number which may be omitted.  
Use train_unbin-VFM.py as follows:
```
python3 train_unbin-VFM.py -vir VIR_FILE -bac BAC_FILE [-cpus THREADS] -model MODEL_NAME
```
where the parameters -vir and -bac are fasta files of phage and bacterium contigs for training.The other parameters are the same as above.  
In order to generate metagenomic bins, binning tools such as [COCACOLA](https://github.com/younglululu/COCACOLA) should be applied to metagenomic contigs.Then run script bin-VFM_predict.py to predict the bins.
