#!/bin/bash

workdir="Secretome_in_Single_Cell_RNAseq_skeletal_muscle"
executable_dir="${workdir}/outcyte_pkg"

cd ${executable_dir}

git clone https://github.com/mariakondili/Zhao_Outcyte ${executable_dir}
#! Attention on the repo : The original on Linlinzhao/outcyte doesn't contain the modifications needed for the workflow

###
###> Create a conda-env for Outcyte
###

conda create -n outcyte
conda activate outcyte # should be in /home/.conda/envs/

##
##> Define directories
##

fasta_dir=${workdir}/Secretome/Outcyte/PepSeq_For_Outcyte
results_dir=${workdir}/Secretome/Outcyte/Results/
mkdir -p ${results_dir}

###
###> RUN
###

#python run_outcyte.py  /path/to/your/fasta/file  /outcyte-sp/outcyte-ups/outcyte

##  mouse-genes Pep-seqs

pepseqs_faps="${fasta_dir}/PepSeq_FAPsMarkers_Oprescu.fasta"
pepseqs_fibro="${fasta_dir}/PepSeq_FibroMarkers_Oprescu.fasta"


## ATTENTION:
##> Have to eliminate "*" from end of certain sequences.fasta,like this (in separ.script):
##for f in  ${pepseqs_mm1[@]}; do
##	new_f=${f%".fasta"}
##	sed 's/*$//g' ${f} > ${new_f}_clean.fasta
##done
##> then :
## rm ${f};
## mv ${new_f}_clean.fasta ${f}


python run_outcyte.py ${pepseqs_faps}  outcyte  ${results_dir}

python run_outcyte.py ${pepseqs_fibro}  outcyte  ${results_dir}


## results saved in "results/" folder, by default
## Note: run_outcyte.py --> modified, by adding a custom results_folder in user-input-arguments
