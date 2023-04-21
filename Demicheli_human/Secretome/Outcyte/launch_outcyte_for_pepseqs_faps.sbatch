#!/bin/bash

workdir="Secretome_in_Single_Cell_RNAseq_skeletal_muscle"
cd ${workdir}
executable_dir="${workdir}/outcyte"
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

# module load python/3.7 -> if you work in HPC-environment
outcyte_dir=${workdir}/DeMicheli_human/Secretome/Outcyte/
results_dir=${outcyte_dir}/Results/
mkdir -p ${results_dir}

###
###> RUN
###

## [cmd]: python run_outcyte.py  /path/to/your/fasta/file  </outcyte-sp/outcyte-ups/outcyte>

pepseqs_1=${outcyte_dir}/PepSeq_For_Outcyte/PepSeq_COL1A1_Fibro_Genes_human.fasta
pepseqs_2=${outcyte_dir}/PepSeq_For_Outcyte/PepSeq_FBN1_Fibro_Genes_human.fasta
pepseqs_3=${outcyte_dir}/PepSeq_For_Outcyte/PepSeq_MYOC_Fibro_Genes_human.fasta

python ${executable_dir}/run_outcyte.py ${pepseqs_1} outcyte  ${results_dir}
python ${executable_dir}/run_outcyte.py ${pepseqs_2} outcyte  ${results_dir}
python ${executable_dir}/run_outcyte.py ${pepseqs_3} outcyte  ${results_dir}


## Note : results saved in "results/" folder, by default
## run_outcyte.py --> modified original script by adding a custom results_folder in user-input-arguments
