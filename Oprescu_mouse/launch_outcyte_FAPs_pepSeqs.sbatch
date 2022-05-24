#!/bin/bash
#SBATCH -A single_cell_skeletal_muscle
#SBATCH --mem 8G
#SBATCH -o slurm.%j.out



executable_dir="/home/Tools_Packages/Local_Install_OutCyte"
workdir="/projects/single_cell_skeletal_muscle/Oprescu"
cd ${workdir}

module load conda
CONDA_ROOT="/shared/software/miniconda"
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate /shared/home/mkondili/.conda/envs/outcyte

module load python/3.7

fasta_dir=${workdir}/Secretome/PepSeq_For_Outcyte
results_dir=${workdir}/Secretome/Outcyte_Results/

mkdir -p ${results_dir}

## TO RUN:python run_outcyte.py  /path/to/your/fasta/file  /outcyte-sp/or/outcyte-ups/or/outcyte


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


python ${executable_dir}/run_outcyte_edMK.py ${pepseqs_faps}  outcyte  ${results_dir}

python ${executable_dir}/run_outcyte_edMK.py ${pepseqs_fibro}  outcyte  ${results_dir}


## results saved in "results/" folder, by default
## Note: run_outcyte.py --> modified to _edMK, by adding a custom results_folder in user-input-arguments
