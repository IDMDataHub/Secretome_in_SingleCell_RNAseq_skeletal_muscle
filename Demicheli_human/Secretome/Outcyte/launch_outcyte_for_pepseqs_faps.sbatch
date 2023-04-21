#!/bin/bash
#SBATCH -A single_cell_skeletal_muscle
#SBATCH -p fast
#SBATCH --cpus-per-task 4
#SBATCH -o slurm.outcyte.%j.out

executable_dir="/home/Tools_Packages/Local_Install_OutCyte"
#git clone https://github.com/mariakondili/Zhao_Outcyte ${executable_dir}

workdir="Secretome_in_Single_Cell_skeletal_muscle/DeMicheli_human/Secretome"
cd ${workdir}

#> Create a conda-env for Outcyte
module load conda
CONDA_ROOT="/shared/software/miniconda"
source ${CONDA_ROOT}/etc/profile.d/conda.sh

conda create -n outcyte
conda activate /home/.conda/envs/outcyte

module load python/3.7

results_dir=${workdir}/Results/
mkdir -p ${results_dir}

## TO RUN:python run_outcyte.py  /path/to/your/fasta/file  /outcyte-sp/or/outcyte-ups/or/outcyte
# pepseqs=${workdir}/PepSeq_For_Outcyte/PepSeq_Fibroblasts_Markers_human.fasta

pepseqs_1=${workdir}/PepSeq_For_Outcyte/PepSeq_COL1A1_Fibro_Genes_human.fasta
pepseqs_2=${workdir}/PepSeq_For_Outcyte/PepSeq_FBN1_Fibro_Genes_human.fasta
pepseqs_3=${workdir}/PepSeq_For_Outcyte/PepSeq_MYOC_Fibro_Genes_human.fasta

#python ${executable_dir}/run_outcyte_edMK.py ${pepseqs}   outcyte  ${results_dir}
python ${executable_dir}/run_outcyte_edMK.py ${pepseqs_1} outcyte  ${results_dir}
python ${executable_dir}/run_outcyte_edMK.py ${pepseqs_2} outcyte  ${results_dir}
python ${executable_dir}/run_outcyte_edMK.py ${pepseqs_3} outcyte  ${results_dir}


## Note : results saved in "results/" folder, by default
## run_outcyte.py --> modified to _edMK, by adding a custom results_folder in user-input-arguments
