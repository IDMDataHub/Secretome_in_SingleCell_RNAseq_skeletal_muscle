# Secretome_in_Single_Cell_skeletal_muscle
In this repository we propose an approach for extracting FAPs (Fibro-adipogenic Progenitors) expression from human & mouse S.C-RNA-seq of muscle from following studies :
  1. [De Micheli et al.,Skeletal Muscle (2020)](https://doi.org/10.1186/s13395-020-00236-3)  &
  2. [Oprescu et al., iScience (2020)](https://doi.org/10.1016/j.isci.2020.100993) \
   in order to study **Virtual Secretome** of FAPs during regeneration and in homeostasis. 
 
Our study is published here : \
[**Negroni et al.,Muscle fibro-adipogenic progenitors from a single-cell perspective:Focus on their "virtual" secretome**]( https://doi.org/10.3389/fcell.2022.952041 )

# Prerequisites 
The scripts were run in R/4.1.1, with Seurat/4.1.1, on CentOS Linux 7(64-bit).

The open-source tools used for the approach :
  - Outcyte : https://github.com/mariakondili/Zhao_Outcyte , a python-based tool, to be run in bash. Further dependencies noted in the repository's README.
  - NicheNet : https://github.com/saeyslab/nichenetr ,as an R-package
    > devtools::install_github("saeyslab/nichenetr") 
    > library(nichenetr)) 
  

The Scripts should be used as follows: 
 > 1. analyse_with_Seurat.Rmd : analyse raw S.C.data published.Extract FAPs cells x FAPs markers table 
 > 2. extract_PepSeqs_for_Outcyte.R : create Fasta files necessary for the Virtual Secretome
 > 3. launch_outcyte_FAPs_pepSeqs.sbatch : run "Outcyte" to obtain the predictions of virtual secretome of each peptide-seq : \
 >    { Intracellular , Extracellular, Signal-peptide, UPS }
 > 4. annotate_outcyte_results.R : \
 >   i. Add "GeneSymbol" column to results of Outcyte,that contain only EntrezID.\
 >   ii. Match the genes labelled as "Signal-peptide"/"UPS" with the expression data for focusing only on the secreted genes.
 > 5. launch_NicheNet_faps.Rmd



# License
The project is licensed under the MIT License.
