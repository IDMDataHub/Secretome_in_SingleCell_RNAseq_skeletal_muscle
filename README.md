# Secretome_in_Single_Cell_skeletal_muscle
Extracted FAPs (Fibro-adipogenic Progenitors) expression from human & mouse S.C-RNA-seq of muscle from studies De Micheli et al.,Skeletal Muscle (2020) & Oprescu et al., iScience (2020)  in order to Study Virtual Secretome of FAPs during regeneration and in homeostasis. 

Scripts used as follows: 
- 1/analyse_with_Seurat.Rmd : analyse raw S.C.data published-> extract FAPs cells x FAPs markers table 
- 2/ extract_PepSeqs_for_Outcyte.R : create Fasta files necessary for the Virtual Secretome
- 3/ launch_outcyte_FAPs_pepSeqs.sbatch : run the tool "Outcyte" to obtain the predictions of virtual secretome : labelling of each peptide-seq : Intracellular , Extracellulare, Signal-peptide, UPS
- 4/ annotate_outcyte_results.R : Add "GeneSymbol" column to results of Outcyte,that contain only EntrezID --> match the genes labelled as "Signal-peptide"/"UPS" to expression data for further analysis and comparisons,only with these secreted genes.
