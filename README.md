# Herbicide_Microbes_PT1
This repo contains scripts for processing of Herb PT1 enzyme and sequence data. 

Collaborators: Josh Harrison, Andrew Kniss, and Linda van Diepen

Descriptions of repo contents:

forModeling_*_otuTables -- Directories containing OTU tables for 16s and ITS data that have been formatted for CNVRG modeling.

post_modeling_ML_estimates -- estimates of Dirichlet and multinomial parameters. Means of PPDs are taken as estimates for each parameter. 

R/CNVRG_modeling_script.R -- perform modeling with CNVRG
R/comparing_HMC_VI_p.R -- diagnostic script to determine similarity in variational inference and Hamiltonian Monte Carlo parameter estimates. 
R/Create_Phyloseq_Object_From_Usearch.Rmd - title says it all!
R/differential_rel_abund.R - does differential relative abundance testing with CNVRG. Outputs are in "Diff_relAbund_summaries"

Diff_relAbund_summaries - contains the output of "differential_rel_abund.R". The files ending in "Summary" provide the number of taxa that differ between treatment groups.