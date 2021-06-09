#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=microbiome
#SBATCH --time=0-5:00:00

module load swset/2018.05 gcc/7.3.0 vsearch/2.9.0

vsearch -sintax /project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/16S/zotus_nonchimeric.fa -db /project/microbiome/ref_db/gg_16s_13.5.fa  -tabbedout /project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/16S/OUT.sintax -strand both -sintax_cutoff 0.8
