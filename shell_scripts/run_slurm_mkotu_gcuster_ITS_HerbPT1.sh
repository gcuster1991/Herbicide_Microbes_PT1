#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=microbiome
#SBATCH --time=0-5:00:00

module load swset/2018.05 gcc/7.3.0 vsearch/2.9.0

WORKDIR=/lscratch/mkotu_$SLURM_JOB_ID
mkdir $WORKDIR

# HerbPt1
# ITS
#1C
#/project/microbiome/data/seq/psomagen_29jan21novaseq1c/tfmergedreads/ITS/HerbPT1
#/project/microbiome/data/seq/psomagen_29jan21novaseq1c/tfmergedreads/ITS/sm18-herb


TFDIR=/project/microbiome/data/seq/psomagen_29jan21novaseq1c/tfmergedreads/ITS/HerbPT1
LIBNAME=1C
cat $TFDIR/*tfmergedreads.fa $TFDIR/joined/*tfmergedreads.fa | sed -E "s/^(>\w+)/\1_$LIBNAME/" >> $WORKDIR/tfmergedreads.fa

TFDIR=/project/microbiome/data/seq/psomagen_29jan21novaseq1c/tfmergedreads/ITS/sm18-herb
LIBNAME=1C
cat $TFDIR/*tfmergedreads.fa $TFDIR/joined/*tfmergedreads.fa | sed -E "s/^(>\w+)/\1_$LIBNAME/" >> $WORKDIR/tfmergedreads.fa

OTUDIR=/project/microbiome/data/seq/combined/gcuster/otu/HerbPt1/ITS
if [ ! -d "$OTUDIR" ]; then
    mkdir $OTUDIR
fi
cd $OTUDIR # this is the otu directory for output

vsearch --derep_fulllength ${WORKDIR}/tfmergedreads.fa --threads 32 --output uniqueSequences.fa --sizein --sizeout;
vsearch --cluster_unoise uniqueSequences.fa --relabel 'otu' --sizein --sizeout --consout zotus.fa --minsize 8 --threads 32 ;
vsearch --uchime3_denovo zotus.fa --nonchimeras zotus_nonchimeric.fa --threads 32;
vsearch --usearch_global ${WORKDIR}/tfmergedreads.fa --db zotus_nonchimeric.fa --otutabout - --id 0.99 --threads 32 | sed 's/^#OTU ID/OTUID/' > otutable
vsearch --search_exact ${WORKDIR}/tfmergedreads.fa --db zotus_nonchimeric.fa --otutabout - --threads 32 | sed 's/^#OTU ID/OTUID/' > otutable.esv
rm -rf $WORKDIR
