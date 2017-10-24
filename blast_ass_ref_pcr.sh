#!/usr/bin/env bash
#$ -cwd -S /bin/bash
source activate tipseqhunterBLAST
# Comment when running on phoenix
SGE_TASK_ID=7

# settings for cluster
# module load python/3.4.3

#folder=${1}
#standard=${2}
#label=${3}
#cd 'input/'${folder}
#F=(${SGE_TASK_ID}.*)
#cd ../..

#python glue.py ${F} ${standard} ${label}


samples_file=${1}
sample=$(sed -n "${SGE_TASK_ID}p" ${samples_file})
python blast_assembly_pcr.py ${sample}
