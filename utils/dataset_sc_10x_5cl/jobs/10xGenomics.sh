#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -p shared,stats
#SBATCH --mem=60000
#SBATCH -J 10x
#SBATCH -o 10x_%j_%N.out
#SBATCH -e 10x_%j_%N.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=debug.pie@gmail.com

module load Anaconda3/5.0.1-fasrc02
source activate my_root

module load GCC/8.2.0-2.31.1
module load R/3.6.3-fasrc01
# module load Python/3.7.2

Rscript --vanilla  ../impute_cluster_de.R  --method ${1} --K0 ${2} --M0 ${3} --test ${4} --ncores 8
