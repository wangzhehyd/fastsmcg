#!/bin/bash
#SBATCH -J rdkit
#SBATCH -o rdkit-%j.log
#SBATCH -e rdkit-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo -n "Start at   " >> rdkit_time.log
date +%T%t%F >> rdkit_time.log

. ~/.conda/setenv.sh && conda activate fastsmcg && python rdkit_gen.py . 0

echo -n "Finish at  " >> rdkit_time.log
date +%T%t%F >> rdkit_time.log
