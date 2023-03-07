#!/bin/bash
#SBATCH -J conformator
#SBATCH -o conformator-%j.log
#SBATCH -e conformator-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo -n "Start at   " >> conformator_time.log
date +%T%t%F >> conformator_time.log

module load conformator
. ~/.conda/setenv.sh && conda activate fastsmcg && python conformator_gen.py . 0

echo -n "Finish at  " >> conformator_time.log
date +%T%t%F >> conformator_time.log
