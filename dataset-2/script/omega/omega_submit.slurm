#!/bin/bash
#SBATCH -J omega
#SBATCH -o omega-%j.log
#SBATCH -e omega-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo -n "Start at   " >> omega_time.log
date +%T%t%F >> omega_time.log

. ~/.conda/setenv.sh && conda activate fastsmcg && python omega_gen.py . 0

echo -n "Finish at  " >> omega_time.log
date +%T%t%F >> omega_time.log
