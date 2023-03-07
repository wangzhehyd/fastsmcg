#!/bin/bash
#SBATCH -J confgenx
#SBATCH -o confgenx-%j.log
#SBATCH -e confgenx-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo -n "Start at   " >> confgenx_time.log
date +%T%t%F >> confgenx_time.log

module load schrodinger && run confgenx_gen.py . 0

echo -n "Finish at  " >> confgenx_time.log
date +%T%t%F >> confgenx_time.log
