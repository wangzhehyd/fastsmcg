#!/bin/bash
#SBATCH -J conformator_001
#SBATCH -o conformator_001-%j.log
#SBATCH -e conformator_001-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=1

echo -n "Start at   " >> conformator_001_time.log
date +%T%t%F >> conformator_001_time.log

ln -s ../../../dataset-1_with-index.smi ./dataset-1.smi
module load conformator && conformator -i dataset-1.smi -o dataset-1-out.sdf -n 1 -q 1

echo -n "Finish at  " >> conformator_001_time.log
date +%T%t%F >> conformator_001_time.log
