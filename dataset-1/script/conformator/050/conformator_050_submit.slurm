#!/bin/bash
#SBATCH -J conformator_050
#SBATCH -o conformator_050-%j.log
#SBATCH -e conformator_050-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=1

echo -n "Start at   " >> conformator_050_time.log
date +%T%t%F >> conformator_050_time.log

ln -s ../../../dataset-1_with-index.smi ./dataset-1.smi
module load conformator && conformator -i dataset-1.smi -o dataset-1-out.sdf -n 50 -q 1

echo -n "Finish at  " >> conformator_050_time.log
date +%T%t%F >> conformator_050_time.log
