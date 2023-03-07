#!/bin/bash
#SBATCH -J omega_001
#SBATCH -o omega_001-%j.log
#SBATCH -e omega_001-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=1

echo -n "Start at   " >> omega_001_time.log
date +%T%t%F >> omega_001_time.log

ln -s ../../../dataset-1_with-index.smi ./dataset-1.smi
module load openeye && oeomega classic -in dataset-1.smi -out dataset-1-out.sdf -maxconfs 1 -warts false -rms 0 -strictStereo false

echo -n "Finish at  " >> omega_001_time.log
date +%T%t%F >> omega_001_time.log
