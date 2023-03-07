#!/bin/bash
#SBATCH -J omega_050
#SBATCH -o omega_050-%j.log
#SBATCH -e omega_050-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=1

echo -n "Start at   " >> omega_050_time.log
date +%T%t%F >> omega_050_time.log

ln -s ../../../dataset-1_with-index.smi ./dataset-1.smi
module load openeye && oeomega classic -in dataset-1.smi -out dataset-1-out.sdf -maxconfs 50 -warts false -rms 0 -strictStereo false

echo -n "Finish at  " >> omega_050_time.log
date +%T%t%F >> omega_050_time.log
