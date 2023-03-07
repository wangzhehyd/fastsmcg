#!/bin/bash
#SBATCH -J run_generate_confs
#SBATCH -o run_generate_confs-%j.log
#SBATCH -e run_generate_confs-%j.err
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
##SBATCH -w xxx # If you want to specify a computing node, you can write its name here and remove the first # 

echo -n "Start at   " >> run_generate_confs_time.log
date +%T%t%F >> run_generate_confs_time.log

bash run_generate_confs.sh

echo -n "Finish at  " >> run_generate_confs_time.log
date +%T%t%F >> run_generate_confs_time.log
