#!/bin/bash
#SBATCH -J run_geodiff_gen
#SBATCH -o run_geodiff_gen-%j.log
#SBATCH -e run_geodiff_gen-%j.err
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -w gpu22 # If you want to specify a computing node, you can write its name here and remove the first # 

echo -n "Start at   " >> run_geodiff_gen_time.log
date +%T%t%F >> run_geodiff_gen_time.log

bash run_geodiff_gen.sh

echo -n "Finish at  " >> run_geodiff_gen_time.log
date +%T%t%F >> run_geodiff_gen_time.log
