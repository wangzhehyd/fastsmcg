#!/bin/bash
#SBATCH -J dmcg_gen
#SBATCH -o dmcg_gen-%j.log
#SBATCH -e dmcg_gen-%j.err
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
##SBATCH -w xxx # If you want to specify a computing node, you can write its name here and remove the first # 

echo -n "Start at   " >> dmcg_gen_time.log
date +%T%t%F >> dmcg_gen_time.log

bash infer.sh

echo -n "Finish at  " >> dmcg_gen_time.log
date +%T%t%F >> dmcg_gen_time.log
