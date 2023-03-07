#!/bin/bash
#SBATCH -J confgf_gen
#SBATCH -o confgf_gen-%j.log
#SBATCH -e confgf_gen-%j.err
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -w gpu06 # If you want to specify a computing node, you can write its name here and remove the first # 

echo -n "Start at   " >> confgf_gen_time.log
date +%T%t%F >> confgf_gen_time.log

. ~/.conda/setenv.sh && conda activate confgf && python -u ./gen.py --config_path ./drugs_default.yml --generator ConfGF --smiles ../dataset-2.csv

echo -n "Finish at  " >> confgf_gen_time.log
date +%T%t%F >> confgf_gen_time.log
