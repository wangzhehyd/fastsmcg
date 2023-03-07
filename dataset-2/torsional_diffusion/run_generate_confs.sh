#!/bin/bash

#===========================================================
# File Name: run_generate_confs.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-09-02 08:10:42
# Last Modified: 2023-03-07 12:22:55
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

. ~/.conda/setenv.sh && conda activate torsional_diffusion && python ./generate_confs.py --test_csv ../dataset-2.csv --inference_steps 20 --model_dir ../ --out ./dataset-2-out.pkl --tqdm --batch_size 128 --no_energy
