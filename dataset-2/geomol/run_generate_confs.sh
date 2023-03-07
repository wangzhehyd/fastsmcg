#!/bin/bash

#===========================================================
# File Name: run_generate_confs.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-08-10 09:35:27
# Last Modified: 2023-03-07 12:21:49
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================


. ~/.conda/setenv.sh && conda activate geomol && python ../generate_confs.py --trained_model_dir ../train --test_csv ../dataset-2.csv --dataset drugs
