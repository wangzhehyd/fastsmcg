#!/bin/bash

#===========================================================
# File Name: run_fastsmcg.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-11-03 23:48:31
# Last Modified: 2023-03-07 10:36:53
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================


model=$1
number=$2

# for ai model
. ~/.conda/setenv.sh && conda activate ${model} && python script/${model}/${model}_gen.py script/${model}/input/${number} 0

# for rdkit
#. ~/.conda/setenv.sh && conda activate fastsmcg && python script/rdkit/rdkit_gen.py script/rdkit/input/${number} 0
