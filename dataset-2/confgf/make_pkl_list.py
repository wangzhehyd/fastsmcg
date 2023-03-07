#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: make_pkl_list.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-10-25 14:28:22
# Last Modified: 2022-10-25 14:29:30
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

input_dir = sys.argv[1]
output_file = sys.argv[2]

pkls = []
pkl_files = [f for f in os.listdir(input_dir) if f.endswith('pkl')]
pkl_sorted_files = sorted(pkl_files, key = lambda x: int(x.split('_')[2]))

for pkl in pkl_sorted_files:
    with open(pkl, 'rb') as f:
        pkl = pickle.load(f)
    pkls.append(pkl)

with open(output_file, 'wb') as f:
    pickle.dump(pkls, f)
