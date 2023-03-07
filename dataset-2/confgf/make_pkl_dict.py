#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: pkl2sdf.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-07-10 21:44:11
# Last Modified: 2023-01-28 15:31:19
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import collections
import pandas as pd
import copy

input_fname = sys.argv[1]
output_fname = sys.argv[2]
csv_file = sys.argv[3]

with open(input_fname, 'rb') as f:
    pkl = pickle.load(f)

test_data = pd.read_csv(csv_file)
#rdkit_smiles = test_data.smiles.values
cs = test_data.corrected_smiles.values

mol_dict = collections.OrderedDict()
writer = Chem.SDWriter(output_fname)

for index, data in enumerate(pkl):
    confs = []
    #smi = rdkit_smiles[index]
    smi = cs[index]

    num_confs = int(data.pos_gen.shape[0] / data.pos.shape[0])
    pos_gen_list = data.pos_gen.tolist()
    coords_list = [pos_gen_list[i:i+data.pos.shape[0]] for i in range(0, data.pos_gen.shape[0], data.pos.shape[0])]

    for conf_id in range(num_confs):
        mol = copy.deepcopy(data.rdmol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom_id in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(atom_id, coords_list[conf_id][atom_id])
        mol.AddConformer(conf, assignId=True)
        confs.append(mol)
    mol_dict[smi] = confs

with open(output_fname, 'wb') as f:
    pickle.dump(mol_dict, f)
