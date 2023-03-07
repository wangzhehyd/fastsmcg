#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: make_data.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-28 20:21:59
# Last Modified: 2023-01-28 22:02:41
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import numpy as np

#deltas = np.arange(0,2.05, 0.05)
deltas = np.arange(0, 2.5, .125)

def get_data(log_file):
    data = []
    with open(log_file, 'r') as f:
        for l in f.readlines():
            d = float(l.split(',')[0].split()[-1])
            data.append(d)
    #print(data)
    return data


ifiles = ['data/evaluate_confgenx_cov-r.log', 'data/evaluate_conformator_cov-r.log', 'data/evaluate_omega_cov-r.log', 'data/evaluate_rdkit_cov-r.log', 
          'data/evaluate_confgf_cov-r.log', 'data/evaluate_dmcg_cov-r.log', 'data/evaluate_geodiff_cov-r.log', 'data/evaluate_geomol_cov-r.log', 'data/evaluate_torsional_diffusion_cov-r.log']

methods = {'data/evaluate_confgenx_cov-r.log':'ConfGenX', 'data/evaluate_conformator_cov-r.log':'Conformator', 'data/evaluate_omega_cov-r.log':'OMEGA', 'data/evaluate_rdkit_cov-r.log':'RDKit ETKDG',
        'data/evaluate_confgf_cov-r.log':'ConfGF', 'data/evaluate_dmcg_cov-r.log':'DMCG', 'data/evaluate_geodiff_cov-r.log':'GeoDiff', 'data/evaluate_geomol_cov-r.log':'GeoMol', 'data/evaluate_torsional_diffusion_cov-r.log':'Torsional Diffusion'}

print('delta,cov,method')
for f in ifiles:
    covs = get_data(f)
    methodn = [methods[f]] * len(covs)
    for delta,cov,method in zip(deltas,covs,methodn):
        print(f'{delta:<.3f},{cov},{method}')
