#!/bin/bash

#===========================================================
# File Name: run_boxplot_amrs.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-03-07 12:49:38
# Last Modified: 2023-03-07 12:49:40
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

python boxplot_arms.py confgenx_arms-r.csv confgenx_arms-r.png ConfGenX
python boxplot_arms.py conformator_arms-r.csv conformator_arms-r.png Conformator
python boxplot_arms.py omega_arms-r.csv omega_arms-r.png OMEGA
python boxplot_arms.py rdkit_arms-r.csv rdkit_arms-r.png 'RDKit ETKDG'
python boxplot_arms.py confgf_arms-r.csv confgf_arms-r.png ConfGF
python boxplot_arms.py dmcg_arms-r.csv dmcg_arms-r.png DMCG
python boxplot_arms.py geodiff_arms-r.csv geodiff_arms-r.png GeoDiff
python boxplot_arms.py geomol_arms-r.csv geomol_arms-r.png GeoMol

python boxplot_arms.py torsional_diffusion_arms-r.csv torsional_diffusion_arms-r.png 'Torsional Diffusion'
