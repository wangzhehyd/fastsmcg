#!/bin/bash

#===========================================================
# File Name: run_all_snsplot_cov-r.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-02-08 13:40:35
# Last Modified: 2023-02-09 19:49:34
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

python all_snsplot_cov-r.py 001/confgenx_all_cov-r.csv 001/conformator_all_cov-r.csv 001/omega_all_cov-r.csv \
      001/rdkit_all_cov-r.csv 001/confgf_all_cov-r.csv 001/geodiff_all_cov-r.csv  001/geomol_all_cov-r.csv 001/torsional_diffusion_all_cov-r.csv all_cov-r_001.pdf

python all_snsplot_cov-r.py 010/confgenx_all_cov-r.csv 010/conformator_all_cov-r.csv 010/omega_all_cov-r.csv \
      010/rdkit_all_cov-r.csv 010/confgf_all_cov-r.csv 010/geodiff_all_cov-r.csv  010/geomol_all_cov-r.csv 010/torsional_diffusion_all_cov-r.csv all_cov-r_010.pdf

python all_snsplot_cov-r.py 050/confgenx_all_cov-r.csv 050/conformator_all_cov-r.csv 050/omega_all_cov-r.csv \
      050/rdkit_all_cov-r.csv 050/confgf_all_cov-r.csv 050/geodiff_all_cov-r.csv  050/geomol_all_cov-r.csv 050/torsional_diffusion_all_cov-r.csv all_cov-r_050.pdf
