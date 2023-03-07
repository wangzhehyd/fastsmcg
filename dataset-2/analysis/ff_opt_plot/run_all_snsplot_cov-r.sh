#!/bin/bash

#===========================================================
# File Name: run_all_snsplot_cov-r.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-02-08 13:40:35
# Last Modified: 2023-02-09 20:16:18
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

python all_snsplot_cov-r.py data/confgenx_all_cov-r.csv data/conformator_all_cov-r.csv data/omega_all_cov-r.csv \
      data/rdkit_all_cov-r.csv data/confgf_all_cov-r.csv data/geodiff_all_cov-r.csv  data/geomol_all_cov-r.csv data/torsional_diffusion_all_cov-r.csv all_cov-r.pdf
