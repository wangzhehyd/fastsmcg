#!/bin/bash

#===========================================================
# File Name: run_heatmap_cov-rb.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-03-07 12:51:13
# Last Modified: 2023-03-07 12:51:41
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

python heatmap_cov-rb.py all_001_covs-r.csv all_001_covs-r.pdf 1
python heatmap_cov-rb.py all_010_covs-r.csv all_010_covs-r.pdf 10
python heatmap_cov-rb.py all_050_covs-r.csv all_050_covs-r.pdf 50
