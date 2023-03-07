#!/bin/bash

#===========================================================
# File Name: geodiff_gen.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-07-22 09:52:34
# Last Modified: 2023-03-07 12:19:55
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

. ~/.conda/setenv.sh && conda activate geodiff && python ./geodiff_gen.py ./2860000.pt --input_file ../dataset-2.csv --out dataset-2-out.pkl
