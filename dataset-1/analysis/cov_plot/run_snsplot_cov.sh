#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: run_snsplot_cov-r.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-29 10:17:26
# Last Modified: 2023-02-02 12:48:04
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

run snsplot_cov-r.py 001/all_cov-r.csv all_001_cov-r.png 1
run snsplot_cov-r.py 010/all_cov-r.csv all_010_cov-r.png 10
run snsplot_cov-r.py 050/all_cov-r.csv all_050_cov-r.png 50

run snsplot_cov-p.py 001/all_cov-p.csv all_001_cov-p.png 1
run snsplot_cov-p.py 010/all_cov-p.csv all_010_cov-p.png 10
run snsplot_cov-p.py 050/all_cov-p.csv all_050_cov-p.png 10

run snsplot_cov-r.py 001/all_cov-r.csv all_001_cov-r.pdf 1
run snsplot_cov-r.py 010/all_cov-r.csv all_010_cov-r.pdf 10
run snsplot_cov-r.py 050/all_cov-r.csv all_050_cov-r.pdf 50

run snsplot_cov-p.py 001/all_cov-p.csv all_001_cov-p.pdf 1
run snsplot_cov-p.py 010/all_cov-p.csv all_010_cov-p.pdf 10
run snsplot_cov-p.py 050/all_cov-p.csv all_050_cov-p.pdf 10
