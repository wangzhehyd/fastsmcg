#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: make_data.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-02-01 15:49:24
# Last Modified: 2023-02-02 19:27:20
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import pandas as pd

ifile = sys.argv[1]
ofile = sys.argv[2]

df = pd.read_csv(ifile, header=None)
df.columns = ['order','method','title','rb','cr']

df.groupby(['order','rb','method'])['cr'].mean().to_csv(ofile)
