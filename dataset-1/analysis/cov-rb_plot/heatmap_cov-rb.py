#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: heatmap_cov-rb.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-03-07 12:51:49
# Last Modified: 2023-03-07 12:51:50
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns

plt.switch_backend('agg')
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
plt.rcParams["font.family"] = "sans-serif"
mpl.rcParams['axes.linewidth'] = 1.0


ifile = sys.argv[1]
ofile = sys.argv[2]
num = sys.argv[3]

df = pd.read_csv(ifile, header=0)

#df = df.pivot("method", "rb", "cr")
df = df.pivot(index="order", columns="rb", values="cr")
#print(df)
yticklabels = ['ConfGenX','Conformator','OMEGA', 'RDKit ETKDG','ConfGF','DMCG','GeoDiff','GeoMol','Torsional Diffusion']


cmap = sns.color_palette("crest", as_cmap=True)
plt.figure(figsize=(9, 4))
ax = sns.heatmap(data=df, yticklabels=yticklabels, annot=True, linewidth=0.5, linecolor='k', cmap=cmap, vmin=0, vmax=100, fmt='.1f', annot_kws={'size':9}, cbar_kws={'label':'COV-R (%)'})

ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

plt.xlabel('Number of rotatable bonds', fontsize=12)
plt.ylabel('Method', fontsize=12)

plt.title('Maximum ensemble size = %s' %num, fontsize=12)

plt.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches = 'tight')
plt.close()
