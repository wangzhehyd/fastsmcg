#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: boxplot_amrs.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-03-07 12:49:32
# Last Modified: 2023-03-07 12:49:33
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
app = sys.argv[3]

df = pd.read_csv(ifile, header=0)

plt.figure(figsize=(9, 3))
sns.boxplot(x='rb',y='arm',hue='es',data=df, palette="Blues", linewidth=1, width=0.6, fliersize=4, **{'flierprops':dict(markerfacecolor='white', markeredgecolor='k', marker='o'), 'capprops':dict(color='k'),'medianprops':dict(color='k'),'whiskerprops':dict(color='k'),'boxprops':dict(edgecolor='k'), })

#plt.xlim(0,17)
plt.ylim(0,6)
#plt.xticks(np.arange(0, +0.5, 0.5))

plt.xlabel('Number of rotatable bonds', fontsize = 12)
plt.ylabel('AMR-R ($\mathrm{\AA}$)', fontsize = 12)

#plt.legend(frameon=True, edgecolor='white', fontsize='small', loc='lower right')
plt.legend(frameon=True, edgecolor='white', fontsize='small', loc='upper left')
plt.title('%s' %app, fontsize = 12)

plt.grid(color='k', linestyle=':', alpha=0.2)
plt.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches = 'tight')
plt.close()
