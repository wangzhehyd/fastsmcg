#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: snsplot_cov-r.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-28 20:42:50
# Last Modified: 2023-02-02 12:47:41
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

plt.figure(figsize=(5, 4))
palette = sns.color_palette("tab10") 

#ax = sns.lineplot(data=df, x="Delta", y="COV", palette=palette, hue="Method", legend='brief')
ax = sns.lineplot(data=df, x="Delta", y="COV", hue="Method", legend='brief')

#legend = ax.legend()
#legend.texts[0].set_text("Whatever else")
handles, labels = ax.get_legend_handles_labels()
if num == '50':
    ax.legend(handles=handles[1:], labels=labels[1:], frameon=True, edgecolor='white', fontsize=9, loc='lower right')
else:
    ax.legend(handles=handles[1:], labels=labels[1:], frameon=True, edgecolor='white', fontsize=9, loc='upper left')


plt.xlim(0,2)
plt.ylim(0,100)
#plt.xticks(np.arange(0, +0.5, 0.5))

plt.xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 12)
plt.ylabel('COV-R (%)', fontsize = 12)

#plt.legend(frameon=True, edgecolor='white', fontsize='small', loc='lower right')
#plt.legend(frameon=True, edgecolor='white', fontsize='small', loc='upper left')
plt.title('Maximum ensemble size = %s' %num)

plt.grid(color='k', linestyle=':', alpha=0.2)
plt.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches = 'tight')
plt.close()
