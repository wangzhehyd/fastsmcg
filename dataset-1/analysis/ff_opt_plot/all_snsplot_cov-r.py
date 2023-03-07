#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: snsplot_cov-r.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-28 20:42:50
# Last Modified: 2023-02-09 19:45:06
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

ifile1 = sys.argv[1]
ifile2 = sys.argv[2]
ifile3 = sys.argv[3]
ifile4 = sys.argv[4]
ifile5 = sys.argv[5]
ifile6 = sys.argv[6]
ifile7 = sys.argv[7]
ifile8 = sys.argv[8]
ofile = sys.argv[9]

df1 = pd.read_csv(ifile1, header=0)
df2 = pd.read_csv(ifile2, header=0)
df3 = pd.read_csv(ifile3, header=0)
df4 = pd.read_csv(ifile4, header=0)
df5 = pd.read_csv(ifile5, header=0)
df6 = pd.read_csv(ifile6, header=0)
df7 = pd.read_csv(ifile7, header=0)
df8 = pd.read_csv(ifile8, header=0)

palette = sns.color_palette("tab10") 

fig, ax = plt.subplots(3, 3, constrained_layout=True, figsize=(10, 8))

f1 = sns.lineplot(data=df1, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[0,0])
handles, labels = f1.get_legend_handles_labels()
f1.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f1.set_xlim(0,2)
f1.set_ylim(0,100)
f1.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f1.set_ylabel('COV-R (%)', fontsize = 10)
f1.grid(color='k', linestyle=':', alpha=0.2)

f2 = sns.lineplot(data=df2, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[0,1])
handles, labels = f2.get_legend_handles_labels()
f2.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f2.set_xlim(0,2)
f2.set_ylim(0,100)
f2.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f2.set_ylabel('COV-R (%)', fontsize = 10)
f2.grid(color='k', linestyle=':', alpha=0.2)


f3 = sns.lineplot(data=df3, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[0,2])
handles, labels = f3.get_legend_handles_labels()
f3.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f3.set_xlim(0,2)
f3.set_ylim(0,100)
f3.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f3.set_ylabel('COV-R (%)', fontsize = 10)
f3.grid(color='k', linestyle=':', alpha=0.2)


f4 = sns.lineplot(data=df4, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[1,0])
handles, labels = f4.get_legend_handles_labels()
f4.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f4.set_xlim(0,2)
f4.set_ylim(0,100)
f4.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f4.set_ylabel('COV-R (%)', fontsize = 10)
f4.grid(color='k', linestyle=':', alpha=0.2)

f5 = sns.lineplot(data=df5, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[1,1])
handles, labels = f5.get_legend_handles_labels()
f5.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f5.set_xlim(0,2)
f5.set_ylim(0,100)
f5.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f5.set_ylabel('COV-R (%)', fontsize = 10)
f5.grid(color='k', linestyle=':', alpha=0.2)

f6 = sns.lineplot(data=df6, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[1,2])
handles, labels = f6.get_legend_handles_labels()
f6.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f6.set_xlim(0,2)
f6.set_ylim(0,100)
f6.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f6.set_ylabel('COV-R (%)', fontsize = 10)
f6.grid(color='k', linestyle=':', alpha=0.2)


f7 = sns.lineplot(data=df7, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[2,0])
handles, labels = f7.get_legend_handles_labels()
f7.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f7.set_xlim(0,2)
f7.set_ylim(0,100)
f7.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f7.set_ylabel('COV-R (%)', fontsize = 10)
f7.grid(color='k', linestyle=':', alpha=0.2)

f8 = sns.lineplot(data=df8, x="Delta", y="COV", hue="Method", legend='brief', ax=ax[2,1])
handles, labels = f8.get_legend_handles_labels()
f8.legend(handles=handles[0:], labels=labels[0:], frameon=True, edgecolor='white', fontsize=8, loc='upper left')
f8.set_xlim(0,2)
f8.set_ylim(0,100)
f8.set_xlabel('Threshold $\mathit{\\delta}$  ($\mathrm{\AA}$)', fontsize = 10)
f8.set_ylabel('COV-R (%)', fontsize = 10)
f8.grid(color='k', linestyle=':', alpha=0.2)

#plt.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches = 'tight')
plt.close()
