#!/bin/bash

#===========================================================
# File Name: make_table_data.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-02-06 11:44:13
# Last Modified: 2023-02-06 11:45:07
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

for i in *-r.log;do echo $i >> table_cov-r.txt && head -n 7 $i | tail -n 1 >> table_cov-r.txt;done
for i in *-p.log;do echo $i >> table_cov-p.txt && head -n 7 $i | tail -n 1 >> table_cov-p.txt;done
