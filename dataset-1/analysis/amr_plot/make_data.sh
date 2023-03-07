#!/bin/bash

#===========================================================
# File Name: make_data.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-27 20:48:30
# Last Modified: 2023-03-07 12:48:14
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================


# confgenx
echo rb,arm,es > confgenx_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_confgenx_001_arms-r.csv >> confgenx_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_confgenx_010_arms-r.csv >> confgenx_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_confgenx_050_arms-r.csv >> confgenx_arms-r.csv

# conformator
echo rb,arm,es > conformator_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_conformator_001_arms-r.csv >> conformator_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_conformator_010_arms-r.csv >> conformator_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_conformator_050_arms-r.csv >> conformator_arms-r.csv


# omega
echo rb,arm,es > omega_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_omega_001_arms-r.csv >> omega_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_omega_010_arms-r.csv >> omega_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_omega_050_arms-r.csv >> omega_arms-r.csv


# rdkit
echo rb,arm,es > rdkit_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_rdkit_001_arms-r.csv >> rdkit_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_rdkit_010_arms-r.csv >> rdkit_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_rdkit_050_arms-r.csv >> rdkit_arms-r.csv


# confgf
echo rb,arm,es > confgf_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_confgf_001_arms-r.csv >> confgf_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_confgf_010_arms-r.csv >> confgf_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_confgf_050_arms-r.csv >> confgf_arms-r.csv

# dmcg
echo rb,arm,es > dmcg_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_dmcg_001_arms-r.csv >> dmcg_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_dmcg_010_arms-r.csv >> dmcg_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_dmcg_050_arms-r.csv >> dmcg_arms-r.csv

# geodiff
echo rb,arm,es > geodiff_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_geodiff_001_arms-r.csv >> geodiff_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_geodiff_010_arms-r.csv >> geodiff_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_geodiff_050_arms-r.csv >> geodiff_arms-r.csv

# geomol
echo rb,arm,es > geomol_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_geomol_001_arms-r.csv >> geomol_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_geomol_010_arms-r.csv >> geomol_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_geomol_050_arms-r.csv >> geomol_arms-r.csv

# torsional_diffusion
echo rb,arm,es > torsional_diffusion_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 1"}' data/evaluate_torsional_diffusion_001_arms-r.csv >> torsional_diffusion_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 10"}' data/evaluate_torsional_diffusion_010_arms-r.csv >> torsional_diffusion_arms-r.csv
awk -F ',' -v OFS=',' '{print $2,$3,"Maximum ensemble size = 50"}' data/evaluate_torsional_diffusion_050_arms-r.csv >> torsional_diffusion_arms-r.csv



#echo rb,arm,es > omega_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 1"}' evaluate_omega_001_arms-r.csv >> omega_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 10"}' evaluate_omega_010_arms-r.csv >> omega_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 50"}' evaluate_omega_050_arms-r.csv >> omega_arms-p.csv

#echo rb,arm,es > torsional_diffusion_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 1"}' evaluate_torsional_diffusion_001_arms-r.csv >> torsional_diffusion_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 10"}' evaluate_torsional_diffusion_010_arms-r.csv >> torsional_diffusion_arms-p.csv
#awk -F ',' -v OFS=',' '{print $2,$4,"Maximum ensemble size = 50"}' evaluate_torsional_diffusion_050_arms-r.csv >> torsional_diffusion_arms-p.csv
