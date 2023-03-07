#!/bin/bash

#===========================================================
# File Name: make_table_data.sh
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2023-01-30 17:09:36
# Last Modified: 2023-01-30 17:23:27
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

cd 001
rm table*.txt
for i in evaluate_confgenx_001_cov-r.log evaluate_conformator_001_cov-r.log evaluate_omega_001_cov-r.log evaluate_rdkit_001_cov-r.log evaluate_confgf_001_cov-r.log evaluate_dmcg_001_cov-r.log evaluate_geodiff_001_cov-r.log evaluate_geomol_001_cov-r.log evaluate_torsional_diffusion_001_cov-r.log;do echo ${i} >> table_cov-r.txt && tail -n 4 ${i} | head -n 1 >> table_cov-r.txt;done
for i in evaluate_confgenx_001_cov-p.log evaluate_conformator_001_cov-p.log evaluate_omega_001_cov-p.log evaluate_rdkit_001_cov-p.log evaluate_confgf_001_cov-p.log evaluate_dmcg_001_cov-p.log evaluate_geodiff_001_cov-p.log evaluate_geomol_001_cov-p.log evaluate_torsional_diffusion_001_cov-p.log;do echo ${i} >> table_cov-p.txt && tail -n 4 ${i} | head -n 1 >> table_cov-p.txt;done
for i in evaluate_confgenx_001_arm-r.log evaluate_conformator_001_arm-r.log evaluate_omega_001_arm-r.log evaluate_rdkit_001_arm-r.log evaluate_confgf_001_arm-r.log evaluate_dmcg_001_arm-r.log evaluate_geodiff_001_arm-r.log evaluate_geomol_001_arm-r.log evaluate_torsional_diffusion_001_arm-r.log;do echo ${i} >> table_arm-r.txt && tail -n 4 ${i} | head -n 1 >> table_arm-r.txt;done
for i in evaluate_confgenx_001_arm-p.log evaluate_conformator_001_arm-p.log evaluate_omega_001_arm-p.log evaluate_rdkit_001_arm-p.log evaluate_confgf_001_arm-p.log evaluate_dmcg_001_arm-p.log evaluate_geodiff_001_arm-p.log evaluate_geomol_001_arm-p.log evaluate_torsional_diffusion_001_arm-p.log;do echo ${i} >> table_arm-p.txt && tail -n 4 ${i} | head -n 1 >> table_arm-p.txt;done

cd ../010
rm table*.txt
for i in evaluate_confgenx_010_cov-r.log evaluate_conformator_010_cov-r.log evaluate_omega_010_cov-r.log evaluate_rdkit_010_cov-r.log evaluate_confgf_010_cov-r.log evaluate_dmcg_010_cov-r.log evaluate_geodiff_010_cov-r.log evaluate_geomol_010_cov-r.log evaluate_torsional_diffusion_010_cov-r.log;do echo ${i} >> table_cov-r.txt && tail -n 4 ${i} | head -n 1 >> table_cov-r.txt;done
for i in evaluate_confgenx_010_cov-p.log evaluate_conformator_010_cov-p.log evaluate_omega_010_cov-p.log evaluate_rdkit_010_cov-p.log evaluate_confgf_010_cov-p.log evaluate_dmcg_010_cov-p.log evaluate_geodiff_010_cov-p.log evaluate_geomol_010_cov-p.log evaluate_torsional_diffusion_010_cov-p.log;do echo ${i} >> table_cov-p.txt && tail -n 4 ${i} | head -n 1 >> table_cov-p.txt;done
for i in evaluate_confgenx_010_arm-r.log evaluate_conformator_010_arm-r.log evaluate_omega_010_arm-r.log evaluate_rdkit_010_arm-r.log evaluate_confgf_010_arm-r.log evaluate_dmcg_010_arm-r.log evaluate_geodiff_010_arm-r.log evaluate_geomol_010_arm-r.log evaluate_torsional_diffusion_010_arm-r.log;do echo ${i} >> table_arm-r.txt && tail -n 4 ${i} | head -n 1 >> table_arm-r.txt;done
for i in evaluate_confgenx_010_arm-p.log evaluate_conformator_010_arm-p.log evaluate_omega_010_arm-p.log evaluate_rdkit_010_arm-p.log evaluate_confgf_010_arm-p.log evaluate_dmcg_010_arm-p.log evaluate_geodiff_010_arm-p.log evaluate_geomol_010_arm-p.log evaluate_torsional_diffusion_010_arm-p.log;do echo ${i} >> table_arm-p.txt && tail -n 4 ${i} | head -n 1 >> table_arm-p.txt;done

cd ../050
rm table*.txt
for i in evaluate_confgenx_050_cov-r.log evaluate_conformator_050_cov-r.log evaluate_omega_050_cov-r.log evaluate_rdkit_050_cov-r.log evaluate_confgf_050_cov-r.log evaluate_dmcg_050_cov-r.log evaluate_geodiff_050_cov-r.log evaluate_geomol_050_cov-r.log evaluate_torsional_diffusion_050_cov-r.log;do echo ${i} >> table_cov-r.txt && tail -n 4 ${i} | head -n 1 >> table_cov-r.txt;done
for i in evaluate_confgenx_050_cov-p.log evaluate_conformator_050_cov-p.log evaluate_omega_050_cov-p.log evaluate_rdkit_050_cov-p.log evaluate_confgf_050_cov-p.log evaluate_dmcg_050_cov-p.log evaluate_geodiff_050_cov-p.log evaluate_geomol_050_cov-p.log evaluate_torsional_diffusion_050_cov-p.log;do echo ${i} >> table_cov-p.txt && tail -n 4 ${i} | head -n 1 >> table_cov-p.txt;done
for i in evaluate_confgenx_050_arm-r.log evaluate_conformator_050_arm-r.log evaluate_omega_050_arm-r.log evaluate_rdkit_050_arm-r.log evaluate_confgf_050_arm-r.log evaluate_dmcg_050_arm-r.log evaluate_geodiff_050_arm-r.log evaluate_geomol_050_arm-r.log evaluate_torsional_diffusion_050_arm-r.log;do echo ${i} >> table_arm-r.txt && tail -n 4 ${i} | head -n 1 >> table_arm-r.txt;done
for i in evaluate_confgenx_050_arm-p.log evaluate_conformator_050_arm-p.log evaluate_omega_050_arm-p.log evaluate_rdkit_050_arm-p.log evaluate_confgf_050_arm-p.log evaluate_dmcg_050_arm-p.log evaluate_geodiff_050_arm-p.log evaluate_geomol_050_arm-p.log evaluate_torsional_diffusion_050_arm-p.log;do echo ${i} >> table_arm-p.txt && tail -n 4 ${i} | head -n 1 >> table_arm-p.txt;done
