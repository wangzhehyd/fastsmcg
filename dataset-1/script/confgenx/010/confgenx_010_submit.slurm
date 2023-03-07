#!/bin/bash
#SBATCH -J confgenx_010
#SBATCH -o confgenx_010-%j.log
#SBATCH -e confgenx_010-%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo -n "Start at   " >> confgenx_010_time.log
date +%T%t%F >> confgenx_010_time.log

ln -s ../../../dataset-1_with-index_2d.sdf ./dataset-1.sdf
module load schrodinger && confgenx -n 10 -m 10 -drop_problematic -LOCAL -WAIT dataset-1.sdf

echo -n "Finish at  " >> confgenx_010_time.log
date +%T%t%F >> confgenx_010_time.log
