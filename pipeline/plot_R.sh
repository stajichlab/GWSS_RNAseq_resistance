#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 2 --out logs/plot_R.log

Rscript Rscripts/kallisto_profile.R
