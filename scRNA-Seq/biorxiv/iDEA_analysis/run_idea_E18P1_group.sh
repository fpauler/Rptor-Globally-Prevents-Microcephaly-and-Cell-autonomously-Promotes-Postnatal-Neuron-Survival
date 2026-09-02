#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --job-name=iDEA
#SBATCH --output=iDEA.E18P1.group.out
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --no-requeue

module load R/4.4.1
Rscript iDEA.analysis.cl.R RDS_files/GO_anno_iDEA.RDS RDS_files/E18_P1_group_deg_list.RDS RDS_files/E18_P1
