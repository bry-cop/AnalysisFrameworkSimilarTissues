#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N WGCNA_MergeModules
#$ -cwd
#$ -j Y
#$ -V
#$ -m be

. /usr/local/Modules/current/init/bash

module load R/3.2.3
## 1. The defined modules input
## 2. The number of CutHeights
## 3:7 The selected CutHeight

Rscript --vanilla 3.4.WGCNA_MergeModules.R ${1} ${2} ${3} ${4} ${5} ${6} ${7}

