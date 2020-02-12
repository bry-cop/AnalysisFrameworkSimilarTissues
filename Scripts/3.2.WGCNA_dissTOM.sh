#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N WGCNA_dissTOM
#$ -cwd
#$ -j Y
#$ -V
#$ -m be

. /usr/local/Modules/current/init/bash

module load R/3.2.3
## 1. The number of selected powers
## 2: The selected powers

Rscript --vanilla 3.2.WGCNA_dissTOM.R ${1} ${2} ${3} ${4} ${5} ${6} ${7}

