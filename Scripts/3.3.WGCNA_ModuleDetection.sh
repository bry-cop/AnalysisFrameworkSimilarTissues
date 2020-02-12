#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N WGCNA_ModuleDetection
#$ -cwd
#$ -j Y
#$ -V
#$ -m be

. /usr/local/Modules/current/init/bash

module load R/3.2.3
## 1. The number of slots
## 2. The dissTOM matrix input
## 3. The number of MinModuleSizes
## 4:6 The selected MinModuleSize
## 7. The number of DeepSplits
## 8:10 The selected DeepSplit

Rscript --vanilla 3.3.WGCNA_ModuleDetection.R ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}

