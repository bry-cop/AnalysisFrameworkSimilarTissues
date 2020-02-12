#!/bin/bash

# Calculating dissimilarity topological matrix
NoDoParallel_1=6

NoPowers=6
powers=$(echo 6 8 10 14 18 22)

FIRST=$(qsub -l h_vmem=20G -pe BWA $NoDoParallel_1 Scripts/3.2.WGCNA_dissTOM.sh $NoPowers $powers)
echo $FIRST

# Detecting modules
NoDoParallel_2=9

NoMinModuleSize=3
MinModuleSizes=$(echo 15 20 30)

NoDeepSplit=3
DeepSplits=$(echo 0 2 4)

for dissTOM in $(ls Outputs/WGCNA_dissTOM/*);
  do qsub -W depend=afterany:$FIRST -l h_vmem=20G -pe BWA $NoDoParallel_2 Scripts/3.3.WGCNA_ModuleDetection.sh $NoDoParallel_2 $dissTOM $NoMinModuleSize $MinModuleSizes $NoDeepSplit $DeepSplits
done

# Merging modules
NoDoParallel_3=5

NoCutHeight=5
CutHeights=$(echo 0.1 0.15 0.2 0.25 0.3)

for Setting in $(ls Outputs/WGCNA_ModuleDetection/*);
  do qsub -hold_jid "WGCNA_ModuleDetection" -l h_vmem=20G -pe BWA $NoDoParallel_3 Scripts/3.4.WGCNA_MergeModules.sh $Setting $NoCutHeight $CutHeights
done

# Coexpressed Pairs
for Setting in $(ls Outputs/WGCNA_MergeModules/All_genes_modules*);
  do ./Scripts/3.5.CoexpressedPairs.R $Setting 
done


