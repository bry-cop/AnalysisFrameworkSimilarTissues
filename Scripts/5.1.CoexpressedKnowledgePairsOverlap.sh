#!/bin/bash
# Overlap results
for database in $(ls Outputs/*_KnowledgePairs.RData)
  do
  OutputName=$(echo $database | sed -e 's/Outputs\///g' | sed 's/_KnowledgePairs.RData//')
  for Setting in $(ls Outputs/CoexpressedPairs/CoexpressedPairs_*);
    do ./Scripts/5.1.CoexpressedKnowledgePairsOverlap.R $database $Setting | awk 'FNR > 1' | sed 's/ /,/g' >> Outputs/${OutputName}_OverlapResults.csv
done


