#!/usr/bin/env Rscript

args <- commandArgs (trailingOnly = TRUE)
# args <-  c("Outputs/WGCNA_ModuleDetection/ModuleDetectionPower6MinModuleSize15deepSplit0.RData", 5, 0.1, 0.15, 0.20,0.25, 0.30)
options (echo = TRUE)

# Preparing environmet and loading required packages 
suppressPackageStartupMessages ({ 
  library (WGCNA)
  library (cluster)
  library (flashClust)
  library (gplots)
  library (foreach)
  library (doParallel) 
  library (data.table) })
registerDoParallel (as.numeric (args [ 2 ]))
options (stringsAsFactors = FALSE)
options (java.parameters = "-Xmx8192m") 
enableWGCNAThreads ()

# Data reading
# Expression data
load ("Outputs/HealthyMusclesDataforWGCNA.RData")
Setting <- gsub (".*ModuleDetection|.RData", "", args [ 1 ])
print ( Setting )
# loading Modules
load (args [ 1 ])
# Merging similar modules
args <- as.numeric (args [ -1 ])
CutHeight <- args [ 2 : (args [ 1 ] + 1) ]
foreach (i = 1: length (CutHeight)) %dopar% 
{
  # Call an automatic merging function
  merge <- WGCNA :: mergeCloseModules (resid_datasetFilt_control, dynamicColors, cutHeight = CutHeight [ i ] ,
                                        verbose = 3, corFnc = bicor) 
  # The merged module colors
  mergedColors <- merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs <- merge$newMEs
  write.csv (mergedMEs, paste0 ("Outputs/WGCNA_MergeModules/ModuleEigengenes_merged_", Setting ,
                                   "CutHeight", CutHeight [ i ], ".csv"), row.names = FALSE)

  # Summary output of network analysis results (after merging)
  module_colors <- unique (mergedColors)
  All_genes_modules <- as.data.frame (cbind (NAMES, mergedColors)) 
  colnames (All_genes_modules) <- NULL
  write.csv (All_genes_modules, paste0 ("Outputs/WGCNA_MergeModules/All_genes_modules_", Setting,
                                        "CutHeight", CutHeight [ i ], ".csv")  ,
             row.names = FALSE)
  # names (colors) of the modules
  nSamples <- nrow (resid_datasetFilt_control)
  modNames <- substring (names (mergedMEs), 3)
  geneModuleMembership <- as.data.frame (bicor (resid_datasetFilt_control, mergedMEs, maxPOutliers = 0.1)) 
  MMPvalue <- as.data.frame (corPvalueStudent (as.matrix (geneModuleMembership), nSamples)) 
  names (geneModuleMembership) <- gsub ("ME", "MM", names (geneModuleMembership))
  names (MMPvalue) <-  gsub ("ME", "p.MM", names (MMPvalue))
  write.csv (merge (geneModuleMembership, MMPvalue, by=0, all = TRUE, sort = FALSE) ,
              paste0 ("Outputs/WGCNA_MergeModules/GeneModule_Membership_MMPvalue_", Setting,
                      "CutHeight", CutHeight [ i ], ".csv"), 
              col.names = TRUE, row.names = FALSE)
  
  # Module membership values and gene significance
  for (j in 3 : length (phenWGCNA) [ 1 ]) 
  {
    # Define interested variable of datTrait
    traitOfInterest <- as.data.frame (phenWGCNA [, j ])
    names (traitOfInterest) <-  colnames (phenWGCNA) [ j ]
    # correlation matrix of each gene and each traits
    if (all (traitOfInterest [, 1 ] %in% c (0, 1))) {
      geneTraitSignificance <- as.data.frame (cor (resid_datasetFilt_control, traitOfInterest)) 
    } else { 
      geneTraitSignificance <- as.data.frame (bicor (resid_datasetFilt_control, traitOfInterest, maxPOutliers = 0.1)) 
    }
    GSPvalue <- as.data.frame (corPvalueStudent (as.matrix (geneTraitSignificance), nSamples)) 
    names (geneTraitSignificance) <- paste0 ("GS.", names (traitOfInterest))
    names (GSPvalue) <- paste0 ("p.GS.", names (traitOfInterest))
    write.csv (merge (geneTraitSignificance,  GSPvalue, by=0, all = TRUE, sort = FALSE), 
                paste0 ("Outputs/WGCNA_MergeModules/", colnames (phenWGCNA)[ j ] ,"_GeneTrait_Significance_GSPvalue_",
                        Setting, "CutHeight",  CutHeight [ i ],
                         ".csv"), col.names = TRUE, row.names = FALSE)
  }
}
  

sessionInfo()
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## locale:
## [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
## [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
## [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base
## 
## other attached packages:
## [1] data.table_1.12.8     doParallel_1.0.15     iterators_1.0.12      foreach_1.4.7         gplots_3.0.1.1
## [6] flashClust_1.01-2     cluster_2.0.7-1       WGCNA_1.68            RSQLite_2.1.4         fastcluster_1.1.25
## [11] dynamicTreeCut_1.63-1
## 
## loaded via a namespace (and not attached):
## [1] Biobase_2.30.0        bit64_0.9-7           splines_3.2.3         gtools_3.8.1          Formula_1.2-3
## [6] assertthat_0.2.1      stats4_3.2.3          latticeExtra_0.6-28   blob_1.2.0            fit.models_0.5-14
## [11] robustbase_0.93-5     impute_1.44.0         pillar_1.4.2          backports_1.1.5       lattice_0.20-38
## [16] glue_1.3.1            digest_0.6.23         RColorBrewer_1.1-2    checkmate_1.9.4       colorspace_1.4-1
## [21] htmltools_0.4.0       preprocessCore_1.32.0 Matrix_1.2-18         pcaPP_1.9-73          pkgconfig_2.0.3
## [26] purrr_0.3.3           GO.db_3.2.2           mvtnorm_1.0-8         scales_1.1.0          gdata_2.18.0
## [31] htmlTable_1.13.3      tibble_2.1.3          IRanges_2.4.8         ggplot2_3.2.1         nnet_7.3-12
## [36] BiocGenerics_0.16.1   lazyeval_0.2.2        survival_2.44-1.1     magrittr_1.5          crayon_1.3.4
## [41] memoise_1.1.0         MASS_7.3-51.4         foreign_0.8-72        tools_3.2.3           lifecycle_0.1.0
## [46] matrixStats_0.55.0    stringr_1.4.0         S4Vectors_0.8.11      munsell_0.5.0         AnnotationDbi_1.32.3
## [51] caTools_1.17.1.3      rlang_0.4.2           grid_3.2.3            rstudioapi_0.10       htmlwidgets_1.5.1
## [56] robust_0.4-18.1       bitops_1.0-6          base64enc_0.1-3       gtable_0.3.0          codetools_0.2-16
## [61] DBI_1.0.0             rrcov_1.4-9           R6_2.4.1              gridExtra_2.3         knitr_1.27
## [66] dplyr_0.8.3           bit_1.1-14            zeallot_0.1.0         Hmisc_4.3-0           KernSmooth_2.23-16
## [71] stringi_1.4.5         Rcpp_1.0.3            vctrs_0.2.0           rpart_4.1-15          acepack_1.4.1
## [76] DEoptimR_1.0-8        tidyselect_0.2.5      xfun_0.12

