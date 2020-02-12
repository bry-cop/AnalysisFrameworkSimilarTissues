#!/usr/bin/env Rscript

args <- commandArgs (trailingOnly = TRUE)
# args <-  c(9, "/exports/humgen/tabbassidaloii/Publication/Outputs/WGCNA_dissTOM/dsfpower6.RData", 3, 15, 20, 30, 3, 0, 2, 4)
options (echo = TRUE)

# Preparing environmet and loading required packages 
suppressPackageStartupMessages ({ 
  library (WGCNA)
  library (cluster)
  library (flashClust)
  library (foreach)
  library (doParallel) })
registerDoParallel (as.numeric (args [ 1 ]))
options (stringsAsFactors = FALSE)
options (java.parameters = "-Xmx8192m") 
enableWGCNAThreads ()

# Data reading
# Expression data
load ("Outputs/HealthyMusclesDataforWGCNA.RData")

# loading dissimilarity matrix
power <- gsub ("[^[:digit:], ]", "", args [ 2 ])
load (args [ 2 ]) 

# Module detection
# Set the minimum module size 
args <- as.numeric (args [ -2 ])
minModuleSize <- args [ 3 : (args [ 2 ] + 2) ]
deepSplit <- args [ 7:(args [ 6 ] + 6) ]
Combin <- expand.grid (minModuleSize, deepSplit)
print (Combin)
foreach (i = 1: nrow (Combin)) %dopar% 
  {
    # Module detection by cutting branches
    dynamicMods <- dynamicTreeCut :: cutreeDynamic (dendro = geneTree, distM = dissTOM, method = "hybrid", 
                                                     deepSplit = Combin [ i, 2 ], pamRespectsDendro = FALSE ,
                                                     minClusterSize = Combin [ i, 1 ])
    # Convert labels to colors for plotting
    dynamicColors <- WGCNA :: labels2colors (dynamicMods)
    
    # Calculate eigengenes (clustring modules based on expression similarities)
    MEList <- WGCNA :: moduleEigengenes (resid_datasetFilt_control, colors = dynamicColors)
    MEs <- MEList$eigengenes
    NAMES <- colnames ( dissTOM )
    save (dynamicColors, geneTree, NAMES, file = paste0 ("Outputs/WGCNA_ModuleDetection/ModuleDetectionPower", power ,
                                                         "MinModuleSize", Combin [ i, 1 ] ,
                                                         "deepSplit", Combin [ i, 2 ], ".RData"))
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
## [1] doParallel_1.0.15     iterators_1.0.12      foreach_1.4.7         flashClust_1.01-2     cluster_2.0.7-1      
## [6] WGCNA_1.68            RSQLite_2.1.4         fastcluster_1.1.25    dynamicTreeCut_1.63-1
## 
## loaded via a namespace (and not attached):
## [1] Biobase_2.30.0        bit64_0.9-7           splines_3.2.3         Formula_1.2-3         assertthat_0.2.1     
## [6] stats4_3.2.3          latticeExtra_0.6-28   blob_1.2.0            fit.models_0.5-14     robustbase_0.93-5    
## [11] impute_1.44.0         pillar_1.4.2          backports_1.1.5       lattice_0.20-38       glue_1.3.1           
## [16] digest_0.6.23         RColorBrewer_1.1-2    checkmate_1.9.4       colorspace_1.4-1      htmltools_0.4.0      
## [21] preprocessCore_1.32.0 Matrix_1.2-18         pcaPP_1.9-73          pkgconfig_2.0.3       purrr_0.3.3          
## [26] GO.db_3.2.2           mvtnorm_1.0-8         scales_1.1.0          htmlTable_1.13.3      tibble_2.1.3         
## [31] IRanges_2.4.8         ggplot2_3.2.1         nnet_7.3-12           BiocGenerics_0.16.1   lazyeval_0.2.2       
## [36] survival_2.44-1.1     magrittr_1.5          crayon_1.3.4          memoise_1.1.0         MASS_7.3-51.4        
## [41] foreign_0.8-72        tools_3.2.3           data.table_1.12.8     lifecycle_0.1.0       matrixStats_0.55.0   
## [46] stringr_1.4.0         S4Vectors_0.8.11      munsell_0.5.0         AnnotationDbi_1.32.3  rlang_0.4.2          
## [51] grid_3.2.3            rstudioapi_0.10       htmlwidgets_1.5.1     robust_0.4-18.1       base64enc_0.1-3      
## [56] gtable_0.3.0          codetools_0.2-16      DBI_1.0.0             rrcov_1.4-9           R6_2.4.1             
## [61] gridExtra_2.3         knitr_1.27            dplyr_0.8.3           bit_1.1-14            zeallot_0.1.0        
## [66] Hmisc_4.3-0           stringi_1.4.5         Rcpp_1.0.3            vctrs_0.2.0           rpart_4.1-15         
## [71] acepack_1.4.1         DEoptimR_1.0-8        tidyselect_0.2.5      xfun_0.12            
