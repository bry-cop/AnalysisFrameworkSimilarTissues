#!/usr/bin/env Rscript

args <- commandArgs (trailingOnly = TRUE)
# args <- "Outputs/WGCNA_MergeModules/All_genes_modules_Power8MinModuleSize15deepSplit2CutHeight0.2.csv" 

# Preparing environmet and loading required packages 
suppressPackageStartupMessages ({ 
  library (dplyr)
})

Genes_Modules <- read.csv (args [ 1 ], header = F)
  
Genes_Modules <- unique (Genes_Modules %>%
                            mutate (V1 = gsub ("___.*", "", Genes_Modules$V1))) 


# remove genes which were not clustered (grey module)
Genes_Modules <- Genes_Modules [ Genes_Modules$V2 != "grey", ]
Modules <- unique (Genes_Modules$V2)

All_pairs <- c ()
for (module in 1 : length (Modules)) { 
  subset <- Genes_Modules [ Genes_Modules$V2 == Modules [ module ], ]
  Pairs <- as.data.frame (t (combn (subset$V1, 2)))
  Pairs$module <- Modules [ module ] 
  All_pairs [[ module ]] <- Pairs
}

All_pairs <- as.data.frame (do.call (rbind, All_pairs))
head (All_pairs)
dim (All_pairs)

write.csv (All_pairs, paste0 ("Outputs/CoexpressedPairs/CoexpressedPairs", gsub (".*All_genes_modules", "", args [ 1 ])), row.names = F)

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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_0.8.3
## 
## loaded via a namespace (and not attached):
## [1] Rcpp_1.0.3        fansi_0.4.0       withr_2.1.2       crayon_1.3.4      assertthat_0.2.1  R6_2.4.1         
## [7] magrittr_1.5      pillar_1.4.2      rlang_0.4.2       cli_2.0.0         rstudioapi_0.10   tools_3.2.3      
## [13] glue_1.3.1        purrr_0.3.3       pkgconfig_2.0.3   sessioninfo_1.1.1 tidyselect_0.2.5  tibble_2.1.3     