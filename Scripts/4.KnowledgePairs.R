#!/usr/bin/env Rscript

# Preparing environmet and loading required packages 
suppressPackageStartupMessages ({
  library (gProfileR)
  library (data.table)
  library (biomaRt)
    })

# Data reading
# Expression data
load ("Outputs/HealthyMusclesDataforWGCNA.RData")
Gene_list <- unique (gsub ("___.*", "", colnames (resid_datasetFilt_control)))

# Gene name to ENSEMBL gene id (including the version which was used for mapping the probes)
human <- useMart ("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://Apr2018.archive.ensembl.org")
annot_ens92 <- getBM (attributes = c ("ensembl_gene_id", "external_gene_name"), mart = human)
annot_ens92_subset <- annot_ens92 [ annot_ens92$external_gene_name %in% Gene_list, ]
human <- useMart ("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://Jul2018.archive.ensembl.org")
annot_ens93 <- getBM (attributes = c ("ensembl_gene_id", "external_gene_name"), mart = human)
annot_ens93_subset <- annot_ens93 [ annot_ens93$external_gene_name %in% Gene_list, ]
Annot_ENS_GeneNames <- merge (annot_ens93_subset, annot_ens92_subset, all = T)
Gene_list <- unique (Annot_ENS_GeneNames$ensembl_gene_id)
length (Gene_list)

# Annotating gene ids using gProfileR
GenePathway_all <- c ()
for (i in seq (0, length (Gene_list), 100)) {
  print (i)
  GS <- as.vector (na.omit (Gene_list [ (i + 1) : (i + 100) ]))
  BG <- unique (as.vector (Gene_list))
  print (length (GS))
  GP <- gprofiler (as.vector (GS), organism = 'hsapiens', significant = F,
                    domain_size = "known", custom_bg = BG, correction_method = "fdr")
  if (dim (GP) [ 1 ] != 0) { GP$query.number <- i }
  if (i == 0) GenePathway_all <- GP else GenePathway_all <- rbind (GenePathway_all, GP)
}
# Adding the gene symbols
Gene_Symbols <- list ()
for (i in 1 : nrow (GenePathway_all)){
  ENS <- GenePathway_all [ i, "intersection" ]
  ENS <- unique (unlist (strsplit (as.character (ENS), split = ",")))
  Gene_Symbol <- c ()
  for (j in ENS) {
    Gene_Symbol [ j ] <- Annot_ENS_GeneNames [ Annot_ENS_GeneNames$ensembl_gene_id == j, 2 ]
    
  }
  Gene_Symbols [[ i ]] <- paste (Gene_Symbol, collapse = ",")
}
Gene_Symbols <- do.call (rbind, Gene_Symbols)
GenePathway_all <- cbind (GenePathway_all, Gene_Symbols)
GenePathway_all [, c (1 : 3, 5, 7, 8, 11, 13) ] <- NULL

GenePathway_all <- aggregate (GenePathway_all [ - c (1, 3, 4, 5) ], by = list (GenePathway_all$term.size, GenePathway_all$term.id, GenePathway_all$ domain, GenePathway_all$term.name), c)
colnames (GenePathway_all) [ 1 : 4 ] <- c ("term.size", "term.id", "domain", "term.name")
x <- GenePathway_all$overlap.size
y <- character (length (x))
for (i in 1 : length (x)) y [ i ] <- as.numeric (sum (x [[ i ]]))
GenePathway_all$overlap.size <- as.numeric (y)
x <- GenePathway_all$intersection
y <- character (length (x)) 
for (i in 1 : length (x)) y [ i ] <- paste0 (unique (x [[ i ]]), collapse = ",")
GenePathway_all$intersection <- y
x <- GenePathway_all$Gene_Symbols
y <- character (length (x))
for (i in 1 : length (x)) y [ i ] <- paste0 (unique (x [[ i ]]), collapse = ",")
GenePathway_all$Gene_Symbols <- y

# Extracting the annotations from Reactome and HPO
GenePathway_REAC <- GenePathway_all [ grep ("REAC:", GenePathway_all$term.id), ]
GenePathway_HPO <- GenePathway_all [ grep ("HP:", GenePathway_all$term.id), ]
# Removing generic terms
GenePathway_REAC <- GenePathway_REAC [ ! GenePathway_REAC$term.name == "Reactome", ]
GenePathway_HPO <- GenePathway_HPO [ ! GenePathway_HPO$term.name %in% c ("All", "Mode of inheritance", "Phenotypic abnormality"), ]

rm (list = setdiff (ls (), c ("GenePathway_REAC", "GenePathway_HPO")))

# Knowledge pairs based on Reactome and HP
for ( database in c ("GenePathway_REAC", "GenePathway_HPO")) {
  database <- get (database)
  Pathways <- unique (database$term.id)
  All_pairs <- c ()
  for (i in 1 : length (Pathways)) {
    genes <- sort (unlist (strsplit (database [ database$term.id == Pathways [ i ], "Gene_Symbols" ], split = ",")))
    if (length (genes) != 1)
      { Pairs <- as.data.frame (t (combn (sort (genes), 2)))
      if (length (All_pairs) == 0) All_pairs <- Pairs else All_pairs <- rbind (All_pairs, Pairs)
    }
  }
  Annotated_genes <- unique ( unlist ( strsplit ( database $ Gene_Symbols , split = "," ) ) ) 
  OutputName <- paste0 ( ifelse (grep ("REAC", database$term.id [ 1 ]), "Outputs/Reactome", "Outputs/HPO"),
                         "_KnowledgePairs.RData")
  save (All_pairs, Annotated_genes, database, file = OutputName)
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] biomaRt_2.37.4    data.table_1.12.8 gProfileR_0.7.0  
## 
## loaded via a namespace (and not attached):
## [1] Rcpp_1.0.3           AnnotationDbi_1.32.3 rstudioapi_0.10      magrittr_1.5         hms_0.5.2           
## [6] progress_1.2.2       IRanges_2.4.8        BiocGenerics_0.16.1  bit_1.1-14           R6_2.4.1            
## [11] rlang_0.4.2          httr_1.4.1           stringr_1.4.0        blob_1.2.0           tools_3.2.3         
## [16] parallel_3.2.3       Biobase_2.30.0       DBI_1.0.0            assertthat_0.2.1     bit64_0.9-7         
## [21] digest_0.6.23        tibble_2.1.3         crayon_1.3.4         S4Vectors_0.8.11     vctrs_0.2.0         
## [26] bitops_1.0-6         RCurl_1.95-4.12      zeallot_0.1.0        memoise_1.1.0        RSQLite_2.1.4       
## [31] stringi_1.4.5        pillar_1.4.2         prettyunits_1.0.2    backports_1.1.5      stats4_3.2.3        
## [36] XML_3.98-1.20        pkgconfig_2.0.3     