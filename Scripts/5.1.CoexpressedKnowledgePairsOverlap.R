#!/usr/bin/env Rscript
library ( data.table )
args <- commandArgs ( trailingOnly = TRUE )
# args <- "Outputs/CoexpressedPairs/CoexpressedPairs_Power8MinModuleSize15deepSplit2CutHeight0.2.csv"

# Data reading
load ( args [ 1 ] )
# length ( Annotated_genes ) ## No. of genes with Reatome annotation in GSE36398 dataset

Coexpressed_pairs <- as.data.frame ( fread ( args [ 2 ]  ) ) [ -3 ]
Annotated_Coexpressed_pairs <- unique ( Coexpressed_pairs [ Coexpressed_pairs $ V1 %in% Annotated_genes & 
                                                              Coexpressed_pairs $ V2 %in% Annotated_genes , ] )
Annotated_Coexpressed_genes <- unique ( c ( Annotated_Coexpressed_pairs $ V1 , Annotated_Coexpressed_pairs $ V2 ) )
N  <- length ( Annotated_Coexpressed_genes ) ## Co-expressed and annotated genes
All_possible_pairs <-  ( length ( Annotated_genes ) * ( length ( Annotated_genes ) - 1 ) ) / 2 


Pairs <- merge ( Annotated_Coexpressed_pairs , All_pairs )
Annotated_Coexpressed_pairs <- Annotated_Coexpressed_pairs [ , c ( 2 , 1 ) ]
colnames ( Annotated_Coexpressed_pairs ) <- colnames ( All_pairs )
Pairs1 <- merge ( Annotated_Coexpressed_pairs , All_pairs )
Pairs <- rbind ( Pairs , Pairs1 )
SM_SP <- dim ( Pairs ) [ 1 ]
SM_nSP <- dim ( Annotated_Coexpressed_pairs ) [ 1 ] - SM_SP
nSM_SP <- dim ( All_pairs ) [ 1 ] - SM_SP
nSM_nSP <- All_possible_pairs - ( SM_SP + SM_nSP + nSM_SP )
table2_2 <- rbind ( c ( SM_SP , nSM_SP ) , c ( SM_nSP , nSM_nSP ) )
Ftest <- fisher.test ( table2_2 )
OutPut <- data.frame ( setting = gsub ( ".*essedPairs_|.csv" , "" , args [ 2 ]  ) ,
                               No_Annotated_Coexpressed_genes = length ( Annotated_Coexpressed_genes ) ,
                               Ftest =  Ftest $ p.value , odds_ratio = Ftest $ estimate )
names(OutPut) <- NULL
rownames(OutPut) <- NULL
print (OutPut)
