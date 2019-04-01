#!/usr/bin/env Rscript

### libraries
library(Biobase)

### load data
load("data/ESetMbolCSFPR2.Rdata")
load("data/NwkDegreePartition.Rdata")

### anonimise data
pheno <- pData(ESetMbolCSFPR2)
pheno.apoe <- pheno[(pheno$D_diag_name=="Probable AD" & pheno$APOE=="E4YES") |
                      (pheno$D_diag_name=="Subjectieve klachten" & 
                         pheno$APOE=="E4NO"), ]
pheno.apoe <- pheno.apoe[, "D_diag_name", drop=FALSE]

metabol <- t(exprs(ESetMbolCSFPR2))
metabol.apoe <- metabol[(pheno$D_diag_name=="Probable AD" & 
                           pheno$APOE=="E4YES") |
                          (pheno$D_diag_name=="Subjectieve klachten" & 
                             pheno$APOE=="E4NO"), ]

feat <- fData(ESetMbolCSFPR2)

rownames(NetworkDegreeClass) <- colnames(metabol.apoe) <- 
  paste0("metabol", 1:nrow(NetworkDegreeClass))

save(pheno.apoe, metabol.apoe, feat, NetworkDegreeClass,
     file="data/metabolomics_alzheimer_data.Rdata")
