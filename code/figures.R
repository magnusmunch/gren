################################## manuscript ##################################
# ---- figures ----
# ---- boxplot_metabolomics_alzheimer_res1_auc ---- 
library(pROC)
library(Biobase)
labels1 <- c(unique(as.character(feat$Platform))[-c(4, 5)], "Oxidative Stress")
res1 <- read.table("results/metabolomics_alzheimer_res1.csv", 
                   stringsAsFactors=FALSE)
load("data/ESetMbolCSFPR2.Rdata")
pheno <- pData(ESetMbolCSFPR2)
pheno.apoe <- pheno[(pheno$D_diag_name=="Probable AD" & pheno$APOE=="E4YES") |
                      (pheno$D_diag_name=="Subjectieve klachten" & 
                         pheno$APOE=="E4NO"), ]
alzheim.apoe <- as.numeric(pheno.apoe$D_diag_name) - 1

auc1 <- sapply(res1[grepl("pred", rownames(res1)), ], function(s) {
  pROC::auc(alzheim.apoe, s)})
psel1 <- colMeans(res1[grepl("psel", rownames(res1)), ])

plot(psel1, auc1)
boxplot(res1[grepl("psel", rownames(res1)), ])



barplot(rbind(fit1.gren1$lambdag$part, fit1.grridge$lambdamults$part), 
        beside=TRUE)
abline(h=1, lty=2)






