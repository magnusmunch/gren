################################## manuscript ##################################
# ---- figures ----
# ---- boxplot_metabolomics_alzheimer_res1_auc ---- 
res1 <- read.table("results/metabolomics_alzheimer_res1.csv", 
                   stringsAsFactors=FALSE)
labels1 <- c(unique(as.character(feat$Platform))[-c(4, 5)], "Oxidative Stress")

