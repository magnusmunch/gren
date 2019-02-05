library(gren)

# load data from the CLL study as provided in MOFAtools
data("CLL_data", package="MOFAtools")

# use methylation data, gene expression data and drug responses as predictors
CLL_data <- CLL_data[1:3]
CLL_data <- lapply(CLL_data,t)
ngr <- sapply(CLL_data,ncol)
CLL_data <- Reduce(cbind, CLL_data)

#only include patient samples profiles in all three omics
CLL_data <- CLL_data[apply(CLL_data,1, function(p) !any(is.na(p))),]
dim(CLL_data)

# prepare design matrix and response
x <- scale(CLL_data[,!grepl("D_002", colnames(CLL_data))])
probs <- rowMeans(CLL_data[,grepl("D_002", colnames(CLL_data))])
y <- as.numeric(probs > mean(probs))
part1 <- rep(1:3, times = ngr-c(5,0,0)) # group annotations to drugs, meth and RNA

fit1.gren1 <- gren(x, y, partitions=list(type=part1), alpha=0.05, psel=TRUE)
fit1.gren2 <- gren(x, y, partitions=list(type=part1), alpha=0.5, psel=TRUE)
fit1.gren3 <- gren(x, y, partitions=list(type=part1), alpha=0.95, psel=TRUE)




