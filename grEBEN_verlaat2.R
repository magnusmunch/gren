##############################  preamble  #############################
# grEBEN on real data                                                 #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 16-10-2017                                                 #
# last edited: 16-10-2017                                             #
#######################################################################

##########################  data description  #########################
# A deep sequencing analysis on small non-coding ribonucleic acid     #
# (miRNAseq) was performed on 56 samples (24 women with high-grade    #
# cervical intraepithelial neoplasia (CIN3) and 32 healthy women) for #
# the purpose of finding relevant screening markers for cervical      #
# cancer screening. The next generation sequencing analysis resulted  #
# in 2,576 transcripts. The data was normalized and pre-processed,    #
# rendering 772 transcripts. More detail description of the data sets #
# and the preprocessing step are available on the supplementary       #
# material of this following publication "Better diagnostic           #
# signatures from RNAseq data through use of auxiliary co-data",      #
# Bioinformatics (2017).                                              #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### loading data
data(dataVerlaat)

### data manipulation
data.mirna <- t(datcenVerlaat)
norm.mirna <- apply(data.mirna, 2, function(x) {(x - mean(x))/sd(x)})
resp.cin3 <- respVerlaat
n <- nrow(norm.mirna)
p <- ncol(norm.mirna)
m <- rep(1, n)

### partitions
annotation <- as.numeric(CpGann)
partitions1 <- list(annotation=CreatePartition(as.factor(annotation)))
partitions2 <- list(annotation=annotation)
partitions3 <- list(pvalues=CreatePartition(pvalFarkas, decreasing=FALSE, mingr=10, ngroup=10))
partitions4 <- list(pvalues=rep(1:length(partitions3$pvalues), 
                                times=lapply(partitions3$pvalues, length))[
                                  order(unlist(partitions3$pvalues))])
partitions5 <- list(annotation=partitions1$annotation,
                    pvalues=partitions3$pvalues)
partitions6 <- list(annotation=partitions2$annotation,
                    pvalues=partitions4$pvalues)

### using annotation as grouping information
# fitting models
enet.pen <- cv.pen(norm.mirna, resp.cin3, unpenalized=NULL, intercept=TRUE, psel=NULL)
fit.ridge <- cv.glmnet(norm.mirna, resp.cin3, family="binomial", alpha=0, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.cin3))
fit.lasso <- cv.glmnet(norm.mirna, resp.cin3, family="binomial", alpha=1, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.cin3))
fit.enet <- glmnet(norm.mirna, resp.cin3, family="binomial", standardize=FALSE,
                   alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                   lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit1.grridge <- grridge(t(norm.mirna), resp.cin3, partitions1)
fit1.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions2,
                     lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                     lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])

fit2.grridge <- grridge(t(norm.mirna), resp.cin3, partitions3, monotone=FALSE)
fit2.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=FALSE, decreasing=FALSE))

fit3.grridge <- grridge(t(norm.mirna), resp.cin3, partitions5, monotone=c(FALSE, TRUE))
fit3.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions6,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=c(FALSE, TRUE), decreasing=c(FALSE, FALSE)))





