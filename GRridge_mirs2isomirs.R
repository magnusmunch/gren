rm(list=ls())

list.of.packages <- c("GRridge","edgeR", "pROC", "glmnet", "reshape", "ggplot2")
lapply(list.of.packages, require, character.only = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  Data preparation                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
load("normIsoProstate.RData")
load("ann2.RData")

resp = factor(c(replicate(25,"Normal"),replicate(36,"PC")))
n = length(resp)
mirs = as.character(ann2$Name)

# mirsData (combine normalized-isomirs)
datMirNorm = aggregate(normIso$datCount, by=list(mirs),sum)
dim(datMirNorm)
datMirTrans = sqrt(datMirNorm[,-1] + 3/8) - sqrt(3/8)
datMirStd = scale(datMirTrans, center = TRUE, scale = TRUE)
rownames(datMirTrans) = rownames(datMirStd) = datMirNorm[,1]

# The number of maximum selected markers
maxsel = c(5,10,15)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                One step approach by GRridge                         #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# one step approach is predictive modeling by using isomiRs-level data
abundance=rowSums(normIso$datCount) 
parAbund = CreatePartition(abundance,ngroup=10,uniform=TRUE,decreasing=TRUE)
Sds = apply(normIso$datTransformed,1,sd)
parSds = CreatePartition(Sds,decreasing=TRUE,uniform=TRUE,ngroup=10)
parIso1 = list(abn=parAbund,  sds=parSds);
ListMonotone = c(TRUE,FALSE)
for(i in maxsel){
  grIsomirs = grridge(normIso$datStd,resp,partitions=parIso1,
                      monotone=ListMonotone, selectionEN = TRUE, maxsel = i)
  grIsomirsCV = grridgeCV(grIsomirs, normIso$datStd, resp)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                             Two-step approach by GRridge                            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Two-step approach is predictive modeling by first using miRs-level data 
# (to select informative miRs)
# and using all isomiRs from the selected miRs in the second step
# note: for model evaluation, i used fully-cross validation (LOOCV)

# 1. FIRST STEP
abnMirs=rowSums(datMirNorm[,-1])
parAbund = CreatePartition(abnMirs,mingr=25,ngroup=10,decreasing=TRUE, uniform=TRUE)
Sds = apply(datMirTrans,1,sd)
parSds = CreatePartition(Sds,decreasing=TRUE,mingr=25,uniform=TRUE,ngroup=10)
parMirs = list(abn=parAbund, sds=parSds)
ListMonotone = c(TRUE,FALSE)
grMirs = grridge(datMirStd,resp,partitions=parMirs,
                 monotone=ListMonotone, selectionEN = TRUE, maxsel = 50)

# 2. SECOND STEP
datIso <- newDatIso(selMirs=names(grMirs$resEN$whichEN), mirs= mirs,
                    datCount=normIso$datCount, datStd=normIso$datStd,
                    datTrns = normIso$datTransformed)
abnIso = rowSums(datIso$datCount[,-1])
parAbundIso = CreatePartition(abnIso,ngroup=3,decreasing=TRUE, uniform=TRUE)
SdsIso = apply(datIso$datTrns[,-1],1,sd)
parSds = CreatePartition(SdsIso,decreasing=TRUE,uniform=TRUE,ngroup=3)
for(i in maxsel){
  grIso <- grridge(datIso$datStd[,-1],resp,partitions=list(abn=parAbundIso,sds=parSds),
                   monotone=c(TRUE,FALSE), selectionEN = TRUE, maxsel = i)
  idSel = grIso$resEN$whichEN
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



newDatIso <- function(selMirs=selMirs, mirs= mirs,datCount=datCount, 
                      datStd=datStd, datTrns=datTrns, trace=TRUE){
  selMirs = as.character(selMirs)  
  newDatCount = newDatStd = newDatTrns = data.frame()
  for(i in 1:length(selMirs)){
    mir_i = selMirs[i]
    idIso = which(mirs==mir_i)
    datCountIso_i = as.matrix(datCount[idIso,])
    datStdIso_i = as.matrix(datStd[idIso,])
    datTrnsIso_i = as.matrix(datTrns[idIso,])
    tempName = colnames(datStdIso_i)
    
    if(length(idIso)==1 & nrow(datStdIso_i)!=1){
      tempName = rownames(datStdIso_i)
      datCountIso_i=as.matrix(t(datCountIso_i))
      datStdIso_i=as.matrix(t(datStdIso_i))
      datTrnsIso_i=as.matrix(t(datTrnsIso_i))
    }
    colnames(datCountIso_i)=colnames(datTrnsIso_i)=colnames(datStdIso_i)=tempName
    rownames(datCountIso_i)=rownames(datTrnsIso_i)=rownames(datStdIso_i)=rownames(datCount)[idIso]
    newDatCount = rbind(newDatCount, data.frame(mirs=replicate(length(idIso),mir_i), 
                                                datCountIso_i))
    newDatStd = rbind(newDatStd, data.frame(mirs=replicate(length(idIso),mir_i), 
                                            datStdIso_i))
    newDatTrns = rbind(newDatTrns, data.frame(mirs=replicate(length(idIso),mir_i), 
                                              datTrnsIso_i))
  }
  if(trace==TRUE)  print(paste("The number of isomirs: ", nrow(newDatCount), " (", length(selMirs), 
                               " mirs)", sep=""))
  newDat = list(datCount=newDatCount, datStd=newDatStd, datTrns=newDatTrns)
  return(newDat)
}
