
setwd("C:\\VUData\\Maarten\\AnalysePlanProject2")
dataraw=read.table("miRseq.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A")) 
cn <- colnames(dataraw)[-(1:2)]
outcomeResp = read.table("Covariaat_response.txt",header=TRUE,row.names=1, sep="\t", na.strings = c("NA","#N/A"))

#    1= Xeloda monotherapie  
#    2= Xeloda + oxaliplatin al dan niet met bevacizumab of cetuximab    
#    3= Xeloda + irinotecan al dan niet met bevacizumab of cetuximab 
#    1= progressive disease; 2= clinical TherBenefit

therapy = factor(sapply(outcomeResp[,4],function(i) {if(i==1) return("XelMono") else {if(i==2) return("XelOxa") else return("XelIri")}}))

#library sizes
dataann <- dataraw[,1:2]
datacounts <- dataraw[,-(1:2)]
rownames(datacounts) <- sapply(1:nrow(datacounts),function(i) paste(rn[i],i))
libs <- colSums(datacounts)
summary(libs)

keep <- which(rowSums(datacounts>0) >= 5) 
length(keep) #2114 mirs left
rn <- as.character(dataann[,1])

rnkeep <- rn[keep]
TumMirs <- read.table("tumor_specific_miRs.txt",header=TRUE, sep="\t", na.strings = c("NA","#N/A"))

whichin <- match(rnkeep,TumMirs[,1])
#logfc <- TumMirs[,2]
whnonNA <- which(!is.na(whichin))



load("mirnorms.Rdata")
dim(mirnormcenTS)
library(GRridge)
library(Iso)

#preparing data
response <- outcomeResp[,7]
exclude <- outcomeResp[,2]
whrespNA <- which(is.na(response))
mirnormcen_resp0 <- mirnormcen[,-whrespNA]
resp0 = factor(sapply(outcomeResp[-whrespNA,7],function(i) {if(i==1) return("TherBenefit") else return("Progr")}))
covs <- outcomeResp[-whrespNA,]
adjth <- as.factor(covs[,3])
thscheme <- as.factor(covs[,4])
age <- covs[,5]
pcrcdiff0<- covs[,6]
whichnadiff <- which(is.na(pcrcdiff0) | pcrcdiff0==1) 
if(length(whichnadiff) > 0) pcrcdiff0[whichnadiff] <- 2
pcrcdiff <- as.factor(pcrcdiff0)
datfr0 <- data.frame(adjth,thscheme,age,pcrcdiff)
whnacov <- which(apply(apply(datfr0,2,is.na),1,max)==1)
if(length(whnacov) > 0) {datfr <- datfr0[-whnacov,];mirnormcen_resp <- mirnormcen_resp0[,-whnacov];
resp <- resp0[-whnacov]}


abundance=rowSums(mirnormcen_resp)
sdscounts <- apply(mirnormcen_resp,1,sd)
partkeep <- list(abund = CreatePartition(abundance,grsize=100,uniform=TRUE),
                 sd=CreatePartition(sdscounts,grsize=100,decreasing=FALSE,uniform=T),
                 TS = CreatePartition(as.character(whnonNA),as.character(1:length(abundance))))
save(TumMirs,whichin, mirnormcen_resp,resp,partkeep,datfr,file="forMagnus.Rdata")
grMaartenAll <- grridge(mirnormcen_resp,resp,partkeep,optl=NULL,unpenal = ~1 + adjth + thscheme + age + pcrcdiff, 
                        niter=1,method="exact",dataunpen = datfr,innfold=10,savepredobj="all",comparelasso=TRUE,optllasso=NULL,
                        compareunpenal=TRUE,selection=TRUE,maxsel=25) 
grMaartenAll.cv <- grridgeCV(grMaartenAll,mirnormcen_resp,resp,outerfold=10) 

save(grMaartenAll,grMaartenAll.cv,file="grMaartenAll.Rdata") #first two contain lasso results, latter two the rest


truth = as.numeric(grMaartenAll.cv[,1])-1
cutoffs <- rev(seq(0,1,by=0.01))

rocridge <- roc(probs=grMaartenAll.cv[,2],true=truth,cutoffs=cutoffs)
auc(rocridge)
rocgrridge <- roc(probs=grMaartenAll.cv[,3],true=truth,cutoffs=cutoffs)
auc(rocgrridge)
rocgrridgesel <- roc(probs=grMaartenAll.cv[,4],true=truth,cutoffs=cutoffs)
auc(rocgrridgesel)
roclasso <- roc(probs=grMaartenAll.cv[,5],true=truth,cutoffs=cutoffs)
auc(roclasso)
rocsimple <- roc(probs=grMaartenAll.cv[,6],true=truth,cutoffs=cutoffs)
auc(rocsimple)