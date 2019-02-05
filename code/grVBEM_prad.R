### libraries
library(edgeR)

### read in data
data <- getTCGA(disease="PRAD", data.type="miRNASeq")
isoDat <- data$dat

### prepare data
# assigning variables
isomirs <- rownames(isoDat)
barcodes <- colnames(isoDat)
PID <- substr(barcodes, 1, 12)
ctrlIndex <- substr(barcodes, 14, 14)=="1"
ctrlPID <- unique(PID[ctrlIndex])
tmrPID <- setdiff(PID[!ctrlIndex], ctrlPID)

# indices of variables to select
setIndex <- ((PID %in% ctrlPID) & ctrlIndex) | ((PID %in% tmrPID) & !ctrlIndex)
setIndex[setIndex & ctrlIndex] <- rev(!duplicated(rev(PID[setIndex & ctrlIndex])))
setIndex[setIndex & !ctrlIndex] <- rev(!duplicated(rev(PID[setIndex & !ctrlIndex])))

# Remove samples from opposite group and any other duplicates
isoDat <- isoDat[, setIndex]
barcodes <- barcodes[setIndex]
PID <- PID[setIndex]
ctrlIndex <- ctrlIndex[setIndex]

# Select subset of data s.t. non-zero samples > 10% (unlikely to be important if less)
isomirIndex <- apply(isoDat, 1, function(x) mean(x > 0)) > 0.1
isoDat <- isoDat[isomirIndex, ]

# normalize data
isoLibSize <- colSums(isoDat)
isoRelLibSize <- isoLibSize/exp(mean(log(isoLibSize)))
isoNF <- edgeR::calcNormFactors(isoDat)*isoRelLibSize 
isoDat <- round(sweep(isoDat, 2, isoNF, "/"))

# transform data
isoDat <- sqrt(isoDat + 3/8) - sqrt(3/8)
str(isoDat)

x <- t(isoDat)
y <- 1 - as.numeric(ctrlIndex)
