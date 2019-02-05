library(sp)
library(gren)
library(pROC)
library(GRridge)


data(dataVerlaat)
x <- apply(t(datcenVerlaat), 2, function(x) {(x - mean(x))/sd(x)})
y <- respVerlaat

part1 <- as.numeric(CpGann)
labels1 <- levels(CpGann)

n <- nrow(x)
p <- ncol(x)

# fit1.gren1 <- gren(x, y, partitions=list(annotation=part1), alpha=0.05)
# fit1.gren2 <- gren(x, y, partitions=list(annotation=part1), alpha=0.5)
# fit1.gren3 <- gren(x, y, partitions=list(annotation=part1), alpha=0.95)

psel <- c(1:20, seq(22, 80, 2))
pred1.gren1 <- pred1.gren2 <- pred1.gren3 <- pred1.enet1 <- pred1.enet2 <- 
  pred1.enet3 <- matrix(0, nrow=n, ncol=length(psel))
set.seed(2018)
nfolds <- n
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=n %% nfolds),
  rep(n %/% nfolds, times=nfolds - (n %% nfolds))))))
for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  
  xtrain <- x[foldid!=k, ]
  xtest <- matrix(x[foldid==k, ], ncol=p)
  ytrain <- y[foldid!=k]
  
  # cv1.gren1 <- gren(xtrain, ytrain, partitions=list(annotation=part1), 
  #                   alpha=0.05, psel=psel, trace=FALSE)
  # cv1.gren2 <- gren(xtrain, ytrain, partitions=list(annotation=part1), 
  #                   alpha=0.5, psel=psel, trace=FALSE)
  cv1.gren3 <- gren(xtrain, ytrain, partitions=list(annotation=part1),
                    alpha=0.95, psel=psel, trace=FALSE)
  
  # pred1.gren1[foldid==k, ] <- predict(cv1.gren1, xtest, type="groupreg",
  #                                     s=cv1.gren1$freq.model$groupreg$lambda)
  # pred1.gren2[foldid==k, ] <- predict(cv1.gren2, xtest, type="groupreg",
  #                                     s=cv1.gren2$freq.model$groupreg$lambda)
  pred1.gren3[foldid==k, ] <- predict(cv1.gren3, xtest, type="groupreg",
                                      s=cv1.gren3$freq.model$groupreg$lambda)
  
  # pred1.enet1[foldid==k, ] <- predict(cv1.gren1, xtest, type="regular",
  #                                     s=cv1.gren1$freq.model$regular$lambda)
  # pred1.enet2[foldid==k, ] <- predict(cv1.gren2, xtest, type="regular",
  #                                     s=cv1.gren2$freq.model$regular$lambda)
  pred1.enet3[foldid==k, ] <- predict(cv1.gren3, xtest, type="regular",
                                      s=cv1.gren3$freq.model$regular$lambda)
  
}
cbind(y, as.numeric(pred1.gren1[, 1] < mean(pred1.gren1[, 1])))

plot(psel, apply(pred1.gren3, 2, function(pred) {pROC::roc(y, pred)$auc}), 
     type="l")
lines(psel, apply(pred1.enet1, 2, function(pred) {pROC::roc(y, pred)$auc}), 
      lty=2)





results1 <- list(idtrain=idtrain,
                 pred=list(enet1=pred1.enet1, enet2=pred1.enet2, 
                           enet3=pred1.enet3, gren1=pred1.gren1, 
                           gren2=pred1.gren2, gren3=pred1.gren3),
                 psel=list(enet1=psel1.enet1, enet2=psel1.enet2, 
                           enet3=psel1.enet3, gren1=psel1.gren1, 
                           gren2=psel1.gren2, gren3=psel1.gren3),
                 auc=list(enet1=auc1.enet1, enet2=auc1.enet2, enet3=auc1.enet3, 
                          gren1=auc1.gren1, gren2=auc1.gren2, gren3=auc1.gren3))











