#Suppose truth = negative linear relation
alpha <- -2; lambda <- 1
xlin <- 1:1000/1000
data <- alpha/lambda*(xlin^(lambda)-1) + rnorm(1000,0,1)
plot(xlin,data)

#estimation
f <- function(par,x,dat) {sum((dat-(par[1]/par[2])*(x^(par[2])-1))^2)}
opt <- optim(c(-1,1/2),f,x=xlin,dat=data)

#fit
trans <- function(par,x){(par[1]/par[2])*(x^(par[2])-1)}
pars <- opt$par
pars
points(xlin,sapply(xlin,trans,par=pars),col="red",type="l")

#Now suppose truth = negative inverse sqrt
alpha <- -2; lambda <- -1/2
xlin <- 1:1000/1000
data <- alpha/lambda*(xlin^(lambda)-1) + rnorm(1000,0,10)
plot(xlin,data)

#estimation
opt <- optim(c(-1,1/2),f,x=xlin,dat=data)
pars <- opt$par
pars
points(xlin,sapply(xlin,trans,par=pars),col="red",type="l")

#Now suppose truth = positive 1/5
alpha <- 10; lambda <- 1/5
xlin <- 1:1000/1000
data <- alpha*(xlin^(lambda)-1) + rnorm(1000,0,1)
plot(xlin,data)

#estimation
opt <- optim(c(1,1/2),f,x=xlin,dat=data)
pars <- opt$par
pars
points(xlin,sapply(xlin,trans,par=pars),col="red",type="l")

#now suppose truth = positive log (limit case for box-cox, lambda -> 0)
xlin <- 1:1000/1000
alpha=1
data <- alpha*log(xlin) + rnorm(1000,0,1)
plot(xlin,data)

#estimation
opt <- optim(c(1,1/2),f,x=xlin,dat=data)
pars <- opt$par
pars
points(xlin,sapply(xlin,trans,par=pars),col="red",type="l")
#conclusion : still fits rather well


xlin <- 1:1000/1000
alpha <- c(-0.5, 1)
data <- ((xlin + alpha[2])^alpha[1] - 1)/alpha[1]
plot(xlin, data, type="l")













