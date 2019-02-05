dy <- function(y, x, beta) {
  exp(-2*max(1 - y*t(x) %*% beta, 0))/
    (exp(-2*max(1 - y*t(x) %*% beta, 0)) + exp(-2*max(1 + y*t(x) %*% beta, 0)))
}


x <- rnorm(10)
beta <- rnorm(10)
dy(-1, x, beta) + dy(1, x, beta)



x <- seq(-5, 5, 0.01)
plot(x, abs(1 - x)^2, type="l")
plot(x, (1 - x)^2, type="l")
plot(x, 1 - 2*x + x^2, type="l")
sum(abs(1 - x)^2!=(1 - x)^2)

