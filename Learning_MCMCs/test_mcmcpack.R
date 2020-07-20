rm(list=ls())
require(MCMCpack)

# x1 = rnorm(1000)
# x2 = rnorm(1000)
# x3 = rnorm(1000)
# 
# y = x1+x2+x3+rnorm(1000)+1
# 
# df.dat = data.frame(y=y, x1=x1, x2=x2, x3=x3)
# 
# reg.out = lm(y~x1+x2+x3, data = df.dat)
# 
# Bayes.reg.out = MCMCregress(y~x1+x2+x3, data = df.dat)
# 
# reg.out


## logistic regression with an improper uniform prior
## X and y are now in the global environment

logitfun <- function(beta){
  eta <- X %*% beta
  p <- 1.0/(1.0+exp(-eta))
  sum( y * log(p) + (1-y)*log(1-p) )
}

x1 <- rnorm(1000)
x2 <- rnorm(1000)
X <- cbind(1,x1,x2)
p <- exp(.5 - x1 + x2)/(1+exp(.5 - x1 + x2))
y <- rbinom(1000, 1, p)

post.samp <- MCMCmetrop1R(logitfun, theta.init=c(0,0,0),
                          thin=1, mcmc=40000, burnin=500,
                          tune=c(1.5, 1.5, 1.5),
                          verbose=TRUE, logfun=TRUE)

raftery.diag(post.samp)
plot(post.samp)
summary(post.samp)
