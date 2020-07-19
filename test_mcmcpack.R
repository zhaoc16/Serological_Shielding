rm(list=ls())
require(MCMCpack)

x1 = rnorm(1000)
x2 = rnorm(1000)
x3 = rnorm(1000)

y = x1+x2+x3+rnorm(1000)+1

df.dat = data.frame(y=y, x1=x1, x2=x2, x3=x3)

reg.out = lm(y~x1+x2+x3, data = df.dat)

Bayes.reg.out = MCMCregress(y~x1+x2+x3, data = df.dat)

reg.out
