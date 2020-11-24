# Simulations-
Normal Distribution, weak law of large numbers, central limit theorem, and Bootstrap. 

### Normal Distributions pnorm gives the distribution function
## Normal distributiion with default creates a nd with mean 0 and sd 1 with a lower tail.
pnorm(1.50)
### This also gives the same results.
pnorm(1.50, mean = 0, sd=1)
### You can change the mean and the standard deviation and the lower.tail 
pnorm(1.50, mean = 2, sd=2, lower.tail = TRUE, log.p = FALSE)
### By returning the upper tail
pnorm(1.90, mean = 5, sd=3, lower.tail = FALSE)


### qnorm gives the quantile function
qnorm(0.025)
### This also gives the same results.
qnorm(0.025, mean = 0, sd=1)
### You can change the mean and the standard deviation and the lower.tail 
qnorm(0.025, mean = 2, sd=1, lower.tail =FALSE)

### rnorm generates random deviates and this is mostly used . it allow us to get random draws from a normal distribution.
### Generating 5 draws from a normal distribution with mean 8 and standard deviation 2 from any set.seed
set.seed(500)
rnorm(5,8,2)

### Generating 10 draws from a normal distribution with mean 8 and standard deviation 2 from any set.seed
set.seed(200)
rnorm(10,8,2)

###The sample function. sample takes a sample of the specified size from the elements of x using either with or without replacement
###Sampling 5 values from 1 to 10 equally likely with replacement.
set.seed(3000)
sample(1:10, size = 5, replace = TRUE )

### Sampling 10 values from a sequence of 2 to 100 with interval 4 without replacement.
set.seed(1000)
sample(seq(2,200, by=4), size = 10, replace = FALSE)
### Another
sample(2:30, size = 5, replace = TRUE)
### Another but being a little bias. 
set.seed(100)
sample(1:5, 3, replace = T, prob = c(0.6, rep(0.1,4)))

### Check this out and why the error?
sample(1:10, size = 11, replace = F)
#### Because it cannot take a sample larger than the population when 'replace = FALSE'

 ## Simulation studies. 
 ### During simulations, it is a good idea to set a seed once at the beginning. 
##Examples. 
library(ggplot2)
set.seed(500)
runif(5)
runif(3)
save.seed<- .Random.seed
runif(4)
runif(3)
runif(10)
runif(5)

### Restoring the state of the Random number generator RNG. 
.Random.seed<-save.seed
runif(5)
runif(3)
runif(10)

### Confidence intervals. 
### Assume mu=15, sigma=7, n=2500, nsim=1200 alpha= 0.05


mu<-15

sigma<-7

n<-2500

n.sim<-1200

alpha<-0.05



F.MakeCI <- function(x, alpha, sigma){

x.bar <- mean(x)

se <- sigma/sqrt(length(x))

z.alpha2 <- qnorm(1-alpha/2)

CI <- c(x.bar - z.alpha2 * se, x.bar + z.alpha2 * se)

return(CI)

}


set.seed(200000)



CIs <- matrix(NA, n.sim, 2)

for (i in 1:n.sim){

X <- rnorm(n, mu, sigma)

CIs[i, ] <- F.MakeCI(X, alpha, sigma)

}


sum((mu > CIs[,1]) & (mu < CIs[,2]))/n.sim


### Another way without a loop. 

set.seed(200000)

sim.data <- replicate(n.sim, rnorm(n, mu, sigma), simplify=FALSE) ### returns a list

CIs <- vapply(sim.data, F.MakeCI, c(lower=0, upper=0), alpha=alpha, sigma=sigma)

sum((mu > CIs["lower",]) & (mu < CIs["upper",]))/n.sim
 

### Another way

average+-z*se

a<-10

s<-4

n<-30

error<-qnorm(0.975) *s/sqrt(n)

a-error

a+error

### Weak law of large numbers
### Assume it follow a poisson distribution. 



epsilon <- 0.01

mu <- 2


all.n <- c(10, 100, 500, 1000, 2500, 5000, 7500, 10000, 20000, 50000, 75000, 100000, 200000, 400000)

n.sim <- 1000



F.error <- function(X, mu){

return(abs(mean(X)-mu))

}

set.seed(39883)

probs <- numeric(length(all.n))

for (i in seq_along(all.n)){

print(i)

n <- all.n[i]



X.data <- replicate(n.sim, rpois(n, mu), simplify=FALSE)

all.error <- vapply(X.data, F.error, c(err=0), mu=mu)



probs[i] <- sum(all.error >= epsilon)/n.sim

}
library(tidyverse)
qplot(all.n, probs, geom = "line")





### Central Limit Theorem. You can change the value of the parameters

n <- 2000

std.X.bar <- NULL



lambda <- 10



for (r in 1:1000){

X = rpois(n, lambda)

std.X.bar[r] <- sqrt(n) * (mean(X)-lambda) / sqrt(lambda)

}



hist(std.X.bar, freq=FALSE)

x <- seq(-5, 10, length=100)

y <- dnorm(x)

lines(x,y, col='red')



pdf("Plot-CLT.pdf", width=8, height=8)

hist(std.X.bar, freq=FALSE)

lines(x,y, col='red')

dev.off()







### Bootstrap

###constructing  a 95% confidence interval for the correlation coefficient of two variables A and B



set.seed(1000)

X.obs <- rnorm(100)

Y.obs <- 0.5*X.obs + rnorm(100, 0, sqrt(0.75))

r.obs <- cor(X.obs, Y.obs)

r.obs





r.boot <- function(data, n.boot){

n <- nrow(data)

r.boot <- rep(NA, n.boot)



for (i in 1:n.boot){

data.boot <- data[sample(1:n, replace=TRUE),]

r.boot[i] <- cor(data.boot[,1], data.boot[,2])

}

return(r.boot)

}





Z.obs <-cbind(X.obs, Y.obs)

r.estimates <- r.boot(Z.obs, 1000)

#### standard error

CI1.obs <- c(r.obs - 1.96*sd(r.estimates), r.obs + 1.96*sd(r.estimates))

CI1.obs

#### quantiles

CI2.obs <- c(2*r.obs - quantile(r.estimates, 0.975), 2*r.obs - quantile(r.estimates, 0.025))



# Generate 1000 samples from the data generating model



n.sim <- 1000   ### 1000 replications

n.boot <- 500   ### 500 bootstrap samples



true.value <- 0.5



CI.mat <- matrix(NA, n.sim, 4)



set.seed(23889)



for (i in 1:n.sim){

X.obs <- rnorm(100)

Y.obs <- 0.5*X.obs + rnorm(100, 0, sqrt(0.75))

r.obs <- cor(X.obs, Y.obs)

r.est <- r.boot(cbind(X.obs, Y.obs), n.boot)

CI1 <- c(r.obs - 1.96*sd(r.est), r.obs + 1.96*sd(r.est))

CI2 <- c(2*r.obs - quantile(r.est, 0.975), 2*r.obs - quantile(r.est, 0.025))

CI.mat[i,] <- c(CI1, CI2)

}



sum((true.value > CI.mat[,1]) & (true.value < CI.mat[,2]))/n.sim

sum((true.value > CI.mat[,3]) & (true.value < CI.mat[,4]))/n.sim







