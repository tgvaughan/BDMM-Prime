library(deSolve)


dpdt <- function(t, p, params) {
    res <- NULL
    for (i in 1:length(params$lambda)) {
        res[i] <- -(params$lambda[i] + params$mu[i] + params$psi[i] + sum(params$M[i,]))*p[i] +
            params$lambda[i]*p[i]^2 + params$mu[i] + (params$M[i,] %*% p)
    }
    list(res)
}

dpgdt <- function(t, pg, params) {
    res <- NULL
    for (i in 1:params$n) {
        res[i] <- -(params$lambda[i] + params$mu[i] + params$psi[i] + sum(params$M[i,]))*pg[i] +
            params$lambda[i]*pg[i]^2 + params$mu[i] + (params$M[i,] %*% pg[1:params$n])
        res[i+params$n] <- -(params$lambda[i] + params$mu[i] + params$psi[i] + sum(params$M[i,]))*pg[i+params$n] +
            params$lambda[i]*pg[i+params$n]*pg[i]
    }
    list(res)
}

params <- list(n=2,
               lambda = c(2,2),
               mu = c(1,1),
               psi = c(0.5,0.5),
               M = matrix(c(0,1,1,0), nrow=2, ncol=2))


## Compute the probability of the following tree:
## ((([type=1]:0.75)[type=2]:0.25),[type=2]:0.5,)[type=2]:0.2;

nsteps <- 101

                                        # Left edge
lpg <- c(0,0,params$psi[1],0)
res <- ode(lpg, seq(0,0.75,length.out=nsteps), dpgdt, params)
lpg <- res[nsteps, -1]
lpg[4] <- params$M[2,1]*lpg[3] 
lpg[3] <- params$M[1,1]*lpg[3]
res <- ode(lpg, seq(0.75,1.0,length.out=nsteps), dpgdt, params)
lpg <- res[nsteps, -1]

print(lpg)

##                                         # Right edge
## p <- c(0,0)
## res <- ode(p, seq(0, 0.5, length.out=nsteps), dpdt, params)
## p <- res[nsteps, -1]
## rpg <- c(p, 0, params$psi[2])
## res <- ode(rpg, seq(0.5,1.0,length.out=nsteps), dpgdt, params)
## rpg <-  res[nsteps, -1]

## print(rpg)

##                                         # Root edge
## rootpg <- c(lpg[1:2], params$lambda*lpg[3:4]*rpg[3:4])
## res <- ode(rootpg, seq(1.0, 1.2, length.out=nsteps), dpgdt, params)
## rootpg <- res[nsteps, -1]

                                        # Tree probability:
## logTreeProb <- log(rootpg[3:4] %*% c(0.5, 0.5))
logTreeProb <- log(lpg[3:4] %*% c(0.5, 0.5))

logTreeProb

