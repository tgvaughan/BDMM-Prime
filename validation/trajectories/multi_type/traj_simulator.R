simBD <- function(lambda, mu, psi, T) {

    t <- 0
    n <- 1
    sample <- FALSE

    idx <- 1
    while (TRUE) {

        a_birth <- n[idx]*lambda
        a_death <- n[idx]*mu
        a_psi <- n[idx]*psi

        a_tot <- a_birth + a_death + a_psi

        if (a_tot > 0)
            tnew <- t[idx] + rexp(1, a_tot)
        else
            tnew <- +Inf

        if (tnew > T)
            break;

        t[idx+1] <- tnew

        u <- runif(1, 0, a_tot)
        if (u < a_birth) {
            ## Birth event
            
            n[idx+1] <- n[idx] + 1
            sample[idx+1] <- FALSE
        } else if ( u - a_birth < a_death) {
            ## Death event
            
            n[idx+1] <- n[idx] - 1
            sample[idx+1] <- FALSE
        } else {
            ## Sampling event
            
            n[idx+1] <- n[idx] - 1
            sample[idx+1] <- TRUE
        }

        idx <- idx + 1
    }

    ## t <- t[1:final_samp_idx]
    ## n <- n[1:final_samp_idx]
        
    return(tibble(t=t, n=n, sample=sample))
}

simBDconditioned <- function(lambda, mu, psi, T, minSamp) {

    while (TRUE) {
        traj <- simBD(lambda, mu, psi, T)
        if (sum(traj$sample) >= minSamp)
            break
    }

    finalSampleIdx <- max(which(traj$sample))

    return(tibble(t=traj$t[1:finalSampleIdx],
                  n=traj$n[1:finalSampleIdx],
                  sample=traj$sample[1:finalSampleIdx]))
}


simBDensemble <- function(lambda, mu, psi, T, minSamp, N=1000) {

    df <- NULL
    for (i in 1:N) {
        thistraj <- simBDconditioned(lambda, mu, psi, T, minSamp)
        df <- bind_rows(df,
                        tibble(traj=i, type=0, time=thistraj$t, I=thistraj$n)) 
    }

    return(df)
}


simSingleTMRCA <- function(traj) {
    first_samp_idx <- min(which(traj$sample))
    last_samp_idx <- max(which(traj$sample))
    idx <- last_samp_idx
    k <- 0

    while (TRUE) {
        if (traj$sample[idx]) {
            k <- k + 1
        } else if (traj$n[idx-1]<traj$n[idx]) {
            pcoal <- k*(k-1)/(traj$n[idx]*(traj$n[idx]-1))

            if (runif(1,0,1) < pcoal)
                k <- k - 1
        }

        if (idx < first_samp_idx && k == 1)
            return(traj$t[last_samp_idx] - traj$t[idx])

        idx <- idx - 1
    }
}

simTMRCAs <- function(lambda, mu, psi, T, minSamp, N=1000) {
    tmrcas <- rep(0, N)
    for (i in 1:N) {
        thistraj <- simBDconditioned(lambda, mu, psi, T, minSamp)
        tmrcas[i] <- simSingleTMRCA(thistraj)
    }

    return(tmrcas)
}
