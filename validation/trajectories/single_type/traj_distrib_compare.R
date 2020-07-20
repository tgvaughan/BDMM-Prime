## Scripts useful for plotting and manipulating trajectory data
## using ggplot and friends.

library(tidyverse)

source("traj_simulator.R")

parseTrajectory <- function(trajStr) {
  values <- apply(str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE), 2, as.numeric)
  time <- values[,1]
  I <- values[,-1]

  res <- list(time=time, I = I)

  return(res)
}

parseEpiInfTrajectory <- function(trajStr, origin) {
    values <- apply(str_split(str_split(trajStr, ",")[[1]], ":", simplify=TRUE)[,c(1,3)], 2, as.numeric)
  time <- origin - values[,1]
  I <- values[,-1]

  res <- list(time=time, I = I)

  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, origins=NULL) {
    df <- NULL

    cat(paste("Loading", filename,"..."))
    df_in <- read_tsv(filename, col_types="ic")

    if (burninFrac>0) {
        N <- dim(df_in)[1]
        df_in <- df_in[-(1:ceiling(burninFrac*N)),]
    }
    
    for (row in 1:(dim(df_in)[1])) {
        trajStr <- df_in[row,2]
        if (length(origins)>0)
            trajStates <- parseEpiInfTrajectory(trajStr, origins[row])
        else
            trajStates <- parseTrajectory(trajStr)
        Idim <- dim(trajStates$I)
        if (length(Idim)==0) {
            ntypes <- 1
            df <- bind_rows(df,
                            tibble(traj=row,
                                   type=0,
                                   time=trajStates$time,
                                   I=trajStates$I))
        } else {
            ntypes <- dim(trajStates$I)[2]
            for (s in 1:ntypes) {
                
                df <- bind_rows(df,
                                tibble(traj=row,
                                       type=s-1,
                                       time=trajStates$time,
                                       I=trajStates$I[,s]))
            }
        }
    }
    
    cat("done.\n")
    
    return(df)
}

gridTrajectories <- function(trajdf, times) {
    df_grid <- NULL

    for (grid_time in times) {

        time_summary <- trajdf %>%
            group_by(traj, type) %>%
            summarize(
                I=I[max(which(time<=grid_time))],
                .groups = "drop_last")

        time_summary$time <- grid_time
        df_grid <- bind_rows(df_grid, time_summary)
    }

    return(df_grid)
}



dftrue <- loadTrajectories("traj_and_tree_simulator.traj", burninFrac=0)

df <- loadTrajectories("traj_inference.traj", burninFrac=0)

dfsim <- simBDensemble(2, 1, 0.5, 5, 2, 1000)

## df <- loadTrajectories("tree_prior_estimates.traj", burninFrac=0)
## origins <- read.table("traj_and_tree_simulator.log", header=T)$origin
## dfepi <- loadTrajectories("epiinf_results.traj", burninFrac=0, origins=origins)

times <- seq(0,5,length.out=51)
df_compare <- bind_rows(gridTrajectories(df, times) %>% mutate(ensemble="filter"),
                        gridTrajectories(dftrue, times) %>% mutate(ensemble="true"),
                        ## gridTrajectories(dfepi, times) %>% mutate(ensemble="epiinf"),
                        gridTrajectories(dfsim, times) %>% mutate(ensemble="R sim"))

df_comb <- bind_rows(dftrue %>% mutate(ensemble="direct"),
                     df %>% mutate(ensemble="filter"),
                     ## dfepi %>% mutate(ensemble="epiinf"),
                     dfsim %>% mutate(ensemble="R sim"))


p <- ggplot(df_compare %>%
       ## filter(ensemble != "R sim", ensemble != "epiinf") %>%
       group_by(time, ensemble) %>%
       summarize(Imed=mean(I), Ilow=quantile(I,0.25), Ihigh=quantile(I,0.75))) +
    geom_ribbon(aes(time, ymin=Ilow, ymax=Ihigh, fill=ensemble, color=ensemble), alpha=0.5) +
    geom_line(aes(time, Imed, colour=ensemble)) +
    ylab("Population size")#  +
    ## scale_y_log10()
p
ggsave("trajectory_comparison.png", p, width=20, height=15, units="cm")

ggplot(df_comb %>%
       filter(traj<=1000)) + #, ensemble != "filter")) +
    geom_line(aes(time, I, colour=ensemble, group=interaction(traj,ensemble)), alpha=0.05) +
    scale_y_log10()




for (i in 1:1000) {
    p <- ggplot(df_comb %>% filter(traj==i, ensemble == "filter")) +
        geom_step(aes(time, I, col=ensemble)) +
        geom_vline(aes(xintercept=max(time)), linetype="dashed")
    print(p)
    print(i)
    x <- readline()
}

ggplotly(p)


ggplot(df_comb %>% group_by(ensemble, traj) %>% summarize(Imed=median(I))) +
    geom_density(aes(Imed, col=ensemble))


## Tree prior estimate comparison

df_prob <- read.table("tree_prior_estimates.log", header=T)
p <- ggplot(df_prob) +
    geom_abline(col="grey") +
    geom_point(aes(-analyticTreeProb, -filterTreeProb), alpha=0.1, col="blue") +
    xlab("Analytic -log(P(tree))") +
    ylab("Importance sampled -log(P(tree))") +
    scale_x_log10() + scale_y_log10()
p
ggsave("prior_estimate_comparison.png", p, width=20, height=15, units="cm")
