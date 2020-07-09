## Scripts useful for plotting and manipulating trajectory data
## using ggplot and friends.

library(tidyverse)

parseTrajectory <- function(trajStr) {
  values <- apply(str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE), 2, as.numeric)
  time <- values[,1]
  I <- values[,-1]

  res <- list(time=time, I = I)

  return(res)
}

loadTrajectories <- function(filenames, burninFrac=0.1) {
    df <- NULL

    for (filename in filenames) {
        cat(paste("Loading", filename,"..."))
        df_in <- read_tsv(filename, col_types="ic")

        if (burninFrac>0) {
            N <- dim(df_in)[1]
            df_in <- df_in[-(1:ceiling(burninFrac*N)),]
        }
        
        for (row in 1:(dim(df_in)[1])) {
            trajStr <- df_in[row,2]
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



df <- loadTrajectories("traj_inference.traj", burninFrac=0)

dftrue <- loadTrajectories("traj_and_tree_simulator.traj", burninFrac=0)

df_compare <- bind_rows(gridTrajectories(df, seq(1,4)) %>% mutate(ensemble="filter"),
                        gridTrajectories(dftrue, seq(1,4)) %>% mutate(ensemble="direct"))

ggplot(df_compare %>% filter(time==3)) +
    geom_freqpoly(aes(I, colour=factor(type), linetype=ensemble)) +
    facet_grid(cols=vars(time)) +
    xlim(c(0, 50))
    ## scale_y_log10() # + scale_x_log10()
