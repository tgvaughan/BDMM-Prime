## Scripts useful for manipulating trajectory data

require(tidyverse)

parseTrajectory <- function(trajStr) {
  values <- apply(str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE), 2, as.numeric)
  time <- values[,1]
  N <- values[,-1]

  res <- list(time=time, N = N)

  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, origins=NULL) {
    df <- NULL

    message("Loading ", filename,"...", appendLF = FALSE)
    df_in <- read_tsv(filename, col_types="ic")

    if (burninFrac>0) {
        n <- dim(df_in)[1]
        df_in <- df_in[-(1:ceiling(burninFrac*n)),]
    }
    
    for (row in 1:(dim(df_in)[1])) {
        trajStr <- df_in[row,2]
        trajStates <- parseTrajectory(trajStr)
        Ndim <- dim(trajStates$N)
        if (length(Ndim)==0) {
            ntypes <- 1
            df <- bind_rows(df,
                            tibble(traj=row,
                                   type=0,
                                   time=trajStates$time,
                                   N=trajStates$N))
        } else {
            ntypes <- dim(trajStates$N)[2]
            for (s in 1:ntypes) {
                
                df <- bind_rows(df,
                                tibble(traj=row,
                                       type=s-1,
                                       time=trajStates$time,
                                       N=trajStates$N[,s]))
            }
        }
    }

    df <- df %>% group_by(traj) %>% mutate(age=max(time)-time)
    
    message("done.")
    
    return(df)
}

gridTrajectories <- function(trajdf, times) {
    df_grid <- NULL

    for (grid_time in times) {
        time_summary <- trajdf %>%
            group_by(traj, type) %>%
            summarize(
                N=N[max(which(time<=grid_time))],
                .groups = "drop_last")

        time_summary$time <- grid_time
        df_grid <- bind_rows(df_grid, time_summary)
    }

    return(df_grid)
}

gridTrajectoriesByAge <- function(trajdf, ages) {
    df_grid <- NULL

    for (grid_age in ages) {
        age_summary <- trajdf %>%
            group_by(traj, type) %>%
            summarize(
                N=N[max(which(age>=grid_age))],
                .groups = "drop_last")

        time_summary$time <- grid_age
        df_grid <- bind_rows(df_grid, age_summary)
    }

    return(df_grid)
}
