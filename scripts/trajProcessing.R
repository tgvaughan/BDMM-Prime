## Scripts useful for manipulating trajectory data

require(tidyverse)

parseTrajectory <- function(trajStr) {
  strValues <- str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE)
  values <- apply(strValues[,-2], 2, as.numeric)
  time <- values[,1]
  src <- values[,2]
  dest <- values[,3]
  mult <- values[,4]
  N <- values[,-(1:4)]
  event<- strValues[,2]

  res <- list(time = time,
              N = N,
              event = event,
              src = src,
              dest = dest,
              mult = mult)

  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, subsample=NA) {
    df <- NULL

    message("Loading ", filename,"...", appendLF = FALSE)
    df_in <- read_tsv(filename, col_types="ic")

    if (burninFrac>0) {
        n <- dim(df_in)[1]
        df_in <- df_in[-(1:ceiling(burninFrac*n)),]
    }

    if (!is.na(subsample)) {
        indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
        df_in <- df_in[indices,]
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
                                   N=trajStates$N,
                                   event=trajStates$event,
                                   src=trajStates$src,
                                   dest=trajStates$dest,
                                   mult=trajStates$mult))
        } else {
            ntypes <- dim(trajStates$N)[2]
            for (s in 1:ntypes) {
                
                df <- bind_rows(df,
                                tibble(traj=row,
                                       type=s-1,
                                       time=trajStates$time,
                                       N=trajStates$N[,s],
                                       event=trajStates$event,
                                       src=trajStates$src,
                                       dest=trajStates$dest,
                                       mult=trajStates$mult))
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

        age_summary$age <- grid_age
        df_grid <- bind_rows(df_grid, age_summary)
    }

    return(df_grid)
}
