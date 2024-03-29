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
    states <- NULL
    events <- NULL

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
            states <- bind_rows(states,
                                tibble(traj=row,
                                       type=0,
                                       time=trajStates$time,
                                       N=trajStates$N))
        } else {
            ntypes <- dim(trajStates$N)[2]
            for (s in 1:ntypes) {
                
                states <- bind_rows(states,
                                    tibble(traj=row,
                                           type=s-1,
                                           time=trajStates$time,
                                           N=trajStates$N[,s]))
            }
        }

        events <- bind_rows(events,
                            tibble(traj=row,
                                   time=trajStates$time,
                                   event=trajStates$event,
                                   src=trajStates$src,
                                   dest=trajStates$dest,
                                   mult=trajStates$mult))
    }

    states <- states %>% group_by(traj) %>% mutate(age=max(time)-time)
    events <- events %>% group_by(traj) %>% mutate(age=max(time)-time)
    
    message("done.")
    
    return(list(states=states, events=events))
}

gridTrajectoriesByTime <- function(trajStates, times) {
    return(trajStates %>%
           group_by(traj, type) %>%
           reframe(N=approx(age, N, times, method="constant", f=1, yleft=0)$y,
                     time=times))

}

gridTrajectoriesByAge <- function(trajStates, ages) {
    return(trajStates %>%
           group_by(traj, type) %>%
           reframe(N=approx(age, N, ages, method="constant", f=0, yright=0)$y,
                     age=ages))
}
