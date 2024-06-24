## Scripts useful for manipulating trajectory data

require(dplyr)

gridTrajectoriesByTime <- function(traj, times) {
traj %>% filter(variable=="N") %>%
    group_by(Sample, type) %>%
    reframe(N=approx(t, value, times, method="constant", f=1, yleft=0)$y,
            t=times)
}

gridTrajectoriesByAge <- function(traj, ages) {
traj %>% filter(variable=="N") %>%
    group_by(Sample, type) %>%
    reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
            age=ages)
}
