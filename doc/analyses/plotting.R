library(tidyverse)

traj <- read_tsv("h3n2.h3n2_2deme.traj") %>% filter(Sample>=100000) %>%
    filter(variable=="N")

p <- ggplot(traj, aes(age,value,col=factor(type),group_by=factor(Sample))) +
    geom_line(alpha=0.2) +
    xlim(0,10) + ylab("Population size")
p

ggsave("trajectories.png", p, width=800, height=400, units="px", dpi=150)

ages <- seq(0,10,length.out=1001)
gridded_traj <- traj %>%
    group_by(Sample, type) %>%
    reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
            age=ages) %>%
    group_by(type,age) %>%
    summarize(low=quantile(N,0.025),
              high=quantile(N,0.975),
              med=quantile(N,0.5))

p <- ggplot(gridded_traj,
       aes(age, N, col=factor(type), fill=factor(type), y=med, ymin=low, ymax=high)) +
    geom_ribbon(alpha=0.5) +
    geom_line() + ylab("Population size")
p

ggsave("trajectories_CI.png", p, width=800, height=400, units="px", dpi=150)

p <- ggplot(gridded_traj %>%
            mutate(Date=ymd("2005/09/01")-age*365) %>%
            mutate(Location=recode(factor(type),
                              "0"="Hong Kong",
                              "1"="New Zealand")),
       aes(Date, N, col=Location, fill=Location, y=med, ymin=low, ymax=high)) +
    geom_ribbon(alpha=0.5) +
    geom_line() + ylab("Population size") +
    scale_x_date(date_breaks="1 year", date_labels="%Y")
p

ggsave("trajectories_CI_years.png", p, width=800, height=400, units="px", dpi=150)
