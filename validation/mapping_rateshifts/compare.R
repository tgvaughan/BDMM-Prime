library(tidyverse)

df_sim <- read_tsv("simulator.log") %>% mutate(tree="Simulated")
df_map <- read_tsv("mapper.log") %>% mutate(tree="Mapped")

df <- bind_rows(df_sim, df_map) %>%
    select(-finalSampleOffset) %>%
    pivot_longer(cols=starts_with("tree."), names_prefix="tree.") %>%
    filter(name != "height") %>%
    filter(name != "treeLength")
    

ggplot(df, aes(col=tree,fill=tree)) +
    ## scale_y_log10()
    ## geom_density(aes(value), position="identity", alpha=0.5) +
    geom_histogram(aes(value), position="identity", alpha=0.5, bins=200) +
    facet_wrap(facets=vars(name), scales="free")
