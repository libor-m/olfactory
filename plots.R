library(tidyverse)

read_tsv("data-results/r-esox.tsv") -> d

d %>% 
  mutate(len1 = end1 - zstart1,
         len2 = end2 - zstart2) %>% 
  ggplot(aes(size2, len1)) +
  geom_point() +
  xlab("query len") + ylab("match len")

ggsave("results/r-esox.png", width = 6, height = 5)
