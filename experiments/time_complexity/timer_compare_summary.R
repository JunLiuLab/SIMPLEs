library(matrixStats)
library(here)
library(tidyverse)

mycnms <- c(1000, seq(2000, 12000, 2000)) / 1000

simple <- c(812, 1366, 2171, 3552, 4937, 6048, 7174)
rmagic <- c(24.4, 32.3, 44.2, 65.8, 71.562, 88.1, 92.246)
saver <- c(1239.5, 2268.3, 2844.1, 4643.690, 11344.58,12967.54, NA)
scimpute <- c(294.342, 841.803, 3708.2, 8041.4, NA, NA, NA)
scrabble <- c(129.3, 987.8, 6162.3, NA, NA, NA,NA)

times <- cbind(simple, rmagic, saver, scimpute, scrabble) / 60.0
times <- as.data.frame(times)
times$cnm <- mycnms
rownames(times) <- mycnms


time_r <- ggplot(times, aes(x = cnm), shape=16) +
  geom_point(aes(y = simple, color = "SIMPLE"), size=5) +
  geom_point(aes(y = saver, color="SAVER"), size=5) +
  geom_point(aes(y = scrabble, color="SCRABBLE"), size=5)+
  geom_point(aes(y = rmagic, color="RMAGIC"), size=5)+
  geom_point(aes(y = scimpute, color="SCIMPUTE"), size=5)+
  scale_x_continuous(name = "cells (k)", breaks = mycnms) +
  geom_line(aes(y=simple)) +
  geom_line(aes(y=saver)) +
  geom_line(aes(y=scrabble)) +
  geom_line(aes(y=rmagic)) +
  geom_line(aes(y=scimpute)) +
  theme(
    axis.text.x = element_text(color = "black", size = 13, family = "Helvetica"),
    axis.text.y = element_text(color = "black", family = "Helvetica", size = 13),
    axis.title.x = element_text(color = "black", size = 20, family = "Helvetica", face = "bold"),
    axis.title.y = element_text(color = "black", family = "Helvetica", size = 20),
    title = element_text(
      family = "Helvetica", color = "black", size = 24,
      face = "bold"
    ),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
  ylab("time (min)")

ggsave("./pdfs/simple_time_r.pdf", device = "pdf", plot = time_r, width = 11, height = 11)



