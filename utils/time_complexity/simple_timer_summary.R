library(matrixStats)
library(here)
library(tidyverse)
load_simple_timer <- function(mycnm) {
  fnm <- paste0("simple_timer_", mycnm, ".rds")
  readRDS(here::here("TabulaMuris", "supplement", "timer", fnm))
}

mycnms <- c(1000, seq(2000, 12000, 2000))

mytimes <- mycnms %>% map(load_simple_timer)

mydf <- rbind(
  c(mytimes[[1]]),
  c(mytimes[[2]]),
  c(mytimes[[3]]),
  c(mytimes[[4]]),
  c(mytimes[[5]]),
  c(mytimes[[6]]),
  c(mytimes[[7]])
)
## proc.time format are hard to handle, use add-hoc codes.
mydf <- data.frame(elapsed = c(812, 1366, 2171, 3552, 4937, 6048, 7174) / 60.0, cnm = mycnms)


time_r <- ggplot(mydf, aes(x = cnm, y = elapsed)) +
  geom_point(shape = 21, fill = "blue", size = 3) +
  scale_x_continuous(name = "Number of cells", breaks = mycnms) +
  geom_smooth(method = "lm", span = 0.3) +
  theme(
    axis.text.x = element_text(color = "black", size = 13, family = "Helvetica"),
    axis.text.y = element_text(color = "black", family = "Helvetica", size = 13),
    axis.title.x = element_text(color = "black", size = 20, family = "Helvetica", face = "bold"),
    axis.title.y = element_text(color = "black", family = "Helvetica", size = 20),
    title = element_text(
      family = "Helvetica", color = "black", size = 24,
      face = "bold"
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  ylab("Time (min)") +
  ggtitle("Time complexity based on cell numbers")

ggsave("./pdfs/simple_time_r.pdf", device = "pdf", plot = time_r, width = 11, height = 11)



## * summary from sacct

# JobID    Elapsed            Eligible               Start                 End     ReqMem     MaxRSS    CPUTime           State
#--------------- ---------- ------------------- ------------------- ------------------- ---------- ---------- ---------- --------------- 
# 58979041.batch   00:15:11 2020-06-07T21:09:13 2020-06-07T21:09:13 2020-06-07T21:24:24    50000Mn   3630276K   02:01:28       COMPLETED
# 58979042.batch   00:23:29 2020-06-07T21:09:22 2020-06-07T21:09:22 2020-06-07T21:32:51    50000Mn   4273976K   03:07:52       COMPLETED
# 58979043.batch   00:36:47 2020-06-07T21:09:31 2020-06-07T21:09:31 2020-06-07T21:46:18    50000Mn   5803312K   04:54:16       COMPLETED
# 58979044.batch   00:59:44 2020-06-07T21:09:31 2020-06-07T21:09:31 2020-06-07T22:09:15    50000Mn   7204120K   07:57:52       COMPLETED
# 58979045.batch   01:22:52 2020-06-07T21:09:35 2020-06-07T21:09:35 2020-06-07T22:32:27    50000Mn   9833056K   11:02:56       COMPLETED
# 58979046.batch   01:41:22 2020-06-07T21:09:45 2020-06-07T21:09:45 2020-06-07T22:51:07    50000Mn  13568264K   13:30:56       COMPLETED
# 58979047.batch   02:00:07 2020-06-07T21:09:59 2020-06-07T21:09:59 2020-06-07T23:10:06    50000Mn  15203420K   16:00:56       COMPLETED

mymem <- data.frame(mem = c(3630, 4273, 5803, 7204, 9833, 13568, 15203) / 1000.0, cnm = mycnms)


mem_r <- ggplot(mymem, aes(x = cnm, y = mem)) +
  geom_point(shape = 21, fill = "blue", size = 3) +
  scale_x_continuous(name = "Number of cells", breaks = mycnms) +
  geom_smooth(method = "lm", span = 0.3) +
  theme(
    axis.text.x = element_text(color = "black", size = 13, family = "Helvetica"),
    axis.text.y = element_text(color = "black", family = "Helvetica", size = 13),
    axis.title.x = element_text(color = "black", size = 20, family = "Helvetica", face = "bold"),
    axis.title.y = element_text(color = "black", family = "Helvetica", size = 20),
    title = element_text(
      family = "Helvetica", color = "black", size = 24,
      face = "bold"
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  ylab("Memory (GB)") +
  ggtitle("Space complexity based on cell numbers")

ggsave("./pdfs/simple_memory_maxrss.pdf", device = "pdf", plot = mem_r, width = 11, height = 11)

