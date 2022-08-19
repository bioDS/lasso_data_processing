#!/usr/bin/Rscript
library(ggplot2)
times <- c()
cores <- 1:32

for (c in cores) {
  c_time = readRDS(sprintf("./time_large_%d.rds", c))
  print(c_time)
  times = c(times, c_time)
}

print(times)

rel_speed = c()
for (c in 1:length(cores)) {
  rel_speed = c(rel_speed, times[1]/times[c])
}

print(rel_speed)

time_df = data.frame(Cores=cores, Time=times, "Relative Speedup"=rel_speed)
print(time_df)

time_plot <- ggplot(time_df[1:16,], aes(x = Cores, y = Time)) +
  geom_line() +
  geom_point(aes(x=Cores, y=Time)) +
  geom_line(data=time_df[16:32,], color='grey') +
  geom_point(aes(x=Cores, y=Time)) +
  geom_point(data=time_df[17:32,], aes(x=Cores, y=Time), color='grey') +
  annotate("text", x=20, y=times[16] + 24, label="Second NUMA Node", color='DarkGrey') +
  annotate("segment", x=20,y=times[16] + 21, xend=24, yend=times[24] + 2, arrow=arrow(), color='DarkGrey') +
  theme_bw() +
  ylab("Time (s)") +
  xlab("Thread(s)")
ggsave("varying_cores_time.pdf", time_plot, width=3, height=3)

speedup_plot <- ggplot(time_df[1:16,], aes(x = Cores, y = Relative.Speedup)) +
  geom_line() +
  geom_point(aes(x=Cores, y=Relative.Speedup)) +
  geom_line(data=time_df[16:32,], color='grey') +
  geom_point(aes(x=Cores, y=Relative.Speedup)) +
  geom_point(data=time_df[17:32,], aes(x=Cores, y=Relative.Speedup), color='grey') +
  annotate("text", x=20, y=times[16] - 7.8, label="Second NUMA Node", color='DarkGrey') +
  annotate("segment", x=20,y=times[16] - 7.5, xend=24, yend=times[24] - 2.5, arrow=arrow(), color='DarkGrey') +
  #geom_line() +
  #geom_point() +
  theme_bw() +
  xlab("Thread(s)") +
  ylab("Relative Speed")
ggsave("varying_cores_speedup.pdf", speedup_plot, width=3, height=3)
