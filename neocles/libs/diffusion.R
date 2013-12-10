args <- commandArgs(trailingOnly = TRUE)
csv <- args[1]
save <- args[2]

library(ggplot2)

setEPS()
postscript(save, horizontal=FALSE)

d <- read.csv(csv, header=TRUE)
data <- data.frame(name=rep(d$name, 4),
                   time=c(rep(10, 5), rep(10e3, 5), rep(10e6, 5), rep(10e9, 5)),
                   distance=c(d$diffusion_in_10_ps, d$diffusion_in_10_ns,
                              d$diffusion_in_10_us, d$diffusion_in_10_ms))

g <- ggplot(data, aes(x=time, y=distance, group=name, color=name)) +
  geom_point() +
  geom_line() +
  scale_x_log10("Time (pS)") +
  scale_y_continuous("Distance (nM)")

print(g)
