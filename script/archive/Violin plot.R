#violin plot
library(tidyverse)
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)
#install.packages("devtools")
library(devtools)
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")
library(reshape2)
df <- iris
dat <- melt(iris)

ggplot(data = dat, mapping = aes(x = Species, y = value, fill = Species)) +
  geom_flat_violin(scale = "count", trim = F) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               position = position_nudge(0.05))


+ 
  geom_dotplot(binaxis = "y", 
               dotsize = 0.5, 
               stackdir = "down", 
               binwidth = 0.1, 
               position = position_nudge(-0.025)) + 
  facet_wrap(~variable, ncol = 2, scales = "free") +
  theme(legend.position = "none") + #theme_bw()
  labs(x = "Species", y = "Length/width (cm)")
