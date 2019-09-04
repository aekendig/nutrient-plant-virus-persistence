## Goal: figure of Experiment 1 virus concentrations from coinfection

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)


# import data
cdat <- read_csv("./output/exp-1-concentration-analysis-coinfection-data.csv")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2


#### edit data ####
cdat <- cdat %>%
  mutate(Nutrient = factor(nutrient, levels = c("low", "N", "P", "N+P")))


#### figure ####
pdf("./output/exp-1-coinfection-correlation-figure.pdf", width = 3, height = 3)
ggplot(cdat, aes(x = rpv_log_conc, y = pav_log_conc, colour = Nutrient)) +
  geom_point(shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.12, 0.85),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.001, "mm")) +
  scale_colour_manual(values = col_pal) +
  xlab("ln(RPV density in coinfection)") +
  ylab("ln(PAV density in coinfection)")
dev.off()