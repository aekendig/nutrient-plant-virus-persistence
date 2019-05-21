## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-analysis.R


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("./output/exp-1-analysis-data.csv")

# load models
load("./output/average-concentration-pav-cl-priors.rda")
load("./output/average-concentration-rpv-cl-priors.rda")


#### edit data ####

# inoculation column
dat <- dat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))


#### figure of raw data ####

# PAV
ggplot(filter(dat, target == "PAV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10),
        strip.text = element_text(color = "black", size = 12),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = c("black", "darkgoldenrod2", "palegreen4", "dodgerblue1"),
                      name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("") +
  ylab(expression(paste(Log[10], "(virus concentration)", sep = "")))

# RPV
ggplot(filter(dat, target == "RPV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10),
        strip.text = element_text(color = "black", size = 12),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = c("black", "darkgoldenrod2", "palegreen4", "dodgerblue1"),
                      name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("Days post inoculation") +
  ylab(expression(paste(Log[10], "(virus concentration)", sep = "")))
