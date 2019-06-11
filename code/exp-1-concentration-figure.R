## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-analysis.R


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(brms)
library(tidyverse)
library(ggridges)
library(cowplot)

# import data
dat <- read_csv("./output/exp-1-analysis-data.csv")

# load models
load("./output/exp-1-analysis-log-informative-rpv.rda")
load("./output/exp-1-analysis-log-informative-pav.rda")


#### edit data ####

# inoculation column
dat <- dat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postr <- posterior_samples(m.li.r)
postp <- posterior_samples(m.li.p)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "ar", "sigma", "lp")

# posterior slopes
sloper <- postr %>%
  transmute(low = int - mean(int),
            high_N = N,
            high_P = P,
            high_NP = N + P + NP,
            low_co =  co,
            N_co = co + co_N,
            P_co = co + co_P,
            NP_co = co + co_N + co_P + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

slopep <- postp %>%
  transmute(low = int - mean(int),
            high_N = N,
            high_P = P,
            high_NP = N + P + NP,
            low_co =  co,
            N_co = co + co_N,
            P_co = co + co_P,
            NP_co = co + co_N + co_P + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

# check treatments
sloper %>% select(treatment, Inoculation, Nutrient) %>% unique()


#### figure of raw data ####

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2

# PAV
plotA <- ggplot(filter(dat, target == "PAV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.42, 0.92),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal,
                      name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("Days post inoculation") +
  ylab(expression(paste(Log[10], "(PAV concentration)", sep = "")))

# RPV
plotB <- ggplot(filter(dat, target == "RPV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation), show.legend = F) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.34, 0.93),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_size_manual(values = c(0.5, 0.5), name = "Inoculation") +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation") +
  xlab("Days post inoculation") +
  ylab(expression(paste(Log[10], "(RPV concentration)", sep = "")))


#### figure of model estimates ####

# annotation text
ann_textp <- data.frame(effect = c(0.52, 0.65), Nutrient = c("N+P", "N+P"), lab = c("single infection", "coinfection"), Inoculation = c("single", "coinfection"))
ann_textr <- data.frame(effect = c(0.66, 0.82), Nutrient = c("N+P", "N+P"), lab = c("single infection", "coinfection"), Inoculation = c("single", "coinfection"))

# PAV
plotC <- ggplot(slopep, aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  facet_wrap(~Inoculation)  +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.8))) +
  xlab("Nutrient effect on PAV concentration") +
  ylab("Density of posterior dist.") +
  geom_text(data = ann_textp, label = c("single infection", "coinfection"), nudge_y = 1.6, size = an_txt)

# RPV
plotD <- ggplot(sloper, aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  facet_wrap(~Inoculation)  +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("Nutrient effect on RPV concentration") +
  ylab("Density of posterior dist.") +
  geom_text(data = ann_textr, label = c("single infection", "coinfection"), nudge_y = 1.7, size = an_txt)


#### combine plots ####

pdf("./output/exp-1-concentration-figure.pdf", width = 6, height = 4)
plot_grid(plotA, plotC, plotB, plotD, labels = c("A", "C", "B", "D"), label_size = lg_txt)
dev.off()
