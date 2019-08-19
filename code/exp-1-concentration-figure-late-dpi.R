## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-concentration-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(ggridges)
library(cowplot)

# import data
pdat <- read_csv("./output/exp-1-concentration-analysis-pav-data.csv")
rdat <- read_csv("./output/exp-1-concentration-analysis-rpv-data.csv")

# load models
load("./output/exp-1-concentration-analysis-log-informative-late-rpv.rda")
load("./output/exp-1-concentration-analysis-log-informative-late-pav.rda")


#### edit data ####

# inoculation column
pdat <- pdat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

rdat <- rdat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postr <- posterior_samples(m.lid.r)
postp <- posterior_samples(m.lid.p)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "ar", "sigma", "lp")

# category average

avgp <- postp %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
   gather(key = "treatment", value = "effect") %>%
   mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
          Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
          Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
          Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

avgr <- postr %>%
  transmute(low = int,
            high_N = int + N,
            high_P = int + P,
            high_NP = int + N + P + NP,
            low_co =  int + co,
            N_co = int + N + co + co_N,
            P_co = int + P + co + co_P,
            NP_co = int + N + P + co + co_N + co_P + NP + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

# percentage increase

slopep <- postp %>%
  transmute(high_N = (exp(N) - 1),
            high_P = (exp(P) - 1),
            high_NP = (exp(N + P + NP) - 1),
            low_co =  (exp(co) - 1),
            N_co = (exp(co + co_N) - 1),
            P_co = (exp(co + co_P) - 1),
            NP_co = (exp(co + co_N + co_P + co_NP) - 1)) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

sloper <- postr %>%
  transmute(high_N = (exp(N) - 1),
            high_P = (exp(P) - 1),
            high_NP = (exp(N + P + NP) - 1),
            low_co =  (exp(co) - 1),
            N_co = (exp(co + co_N) - 1),
            P_co = (exp(co + co_P) - 1),
            NP_co = (exp(co + co_N + co_P + co_NP) - 1)) %>%
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

# PAV (concentration over time)
plotA <- ggplot(pdat, aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
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
  ylab("ln(PAV concentration)")

# RPV (concentration over time)
plotB <- ggplot(rdat, aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation), show.legend = F) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
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
  ylab("ln(RPV concentration)")


#### figure of category averages ####

plotC <- avgp %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21)) +
  xlab("Nutrient") +
  ylab("Est. ln(PAV concentration)")

plotD <- avgr %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21)) +
  xlab("Nutrient") +
  ylab("Est. ln(RPV concentration)")


#### combine plots ####

# combine
plot <- plot_grid(plotA, plotC, plotB, plotD, 
                  labels = c("A", "B", "C", "D"), 
                  label_size = lg_txt, 
                  rel_widths = c(1, 0.55, 1, 0.55),
                  label_x = c(0, 0, 0, 0))

# print
pdf("./output/exp-1-concentration-figure-late-dpi.pdf", width = 5, height = 4)
plot
dev.off()


#### numbers for text ####

# model summaries
summary(m.lid.r)
summary(m.lid.p)

# mean values in proportion change
sloper %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()

slopep %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()

