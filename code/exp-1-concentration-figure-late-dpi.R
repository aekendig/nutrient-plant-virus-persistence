## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-concentration-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(ggridges)
library(cowplot)

# load models
load("./output/exp-1-concentration-analysis-log-informative-rpv.rda")
load("./output/exp-1-concentration-analysis-log-informative-pav.rda")
load("./output/exp-1-concentration-analysis-log-informative-late-rpv.rda")
load("./output/exp-1-concentration-analysis-log-informative-late-pav.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2


#### edit data ####

# posterior samples
postr <- posterior_samples(m.li.r)
postp <- posterior_samples(m.li.p)
postrd <- posterior_samples(m.lid.r)
postpd <- posterior_samples(m.lid.p)

# rename columns
colnames(postr) <- colnames(postp) <- colnames(postrd) <- colnames(postpd) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "ar", "sigma", "lp")

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

avgpd <- postpd %>%
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

avgrd <- postrd %>%
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


#### figure of category averages ####

plotA <- avgp %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = lg_txt, hjust = 0.5),
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
  ylab("Est. ln(PAV concentration)") +
  ggtitle("Full dataset")

plotB <- avgr %>%
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

plotC <- avgpd %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "black", size = lg_txt),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = lg_txt, hjust = 0.5),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  xlab("Nutrient") +
  ylab("Est. ln(PAV concentration)") +
  ggtitle("Later dpi dataset")

plotD <- avgrd %>%
  mutate(Inoculation = recode(Inoculation, coinfection = "co.", single = "sing.")) %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
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
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), name = "Inoc.") +
  xlab("Nutrient") +
  ylab("Est. ln(RPV concentration)")


#### combine plots ####

# combine
plot <- plot_grid(plotA, plotC, plotB, plotD, 
                  labels = c("A", "B", "C", "D"), 
                  label_size = lg_txt, 
                  rel_widths = c(0.7, 1, 0.7, 1),
                  label_x = c(0, 0, 0, 0))

# print
pdf("./output/exp-1-concentration-figure-late-dpi.pdf", width = 4.75, height = 4)
plot
dev.off()


