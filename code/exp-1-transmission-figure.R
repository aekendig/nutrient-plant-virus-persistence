## Goal: figure of Experiment 1 transmission and model estimates from exp-1-transmission-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(ggridges)
library(cowplot)

# import data
datp <- read_csv("./output/exp-1-transmission-analysis-pav-data.csv")
datr <- read_csv("./output/exp-1-transmission-analysis-rpv-data.csv")

# load models
load("./output/exp-1-transmission-pav-up-concentration-informative-priors.rda")
load("./output/exp-1-transmission-rpv-up-concentration-informative-priors.rda")


#### edit data ####

# inoculation column
datp <- datp %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

datr <- datr %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postr <- posterior_samples(mruci)
postp <- posterior_samples(mpuci)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "conc", "co", "N", "P", "N_t", "P_t", "NP", "NP_t", "co_N", "co_P", "co_N_t", "co_P_t", "co_NP", "co_NP_t", "sd_round", "sd_time", "round_1", "round_2", "round_3", "round_4", "time_1", "time_2", "time_3", "time_4", "time_5", "time_6", "time_7", "time_8", "lp")

# percentage increase
sloper <- postr %>%
  transmute(low_co =  exp(int + co)/(1 + exp(int + co)) - exp(int)/(1 + exp(int)),
            high_N = exp(int + N)/(1 + exp(int + N)) - exp(int)/(1 + exp(int)),
            high_P = exp(int + P)/(1 + exp(int + P)) - exp(int)/(1 + exp(int)),
            high_NP = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)) - exp(int)/(1 + exp(int)),
            high_N_t = exp(int + N_t)/(1 + exp(int + N_t)) - exp(int)/(1 + exp(int)),
            high_P_t = exp(int + P_t)/(1 + exp(int + P_t)) - exp(int)/(1 + exp(int)),
            high_NP_t = exp(int + N_t + P_t + NP_t)/(1 + exp(int + N_t + P_t + NP_t)) - exp(int)/(1 + exp(int)),
            N_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)) - exp(int + N)/(1 + exp(int + N)),
            P_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)) - exp(int + P)/(1 + exp(int + P)),
            NP_co = exp(int + N + P + NP + co + co_N + co_P + co_NP)/(1 + exp(int + N + P + NP + co + co_N + co_P + co_NP)) - exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            N_co_t = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)) - exp(int + N_t)/(1 + exp(int + N_t)),
            P_co_t = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)) - exp(int + P_t)/(1 + exp(int + P_t)),
            NP_co_t = exp(int + N_t + P_t + NP_t + co + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t + P_t + NP_t + co + co_N_t + co_P_t + co_NP_t)) - exp(int + N_t + P_t + NP_t)/(1 + exp(int + N_t + P_t + NP_t))
            ) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Transmission = ifelse(grepl("t", treatment, fixed = T), "receiving", "source"),
         Transmission = factor(Transmission, levels = c("source", "receiving")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P", high_N_t = "N", high_P_t = "P", high_NP_t = "N+P", N_co_t = "N", P_co_t = "P", NP_co_t = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

slopep <- postp %>%
  transmute(low_co =  exp(int + co)/(1 + exp(int + co)) - exp(int)/(1 + exp(int)),
            high_N = exp(int + N)/(1 + exp(int + N)) - exp(int)/(1 + exp(int)),
            high_P = exp(int + P)/(1 + exp(int + P)) - exp(int)/(1 + exp(int)),
            high_NP = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)) - exp(int)/(1 + exp(int)),
            high_N_t = exp(int + N_t)/(1 + exp(int + N_t)) - exp(int)/(1 + exp(int)),
            high_P_t = exp(int + P_t)/(1 + exp(int + P_t)) - exp(int)/(1 + exp(int)),
            high_NP_t = exp(int + N_t + P_t + NP_t)/(1 + exp(int + N_t + P_t + NP_t)) - exp(int)/(1 + exp(int)),
            N_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)) - exp(int + N)/(1 + exp(int + N)),
            P_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)) - exp(int + P)/(1 + exp(int + P)),
            NP_co = exp(int + N + P + NP + co + co_N + co_P + co_NP)/(1 + exp(int + N + P + NP + co + co_N + co_P + co_NP)) - exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            N_co_t = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)) - exp(int + N_t)/(1 + exp(int + N_t)),
            P_co_t = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)) - exp(int + P_t)/(1 + exp(int + P_t)),
            NP_co_t = exp(int + N_t + P_t + NP_t + co + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t + P_t + NP_t + co + co_N_t + co_P_t + co_NP_t)) - exp(int + N_t + P_t + NP_t)/(1 + exp(int + N_t + P_t + NP_t))
  ) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Transmission = ifelse(grepl("t", treatment, fixed = T), "receiving", "source"),
         Transmission = factor(Transmission, levels = c("source", "receiving")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P", high_N_t = "N", high_P_t = "P", high_NP_t = "N+P", N_co_t = "N", P_co_t = "P", NP_co_t = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))%>%
  as_tibble()

# check treatments
sloper %>% select(treatment, Inoculation, Nutrient, Transmission) %>% unique()

# fitted values for concentration
datp_pred <- datp %>%
  select(conc, conc_s) %>%
  mutate(co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 0,
         high_P_t = 0) 

datp_pred <- datp_pred %>%
  cbind(fitted(mpuci, newdata = datp_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datr_pred <- datr %>%
  select(conc, conc_s) %>%
  mutate(co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 0,
         high_P_t = 0) 

datr_pred <- datr_pred %>%
  cbind(fitted(mruci, newdata = datr_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)


#### concentration figure ####

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 1

# PAV concentration
pconcp <- ggplot(datp, aes(x = conc)) +
    geom_point(aes(y = t_up), position = position_jitter(width = 0, height = 0.03), alpha = 0.5) +
    geom_ribbon(data = datp_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.5, size = 0.5) +
    geom_line(data = datp_pred, aes(y = pred), color = "black") +
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
  xlab("PAV concentration") +
  ylab("PAV transmission")

# RPV concentration
pconcr <- ggplot(datr, aes(x = conc)) +
    geom_point(aes(y = t_up), position = position_jitter(width = 0, height = 0.03), alpha = 0.5) +
    geom_ribbon(data = datr_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.5, size = 0.5) +
    geom_line(data = datr_pred, aes(y = pred), color = "black") +
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
    xlab("RPV concentration") +
    ylab("RPV transmission")

pconc <- plot_grid(pconcp, pconcr, 
                  labels = c("A", "B"), 
                  label_size = lg_txt)

pdf("./output/exp-1-concentration-transmission-figure.pdf", width = 4, height = 2)
pconc
dev.off()


#### figure of model estimates ####

# PAV single source nutrients
plotA <- ggplot(filter(slopep, Inoculation == "single" & Transmission == "source"), aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  stat_density_ridges(data = filter(slopep, Inoculation == "coinfection" & Transmission == "source"), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.5, 0.1),
        legend.direction = "horizontal",
        legend.key.size = unit(0.1, units = "in"),
        legend.box.spacing = unit(0.05, units = "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        plot.background = element_rect(fill = "gray")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  ylab("Density of posterior dist.") +
  xlim(-0.25, 0.5) +
  ggtitle("")

# PAV coinfection, source nutrients
plotB <- ggplot(filter(slopep, Inoculation == "coinfection" & Transmission == "source"), 
       aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        plot.background = element_rect(fill = "gray")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  xlim(-0.75, 0.4) +
  ggtitle("")

# PAV single infection, receiving nutrients
plotC <- ggplot(filter(slopep, Inoculation == "single" & Transmission == "receiving"), 
       aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  stat_density_ridges(data = filter(slopep, (Inoculation == "coinfection" & Transmission == "receiving") | (Inoculation == "coinfection" & Transmission == "source" & Nutrient == "low")), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  ylab("Density of posterior dist.\n(receiving nutrients)") +
  xlim(-0.6, 0.4) +
  ggtitle("")

# PAV coinfection, receiving nutrients
plotD <- ggplot(filter(slopep, (Inoculation == "coinfection" & Transmission == "receiving") | (Inoculation == "coinfection" & Transmission == "source" & Nutrient == "low")), 
       aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  xlim(-0.5, 0.7) +
  ggtitle("")

# RPV single source nutrients
plotE <- ggplot(filter(sloper, Inoculation == "single" & Transmission == "source"), aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  stat_density_ridges(data = filter(sloper, Inoculation == "coinfection" & Transmission == "source"), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.5, 0.13),
        legend.direction = "horizontal",
        legend.key.size = unit(0.15, units = "in"),
        legend.box.spacing = unit(0.1, units = "in"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        plot.background = element_rect(fill = "gray")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  ylab("Density of posterior dist.") +
  xlim(-0.5, 0.3)

# RPV coinfection, source nutrients
plotF <- ggplot(filter(sloper, Inoculation == "coinfection" & Transmission == "source"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        plot.background = element_rect(fill = "gray")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  xlim(-0.15, 0.7)

# RPV single infection, receiving nutrients
plotG <- ggplot(filter(sloper, Inoculation == "single" & Transmission == "receiving"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  stat_density_ridges(data = filter(sloper, (Inoculation == "coinfection" & Transmission == "receiving") | (Inoculation == "coinfection" & Transmission == "source" & Nutrient == "low")), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  xlim(-0.6, 0.25)

# RPV coinfection, receiving nutrients
plotH <- ggplot(filter(sloper, (Inoculation == "coinfection" & Transmission == "receiving") | (Inoculation == "coinfection" & Transmission == "source" & Nutrient == "low")), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = an_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  xlim(-0.05, 0.85) 

#### combine density plots ####

# combine
peff <- plot_grid(plotA, plotB, plotC, plotD, plotE, plotF, plotG, plotH,
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                  label_size = lg_txt,
                  rel_widths = c(1, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8),
                  rel_heights = c(rep(1, 4), rep(0.8, 4)),
                  nrow = 2,
                  label_x = c(0, -0.02, -0.02, -0.02, 0, -0.02, -0.03, -0.02),
                  label_y = c(rep(0.94, 4), rep(1, 4)))

# save with new x-axis labels
pdf("./output/exp-1-transmission-figure.pdf", width = 6, height = 4)
ggdraw(peff) +
  draw_line(x = c(0.002, 0.527), y = c(0.5, 0.5), color = "gray") +
  draw_line(x = c(0.294, 0.294), y = c(0.003, 0.997), color = "gray") +
  draw_label(label = "Proportion change in RPV concentration", x = 0.5, y = 0.04, size = lg_txt) +
  draw_label(label = "Proportion change in PAV concentration", x = 0.5, y = 0.54, size = lg_txt) +
  draw_label(label = "Source nutrients", x = 0.3, y = 0.98, size = lg_txt) +
  draw_label(label = "Receiving nutrients", x = 0.75, y = 0.98, size = lg_txt)
dev.off()


#### numbers for text ####

# model summaries
summary(m.li.r)
summary(m.li.p)

# mean values in proportion change
sloper %>%
  group_by(Inoculation, Nutrient) %>%
  summarise(mean_prop = mean(effect),
            se_prop = sd(effect) / sqrt(length(effect)),
            l_prop = quantile(effect, probs = 0.025),
            u_prop = quantile(effect, probs = 0.975))

slopep %>%
  group_by(Inoculation, Nutrient) %>%
  summarise(mean_prop = mean(effect),
            se_prop = sd(effect) / sqrt(length(effect)),
            l_prop = quantile(effect, probs = 0.025),
            u_prop = quantile(effect, probs = 0.975))
