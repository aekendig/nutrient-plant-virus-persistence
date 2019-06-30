## Goal: figure of Experiment 1 transmission and model estimates from exp-1-transmission-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(ggridges)
library(cowplot)

# import data
dat <- read_csv("./output/exp-1-transmission-analysis-data.csv")

# load models
load("./output/exp-1-transmission-pav-up-concentration-informative-priors.rda")
load("./output/exp-1-transmission-rpv-up-concentration-informative-priors.rda")


#### edit data ####

# inoculation column
dat <- dat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# split by virus
datp <- filter(dat, inoc != "RPV")
datr <- filter(dat, inoc != "PAV")

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
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

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
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

# check treatments
sloper %>% select(treatment, Inoculation, Nutrient, Transmission) %>% unique()

# fitted values for concentration
datp_pred <- datp %>%
  select(PAV_conc.mg, PAV_conc_s) %>%
  mutate(co = 1,
         high_N = 1,
         high_P = 1,
         high_N_t = 1,
         high_P_t = 1)

datp <- datp %>%
  mutate(pred = fitted(mpuci, newdata = datp_pred, re_formula = NA)[,1],
         lower = fitted(mpuci, newdata = datp_pred, re_formula = NA)[,3],
         upper = fitted(mpuci, newdata = datp_pred, re_formula = NA)[,4])

#### concentration figure ####

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2

# PAV
  ggplot(datp, aes(x = PAV_conc.mg, y = PAV_t_up, colour = nutrient, shape = inoculation)) +
    geom_point(position = position_jitter(width = 0, height = 0.03)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.5) +
    geom_line(aes(y = pred), color = "black") +
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
  ylab("ln(RPV concentration)")


#### figure of model estimates ####

# PAV single
plotC <- ggplot(filter(slopep, Inoculation == "single"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
  stat_density_ridges(data = filter(slopep, Inoculation == "coinfection"), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "solid") +
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
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  ylab("Density of posterior dist.")  +
  xlim(-0.7, 1.5)

# PAV coinfection
plotD <- ggplot(filter(slopep, Inoculation == "coinfection"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
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
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("") +
  ylab("Density of posterior dist.") +
  xlim(-0.9, 3.3)

# RPV single
plotE <- ggplot(filter(sloper, Inoculation == "single"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
  stat_density_ridges(data = filter(sloper, Inoculation == "coinfection"), alpha = 0, color = "white") +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "solid") +
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
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.2))) +
  xlab("") +
  ylab("Density of posterior dist.") +
  xlim(-0.7, 3.5)

# RPV coinfection
plotF <- ggplot(filter(sloper, Inoculation == "coinfection"), 
                aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
  geom_vline(xintercept = 0, color= "gray", size = 0.3) +
  stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black", size = sm_txt),
        axis.text.y = element_blank(),
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
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.2))) +
  xlab("") +
  ylab("Density of posterior dist.") +
  xlim(-1, 4)

#### combine plots ####

# combine
plot <- plot_grid(plotA, plotC, plotD, plotB, plotE, plotF, 
                  labels = c("A", "B", "C", "D", "E", "F"), 
                  label_size = lg_txt, 
                  rel_widths = c(1, 0.55, 0.45, 1, 0.55, 0.45),
                  label_x = c(0, 0.15, -0.05, 0, 0.15, -0.05))

# save with new x-axis labels
pdf("./output/exp-1-concentration-figure.pdf", width = 6, height = 4)
ggdraw(plot) +
  draw_label(label = "Proportion change in RPV concentration", x = 0.77, y = 0.04, size = lg_txt) +
  draw_label(label = "Proportion change in PAV concentration", x = 0.77, y = 0.54, size = lg_txt)
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
