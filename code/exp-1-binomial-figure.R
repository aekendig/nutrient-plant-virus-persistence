## Goal: figure of Experiment 1 inoculation and model estimates from exp-1-binomial-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(cowplot)
library(brms)
library(sjPlot)

# import data
datp <- read_csv("./output/exp-1-binomial-analysis-pav-data.csv")
datr <- read_csv("./output/exp-1-binomial-analysis-rpv-data.csv")

# load models
load("./output/exp-1-binomial-analysis-uninformative-pav.rda")
load("./output/exp-1-binomial-analysis-uninformative-rpv.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 1


#### print model summaries ####

tab_model(m.bu.p)
summary(m.bu.p)
prior_summary(m.bu.p)

tab_model(m.bu.r)
summary(m.bu.r)
prior_summary(m.bu.r)

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
postr <- posterior_samples(m.bu.r)
postp <- posterior_samples(m.bu.p)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "sd_time", "time_1", "time_2", "time_3", "time_4", "time_5", "time_6", "time_7", "time_8", "lp")

# category average

combp <- postp %>%
  transmute(low = exp(int)/(1 + exp(int)),
            n = exp(int + N)/(1 + exp(int + N)),
            p = exp(int + P)/(1 + exp(int + P)),
            coinf = exp(int + co)/(1 + exp(int + co)),
            n_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            p_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            b = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)))

avgp <- combp %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) %in% c("l", "c") ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combr <- postr %>%
  transmute(low = exp(int)/(1 + exp(int)),
            n = exp(int + N)/(1 + exp(int + N)),
            p = exp(int + P)/(1 + exp(int + P)),
            coinf = exp(int + co)/(1 + exp(int + co)),
            n_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            p_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            b = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)))

avgr <- combr %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) %in% c("l", "c") ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()


#### raw data figures ####

plotA <- datp %>%
  ggplot(aes(x = dpi, y = present, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), alpha = 0.5, aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.5, 0.95),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("PAV infection")

plotB <- datr %>%
  ggplot(aes(x = dpi, y = present, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), alpha = 0.5, aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.spacing.y = unit(-0.1, "mm"), 
        legend.key = element_rect(size = 0.5, color = "white"),
        legend.key.size = unit(0.7, 'lines')) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal, name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation") +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("RPV infection")


#### figure of category averages ####

plotC <- avgp %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  xlab("Nutrient") +
  ylab("Est. PAV infection")

plotD <- avgr %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  median_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = sm_txt),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  xlab("Nutrient") +
  ylab("Est. RPV infection")

#### combine plots ####

# extract legend
legB <- get_legend(plotB)

# combine plots
plots <- align_plots(plotA, plotB + theme(legend.position = "none"), plotC, plotD, align = 'v', axis = 'l')

# combine top row
top_row <- cowplot::plot_grid(plots[[1]], plots[[2]],
                              labels = c("A", "B"), 
                              label_size = lg_txt, 
                              nrow = 1)
# combine bottom row
bottom_row <- cowplot::plot_grid(plots[[3]], plots[[4]], legB, 
                                 labels = c("C", "D"), 
                                 label_size = lg_txt, 
                                 rel_widths = c(1, 1, 0.3),
                                 nrow = 1,
                                 align = "h",
                                 axis = "t")

# combine all
plot <- cowplot::plot_grid(top_row, bottom_row, ncol = 1)

# print
pdf("./output/exp-1-binomial-figure.pdf", width = 6, height = 4)
plot
dev.off()
