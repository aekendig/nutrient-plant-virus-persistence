## Goal: figure to compare models with full and truncated datasets from concentration_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4

# load models
load("./output/concentration_analysis_informative_pav.rda")
load("./output/concentration_analysis_informative_rpv.rda")
load("./output/concentration_analysis_informative_late_pav.rda")
load("./output/concentration_analysis_informative_late_rpv.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


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
  mean_hdi() %>%
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
  ylab("Est. ln(PAV density)") +
  ggtitle("Full dataset")

plotC <- avgr %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
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
  ylab("Est. ln(RPV density)")+
  ylim(6.9, 9.4)

plotB <- avgpd %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        plot.title = element_text(color = "black", size = lg_txt, hjust = 0.5),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(-0.1, "mm"), 
        legend.key = element_rect(size = 0.5, color = "white"),
        legend.key.size = unit(0.7, 'lines')) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  xlab("Nutrient") +
  ylab("Est. ln(PAV density)") +
  ggtitle("Later dpi dataset")

plotD <- avgrd %>%
  mutate(Inoculation = recode(Inoculation, coinfection = "co")) %>%
  group_by(treatment, Nutrient, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect,  color = Nutrient)) +
  geom_point(aes(shape = Inoculation), position = position_dodge(0.3)) +
  geom_pointinterval(aes(shape = Inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = lg_txt),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(-0.1, "mm"), 
        legend.key = element_rect(size = 0.5, color = "white"),
        legend.key.size = unit(0.7, 'lines')) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  xlab("Nutrient") +
  ylab("Est. ln(RPV density)")+
  ylim(6.9, 9.4)


#### combine plots ####

# extract legends
legB <- get_legend(plotB)
legD <- get_legend(plotD)

# combine
top_row <- cowplot::plot_grid(plotA, plotB + theme(legend.position = "none"), legB, 
                  labels = c("a", "b"), 
                  nrow = 1,
                  label_size = lg_txt, 
                  rel_widths = c(1, 1, 0.3),
                  label_x = c(0, -0.03))

bottom_row <- cowplot::plot_grid(plotC, plotD + theme(legend.position = "none"), legD, 
                              labels = c("c", "d"), 
                              nrow = 1,
                              label_size = lg_txt, 
                              rel_widths = c(1, 1, 0.3),
                              label_x = c(0, -0.03))

plot <- cowplot::plot_grid(top_row, bottom_row, ncol = 1)

# print
pdf("./output/figure_S3_truncated_dataset.pdf", width = 6, height = 4)
plot
dev.off()


