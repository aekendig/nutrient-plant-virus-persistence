## Goal: summary results and figure of transmission from transmission_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(tidybayes)  # version used: 1.0.4
library(sjPlot)  # version used: 2.7.0

# import data
datp <- read_csv("./output/transmission_analysis_pav_data.csv")
datr <- read_csv("./output/transmission_analysis_rpv_data.csv")

# load models
load("./output/transmission_pav_up_concentration_all_informative.rda")
load("./output/transmission_rpv_up_concentration_all_informative.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


#### print model summaries ####

tab_model(mpcai)
summary(mpcai)
prior_summary(mpcai)

tab_model(mrcai)
summary(mrcai)
prior_summary(mrcai)


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
postr <- posterior_samples(mrcai)
postp <- posterior_samples(mpcai)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "conc", "N", "P", "N_t", "P_t", "co_conc", "NP", "NP_t", "co_N", "co_P", "co_N_t", "co_P_t", "conc_N", "conc_P", "conc_N_t", "conc_P_t", "co_NP", "co_NP_t", "conc_NP", "conc_NP_t", "co_conc_N", "co_conc_P", "co_conc_N_t", "co_conc_P_t", "co_conc_NP", "co_conc_NP_t", "sd_round", "sd_time", "round_1", "round_2", "round_3", "round_4", "time_1", "time_2", "time_3", "time_4", "time_5", "time_6", "time_7", "time_8", "lp")

# category average

combp <- postp %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgp <- combp %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

combr <- postr %>%
  transmute(l_l = exp(int)/(1 + exp(int)),
            l_l_co = exp(int + co)/(1 + exp(int + co)),
            l_n = exp(int + N_t)/(1 + exp(int + N_t)),
            l_n_co = exp(int + N_t + co + co_N_t)/(1 + exp(int + N_t + co + co_N_t)),
            l_p = exp(int + P_t)/(1 + exp(int + P_t)),
            l_p_co = exp(int + P_t + co + co_P_t)/(1 + exp(int + P_t + co + co_P_t)),
            l_b = exp(int + N_t+ P_t + NP_t)/(1 + exp(int + N_t+ P_t + NP_t)),
            l_b_co = exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N_t+ P_t + co + NP_t + co_N_t + co_P_t + co_NP_t)),
            n_l = exp(int + N)/(1 + exp(int + N)),
            n_l_co = exp(int + N + co + co_N)/(1 + exp(int + N + co + co_N)),
            n_n = exp(int + N + N_t)/(1 + exp(int + N + N_t)),
            n_n_co = exp(int + N + N_t + co + co_N + co_N_t)/(1 + exp(int + N + N_t + co + co_N + co_N_t)),
            n_p = exp(int + N + P_t)/(1 + exp(int + N + P_t)),
            n_p_co = exp(int + N + P_t + co + co_N + co_P_t)/(1 + exp(int + N + P_t + co + co_N + co_P_t)),
            n_b = exp(int + N + N_t+ P_t + NP_t)/(1 + exp(int + N + N_t+ P_t + NP_t)),
            n_b_co = exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + N_t+ P_t + co + co_N + NP_t + co_N_t + co_P_t + co_NP_t)),
            p_l = exp(int + P)/(1 + exp(int + P)),
            p_l_co = exp(int + P + co + co_P)/(1 + exp(int + P + co + co_P)),
            p_n = exp(int + P + N_t)/(1 + exp(int + P + N_t)),
            p_n_co = exp(int + P + N_t + co + co_P + co_N_t)/(1 + exp(int + P + N_t + co + co_P + co_N_t)),
            p_p = exp(int + P + P_t)/(1 + exp(int + P + P_t)),
            p_p_co = exp(int + P + P_t + co + co_P + co_P_t)/(1 + exp(int + P + P_t + co + co_P + co_P_t)),
            p_b = exp(int + P + N_t+ P_t + NP_t)/(1 + exp(int + P + N_t+ P_t + NP_t)),
            p_b_co = exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + P + N_t+ P_t + co + co_P + NP_t + co_N_t + co_P_t + co_NP_t)),
            b_l = exp(int + N + P + NP)/(1 + exp(int + N + P + NP)),
            b_l_co = exp(int + N + P + co + NP + co_N + co_P + co_NP)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP)),
            b_n = exp(int + N + P + NP + N_t)/(1 + exp(int + N + P + NP + N_t)),
            b_n_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t + co_N_t)),
            b_p = exp(int + N + P + NP + P_t)/(1 + exp(int + N + P + NP + P_t)),
            b_p_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + P_t + co_P_t)),
            b_b = exp(int + N + P + NP + N_t+ P_t + NP_t)/(1 + exp(int + N + P + NP + N_t+ P_t + NP_t)),
            b_b_co = exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)/(1 + exp(int + N + P + co + NP + co_N + co_P + co_NP + N_t+ P_t + NP_t + co_N_t + co_P_t + co_NP_t)))

avgr <- combr %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()


#### raw data figures ####

plotA <- datp %>%
  ggplot(aes(x = dpi, y = t_up, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), alpha = 0.5, aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("PAV transmission")

plotB <- datr %>%
  ggplot(aes(x = dpi, y = t_up, color = nutrient)) +
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
  scale_colour_manual(values = col_pal,
                      name = "Source plant\nnutrient") +
  scale_shape_manual(values = c(19, 21), name = "Source plant\ninfection") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Source plant\ninfection") +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("RPV transmission")

plotA_S <- datp %>%
  ggplot(aes(x = dpi, y = t_up, color = nutrient_t)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), alpha = 0.5, aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_text(color = "black", size = lg_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("PAV transmission")

plotB_S <- datr %>%
  ggplot(aes(x = dpi, y = t_up, color = nutrient_t)) +
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
        legend.key.size = unit(0.7, 'lines'),
        legend.box = "horizontal",
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal,
                      name = "Recipient plant\nnutrient") +
  scale_shape_manual(values = c(19, 21), name = "Source plant\ninfection") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Source plant\ninfection") +
  ylim(0, 1.07) + 
  xlab("Days post inoculation") +
  ylab("RPV transmission")


#### figure of category averages ####

plotC <- avgp %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  xlab("Recipient plant nutrient") +
  ylab("Est. PAV transmission")

plotD <- avgr %>%
  group_by(treatment, Nutrient, Nutrient_t, Inoculation) %>%
  mean_hdi() %>%
  ggplot(aes(x = Nutrient, y = effect)) +
  geom_pointinterval(aes(shape = Inoculation,  color = Nutrient), fatten_point = 2.5, size_range = c(0.3, 0.4), position = position_dodge(0.5), fill = "white") +
  facet_wrap(~Nutrient_t, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text.y = element_text(color = "black", size = sm_txt),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_y_continuous(breaks = c(0.0, 0.3, 0.6, 0.9)) +
  xlab("Recipient plant nutrient") +
  ylab("Est. RPV transmission")


#### combine plots ####

# extract legend
legB <- get_legend(plotB)
legB_S <- get_legend(plotB_S)

# combine plots
plots <- align_plots(plotA, plotB + theme(legend.position = "none"), plotC, plotD, align = 'v', axis = 'l')

plots_S <- align_plots(plotA_S, plotB_S + theme(legend.position = "none"), plotC, plotD, align = 'v', axis = 'l')

# combine top row
top_row <- cowplot::plot_grid(plots[[1]], plots[[2]],
                              labels = c("a", "b"), 
                              label_size = lg_txt, 
                              nrow = 1)

top_row_S <- cowplot::plot_grid(plots_S[[1]], plots_S[[2]],
                              labels = c("a", "b"), 
                              label_size = lg_txt, 
                              nrow = 1)

# combine bottom row
bottom_row <- cowplot::plot_grid(plots[[3]], plots[[4]], legB, 
                                 labels = c("c", "d"), 
                                 label_size = lg_txt, 
                                 rel_widths = c(1, 1, 0.3),
                                 nrow = 1,
                                 align = "h",
                                 axis = "t")


# combine all
plot <- cowplot::plot_grid(top_row, bottom_row, ncol = 1)

plot_S <- cowplot::plot_grid(top_row_S, legB_S, 
                             rel_heights = c(1, 0.1),
                             ncol = 1)

# print
pdf("./output/figure_3_transmission.pdf", width = 6, height = 4)
plot
dev.off()

pdf("./output/figure_S3_transmission.pdf", width = 6, height = 2.5)
plot_S
dev.off()


#### numbers for text ####

# model summaries
summary(mpuci)
summary(mruci)

# percentage change due to coinfection across all nutrient combinations
percp_co <- combp %>% 
  transmute(l.l.co = l_l_co - l_l,
            l.n.co = l_n_co - l_n,
            l.p.co = l_p_co - l_p,
            l.b.co = l_b_co - l_b,
            n.l.co = n_l_co - n_l,
            n.n.co = n_n_co - n_n,
            n.p.co = n_p_co - n_p,
            n.b.co = n_b_co - n_b,
            p.l.co = p_l_co - p_l,
            p.n.co = p_n_co - p_n,
            p.p.co = p_p_co - p_p,
            p.b.co = p_b_co - p_b,
            b.l.co = b_l_co - b_l,
            b.n.co = b_n_co - b_n,
            b.p.co = b_p_co - b_p,
            b.b.co = b_b_co - b_b) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

percr_co <- combr %>% 
  transmute(l.l.co = l_l_co - l_l,
            l.n.co = l_n_co - l_n,
            l.p.co = l_p_co - l_p,
            l.b.co = l_b_co - l_b,
            n.l.co = n_l_co - n_l,
            n.n.co = n_n_co - n_n,
            n.p.co = n_p_co - n_p,
            n.b.co = n_b_co - n_b,
            p.l.co = p_l_co - p_l,
            p.n.co = p_n_co - p_n,
            p.p.co = p_p_co - p_p,
            p.b.co = p_b_co - p_b,
            b.l.co = b_l_co - b_l,
            b.n.co = b_n_co - b_n,
            b.p.co = b_p_co - b_p,
            b.b.co = b_b_co - b_b) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

# percentage change due to nutrients across infection treatments
percp_nut <- combp %>% 
  transmute(n.l.co = n_l_co - l_l_co,
            n.l = n_l - l_l,
            p.l.co = p_l_co - l_l_co,
            p.l = p_l - l_l,
            b.l.co = b_l_co - l_l_co,
            b.l = b_l - l_l,
            l.n.co = l_n_co - l_l_co,
            l.n = l_n - l_l,
            l.p.co = l_p_co - l_l_co,
            l.p = l_p - l_l,
            l.b.co = l_b_co - l_l_co,
            l.b = l_b - l_l) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Coinoculation = case_when(substr(treatment, 5, 6) == "co" ~ "co",
                                   TRUE ~ "single"),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()

percr_nut <- combr %>% 
  transmute(n.l.co = n_l_co - l_l_co,
            n.l = n_l - l_l,
            p.l.co = p_l_co - l_l_co,
            p.l = p_l - l_l,
            b.l.co = b_l_co - l_l_co,
            b.l = b_l - l_l,
            l.n.co = l_n_co - l_l_co,
            l.n = l_n - l_l,
            l.p.co = l_p_co - l_l_co,
            l.p = l_p - l_l,
            l.b.co = l_b_co - l_l_co,
            l.b = l_b - l_l) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Coinoculation = case_when(substr(treatment, 5, 6) == "co" ~ "co",
                                   TRUE ~ "single"),
         Nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")),
         Nutrient_t = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         Nutrient_t = factor(Nutrient_t, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble()


# mean values in proportion change, comparing co-inoculated with singly-inoculated plants
percp_co %>%
  group_by(Nutrient, Nutrient_t, treatment) %>%
  mean_hdi() %>%
  data.frame()

percr_co %>%
  group_by(Nutrient, Nutrient_t, treatment) %>%
  mean_hdi() %>%
  data.frame()

# mean values in proportion change, comparing nutrient treatments
percp_nut %>%
  group_by(Coinoculation, Nutrient, Nutrient_t, treatment) %>%
  mean_hdi() %>%
  data.frame()

percr_nut %>%
  group_by(Coinoculation, Nutrient, Nutrient_t, treatment) %>%
  mean_hdi() %>%
  data.frame()

# transmission rates for each treatment
transp <- combp %>% 
  gather(key = "treatment", value = "transmission") %>%
  mutate(infection = case_when(substr(treatment, 5, 6) == "co" ~ "co",
                                   TRUE ~ "single"),
         source_nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                              substr(treatment, 1, 1) == "n" ~ "N",
                              substr(treatment, 1, 1) == "p" ~ "P",
                              substr(treatment, 1, 1) == "b" ~ "N+P"),
         source_nutrient = factor(source_nutrient, levels = c("low", "N", "P", "N+P")),
         recipient_nutrient = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                substr(treatment, 3, 3) == "n" ~ "N",
                                substr(treatment, 3, 3) == "p" ~ "P",
                                substr(treatment, 3, 3) == "b" ~ "N+P"),
         recipient_nutrient = factor(recipient_nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble() %>%
  group_by(infection, source_nutrient, recipient_nutrient, treatment) %>%
  mean_hdi() %>%
  mutate(virus = "PAV")

transr <- combr %>% 
  gather(key = "treatment", value = "transmission") %>%
  mutate(infection = case_when(substr(treatment, 5, 6) == "co" ~ "co",
                               TRUE ~ "single"),
         source_nutrient = case_when(substr(treatment, 1, 1) == "l" ~ "low",
                                     substr(treatment, 1, 1) == "n" ~ "N",
                                     substr(treatment, 1, 1) == "p" ~ "P",
                                     substr(treatment, 1, 1) == "b" ~ "N+P"),
         source_nutrient = factor(source_nutrient, levels = c("low", "N", "P", "N+P")),
         recipient_nutrient = case_when(substr(treatment, 3, 3) == "l" ~ "low",
                                        substr(treatment, 3, 3) == "n" ~ "N",
                                        substr(treatment, 3, 3) == "p" ~ "P",
                                        substr(treatment, 3, 3) == "b" ~ "N+P"),
         recipient_nutrient = factor(recipient_nutrient, levels = c("low", "N", "P", "N+P"))) %>%
  as_tibble() %>%
  group_by(infection, source_nutrient, recipient_nutrient, treatment) %>%
  mean_hdi() %>%
  mutate(virus = "RPV")

# combine
trans_dat <- full_join(transp, transr)

# save
write_csv(trans_dat, "./data/model_transmission_parameters.csv")
