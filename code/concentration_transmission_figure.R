## Goal: figure of relationship between concentration and transmission from transmission_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(sjPlot)  # version used: 2.7.0

# import data
datp <- read_csv("./output/transmission_analysis_pav_data.csv")
datr <- read_csv("./output/transmission_analysis_rpv_data.csv")

# load models
load("./output/transmission_pav_up_concentration_all_informative.rda")
load("./output/transmission_rpv_up_concentration_all_informative.rda")


#### edit data ####

# inoculation and nutrient columns
datp <- datp %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"),
         log_conc = log(conc))

datr <- datr %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"),
         log_conc = log(conc))

# fitted values for concentration
datp_pred <- datp %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(min_conc = min(log_conc),
            max_conc = max(log_conc)) %>%
  expand_grid(i = 1:1000) %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  mutate(log_conc = seq(unique(min_conc), unique(max_conc), length.out = 1000),
         conc = exp(log_conc),
         conc_s = (conc - mean(datp$conc)) / sd(datp$conc)) %>%
  ungroup() %>%
  select(-c(min_conc, max_conc, i))

datp_pred <- datp_pred %>%
  cbind(fitted(mpcai, newdata = datp_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datr_pred <- datr %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  summarize(min_conc = min(log_conc),
            max_conc = max(log_conc)) %>%
  expand_grid(i = 1:1000) %>%
  group_by(inoculation, co, nutrient, nutrient_t, high_N, high_P, high_N_t, high_P_t) %>%
  mutate(log_conc = seq(unique(min_conc), unique(max_conc), length.out = 1000),
         conc = exp(log_conc),
         conc_s = (conc - mean(datp$conc)) / sd(datp$conc)) %>%
  ungroup() %>%
  select(-c(min_conc, max_conc, i))

datr_pred <- datr_pred %>%
  cbind(fitted(mrcai, newdata = datr_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

# check that concentration scaling is consistent
datp %>%
  select(conc, conc_s) %>%
  mutate(type = "data") %>%
  full_join(datp_pred %>%
              select(conc, conc_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc, y = conc_s, color = type)) +
  geom_line()

datr %>%
  select(conc, conc_s) %>%
  mutate(type = "data") %>%
  full_join(datr_pred %>%
              select(conc, conc_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc, y = conc_s, color = type)) +
  geom_line()


#### concentration-transmission figure ####


# palettes
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")
line_pal = c("solid", "dashed")

# text sizes
sm_txt = 6
lg_txt = 8

# PAV concentration
plotp <- ggplot(datp, aes(x = log_conc)) +
  geom_point(aes(y = t_up, color = nutrient, fill = interaction(nutrient, inoculation)), shape = 21, alpha = 0.3, size = 0.75) +
  geom_line(data = datp_pred, aes(y = pred, color = nutrient, linetype = inoculation), size = 0.6) +
  facet_wrap(~nutrient_t, nrow = 1) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        plot.title = element_text(color = "black", size = lg_txt, hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  scale_linetype_manual(values = line_pal, guide = F) +
  xlab("ln(PAV density)") +
  ylab("PAV transmission") +
  ggtitle("Recipient plant nutrient") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# RPV concentration
plotr <- ggplot(datr, aes(x = log_conc)) +
  geom_point(aes(y = t_up, color = nutrient, fill = interaction(nutrient, inoculation)), shape = 21, alpha = 0.3, size = 0.75) +
  geom_line(data = datr_pred, aes(y = pred, color = nutrient, linetype = inoculation), size = 0.6) +
  facet_wrap(~nutrient_t, nrow = 1) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  scale_linetype_manual(values = line_pal, guide = F) +
  xlab("ln(RPV density)") +
  ylab("RPV transmission") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# legends
leg <- ggplot(datr, aes(x = log_conc, y = t_up)) +
  geom_point(aes(color = nutrient, fill = inoculation), shape = 21) +
  geom_line(aes(color = nutrient, linetype = inoculation)) +
  theme(legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(1.5, "lines"),
        legend.spacing.x = unit(0.001, "mm"),
        legend.box = "horizontal",
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, name = "Source plant\nnutrient") +
  scale_linetype_manual(values = line_pal, name = "Source plant\ninfection") +
  scale_fill_manual(values = c(col_pal[1], "white"), name = "Source plant\ninfection") +
  guides(color = guide_legend(override.aes = list(shape = 19)))

plotl <- get_legend(leg)

# combine plots
concplots <- cowplot::plot_grid(plotp, plotr,
                                labels = c("a", "b"),
                                label_size = lg_txt,
                                nrow = 2,
                                hjust = -1,
                                vjust = 2.5)

pconc <- cowplot::plot_grid(concplots, plotl, 
                            rel_heights = c(1, 0.08), 
                            nrow = 2)


pdf("./output/figure_4_concentration_transmission.pdf", height = 4, width = 6)
pconc
dev.off()


#### values for text ####

# PAV: change with low nutrients, 

# posterior samples
postr <- posterior_samples(mrcai)
postp <- posterior_samples(mpcai)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "conc", "N", "P", "N_t", "P_t", "co_conc", "NP", "NP_t", "co_N", "co_P", "co_N_t", "co_P_t", "conc_N", "conc_P", "conc_N_t", "conc_P_t", "co_NP", "co_NP_t", "conc_NP", "conc_NP_t", "co_conc_N", "co_conc_P", "co_conc_N_t", "co_conc_P_t", "co_conc_NP", "co_conc_NP_t", "sd_round", "sd_time", "round_1", "round_2", "round_3", "round_4", "time_1", "time_2", "time_3", "time_4", "time_5", "time_6", "time_7", "time_8", "lp")

# category average

combp <- postp %>%
  transmute(l_l = exp(conc)/(1 + exp(conc)),
            l_l_co = exp(conc + co_conc)/(1 + exp(conc + co_conc)),
            l_n = exp(conc + conc_N_t)/(1 + exp(conc + conc_N_t)),
            l_n_co = exp(conc + conc_N_t + co_conc + co_conc_N_t)/(1 + exp(conc + conc_N_t + co_conc + co_conc_N_t)))

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
  as_tibble() %>%
  group_by(Inoculation, Nutrient, Nutrient_t, treatment) %>%
  mean_hdi() %>%
  data.frame()

avgp

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
