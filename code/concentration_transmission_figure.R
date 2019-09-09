## Goal: figure of relationship between concentration and transmission from transmission_analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4

# import data
datp <- read_csv("./output/transmission_analysis_pav_data.csv")
datr <- read_csv("./output/transmission_analysis_rpv_data.csv")

# load models
load("./output/transmission_pav_up_concentration_informative.rda")
load("./output/transmission_rpv_up_concentration_informative.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


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


#### concentration-transmission figure ####

# PAV concentration
pconcp <- ggplot(datp, aes(x = conc)) +
  geom_point(aes(y = t_up, color = nutrient, shape = inoculation), position = position_jitter(width = 0, height = 0.03), alpha = 0.5) +
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
                      name = "Nutrient", guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("PAV density") +
  ylab("PAV transmission")

# RPV concentration
pconcr <- ggplot(datr, aes(x = conc)) +
  geom_point(aes(y = t_up, color = nutrient, shape = inoculation), position = position_jitter(width = 0, height = 0.03), alpha = 0.5) +
  geom_ribbon(data = datr_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.5, size = 0.5) +
  geom_line(data = datr_pred, aes(y = pred), color = "black") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.001, "mm")) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal,
                      name = "Source plant\n nutrient") +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  xlab("RPV density") +
  ylab("RPV transmission")

pconcleg <- get_legend(pconcr)

pconc <- plot_grid(pconcp, pconcr + theme(legend.position = "none"), pconcleg, 
                            labels = c("A", "B"), 
                            rel_widths = c(1, 1, 0.3),
                            label_size = lg_txt,
                            nrow = 1)

pdf("./output/figure_4_concentration_transmission.pdf", width = 6, height = 2.5)
pconc
dev.off()
