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
datpco <- read_csv("./output/transmission_analysis_coinfected_pav_data.csv")
datrco <- read_csv("./output/transmission_analysis_coinfected_rpv_data.csv")

# load models
load("./output/transmission_pav_up_concentration_informative.rda")
load("./output/transmission_rpv_up_concentration_informative.rda")
load("./output/transmission_coinfected_pav_up_concentration_informative.rda")
load("./output/transmission_coinfected_rpv_up_concentration_informative.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


#### edit data ####

# inoculation and nutrient columns
datp <- datp %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"))

datr <- datr %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"))

datpco <- datpco %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"))

datrco <- datrco %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         nutrient_t = fct_relevel(nutrient_t, "low", "N", "P"))

# average values by treatment
datp_avg <- datp %>%
  group_by(nutrient, nutrient_t, inoculation) %>%
  summarise(conc_avg = mean(conc),
            conc_se = sd(conc)/sqrt(length(conc)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up))) %>%
  mutate(treatment = paste(nutrient, nutrient_t, inoculation, sep = ".")) %>%
  ungroup()

datr_avg <- datr %>%
  group_by(nutrient, nutrient_t, inoculation) %>%
  summarise(conc_avg = mean(conc),
            conc_se = sd(conc)/sqrt(length(conc)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up))) %>%
  mutate(treatment = paste(nutrient, nutrient_t, inoculation, sep = ".")) %>%
  ungroup()

datpco_avg <- datpco %>%
  group_by(nutrient, nutrient_t, inoculation) %>%
  summarise(conc_avg = mean(conc_r),
            conc_se = sd(conc_r)/sqrt(length(conc_r)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up))) %>%
  mutate(treatment = paste(nutrient, nutrient_t, inoculation, sep = ".")) %>%
  ungroup()

datrco_avg <- datrco %>%
  group_by(nutrient, nutrient_t, inoculation) %>%
  summarise(conc_avg = mean(conc_p),
            conc_se = sd(conc_p)/sqrt(length(conc_p)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up))) %>%
  mutate(treatment = paste(nutrient, nutrient_t, inoculation, sep = ".")) %>%
  ungroup()

# extreme concentration values
datp_avg %>%
  mutate(top = conc_avg + conc_se,
         bottom = conc_avg - conc_se) %>%
  summarise(min = min(bottom),
            max = max(top))

datr_avg %>%
  mutate(top = conc_avg + conc_se,
         bottom = conc_avg - conc_se) %>%
  summarise(min = min(bottom),
            max = max(top))

datpco_avg %>%
  mutate(top = conc_avg + conc_se,
         bottom = conc_avg - conc_se) %>%
  summarise(min = min(bottom),
            max = max(top))

datrco_avg %>%
  mutate(top = conc_avg + conc_se,
         bottom = conc_avg - conc_se) %>%
  summarise(min = min(bottom),
            max = max(top))

# fitted values for concentration
datp_pred <- tibble(conc = seq(0, max(datp$conc), length.out = 1000)) %>%
  mutate(conc_s = (conc - mean(datp$conc)) / sd(datp$conc),
         co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 1,
         high_P_t = 0)

datp_pred <- datp_pred %>%
  cbind(fitted(mpuci, newdata = datp_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datr_pred <- tibble(conc = seq(0, max(datr$conc), length.out = 1000)) %>%
  mutate(conc_s = (conc - mean(datr$conc)) / sd(datr$conc),
         co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 1,
         high_P_t = 0)

datr_pred <- datr_pred %>%
  cbind(fitted(mruci, newdata = datr_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datpco_pred <- tibble(conc_r = seq(0, max(datpco$conc_r), length.out = 1000)) %>%
  mutate(conc_r_s = (conc_r - mean(datpco$conc_r)) / sd(datpco$conc_r),
         conc_s = 0,
         co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 1,
         high_P_t = 0)

datpco_pred <- datpco_pred %>%
  cbind(fitted(mpco, newdata = datpco_pred, re_formula = NA, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5)

datrco_pred <- tibble(conc_p = seq(0, max(datrco$conc_p), length.out = 1000)) %>%
  mutate(conc_p_s = (conc_p - mean(datrco$conc_p)) / sd(datrco$conc_p),
         conc_s = 0,
         co = 0,
         high_N = 0,
         high_P = 0,
         high_N_t = 1,
         high_P_t = 0)

datrco_pred <- datrco_pred %>%
  cbind(fitted(mrco, newdata = datrco_pred, re_formula = NA, nsamples = 100)) %>%
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

datpco %>%
  select(conc_r, conc_r_s) %>%
  mutate(type = "data") %>%
  full_join(datpco_pred %>%
              select(conc_r, conc_r_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc_r, y = conc_r_s, color = type)) +
  geom_line()

datrco %>%
  select(conc_p, conc_p_s) %>%
  mutate(type = "data") %>%
  full_join(datrco_pred %>%
              select(conc_p, conc_p_s) %>%
              mutate(type = "pred")) %>%
  ggplot(aes(x = conc_p, y = conc_p_s, color = type)) +
  geom_line()


#### concentration-transmission figure ####

# PAV concentration
plotp_raw <- ggplot(datp, aes(x = conc)) +
  geom_point(aes(y = t_up, color = nutrient, shape = nutrient_t, fill = interaction(nutrient, inoculation))) +
  geom_ribbon(data = datp_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = datp_pred, aes(y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  xlab("PAV density") +
  ylab("PAV transmission") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

plotp_avg <- ggplot(datp_avg) +
  geom_point(aes(x = conc_avg, y = t_up_avg, color = nutrient, shape = nutrient_t, fill = interaction(nutrient, inoculation))) +
  geom_errorbar(aes(x = conc_avg, ymin = t_up_avg - t_up_se, ymax = t_up_avg + t_up_se, color = nutrient, group = treatment), width = 0) +
  geom_errorbarh(aes(y = t_up_avg, xmin = conc_avg - conc_se, xmax = conc_avg + conc_se, color = nutrient, group = treatment), height = 0) +
  geom_ribbon(data = filter(datp_pred, conc < 950 & conc > 50), aes(x = conc, ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = filter(datp_pred, conc < 950 & conc > 50), aes(x = conc, y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# RPV concentration
plotr_raw <- ggplot(datr) +
  geom_point(aes(x = conc, y = t_up, color = nutrient, shape = nutrient_t, fill = interaction(nutrient, inoculation))) +
  geom_ribbon(data = datr_pred, aes(x = conc, ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = datr_pred, aes(x = conc, y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  xlab("RPV density") +
  ylab("RPV transmission") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

plotr_avg <- ggplot(datr_avg) +
  geom_point(aes(x = conc_avg, y = t_up_avg, color = nutrient, shape = nutrient_t, fill = interaction(nutrient, inoculation))) +
  geom_errorbar(aes(x = conc_avg, ymin = t_up_avg - t_up_se, ymax = t_up_avg + t_up_se, color = nutrient, group = treatment), width = 0) +
  geom_errorbarh(aes(y = t_up_avg, xmin = conc_avg - conc_se, xmax = conc_avg + conc_se, color = nutrient, group = treatment), height = 0) +
  geom_ribbon(data = filter(datr_pred, conc > 2850 & conc < 29650), aes(x = conc, ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = filter(datr_pred, conc > 2850 & conc < 29650), aes(x = conc, y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_fill_manual(values = c(col_pal, rep("white", 4)), guide = F) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# coinfected PAV
plotpco_raw <- ggplot(datpco, aes(x = conc_r)) +
  geom_point(aes(y = t_up, color = nutrient, shape = nutrient_t), fill = "white") +
 geom_ribbon(data = datpco_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
 geom_line(data = datpco_pred, aes(y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  xlab("RPV density") +
  ylab("PAV transmission") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

plotpco_avg <- ggplot(datpco_avg) +
  geom_point(aes(x = conc_avg, y = t_up_avg, color = nutrient, shape = nutrient_t), fill = "white") +
  geom_errorbar(aes(x = conc_avg, ymin = t_up_avg - t_up_se, ymax = t_up_avg + t_up_se, color = nutrient, group = treatment), width = 0) +
  geom_errorbarh(aes(y = t_up_avg, xmin = conc_avg - conc_se, xmax = conc_avg + conc_se, color = nutrient, group = treatment), height = 0) +
  geom_ribbon(data = filter(datpco_pred, conc_r < 61150 & conc_r > 9500), aes(x = conc_r, ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = filter(datpco_pred, conc_r < 61150 & conc_r > 9500), aes(x = conc_r, y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# coinfected RPV
plotrco_raw <- ggplot(datrco, aes(x = conc_p)) +
  geom_point(aes(y = t_up, color = nutrient, shape = nutrient_t), fill = "white") +
  geom_ribbon(data = datrco_pred, aes(ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = datrco_pred, aes(y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  xlab("PAV density") +
  ylab("RPV transmission") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

plotrco_avg <- ggplot(datrco_avg) +
  geom_point(aes(x = conc_avg, y = t_up_avg, color = nutrient, shape = nutrient_t), fill = "white") +
  geom_errorbar(aes(x = conc_avg, ymin = t_up_avg - t_up_se, ymax = t_up_avg + t_up_se, color = nutrient, group = treatment), width = 0) +
  geom_errorbarh(aes(y = t_up_avg, xmin = conc_avg - conc_se, xmax = conc_avg + conc_se, color = nutrient, group = treatment), height = 0) +
  geom_ribbon(data = filter(datrco_pred, conc_p < 1100 & conc_p > 50), aes(x = conc_p, ymin = lower, ymax = upper), fill = "black", color = NA, alpha = 0.3, size = 0.5) +
  geom_line(data = filter(datrco_pred, conc_p < 1100 & conc_p > 50), aes(x = conc_p, y = pred), color = "gray45", linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(21, 22, 23, 24), guide = F) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))

# legends
leg <- ggplot(filter(datr)) +
  geom_point(aes(x = conc, y = t_up, color = nutrient, shape = nutrient_t, fill = inoculation)) +
  theme(legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.background = element_blank(),
        legend.key = element_rect(color = "white", size = 0.5),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0.001, "mm"),
        legend.box = "horizontal",
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_manual(values = col_pal, name = "Source plant\nnutrient") +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Recipient plant\nnutrient") +
  scale_fill_manual(values = c(col_pal[1], "white"), name = "Source plant\ninfection") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

plotl <- get_legend(leg)

# combine plots
plotp <- plotp_raw + annotation_custom(ggplotGrob(plotp_avg), xmin=1200, xmax=3200, ymin=0.03, ymax=0.6)

plotr <- plotr_raw + annotation_custom(ggplotGrob(plotr_avg), xmin=60000, xmax=160000, ymin=0.03, ymax=0.6)

plotpco <- plotpco_raw + annotation_custom(ggplotGrob(plotpco_avg), xmin=60000, xmax=160000, ymin=0.03, ymax=0.6)

plotrco <- plotrco_raw + annotation_custom(ggplotGrob(plotrco_avg), xmin=1200, xmax=3200, ymin=0.03, ymax=0.6)

concplots <- cowplot::plot_grid(plotp, plotpco, plotrco, plotr,
                                labels = c("a", "b", "c", "d"),
                                label_size = lg_txt,
                                nrow = 2)

pconc <- cowplot::plot_grid(concplots, plotl, 
                            rel_heights = c(1, 0.08), 
                            nrow = 2)


pdf("./output/figure_4_concentration_transmission.pdf", height = 6, width = 6)
pconc
dev.off()


#### values for text ####

# range of RPV transmission rates
datr_pred %>%
  filter(conc %in% c(min(datr_pred$conc), max(datr_pred$conc))) %>%
  select(conc, pred)

# observed RPV transmission for coinfection
datr %>%
  group_by(co) %>%
  summarise(conc_avg = mean(conc),
            conc_se = sd(conc)/sqrt(length(conc)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up)))

datp %>%
  group_by(co) %>%
  summarise(conc_avg = mean(conc),
            conc_se = sd(conc)/sqrt(length(conc)),
            t_up_avg = mean(t_up),
            t_up_se = sd(t_up)/sqrt(length(t_up)))
