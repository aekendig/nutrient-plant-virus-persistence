## Goal: summary results and figure of virus concentration and model estimates from concentration_analysis.R

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
pdat <- read_csv("./output/concentration_analysis_pav_data.csv")
rdat <- read_csv("./output/concentration_analysis_rpv_data.csv")

# load models
load("./output/concentration_analysis_informative_pav.rda")
load("./output/concentration_analysis_informative_rpv.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8


#### print model summaries ####

tab_model(m.li.p)
summary(m.li.p)
prior_summary(m.li.p)

tab_model(m.li.r)
summary(m.li.r)
prior_summary(m.li.r)


#### edit data ####

# inoculation column, nutrient column, predicted values (tried including these in the time series)
pdat <- pdat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

rdat <- rdat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postr <- posterior_samples(m.li.r)
postp <- posterior_samples(m.li.p)

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


#### figures of raw data ####
 
# PAV (concentration over time)
plotA <- ggplot(pdat, aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(aes(linetype = inoculation), fun.y = "mean", geom = "line", position = position_dodge(0.6)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        #legend.position = c(0.42, 0.92),
        legend.background = element_blank(),
        legend.key = element_blank(),
        #legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  scale_size_manual(values = c(0.5, 0.5), guide = F) +
  scale_colour_manual(values = col_pal,
                      name = "Nutrient", guide = F) +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("Days post inoculation") +
  ylab("ln(PAV density)")

# RPV (concentration over time)
plotB <- ggplot(rdat, aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation), show.legend = F) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        #legend.position = c(0.34, 0.93),
        legend.background = element_blank(),
        legend.key = element_blank(),
        #legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_size_manual(values = c(0.5, 0.5), name = "Inoculation") +
  scale_colour_manual(values = col_pal, name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation") +
  xlab("Days post inoculation") +
  ylab("ln(RPV density)")


#### figures of category averages ####

plotC <- avgp %>%
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
  ylab("Est. ln(PAV density)")

plotD <- avgr %>%
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
  ylab("Est. ln(RPV density)")


#### combine plots ####

# extract legend
legB <- get_legend(plotB + theme(legend.spacing.y = unit(-0.1, "mm"), 
                                 legend.key = element_rect(size = 0.5, color = "white"),
                                 legend.key.size = unit(0.7, 'lines')))

# combine plots
plots <- align_plots(plotA, plotB + theme(legend.position = "none"), plotC, plotD, align = 'v', axis = 'l')

# combine top row
top_row <- cowplot::plot_grid(plots[[1]], plots[[2]],
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

# print
pdf("./output/figure_1_concentration.pdf", width = 6, height = 4)
plot
dev.off()


#### numbers for text ####

# model summaries
summary(m.li.p)
summary(m.li.r)

# mean values in proportion change
slopep %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()

sloper %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()



