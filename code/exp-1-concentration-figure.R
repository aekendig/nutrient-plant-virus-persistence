## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-concentration-analysis.R

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(sjPlot)

# import data
pdat <- read_csv("./output/exp-1-concentration-analysis-pav-data.csv")
rdat <- read_csv("./output/exp-1-concentration-analysis-rpv-data.csv")

# load models
load("./output/exp-1-concentration-analysis-log-informative-rpv.rda")
load("./output/exp-1-concentration-analysis-log-informative-pav.rda")

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2


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
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         pred = predict(m.li.p)[, 1])

rdat <- rdat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         pred = predict(m.li.r)[, 1])

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


#### figure of raw data ####
 
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


#### figure of category averages ####

plotC <- avgp %>%
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
  ylab("Est. ln(PAV density)")

plotD <- avgr %>%
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
  ylab("Est. ln(RPV density)")


#### figure of differences between treatments ####
# decided not to use this one because it's less intuitive than above

# PAV single
# plotC <- ggplot(filter(slopep, Inoculation == "single"), 
#        aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
#   stat_density_ridges(data = filter(slopep, Inoculation == "coinfection"), alpha = 0, color = "white") +
#   geom_vline(xintercept = 0, color= "gray", size = 0.3) +
#   stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "solid") +
#   theme_bw() +
#   theme(axis.title = element_text(color = "black", size = lg_txt),
#         axis.text = element_text(color = "black", size = sm_txt),
#         strip.text = element_blank(),
#         legend.title = element_text(color = "black", size = sm_txt),
#         legend.text = element_text(color = "black", size = sm_txt),
#         legend.position = c(0.25, 0.92),
#         legend.direction = "horizontal",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing.x = unit(0, "lines")) +
#   scale_fill_manual(values = col_pal, guide = F) +
#   scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
#   xlab("") +
#   ylab("Density of posterior dist.")  +
#   xlim(-0.7, 1.5)
# 
# # PAV coinfection
# plotD <- ggplot(filter(slopep, Inoculation == "coinfection"), 
#        aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
#   geom_vline(xintercept = 0, color= "gray", size = 0.3) +
#   stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "dashed") +
#   theme_bw() +
#   theme(axis.title.x = element_text(color = "black", size = lg_txt),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(color = "black", size = sm_txt),
#         axis.text.y = element_blank(),
#         strip.text = element_blank(),
#         legend.title = element_text(color = "black", size = sm_txt),
#         legend.text = element_text(color = "black", size = sm_txt),
#         legend.position = c(0.25, 0.92),
#         legend.direction = "horizontal",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing.x = unit(0, "lines")) +
#   scale_fill_manual(values = col_pal, guide = F) +
#   scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
#   xlab("") +
#   ylab("Density of posterior dist.") +
#   xlim(-0.9, 3.3)
# 
# # RPV single
# plotE <- ggplot(filter(sloper, Inoculation == "single"), 
#        aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
#   stat_density_ridges(data = filter(sloper, Inoculation == "coinfection"), alpha = 0, color = "white") +
#   geom_vline(xintercept = 0, color= "gray", size = 0.3) +
#   stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "solid") +
#   theme_bw() +
#   theme(axis.title = element_text(color = "black", size = lg_txt),
#         axis.text = element_text(color = "black", size = sm_txt),
#         strip.text = element_blank(),
#         legend.title = element_text(color = "black", size = sm_txt),
#         legend.text = element_text(color = "black", size = sm_txt),
#         legend.position = c(0.25, 0.92),
#         legend.direction = "horizontal",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing.x = unit(0, "lines")) +
#   scale_fill_manual(values = col_pal, guide = F) +
#   scale_y_discrete(expand = expand_scale(add = c(0.2, 1.2))) +
#   xlab("") +
#   ylab("Density of posterior dist.") +
#   xlim(-0.7, 3.5)
# 
# # RPV coinfection
# plotF <- ggplot(filter(sloper, Inoculation == "coinfection"), 
#        aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient)) +
#   geom_vline(xintercept = 0, color= "gray", size = 0.3) +
#   stat_density_ridges(alpha = 0.7, rel_min_height = 0.005, quantile_lines = T, quantiles = c(0.025, 0.5, 0.975), linetype = "dashed") +
#   theme_bw() +
#   theme(axis.title.x = element_text(color = "black", size = lg_txt),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(color = "black", size = sm_txt),
#         axis.text.y = element_blank(),
#         strip.text = element_blank(),
#         legend.title = element_text(color = "black", size = sm_txt),
#         legend.text = element_text(color = "black", size = sm_txt),
#         legend.position = c(0.25, 0.92),
#         legend.direction = "horizontal",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing.x = unit(0, "lines")) +
#   scale_fill_manual(values = col_pal, guide = F) +
#   scale_y_discrete(expand = expand_scale(add = c(0.2, 1.2))) +
#   xlab("") +
#   ylab("Density of posterior dist.") +
#   xlim(-1, 4)


#### combine plots ####

# extract legend
legB <- get_legend(plotB + theme(legend.spacing.y = unit(-0.1, "mm"), 
                                 legend.key = element_rect(size = 0.5, color = "white"),
                                 legend.key.size = unit(0.7, 'lines')))

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
pdf("./output/exp-1-concentration-figure.pdf", width = 6, height = 4)
plot
dev.off()


# plot <- cowplot::plot_grid(plotA, plotC, plotD, plotB, plotE, plotF, 
#                   labels = c("A", "B", "C", "D", "E", "F"), 
#                   label_size = lg_txt, 
#                   rel_widths = c(1, 0.55, 0.45, 1, 0.55, 0.45),
#                   label_x = c(0, 0.15, -0.05, 0, 0.15, -0.05))

# save with new x-axis labels
# pdf("./output/exp-1-concentration-figure.pdf", width = 6, height = 4)
# ggdraw(plot) +
#   draw_label(label = "Proportion change in RPV concentration", x = 0.77, y = 0.04, size = lg_txt) +
#   draw_label(label = "Proportion change in PAV concentration", x = 0.77, y = 0.54, size = lg_txt)
# dev.off()

#### numbers for text ####

# model summaries
summary(m.li.r)
summary(m.li.p)

# mean values in proportion change
sloper %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()

slopep %>%
  group_by(treatment, Inoculation, Nutrient) %>%
  mean_hdi()

