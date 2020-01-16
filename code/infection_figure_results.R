## Goal: summary results and figure of infection model estimates from infection_analysis.R

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
pdat <- read_csv("./output/infection_analysis_pav_data.csv")
rdat <- read_csv("./output/infection_analysis_rpv_data.csv")

# load models
load("./output/infection_analysis_uninformative_pav.rda")
load("./output/infection_analysis_uninformative_rpv.rda")


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
postp <- posterior_samples(m.bu.p)
postr <- posterior_samples(m.bu.r)

# treatments
trtp <- pdat %>%
  select(inoculation, nutrient, co, high_N, high_P) %>%
  unique()

trtr <- rdat %>%
  select(inoculation, nutrient, co, high_N, high_P) %>%
  unique()

# merge treatments and posterior samples
combp <- merge(trtp, postp, all = T)
combr <- merge(trtr, postr, all = T)

# category average
prevfun <- function(dat){
  dat2 <- dat %>%
    mutate(exp_val = exp(b_Intercept +
             b_high_N*high_N +
             b_high_P*high_P +
             b_co*co +
             `b_high_N:high_P`*high_N*high_P +
             `b_co:high_N`*co*high_N +
             `b_co:high_P`*co*high_P +
             `b_co:high_N:high_P`*co*high_N*high_P),
           prev = exp_val / (1 + exp_val))
  
  return(dat2)
} 

avgp <- prevfun(combp) %>%
  select(nutrient, inoculation, prev)  %>% 
  group_by(nutrient, inoculation) %>%
  mean_hdi() %>%
  ungroup()

avgr <- prevfun(combr) %>%
  select(nutrient, inoculation, prev)  %>% 
  group_by(nutrient, inoculation) %>%
  mean_hdi() %>%
  ungroup() 

# percentage change
percfun <- function(dat){
  
  dat2 <- prevfun(dat) %>%
    select(inoculation, nutrient, prev) %>%
    arrange(inoculation, nutrient) %>%
    mutate(treatment = paste(inoculation, nutrient, sep = "_"),
           count = rep(1:nrow(postp), nrow(trtp))) %>%
    select(-c(inoculation, nutrient)) %>%
    spread(key = treatment, value = prev) %>%
    transmute(high_N = single_N - single_low,
              high_P = single_P - single_low,
              high_NP = `single_N+P` - single_low,
              low_co = co_low - single_low,
              N_co = co_N - single_N,
              P_co = co_P - single_P,
              NP_co = `co_N+P` - `single_N+P`)%>%
    gather(key = "treatment", value = "perc") %>%
    mutate(inoculation = ifelse(grepl("co", treatment, fixed = T), "co", "single"),
           inoculation = factor(inoculation, levels = c("single", "co")),
           nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
           nutrient = factor(nutrient, levels = c("low", "N", "P", "N+P"))) %>%
    select(-treatment) %>%
    as_tibble() %>%
    group_by(inoculation, nutrient) %>%
    mean_hdi() %>%
    ungroup() 
  
  return(dat2)
} 

percp <- percfun(combp)
percr <- percfun(combr)



#### figure settings ####

# palettes
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")
line_pal = c("solid", "dashed")
shape_pal = c(19, 21)

# text sizes
sm_txt = 6
lg_txt = 8

# base figure
base_theme <- theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        plot.title = element_text(color = "black", size = lg_txt, hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(color = "black", size = sm_txt),
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "none")


#### figures of raw data ####

plotA <- ggplot(pdat, aes(x = dpi, y = present, color = nutrient)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation), fill = "white") +
  stat_summary(aes(linetype = inoculation), fun.y = "mean", geom = "line", position = position_dodge(0.6)) +
  base_theme +
  scale_size_manual(values = c(0.5, 0.5)) +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  scale_linetype_manual(values = line_pal) +
  xlab("Days post inoculation") +
  ylab("PAV infection prevalence")

plotB <- plotA %+% rdat %+%
  ylab("RPV infection prevalence")


#### figures of category averages ####

plotC <- ggplot(avgp, aes(x = nutrient, y = prev,  color = nutrient)) +
  geom_pointinterval(aes(shape = inoculation), fatten_point = 2.5, size_range = c(0.4, 0.6), position = position_dodge(0.3), fill = "white", show.legend = F) +
  base_theme +
  scale_colour_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal) +
  xlab("Nutrient") +
  ylab("Est. PAV infection prevalence")

plotD <- plotC %+% avgr %+%
  ylab("Est. RPV infection prevalence")



#### legend ####

leg <- get_legend(plotA %+% 
                    theme(legend.position = "right",
                          legend.spacing.y = unit(-0.1, "mm"),
                          legend.key = element_rect(size = 0.5, color = "white"),
                          legend.key.size = unit(0.7, 'lines'),
                          legend.title = element_text(color = "black", size = lg_txt),
                          legend.text = element_text(color = "black", size = sm_txt),
                          legend.background = element_blank(),
                          legend.key.width = unit(1, "cm")) %+%
                    scale_size_manual(values = c(0.5, 0.5), name = "Inoculation") +
                    scale_colour_manual(values = col_pal, name = "Nutrient") +
                    scale_shape_manual(values = c(19, 21), name = "Inoculation") +
                    scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation"))


#### combine plots ####

# combine plots
plots <- align_plots(plotA, plotB, plotC, plotD, align = 'v', axis = 'l')

# combine top row
top_row <- cowplot::plot_grid(plots[[1]], plots[[2]],
                              labels = c("a", "b"), 
                              label_size = lg_txt, 
                              nrow = 1)
# combine bottom row
bottom_row <- cowplot::plot_grid(plots[[3]], plots[[4]], leg, 
                                 labels = c("c", "d"), 
                                 label_size = lg_txt, 
                                 rel_widths = c(1, 1, 0.3),
                                 nrow = 1,
                                 align = "h",
                                 axis = "t")

# combine all
plot <- cowplot::plot_grid(top_row, bottom_row, ncol = 1)

# print
pdf("./output/figure_1_infection.pdf", width = 6, height = 4)
plot
dev.off()


#### print model summaries ####

tab_model(m.bu.p)
summary(m.bu.p)
prior_summary(m.bu.p)
percp

tab_model(m.bu.r)
summary(m.bu.r)
prior_summary(m.bu.r)
percr