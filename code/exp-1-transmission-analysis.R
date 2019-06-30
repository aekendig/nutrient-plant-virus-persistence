## Goal: Integrate transmisison data and qPCR data from experiment 1 and analyze treatment effects

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(brms)
library(tidyverse)
library(bayesplot)
library(cowplot)

# import data
qdatp <- read_csv("./output/exp-1-qPCR-analysis-pav-data.csv")
qdatr <- read_csv("./output/exp-1-qPCR-analysis-rpv-data.csv")
tdat <- read_csv("./data/transmission-data.csv") # came from MeanqPCR_VCEBaseDataset_Trans_ManUp_123016.csv (unnecessary columns and rows removed)
sdat <- read_csv("./data/sample-exp-molc-data.csv")

# import models
load("./output/lacroix-transmission-pav.rda")
load("./output/lacroix-transmission-rpv.rda")

# functions
plot_intervals <- function(data) {
  ggplot(data, aes(y = parameter, yend = parameter)) + 
    geom_vline(xintercept = 0) +
    geom_segment(aes(x = ll, xend = hh), size = 1) + 
    geom_segment(aes(x = l, xend = h), size = 2) +
    geom_point(aes(x = m), size = 3, color = "red") +
    theme_bw() +
    facet_wrap(~model)
}
# modified from Tristan Mahr (https://www.tjmahr.com/ggplot2-how-to-do-nothing/)


#### edit qPCR data ####

# add aphid mass to dataset and calculate new concentration metrics
sdat2 <- sdat %>%
  filter(material == "shoot") %>%
  select(time, inoc, nutrient, round, replicate, aphid_mass.g) %>%
  mutate(aphid_mass.mg = aphid_mass.g * 1000,
         aphid_mass.mg = replace_na(aphid_mass.mg, 0)) %>%
  select(-aphid_mass.g)

qdatp2 <- qdatp %>%
  left_join(sdat2) %>%
  mutate(quant_t = conc * aphid_mass.mg)

qdatr2 <- qdatr %>%
  left_join(sdat2) %>%
  mutate(quant_t = conc * aphid_mass.mg)

# replicates were combined in round 1 for transmission trial - average their concentration values and combine the amount of tissue provided

# subset data for round 1 plants and average over replicates

qdatpR1 <- qdatp2 %>%
  filter(round == 1) %>%
  group_by(time, inoc, nutrient, round) %>%
  mutate(total_aphid_mass.mg = sum(aphid_mass.mg),
         conc2 = conc * (aphid_mass.mg / total_aphid_mass.mg)) %>%
  summarise(aphid_mass.mg = unique(total_aphid_mass.mg),
            conc = sum(conc2),
            quant_t = sum(quant_t)) %>%
  ungroup() %>%
  mutate(replicate = 1)

qdatrR1 <- qdatr2 %>%
  filter(round == 1) %>%
  group_by(time, inoc, nutrient, round) %>%
  mutate(total_aphid_mass.mg = sum(aphid_mass.mg),
         conc2 = conc * (aphid_mass.mg / total_aphid_mass.mg)) %>%
  summarise(aphid_mass.mg = unique(total_aphid_mass.mg),
            conc = sum(conc2),
            quant_t = sum(quant_t)) %>%
  ungroup() %>%
  mutate(replicate = 1)

# recombine with rest of data
qdatp3 <- qdatp2 %>%
  filter(round != 1) %>%
  full_join(qdatpR1)

qdatr3 <- qdatr2 %>%
  filter(round != 1) %>%
  full_join(qdatrR1)


#### edit transmission data ####

# examine data
unique(tdat$material)
tdat %>% filter(is.na(material)) %>% data.frame() # remove these
unique(tdat$PAV_t)
unique(tdat$RPV_t)
unique(tdat$inoc)
unique(tdat$nutrient)

# remove root samples and unknown samples, create ID column
tdat2 <- tdat %>%
  filter(material == "shoot") %>%
  mutate(ID = paste(round, time, inoc, nutrient, replicate, nutrient_t, sep = "."))

# check duplicates
tdat2 %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  data.frame()
# 1.2.coinfection.N+P.2.N+P: different extraction dates
# 1.8.RPV.N.2.low: take out the one with PAV_t = 0.5 (checked gel)
# 2.4.PAV.P.3.N: take out the one with RPV_t = 0 (checked gel)
# 2.7.RPV.N.3.N: take out the one with PAV_t = 0 (checked gel)

tdat3 <- tdat2 %>%
  filter(!(ID == "1.2.coinfection.N+P.2.N+P" & ext_date_t == "3/23/16")) %>%
  filter(!(ID == "1.8.RPV.N.2.low" & PAV_t == 0.5)) %>%
  filter(!(ID == "2.4.PAV.P.3.N" & RPV_t == 0)) %>%
  filter(!(ID == "2.7.RPV.N.3.N" & PAV_t == 0))


#### combine qPCR and transmission data ####

# only keep overlapping data
datp <- inner_join(qdatp3, tdat3)
datr <- inner_join(qdatr3, tdat3)

# make N and P columns for transmission
# coinfection column for source plant
# transmission columns
# 0 concentration and quantity if missing
# scale concentration and quantity
datp <- datp %>%
  mutate(high_N = ifelse(nutrient %in% c("N", "N+P"), 1, 0),
         high_P = ifelse(nutrient %in% c("P", "N+P"), 1, 0),
         high_N_t = ifelse(nutrient_t %in% c("N", "N+P"), 1, 0),
         high_P_t = ifelse(nutrient_t %in% c("P", "N+P"), 1, 0),
         co = ifelse(inoc == "coinfection", 1, 0),
         t_up = ceiling(PAV_t),
         t_dn = floor(PAV_t),
         conc_s = scale(conc),
         quant_s = scale(quant_t))

datr <- datr %>%
  mutate(high_N = ifelse(nutrient %in% c("N", "N+P"), 1, 0),
         high_P = ifelse(nutrient %in% c("P", "N+P"), 1, 0),
         high_N_t = ifelse(nutrient_t %in% c("N", "N+P"), 1, 0),
         high_P_t = ifelse(nutrient_t %in% c("P", "N+P"), 1, 0),
         co = ifelse(inoc == "coinfection", 1, 0),
         t_up = ceiling(PAV_t),
         t_dn = floor(PAV_t),
         conc_s = scale(conc),
         quant_s = scale(quant_t))

# check for NA's
datp %>% 
  select(t_up, t_dn, conc_s, quant_s, high_N, high_P, high_N_t, high_P_t, co) %>%
  is.na() %>%
  colSums()

datr %>% 
  select(t_up, t_dn, conc_s, quant_s, high_N, high_P, high_N_t, high_P_t, co) %>%
  is.na() %>%
  colSums()


#### visualize ####

datp %>%
  ggplot(aes(x = round, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # 3 lower for single

datr %>%
  ggplot(aes(x = round, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # 1 higher for coinfection

datp %>%
  ggplot(aes(x = as.factor(time), y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # increases with time

datp %>%
  ggplot(aes(x = as.factor(time), y = conc_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # so does conc

cor.test(datp$conc_s, datp$time) # 0.22

datp %>%
  ggplot(aes(x = as.factor(time), y = quant_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # and quantity

cor.test(datp$quant_s, datp$time) #0.25

datr %>%
  ggplot(aes(x = as.factor(time), y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # increases with time

datr %>%
  ggplot(aes(x = as.factor(time), y = conc_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) # so does conc

cor.test(datr$conc_s, datr$time) # 0.17

datr %>%
  ggplot(aes(x = as.factor(time), y = quant_s)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_wrap(~inoc) 

cor.test(datr$conc_s, datr$time) # 0.17

datp %>%
  ggplot(aes(x = nutrient, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_grid(nutrient_t ~ inoc)

datr %>%
  ggplot(aes(x = nutrient, y = t_up)) +
  stat_summary(geom = "point", size = 2, fun.y = "mean") +
  stat_summary(geom = "errorbar", width = 0.1, fun.data = "mean_cl_boot") +
  facet_grid(nutrient_t ~ inoc)


#### replicates ####

datp %>%
  group_by(co, nutrient, nutrient_t) %>%
  summarize(n()) %>%
  data.frame()

datr %>%
  group_by(co, nutrient, nutrient_t) %>%
  summarize(n()) %>%
  data.frame()


#### statistical models ####

# can't do autoregressive models with bernoulli

# mpuc <- brm(data = datp, family = bernoulli,
#            t_up ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time),
#            prior = c(prior(normal(0, 10), class = Intercept),
#                      prior(normal(0, 10), class = b),
#                      prior(cauchy(0, 1), class = sd)),
#            iter = 6000, warmup = 1000, chains = 3, cores = 2,
#            control = list(adapt_delta = 0.99))
# summary(mpuc)
# save(mpuc, file = "./output/exp-1-transmission-pav-up-concentration.rda")
load("./output/exp-1-transmission-pav-up-concentration.rda")

# mpuq <- update(mpuc, formula. = t_up ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp)
# summary(mpuq)
# save(mpuq, file = "./output/exp-1-transmission-pav-up-quantity.rda")
load("./output/exp-1-transmission-pav-up-quantity.rda")

# mruc <- update(mpuc, formula. = t_up ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) +(1|time), newdata = datr)
# summary(mruc)
# save(mruc, file = "./output/exp-1-transmission-rpv-up-concentration.rda")
load("./output/exp-1-transmission-rpv-up-concentration.rda")

# mruq <- update(mruc, formula. = t_up ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr)
# summary(mruq)
# save(mruq, file = "./output/exp-1-transmission-rpv-up-quantity.rda")
load("./output/exp-1-transmission-rpv-up-quantity.rda")

# mpdc <- update(mpuc, formula. = t_dn ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp, control = list(adapt_delta = 0.99999))
# summary(mpdc)
# save(mpdc, file = "./output/exp-1-transmission-pav-down-concentration.rda")
load("./output/exp-1-transmission-pav-down-concentration.rda")

# mpdq <- update(mpdc, formula. = t_dn ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datp, control = list(adapt_delta = 0.9999999)) 
# summary(mpdq)
# save(mpdq, file = "./output/exp-1-transmission-pav-down-quantity.rda")
load("./output/exp-1-transmission-pav-down-quantity.rda")

# mrdc <- update(mruc, formula. = t_dn ~ conc_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr)
# summary(mrdc)
# save(mrdc, file = "./output/exp-1-transmission-rpv-down-concentration.rda")
load("./output/exp-1-transmission-rpv-down-concentration.rda")

# mrdq <- update(mruc, formula. = t_dn ~ quant_s + co * (high_N * high_P + high_N_t * high_P_t) + (1|round) + (1|time), newdata = datr)
# summary(mrdq)
# save(mrdq, file = "./output/exp-1-transmission-rpv-down-quantity.rda")
load("./output/exp-1-transmission-rpv-down-quantity.rda")


#### check models ####

# convergence of chains
plot(mpuc)
plot(mpuq)
plot(mpdc)

plot(mruc)
plot(mruq)
plot(mrdc)
plot(mrdq)

# posterior predictive check
pp_check(mpuc, nsamples = 50)
pp_check(mpuq, nsamples = 50)
pp_check(mpdc, nsamples = 50)
pp_check(mpdq, nsamples = 50)

pp_check(mruc, nsamples = 50)
pp_check(mruq, nsamples = 50)
pp_check(mrdc, nsamples = 50)
pp_check(mrdq, nsamples = 50)

# compare concentration and quantity
loo_mpu <- list(loo(mpuc), loo(mpuq))
loo_mpu
loo::loo_compare(loo_mpu) # pretty much equal

loo_mpd <- list(loo(mpdc), loo(mpdq))
loo_mpd
loo::loo_compare(loo_mpd) # equal

loo_mru <- list(loo(mruc), loo(mruq))
loo_mru
loo::loo_compare(loo_mru) # equal

loo_mrd <- list(loo(mrdc), loo(mrdq))
loo_mrd
loo::loo_compare(loo_mrd) # equal

# compare estimates for up and down
params = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_high_N_t", "b_high_P_t", "b_high_N:high_P", "b_high_N_t:high_P_t", "b_co:high_N", "b_co:high_P", "b_co:high_N_t", "b_co:high_P_t", "b_co:high_N:high_P", "b_co:high_N_t:high_P_t")

apuc <- as.array(mpuc) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "PAV up conc uninformative")

apdc <- as.array(mpdc) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "PAV down conc uninformative")

plot_intervals(full_join(apuc, apdc)) 
# rounding up leads to larger effect sizes

aruc <- as.array(mruc) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "RPV up conc uninformative")

ardc <- as.array(mrdc) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "RPV down conc uninformative")

plot_intervals(full_join(aruc, ardc))
# very similar, up a little more precise


#### models with informative priors ####

# compare estimates for prior model and model without priors
summary(mp)
aplc <- (as.array(mp)) %>%
  mcmc_intervals_data(pars = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_co:high_N", "b_high_N:high_P")) %>%
  mutate(model = "PAV Lacroix")
plot_intervals(full_join(apuc, aplc))
# some estimates in the opposite direction, but error bars generally larger

summary(mr)
arlc <- (as.array(mr)) %>%
  mcmc_intervals_data(pars = c("b_Intercept", "b_conc_s", "b_co", "b_high_N", "b_high_P", "b_co:high_N", "b_co:high_P", "b_high_N:high_P", "b_co:high_N:high_P")) %>%
  mutate(model = "RPV Lacroix")
plot_intervals(full_join(aruc, arlc))
# estimates closer to zero

# PAV model
summary(mp)

mpuci <- update(mpuc,
                prior = c(prior(normal(1.53, 0.60), class = Intercept),
                          prior(normal(-0.26, 0.28), class = b, coef = conc_s),
                          prior(normal(6.05, 4.17), class = b, coef = co),
                          prior(normal(0.69, 0.94), class = b, coef = high_N),
                          prior(normal(0.63, 1.42), class = b, coef = high_P),
                          prior(normal(-6.51, 4.23), class = b, coef = co:high_N),
                          prior(normal(0, 10), class = b, coef = co:high_P),
                          prior(normal(-0.20, 1.63), class = b, coef = high_N:high_P),
                          prior(normal(0, 10), class = b, coef = co:high_N:high_P),
                          prior(normal(0, 10), class = b, coef = high_N_t),
                          prior(normal(0, 10), class = b, coef = high_P_t),
                          prior(normal(0, 10), class = b, coef = co:high_N_t),
                          prior(normal(0, 10), class = b, coef = co:high_P_t),
                          prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
                          prior(normal(0, 10), class = b, coef = co:high_N_t:high_P_t),
                          prior(cauchy(0, 1), class = sd)))

plot(mpuci)

save(mpuci, file = "./output/exp-1-transmission-pav-up-concentration-informative-priors.rda")

apuci <- as.array(mpuci) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "PAV up conc informative")

plot_intervals(full_join(apuc, apuci))

# RPV model
summary(mr)

mruci <- update(mruc,
                prior = c(prior(normal(0.54, 0.30), class = Intercept),
                          prior(normal(0.19, 0.13), class = b, coef = conc_s),
                          prior(normal(-0.19, 0.44), class = b, coef = co),
                          prior(normal(0.62, 0.42), class = b, coef = high_N),
                          prior(normal(0.13, 0.71), class = b, coef = high_P),
                          prior(normal(0.33, 0.60), class = b, coef = co:high_N),
                          prior(normal(-0.10, 1.12), class = b, coef = co:high_P),
                          prior(normal(-0.01, 0.86), class = b, coef = high_N:high_P),
                          prior(normal(-0.60, 1.28), class = b, coef = co:high_N:high_P),
                          prior(normal(0, 10), class = b, coef = high_N_t),
                          prior(normal(0, 10), class = b, coef = high_P_t),
                          prior(normal(0, 10), class = b, coef = co:high_N_t),
                          prior(normal(0, 10), class = b, coef = co:high_P_t),
                          prior(normal(0, 10), class = b, coef = high_N_t:high_P_t),
                          prior(normal(0, 10), class = b, coef = co:high_N_t:high_P_t),
                          prior(cauchy(0, 1), class = sd)))

plot(mruci)

save(mruci, file = "./output/exp-1-transmission-rpv-up-concentration-informative-priors.rda")

aruci <- as.array(mruci) %>%
  mcmc_intervals_data(pars = params) %>%
  mutate(model = "RPV up conc informative")

plot_intervals(full_join(aruc, aruci))


#### save data for plotting ####

# save file
write_csv(datp, "./output/exp-1-transmission-analysis-pav-data.csv")
write_csv(datr, "./output/exp-1-transmission-analysis-rpv-data.csv")

