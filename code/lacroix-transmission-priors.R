## Goal: estimate priors for exp-1-transmission-analysis based on C. Lacroix's data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
pd <- read_csv("./output/lacroix-concentration-pav-data.csv") # from lacroix-concentration-priors.R
rd <- read_csv("./output/lacroix-concentration-rpv-data.csv") # from lacroix-concentration-priors.R
td <- read_csv("./data/NVE_TransmissionFromS1ToS2_ForModellingProject_12-4-14.csv")
# data collected by C. Lacroix


#### edit data ####

# edit transmission data
tdr <- td %>%
  filter(Group.S1.x %in% c("RPV", "PAVRPV") & InfectionSatus.RPV.S1 == 1) %>%
  left_join(rd) %>%
  select(LabelCode.S1, LabelCode.S2, Treatment.S1.x, Group.S1.x, Individual.S1.x, ID, InfectionSatus.RPV.S1, InfectionSatus.RPV.S2, conc) %>%
  rename(RPV_t = InfectionSatus.RPV.S2, nutrient = Treatment.S1.x) %>%
  mutate(co = ifelse(Group.S1.x == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))
  
# check that ID's match
filter(tdr, Individual.S1.x != ID)

tdp <- td %>%
  filter(Group.S1.x %in% c("PAV", "PAVRPV") & InfectionSatus.PAV.S1 == 1) %>%
  left_join(pd) %>%
  select(LabelCode.S1, LabelCode.S2, Treatment.S1.x, Group.S1.x, Individual.S1.x, ID, InfectionSatus.PAV.S1, InfectionSatus.PAV.S2, conc) %>%
  rename(PAV_t = InfectionSatus.PAV.S2, nutrient = Treatment.S1.x) %>%
  mutate(co = ifelse(Group.S1.x == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))

# check that ID's match
filter(tdp, Individual.S1.x != ID)

# concentration datasets
tdpc <- tdp %>%
  filter(!is.na(conc)) %>%
  mutate(conc_s = scale(conc))

tdrc <- tdr %>%
  filter(!is.na(conc)) %>%
  mutate(conc_s = scale(conc))


#### visualize data ####

tdrc %>%
  ggplot(aes(x = conc_s, y = RPV_t)) +
  geom_point() +
  facet_wrap(~nutrient)

tdpc %>%
  ggplot(aes(x = conc_s, y = PAV_t)) +
  geom_point() +
  facet_wrap(~nutrient)

tdr %>%
  ggplot(aes(x = nutrient, y = RPV_t, colour = as.factor(co))) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position = position_dodge(0.4), width = 0.1)

tdrc %>%
  ggplot(aes(x = nutrient, y = RPV_t, colour = as.factor(co))) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position = position_dodge(0.4), width = 0.1)
# similar relationships, higher values

tdp %>%
  ggplot(aes(x = nutrient, y = PAV_t, colour = as.factor(co))) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position = position_dodge(0.4), width = 0.1)

tdpc %>%
  ggplot(aes(x = nutrient, y = PAV_t, colour = as.factor(co))) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position = position_dodge(0.4), width = 0.1) # differences in averages, missing some data

tdpc %>% 
  group_by(co, nutrient) %>%
  summarize(n()) # no co P, low others

tdrc %>%
  group_by(co, nutrient) %>%
  summarize(n()) 


#### analyze treatments ####

# PAV model
mp <- brm(data = tdpc, family = bernoulli(),
              PAV_t ~ conc_s + co * high_N + high_P + high_N:high_P,
              prior <- c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)
plot(mp)
summary(mp)
save(mp, file = "./output/lacroix-transmission-pav.rda")

# RPV model 
mr <- brm(data = tdrc, family = bernoulli(),
          RPV_t ~ conc_s + co * high_N * high_P,
          prior <- c(prior(normal(0, 10), class = Intercept),
                     prior(normal(0, 10), class = b)),
          iter = 6000, warmup = 1000, chains = 3, cores = 2)
plot(mr)
summary(mr)
save(mr, file = "./output/lacroix-transmission-rpv.rda")
