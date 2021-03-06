## Goal: estimate priors for transmission analysis based on C. Lacroix's data

# This script is for illustration only. You will not be able to run it without access to the directories lacroix_data and lacroix_output


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(brms) # version used: 2.7.0

# import data
pd <- read_csv("./output/lacroix_output/lacroix_concentration_pav_data.csv") # from lacroix_concentration_priors.R
rd <- read_csv("./output/lacroix_output/lacroix_concentration_rpv_data.csv") # from lacroix_concentration_priors.R
td <- read_csv("./data/lacroix_data/NVE_TransmissionFromS1ToS2_ForModellingProject_12-4-14.csv")
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
# too few points to estimate concentration effect within P

tdpc %>%
  ggplot(aes(x = conc_s, y = PAV_t)) +
  geom_point() +
  facet_wrap(~nutrient)
# too few points to estimate concentration effect within P, maybe N too

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
save(mp, file = "./output/lacroix_output/lacroix_transmission_pav.rda")

# RPV model 
mr <- brm(data = tdrc, family = bernoulli(),
          RPV_t ~ conc_s + co * high_N * high_P,
          prior <- c(prior(normal(0, 10), class = Intercept),
                     prior(normal(0, 10), class = b)),
          iter = 6000, warmup = 1000, chains = 3, cores = 2)
plot(mr)
summary(mr)
save(mr, file = "./output/lacroix_output/lacroix_transmission_rpv.rda")
