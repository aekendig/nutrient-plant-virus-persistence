## Goal: estimate priors for exp-1-concentration-analysis based on C. Lacroix's data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
pd <- read_csv("./data/SamplesVirusTiterData_08.14.16_PAV_bis.csv") 
rd <- read_csv("./data/SamplesVirusTiterData_08.14.16_RPV_bis.csv") 
st <- read_csv("./data/lacroix-qPCR-std-limits.csv")
ld <- read_csv("./data/NVE_InoculationOfS1_ForModellingProject_12-4-14.csv")
# quantities estimated and data compiled by C. Lacroix
# AK compiled st data from CL's qPCR Excel sheets


#### edit data ####

# check for infection in controls and healthy oats
pd %>%
  filter(SampleName %in% c("NTC", "NAC", "HO")) %>%
  select(SampleName, Ct, Quantity_Adj) %>%
  data.frame() # NAC's

rd %>%
  filter(SampleName %in% c("NTC", "NAC", "HO")) %>%
  select(SampleName, Ct, Quantity_Adj) %>%
  data.frame() # NAC's

# check for high values
pd %>%
  full_join(st) %>%
  filter(Quantity_Adj > PAVmax) # 0

rd %>%
  full_join(st) %>%
  filter(Quantity_Adj > RPVmax) # 0

# remove controls, healthy oats, and RPV-only samples
# summarise by sample
# convert quantity to per ul
# replace values below the standard curve with 0
# add in sample info
pd2 <- pd %>%
  mutate(Ct = as.numeric(Ct)) %>%
  filter(grepl("PAV", SampleName, fixed = T) == T) %>%
  left_join(st) %>%
  group_by(Run, SampleName, TargetName, PAVmin) %>%
  summarise(tech_cycle = mean(Ct, na.rm = T),
            quant = mean(Quantity_Adj, na.rm = T) / 2.5) %>%
  mutate(quant_adj = ifelse(quant < PAVmin | is.na(quant), 0, quant),
         nutrient = case_when(grepl("Ctrl", SampleName, fixed = T) == T ~ "Ctrl",
                              grepl("NPPAV", SampleName, fixed = T) == T | grepl("PAVNP", SampleName, fixed = T) == T | grepl("PAVRPVNP", SampleName, fixed = T) == T ~ "NP",
                              grepl("NPAV", SampleName, fixed = T) == T | grepl("PAVN", SampleName, fixed = T) == T | grepl("PAVRPVN", SampleName, fixed = T) == T ~ "N", 
                              grepl("PPAV", SampleName, fixed = T) == T | grepl("PAVP", SampleName, fixed = T) == T | grepl("PAVRPVP", SampleName, fixed = T) == T ~ "P"),
         inoc = case_when(grepl("RPV", SampleName, fixed = T) == T ~ "PAVRPV",
                          TRUE ~ "PAV"),
         ID = gsub("[^0-9.]", "", SampleName) %>% as.numeric(),
         ID = ifelse(ID == 473472, 473, ID),
         LabelCode.S1 = paste(nutrient, inoc, ID, sep = " ")) %>%
  left_join(select(ld, c(LabelCode.S1, Extraction.WetWeighTaken.mg.S1)))

# check treatments
pd2 %>%
  select(SampleName, nutrient, inoc, ID, LabelCode.S1) %>%
  data.frame()

sum(is.na(pd2$Extraction.WetWeighTaken.mg.S1))

# same as above for RPV
rd2 <- rd %>%
  mutate(Ct = as.numeric(Ct)) %>%
  filter(grepl("RPV", SampleName, fixed = T) == T) %>%
  left_join(st) %>%
  group_by(Run, SampleName, TargetName, RPVmin) %>%
  summarise(tech_cycle = mean(Ct, na.rm = T),
            quant = mean(Quantity_Adj, na.rm = T) / 2.5) %>%
  mutate(quant_adj = ifelse(quant < RPVmin | is.na(quant), 0, quant),
         nutrient = case_when(grepl("Ctrl", SampleName, fixed = T) == T ~ "Ctrl",
                              grepl("NPRPV", SampleName, fixed = T) == T | grepl("RPVNP", SampleName, fixed = T) == T | grepl("PAVRPVNP", SampleName, fixed = T) == T | grepl("NPPAV", SampleName, fixed = T) == T ~ "NP",
                              grepl("NRPV", SampleName, fixed = T) == T | grepl("RPVN", SampleName, fixed = T) == T | grepl("PAVRPVN", SampleName, fixed = T) == T | grepl("NPAV", SampleName, fixed = T) == T ~ "N", 
                              grepl("PRPV", SampleName, fixed = T) == T | grepl("RPVP", SampleName, fixed = T) == T | grepl("PAVRPVP", SampleName, fixed = T) == T | grepl("PPAV", SampleName, fixed = T) == T ~ "P"),
         inoc = case_when(grepl("PAV", SampleName, fixed = T) == T ~ "PAVRPV",
                          TRUE ~ "RPV"),
         ID = gsub("[^0-9.]", "", SampleName) %>% as.numeric(),
         LabelCode.S1 = paste(nutrient, inoc, ID, sep = " ")) %>%
  left_join(select(ld, c(LabelCode.S1, Extraction.WetWeighTaken.mg.S1)))

# check treatments
rd2 %>%
  select(SampleName, nutrient, inoc, LabelCode.S1) %>%
  data.frame()

sum(is.na(rd2$Extraction.WetWeighTaken.mg.S1))

# modify data for models
p.d <- pd2 %>%
  filter(quant_adj > 0) %>%
  mutate(conc = quant_adj / Extraction.WetWeighTaken.mg.S1,
         log_conc = log(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))
write_csv(p.d, "./output/lacroix-concentration-pav-data.csv")

r.d <- rd2 %>%
  filter(quant_adj > 0) %>%
  mutate(conc = quant_adj / Extraction.WetWeighTaken.mg.S1,
         log_conc = log(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))
write_csv(r.d, "./output/lacroix-concentration-rpv-data.csv")


#### visualize data ####

p.d %>%
  ggplot(aes(x = nutrient, y = log_conc, colour = inoc)) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.4), width = 0.1)
# can compare N and control for single infection
# can compare coinfection and single infection for N

r.d %>%
  ggplot(aes(x = nutrient, y = log_conc, colour = inoc)) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.4), width = 0.1)
# can compare all nutrients for single infection
# can compare coinfection and single infection for ctrl, N, N+P


#### analyze treatments ####

# PAV N model
m.n.p <- brm(data = filter(p.d, inoc == "PAV" & high_P == 0), family = gaussian,
              log_conc ~ high_N,
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.n.p)
summary(m.n.p)
save(m.n.p, file = "./output/lacroix-concentration-n-pav.rda")

# PAV co model
m.c.p <- brm(data = filter(p.d, nutrient == "N"), family = gaussian,
             log_conc ~ co,
             prior <- c(prior(normal(0, 100), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.c.p)
summary(m.c.p)
save(m.c.p, file = "./output/lacroix-concentration-co-pav.rda")

# RPV nutrient model
m.n.r <- brm(data = filter(r.d, inoc == "RPV"), family = gaussian,
             log_conc ~ high_N * high_P,
             prior <- c(prior(normal(0, 100), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.n.r)
summary(m.n.r)
save(m.n.r, file = "./output/lacroix-concentration-n-rpv.rda")

# RPV coinfection model
m.c.r <- brm(data = filter(r.d, high_P == 0), family = gaussian,
             log_conc ~ co * high_N,
             prior <- c(prior(normal(0, 100), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.c.r)
summary(m.c.r)
save(m.c.r, file = "./output/lacroix-concentration-co-rpv.rda")
