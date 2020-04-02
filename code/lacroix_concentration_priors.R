## Goal: estimate priors for concentration analysis based on C. Lacroix's data

# This script is for illustration only. You will not be able to run it without access to the directories lacroix_data and lacroix_output


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(brms) # version used: 2.7.0

# import data
pd <- read_csv("./data/lacroix_data/SamplesVirusTiterData_08.14.16_PAV_bis.csv") 
rd <- read_csv("./data/lacroix_data/SamplesVirusTiterData_08.14.16_RPV_bis.csv") 
st <- read_csv("./data/lacroix_data/lacroix_qPCR_std_limits.csv")
ld <- read_csv("./data/lacroix_data/NVE_InoculationOfS1_ForModellingProject_12-4-14.csv")
# quantities estimated and data compiled by C. Lacroix
# AK compiled st data from CL's qPCR Excel sheets


#### edit data ####

# check for infection in controls and healthy oats
pd %>%
  filter(SampleName %in% c("NTC", "NAC", "HO")) %>%
  select(SampleName, Ct, Quantity_Adj) %>%
  data.frame() # NAC's (contain enzyme, not used in main experiment)

rd %>%
  filter(SampleName %in% c("NTC", "NAC", "HO")) %>%
  select(SampleName, Ct, Quantity_Adj) %>%
  data.frame() # NAC's (contain enzyme, not used in main experiment)

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
            quant = mean(Quantity_Adj, na.rm = T)) %>%
  mutate(quant_adj = ifelse(quant < PAVmin | is.na(quant), 0, quant),
         quant_ul = quant_adj / 2.5,
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
            quant = mean(Quantity_Adj, na.rm = T)) %>%
  mutate(quant_adj = ifelse(quant < RPVmin | is.na(quant), 0, quant),
         quant_ul = quant_adj  / 2.5,
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
# select values with quantities 
# calculate concentration and add columns
p.d <- pd2 %>%
  filter(quant_ul > 0) %>%
  mutate(conc = quant_ul * 50 / Extraction.WetWeighTaken.mg.S1,
         log_conc = log(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))
write_csv(p.d, "./output/lacroix_output/lacroix_concentration_pav_data.csv")

r.d <- rd2 %>%
  filter(quant_ul > 0) %>%
  mutate(conc = quant_ul * 50 / Extraction.WetWeighTaken.mg.S1,
         log_conc = log(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))
write_csv(r.d, "./output/lacroix_output/lacroix_concentration_rpv_data.csv")


#### visualize ####

p.d %>%
  ggplot(aes(x = nutrient, y = log_conc, colour = inoc)) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.4), width = 0.1)
# can compare all nutrients for single infection
# can compare coinfection and single infection for N

p.d %>%
  group_by(inoc, nutrient) %>%
  summarise(n = n())

r.d %>%
  ggplot(aes(x = nutrient, y = log_conc, colour = inoc)) +
  stat_summary(fun.y = mean, geom = "point", position = position_dodge(0.4)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.4), width = 0.1)
# can compare all nutrients for single infection
# can compare coinfection and single infection for all nutrients

r.d %>%
  group_by(inoc, nutrient) %>%
  summarise(n = n())


#### analyze treatments ####

# subset PAV data
p.d.mod <- p.d %>%
  filter(!(inoc == "PAVRPV" & nutrient %in% c("Ctrl", "NP")))

# PAV model
m.p <- brm(data = p.d.mod, family = gaussian,
              log_conc ~ high_N * high_P + co,
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.p)
summary(m.p)
save(m.p, file = "./output/lacroix_output/lacroix_concentration_pav.rda")

# RPV model
m.r <- update(m.p, log_conc ~ co * high_N * high_P, newdata = r.d)
plot(m.r)
summary(m.r)
save(m.r, file = "./output/lacroix_output/lacroix_concentration_rpv.rda")
