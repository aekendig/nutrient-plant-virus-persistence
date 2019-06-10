## Goal: estimate priors for exp-1-analysis based on C. Lacroix's data


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
         log_conc = log10(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))

r.d <- rd2 %>%
  filter(quant_adj > 0) %>%
  mutate(conc = quant_adj / Extraction.WetWeighTaken.mg.S1,
         log_conc = log10(conc),
         co = ifelse(inoc == "PAVRPV", 1, 0),
         high_N = ifelse(nutrient == "N" | nutrient == "NP", 1, 0),
         high_P = ifelse(nutrient == "P" | nutrient == "NP", 1, 0))


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
             conc ~ high_N,
             prior <- c(prior(normal(0, 1e3), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.n.p)
summary(m.n.p)
save(m.n.p, file = "./output/lacroix-n-pav.rda")

# PAV co model
m.c.p <- brm(data = filter(p.d, nutrient == "N"), family = gaussian,
             conc ~ co,
             prior <- c(prior(normal(0, 1e3), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.c.p)
summary(m.c.p)
save(m.c.p, file = "./output/lacroix-co-pav.rda")

# RPV nutrient model
m.n.r <- brm(data = filter(r.d, inoc == "RPV"), family = gaussian,
             conc ~ high_N * high_P,
             prior <- c(prior(normal(0, 1e4), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.n.r)
summary(m.n.r)
save(m.n.r, file = "./output/lacroix-n-rpv.rda")

# RPV coinfection model
m.c.r <- brm(data = filter(r.d, high_P == 0), family = gaussian,
             conc ~ co * high_N,
             prior <- c(prior(normal(0, 1e4), class = Intercept),
                        prior(normal(0, 10), class = b)),
             iter = 6000, warmup = 1000, chains = 3, cores = 2,
             control = list(adapt_delta = 0.99))
plot(m.c.r)
summary(m.c.r)
save(m.c.r, file = "./output/lacroix-co-rpv.rda")

## Goal: evaluate treatment effects on virus concentration in Experiment 1


#### set up ####

# import data
source("./code/exp-1-qPCR-raw-data-processing.R")
# clears environment
# loads tidyverse
# sets working directory to data folder

# load libraries
library(brms)
library(tidybayes)
library(bayesplot)

# clear all except dataset
rm(list = setdiff(ls(), c("dat")))

# load models for priors
load("./output/lacroix-n-pav.rda")
load("./output/lacroix-co-pav.rda")
load("./output/lacroix-n-rpv.rda")
load("./output/lacroix-co-rpv.rda")


#### edit data ####

# days post inoculation
dpi <- tibble(
  time = 1:8,
  dpi = c(5, 8, 12, 16, 19, 22, 26, 29)
)

# remove samples:
# poor standard curve efficiency
# quantities below standard curve, but greater than 1e3 (not sure if these should be zeros or not; standards removed if contamination had higher concentration)
# multiple qPCR tests of the same sample and the sample wasn't detected in one or had the higher variance in detected in multiple
# specific cases: low volume, known contamination, mis-labelling
dat <- dat %>%
  filter(remove == 0 & material == "shoot") %>%
  mutate(quant_adj = case_when(target == "PAV" & quant_adj < PAVmin ~ 0,
                               target == "RPV" & quant_adj < RPVmin ~ 0,
                               is.na(quant_adj) ~ 0,
                               TRUE ~ quant_adj),
         quant_zero = case_when(quant_adj == 0 ~ 1,
                                TRUE ~ 0)) %>%
  full_join(dpi)


# remove values above standard curve
dat <- dat %>%
  filter((target == "RPV" & quant_adj <= RPVmax) | (target == "PAV" & quant_adj <= PAVmax))


#### average technical replicates ####

dat2 <- dat %>%
  group_by(target, dpi, time, inoc, high_N, high_P, nutrient, round, replicate, sample, shoot_mass.g, root_mass.g, leaf_area.mm2, leaves, mass_ext.mg, PAVmin, PAVint, PAVslope, RPVmin, RPVint, RPVslope, q_group) %>%
  summarise(tech_cycle = mean(cycle, na.rm = T)) %>%
  mutate(quant = case_when(
    target == "PAV" ~ 10 ^ ((tech_cycle - PAVint) / PAVslope),
    target == "RPV" ~ 10 ^ ((tech_cycle - RPVint) / RPVslope)),
    quant_adj = case_when(target == "PAV" & quant < PAVmin ~ 0,
                          target == "RPV" & quant < RPVmin ~ 0,
                          is.na(quant) ~ 0,
                          TRUE ~ quant),
    quant_zero = case_when(quant_adj == 0 ~ 1,
                           TRUE ~ 0))



#### overall average titer for successful infections ####

# edit data
d.at <- dat2 %>%
  filter(quant_zero == 0 &
           inoc %in% c("PAV", "coinfection", "RPV")) %>%
  mutate(co = ifelse(inoc == "coinfection", 1, 0),
         log_quant = log10(quant_adj),
         exp_round = round,
         conc = quant_adj/mass_ext.mg,
         log_conc = log10(conc))



## PAV model ##

# data
d.at.p <- d.at %>%
  filter(target == "PAV" & inoc != "RPV")

# concentration model
m.ac.p <- brm(data = d.at.p, family = gaussian,
              conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(450, 100), class = Intercept),
                         prior(normal(0, 100), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# save model
save(m.ac.p, file = "./output/average-concentration-pav.rda")

# summary
summary(m.ac.p) # coinfection increases concentration
plot(m.ac.p)

# simulate data based on posterior and compare to observed
s.ac.p <- posterior_predict(m.ac.p)
color_scheme_set("brightblue")
ppc_dens_overlay(d.at.p$conc, s.ac.p[1:100,])
# including negative numbers


# with more specific priors
# use the estimates for coinfection and N addition from CL's models for the mean, multiply the standard error from her models by 10
m.acp.p <- update(m.ac.p, 
                  prior = c(prior(normal(1.34, 9.97), class = b, coef = co),
                            prior(normal(-5.95, 9.47), class = b, coef = high_N),
                            prior(normal(0, 10), class = b, coef = high_P),
                            prior(normal(0, 10), class = b, coef = co:high_N),
                            prior(normal(0, 10), class = b, coef = co:high_P),
                            prior(normal(0, 10), class = b, coef = high_N:high_P),
                            prior(normal(0, 10), class = b, coef = co:high_N:high_P),
                            prior(normal(0, 100), class = Intercept)))

# save model
save(m.acp.p, file = "./output/average-concentration-pav-cl-priors.rda")

# inspect model
summary(m.acp.p) # more specific priors hardly changed output
plot(m.acp.p)


## RPV model ##

# data
d.at.r <- d.at %>%
  filter(target == "RPV" & inoc != "PAV")

# concentration model
m.ac.r <- brm(data = d.at.r, family = gaussian,
              log_conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# save model
save(m.ac.r, file = "./output/average-concentration-rpv.rda")

# inspect model
summary(m.ac.r) # no significant effect
marginal_effects(m.ac.r)
plot(m.ac.r)

# with more specific priors
# use the estimates for coinfection and N addition from CL's models for the mean, multiply the standard error from her models by 10
m.acp.r <- update(m.ac.r, 
                  prior = c(prior(normal(-0.15, 3.0), class = b, coef = co),
                            prior(normal(-0.29, 10.03), class = b, coef = high_N),
                            prior(normal(0.03, 9.94), class = b, coef = high_P),
                            prior(normal(-0.37, 10.07), class = b, coef = co:high_N),
                            prior(normal(0, 10), class = b, coef = co:high_P),
                            prior(normal(-0.05, 9.90), class = b, coef = high_N:high_P),
                            prior(normal(0, 10), class = b, coef = co:high_N:high_P),
                            prior(normal(0, 100), class = Intercept)))

# save model
save(m.acp.r, file = "./output/average-concentration-rpv-cl-priors.rda")

# inspect model
summary(m.acp.r) # similar to model without specific priors
plot(m.acp.r)


#### save dat afor plotting ####

# merge PAV and RPV
d.out = full_join(d.at.p, d.at.r)

# save file
write_csv(d.out, "./output/exp-1-analysis-data.csv")


#### models not used ####

# PAV titer model
m.at.p <- brm(data = d.at.p, family = gaussian,
              log_quant ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

save(m.at.p, file = "./output/average-titer-pav.rda")

# inspect model
summary(m.at.p) # coinfection increases titer
marginal_effects(m.at.p)
plot(m.at.p)

# PAV concentration model withouth 3-way interaction
m.ac.p.1 <- update(m.ac.p, formula. = ~ . - co:high_N:high_P)

# save model
save(m.ac.p.1, file = "./output/average-concentration-pav-simplified-1.rda")

# inspect model
summary(m.ac.p.1) # coinfection increases concentration, but the effect goes away with high N

# concentration model without 1 2-way interaction
m.ac.p.2 <- update(m.ac.p.1, formula. = ~ . - co:high_N)
m.ac.p.3 <- update(m.ac.p.1, formula. = ~ . - co:high_P)
m.ac.p.4 <- update(m.ac.p.1, formula. = ~ . - high_N:high_P)

# concentration model with only 1 2-way interaction
m.ac.p.5 <- update(m.ac.p.2, formula. = ~ . - co:high_P) # N:P
m.ac.p.6 <- update(m.ac.p.3, formula. = ~ . - high_N:high_P) # co:N
m.ac.p.7 <- update(m.ac.p.4, formula. = ~ . - co:high_N) # co:P

# concentration model with only main effects
m.ac.p.8 <- update(m.ac.p.5, formula. = ~ . - high_N:high_P)

# save models
save(m.ac.p.2, file = "./output/average-concentration-pav-simplified-2.rda")
save(m.ac.p.3, file = "./output/average-concentration-pav-simplified-3.rda")
save(m.ac.p.4, file = "./output/average-concentration-pav-simplified-4.rda")
save(m.ac.p.5, file = "./output/average-concentration-pav-simplified-5.rda")
save(m.ac.p.6, file = "./output/average-concentration-pav-simplified-6.rda")
save(m.ac.p.7, file = "./output/average-concentration-pav-simplified-7.rda")
save(m.ac.p.8, file = "./output/average-concentration-pav-simplified-8.rda")

# loo fits
l.ac.p = loo(m.ac.p, save_psis = T)
l.ac.p.1 = loo(m.ac.p.1, save_psis = T)
l.ac.p.2 = loo(m.ac.p.2, save_psis = T)
l.ac.p.3 = loo(m.ac.p.3, save_psis = T)
l.ac.p.4 = loo(m.ac.p.4, save_psis = T)
l.ac.p.5 = loo(m.ac.p.5, save_psis = T)
l.ac.p.6 = loo(m.ac.p.6, save_psis = T)
l.ac.p.7 = loo(m.ac.p.7, save_psis = T)
l.ac.p.8 = loo(m.ac.p.8, save_psis = T)

# compare models
compare(l.ac.p, l.ac.p.1, l.ac.p.2, l.ac.p.3, l.ac.p.4, l.ac.p.5, l.ac.p.6, l.ac.p.7, l.ac.p.8) 
# elpd_diff + se_diff demonstrates whether the elpd_diff overlaps with 0
# many of them do, except 1, the full model, 5, and 2
# interpretation: don't need 3-way interaction, better with co:N, having co:P, N:P, or neither are all about the same (co:P slightly preferred)

# inspect best fitting models
summary(m.ac.p.4)
l.ac.p.4 # no high k values
summary(m.ac.p.6)
l.ac.p.6 # no high k values
summary(m.ac.p.1)
l.ac.p.1 # no high k values

# RPV titer model
m.at.r <- brm(data = d.at.r, family = gaussian,
              log_quant ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# save model
save(m.at.r, file = "./output/average-titer-rpv.rda")

# inspect model
summary(m.at.r) # no significant effect
marginal_effects(m.at.r)
plot(m.at.r)

## Goal: figure of Experiment 1 virus concentration and model estimates from exp-1-analysis.R


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(brms)
library(tidyverse)
library(ggridges)
library(cowplot)

# import data
dat <- read_csv("./output/exp-1-analysis-data.csv")

# load models
load("./output/average-concentration-pav-cl-priors.rda")
load("./output/average-concentration-rpv-cl-priors.rda")


#### edit data ####

# inoculation column
dat <- dat %>%
  mutate(inoculation = ifelse(co == 0, "single", "co"),
         inoculation = fct_relevel(inoculation, "single"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"))

# posterior samples
postr <- posterior_samples(m.acp.r)
postp <- posterior_samples(m.acp.p)

# rename columns
colnames(postr) <- colnames(postp) <- c("int", "co", "N", "P", "co_N", "co_P", "NP", "co_NP", "ar", "sigma", "lp")

# posterior slopes
sloper <- postr %>%
  transmute(low = int - mean(int),
            high_N = N,
            high_P = P,
            high_NP = N + P + NP,
            low_co =  co,
            N_co = co + co_N,
            P_co = co + co_P,
            NP_co = co + co_N + co_P + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

slopep <- postp %>%
  transmute(low = int - mean(int),
            high_N = N,
            high_P = P,
            high_NP = N + P + NP,
            low_co =  co,
            N_co = co + co_N,
            P_co = co + co_P,
            NP_co = co + co_N + co_P + co_NP) %>%
  gather(key = "treatment", value = "effect") %>%
  mutate(Inoculation = ifelse(grepl("co", treatment, fixed = T), "coinfection", "single"),
         Inoculation = factor(Inoculation, levels = c("single", "coinfection")),
         Nutrient = recode(treatment, high_N = "N", high_P = "P", high_NP = "N+P", low_co = "low", N_co = "N", P_co = "P", NP_co = "N+P"),
         Nutrient = factor(Nutrient, levels = c("low", "N", "P", "N+P")))

# check treatments
sloper %>% select(treatment, Inoculation, Nutrient) %>% unique()


#### figure of raw data ####

# color palette
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")

# text sizes
sm_txt = 6
lg_txt = 8
an_txt = 2

# PAV
plotA <- ggplot(filter(dat, target == "PAV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation)) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
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
                      name = "Nutrient") +
  scale_shape_manual(values = c(19, 21), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  xlab("Days post inoculation") +
  ylab(expression(paste(Log[10], "(PAV concentration)", sep = "")))

# RPV
plotB <- ggplot(filter(dat, target == "RPV"), aes(x = dpi, y = log_conc, colour = nutrient)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, position = position_dodge(0.6), aes(size = inoculation), show.legend = F) +
  stat_summary(fun.y = "mean", geom = "point", size = 1.5, position = position_dodge(0.6), aes(shape = inoculation)) +
  stat_summary(fun.y = "mean", geom = "line", position = position_dodge(0.6), aes(linetype = inoculation)) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.34, 0.93),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm")) +
  scale_size_manual(values = c(0.5, 0.5), name = "Inoculation") +
  scale_colour_manual(values = col_pal, guide = F) +
  scale_shape_manual(values = c(19, 21), name = "Inoculation") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Inoculation") +
  xlab("Days post inoculation") +
  ylab(expression(paste(Log[10], "(RPV concentration)", sep = "")))


#### figure of model estimates ####

# annotation text
ann_textp <- data.frame(effect = c(0.64, 0.75), Nutrient = c("N+P", "N+P"), lab = c("single infection", "coinfection"), Inoculation = c("single", "coinfection"))
ann_textr <- data.frame(effect = c(0.8, 0.96), Nutrient = c("N+P", "N+P"), lab = c("single infection", "coinfection"), Inoculation = c("single", "coinfection"))

# PAV
plotC <- ggplot(slopep, aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  facet_wrap(~Inoculation)  +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.8))) +
  xlab("Nutrient effect on PAV concentration") +
  ylab("Density of posterior dist.") +
  geom_text(data = ann_textp, label = c("single infection", "coinfection"), nudge_y = 1.6, size = an_txt)

# RPV
plotD <- ggplot(sloper, aes(x = effect, y = Nutrient, group = Nutrient, fill = Nutrient, linetype = Inoculation)) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.3) +
  facet_wrap(~Inoculation)  +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        strip.text = element_blank(),
        legend.title = element_text(color = "black", size = sm_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.position = c(0.25, 0.92),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_manual(values = col_pal, guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = F) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1.9))) +
  xlab("Nutrient effect on RPV concentration") +
  ylab("Density of posterior dist.") +
  geom_text(data = ann_textr, label = c("single infection", "coinfection"), nudge_y = 1.7, size = an_txt)


#### combine plots ####

pdf("./output/exp-1-concentration-figure.pdf", width = 6, height = 4)
plot_grid(plotA, plotC, plotB, plotD, labels = c("A", "C", "B", "D"), label_size = lg_txt)
dev.off()
