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


#### examine distribution ####

# histogram of values
dat %>%
  ggplot() +
  geom_histogram(aes(x = quant_adj, fill = quant_zero)) 

# one value is very high
filter(dat, quant_adj > 1e10) %>% data.frame() # a PAV value much higher than max

# look at all values above max
filter(dat, (target == "RPV" & quant_adj > RPVmax) | (target == "PAV" & quant_adj > PAVmax)) # only 3, remove these too

# remove values above standard curve
dat <- dat %>%
  filter((target == "RPV" & quant_adj <= RPVmax) | (target == "PAV" & quant_adj <= PAVmax))

# re-do histogram
dat %>%
  ggplot() +
  geom_histogram(aes(x = quant_adj, fill = quant_zero))

# log-tranformed values
dat %>%
  filter(quant_zero == 0) %>%
  ggplot() +
  geom_histogram(aes(x = log10(quant_adj)))


#### visualize sources of variation ####

# look at qPCR wells
dat %>%
  ggplot(aes(x = well, y = cycle, colour = target)) +
  geom_point(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_boot")
# not a strong trend due to wells and the ones that stick out have fewer replicates

# look at PAV and RPV concentrations in same well
dat %>%
  select(q_group, well, target, cycle) %>%
  spread(target, cycle) %>%
  group_by(well) %>%
  summarise(mPAV = mean(PAV, na.rm = T), seP = sd(PAV, na.rm = T)/sqrt(length(!is.na(PAV))), mRPV = mean(RPV, na.rm = T), seR = sd(RPV, na.rm = T)/sqrt(length(!is.na(RPV)))) %>%
  ggplot(aes(x = mPAV, y = mRPV)) +
  geom_point() +
  geom_errorbar(aes(ymin = mRPV - seR, ymax = mRPV + seR)) +
  geom_errorbarh(aes(xmin = mPAV - seP, xmax = mPAV + seP)) +
  geom_smooth(method = "lm")
# wells that produce higher PAV values don't necessarily produce high RPV values and vice versa

dat %>%
  select(q_group, well, target, cycle) %>%
  spread(target, cycle) %>%
  group_by(well) %>%
  summarise(mPAV = mean(PAV, na.rm = T), seP = sd(PAV, na.rm = T)/sqrt(length(!is.na(PAV))), mRPV = mean(RPV, na.rm = T), seR = sd(RPV, na.rm = T)/sqrt(length(!is.na(RPV)))) %$% 
  cor.test(mRPV, mPAV)
# not a strong correlation - well does not seem to affect concentration

# correlation between well's distance from mean and number of reps
dat %>%
  mutate(mean_cycle = mean(cycle, na.rm = T)) %>%
  group_by(well, target, mean_cycle) %>%
  summarise(well_mean = mean(cycle, na.rm = T),
            well_reps = length(well)) %>%
  mutate(well_dist = well_mean - mean_cycle) %>%
  ggplot(aes(x = well_reps, y = well_dist)) +
  geom_point() 
# all are within about 5 cycles of the mean and those outside have very low replication

# look at standard deviation among technical replicates
dat %>%
  filter(inoc != "healthy") %>%
  group_by(q_group, target, dpi, inoc, nutrient, sample) %>%
  summarise(tech_sd = sd(cycle, na.rm = T),
            zero = ifelse(sum(quant_zero) > 0, 1, 0)) %>%
  ggplot(aes(x = dpi, y = tech_sd, colour = inoc, shape = as.factor(zero))) + 
  geom_point() +
  facet_grid(target ~ nutrient, scales = "free")
# all without zero's are less than 3 and there isn't a strong relationship to experimental treatment

# look at mean and standard deviation of technical replicates without zero's
dat %>%
  filter(quant_zero == 0 & inoc != "healthy") %>%
  group_by(q_group, target, dpi, inoc, nutrient, sample) %>%
  summarise(tech_sd = sd(cycle),
            tech_mean = mean(cycle)) %>%
  ggplot(aes(x = tech_mean, y = tech_sd, colour = target)) + 
  geom_point()
# higher sd samples generally have higher mean values

# look at experimental rounds
dat %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  geom_point(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5, position = position_dodge(0.3)) +
    facet_wrap(~target)
# rounds are pretty similar, replicates within rounds are not necessarily more similar than across rounds (except 6 and 7 for RPV)

# look at rounds by treatment
dat %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.3), alpha = 0.7) +
  geom_point(alpha = 0.3) +
  facet_grid(nutrient~target)

# look at round by day
dat %>%
  ggplot(aes(x = round, y = cycle, colour = as.factor(replicate))) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(0.3), alpha = 0.7) +
  geom_point(alpha = 0.3) +
  facet_grid(target~dpi)
# not all days are in all groups, may be causing group differences


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

# check for same sample in multiple qPCR groups
dups <-dat2 %>%
  group_by(target, sample) %>%
  mutate(dup = duplicated(sample)) %>%
  filter(dup == T) %>%
  select(sample, target)

dups %>%
  left_join(dat2) %>%
  select(sample, target, quant_adj) %>%
  data.frame() 
# all are zero's - combine samples if needed


#### visualize treatment effects ####

# look at mean values
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1), color = "blue", alpha = 0.5) + 
  facet_grid(target ~ inoc, scales = "free") # PAV growth is delayed by coinfection, RPV is enhanced (but more variable)

# look at mean values by nutrient - PAV
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # the peak in coinfection is driven by low nutrients, but it is driven by N when PAV is alone

# look at PAV in coinfection
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV" & inoc == "coinfection"), color = "blue", alpha = 0.5) + 
  facet_wrap(~nutrient, nrow = 2, scales = "free") # peaks in the middle with high nutrients (like in single infection), especially P, peaks later with lower

# look at mean values by nutrient - log-transformed PAV
dat2 %>%
  filter(quant_zero == 0 & target == "PAV" & !(inoc %in% c("healthy", "RPV"))) %>%
  ggplot(aes(x = dpi, y = log10(quant_adj))) + 
  stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")

# look at mean values by nutrient - RPV
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "RPV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # it looks like there are multiple peaks in the temporal dynamics and they occur at different times depending on the nutrient and inoculation treatment, highest peaks with P addition

# look at mean values by nutrient - log-transformed RPV
dat2 %>%
  filter(quant_zero == 0 & target == "RPV" & !(inoc %in% c("healthy", "PAV"))) %>%
  ggplot(aes(x = dpi, y = log10(quant_adj))) + 
  stat_summary(fun.data = "mean_cl_boot") +
  stat_summary(fun.data = "mean_cl_boot", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")


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

# concentration values
d.at %>%
  ggplot(aes(x = conc)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

d.at %>%
  ggplot(aes(x = log_conc)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

# replicates within same round and qPCR group
d.at %>%
  group_by(target, exp_round, dpi, inoc, high_N, high_P,  q_groups) %>%
  summarise(reps = length(unique(sample))) %>%
  filter(reps >1) # 19


## PAV model ##

# data
d.at.p <- d.at %>%
  filter(target == "PAV" & inoc != "RPV")

# concentration model
m.ac.p <- brm(data = d.at.p, family = gaussian,
              log_conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 100), class = Intercept),
                         prior(normal(0, 10), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# save model
save(m.ac.p, file = "./output/average-concentration-pav.rda")

# inspect model
summary(m.ac.p) # coinfection increases concentration
marginal_effects(m.ac.p)
plot(m.ac.p)

# with more specific priors
# use the estimates for coinfection and N addition from CL's models for the mean, multiply the standard error from her models by 10
m.acp.p <- update(m.ac.p, 
                  prior = c(prior(normal(0.23, 5.9), class = b, coef = co),
                            prior(normal(-0.22, 2.3), class = b, coef = high_N),
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
                            prior(normal(-0.09, 2.0), class = b, coef = high_N),
                            prior(normal(0.36, 3.3), class = b, coef = high_P),
                            prior(normal(-0.24, 3.8), class = b, coef = co:high_N),
                            prior(normal(0, 10), class = b, coef = co:high_P),
                            prior(normal(-0.13, 4.0), class = b, coef = high_N:high_P),
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