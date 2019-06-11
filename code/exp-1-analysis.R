## Goal: evaluate treatment effects on virus concentration in Experiment 1


#### set up ####

# import data
source("./code/exp-1-qPCR-raw-data-processing.R")
# clears environment
# loads tidyverse
# sets working directory to data folder

# load libraries
library(brms)
library(tidyverse)
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
wcdat <- dat %>%
  select(q_group, well, target, cycle) %>%
  spread(target, cycle) %>%
  group_by(well) %>%
  summarise(mPAV = mean(PAV, na.rm = T), seP = sd(PAV, na.rm = T)/sqrt(length(!is.na(PAV))), mRPV = mean(RPV, na.rm = T), seR = sd(RPV, na.rm = T)/sqrt(length(!is.na(RPV))))

ggplot(wcdat, aes(x = mPAV, y = mRPV)) +
  geom_point() +
  geom_errorbar(aes(ymin = mRPV - seR, ymax = mRPV + seR)) +
  geom_errorbarh(aes(xmin = mPAV - seP, xmax = mPAV + seP)) +
  geom_smooth(method = "lm")
# wells that produce higher PAV values don't necessarily produce high RPV values and vice versa

cor.test(wcdat$mRPV, wcdat$mPAV)
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
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_se") +
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1), color = "blue", alpha = 0.5) + 
  facet_grid(target ~ inoc, scales = "free") # PAV growth is delayed by coinfection, RPV is enhanced (but more variable)

# look at mean values by nutrient - PAV
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_se") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # the peak in coinfection is driven by low nutrients, but it is driven by N when PAV is alone

# look at PAV in coinfection
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_se") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV" & inoc == "coinfection"), color = "blue", alpha = 0.5) + 
  facet_wrap(~nutrient, nrow = 2, scales = "free") # peaks in the middle with high nutrients (like in single infection), especially P, peaks later with lower

# look at mean values by nutrient - log-transformed PAV
dat2 %>%
  filter(quant_zero == 0 & target == "PAV" & !(inoc %in% c("healthy", "RPV"))) %>%
  ggplot(aes(x = dpi, y = log10(quant_adj))) + 
  stat_summary(fun.data = "mean_se") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")

# look at mean values by nutrient - RPV
dat2 %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_se") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_se", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "RPV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # it looks like there are multiple peaks in the temporal dynamics and they occur at different times depending on the nutrient and inoculation treatment, highest peaks with P addition

# look at mean values by nutrient - log-transformed RPV
dat2 %>%
  filter(quant_zero == 0 & target == "RPV" & !(inoc %in% c("healthy", "PAV"))) %>%
  ggplot(aes(x = dpi, y = log10(quant_adj))) + 
  stat_summary(fun.data = "mean_se") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  facet_grid(inoc ~ nutrient, scales = "free")


#### format data for models ####

# edit data
d.at <- dat2 %>%
  filter(quant_zero == 0 &
           inoc %in% c("PAV", "coinfection", "RPV")) %>%
  mutate(co = ifelse(inoc == "coinfection", 1, 0),
         log_quant = log10(quant_adj),
         exp_round = round,
         conc = quant_adj/mass_ext.mg,
         log_conc = log10(conc),
         quant_rd = round(quant_adj))

# concentration values
d.at %>%
  ggplot(aes(x = quant_rd)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

d.at %>%
  ggplot(aes(x = log_conc)) +
  geom_histogram() + 
  facet_wrap(~target, scales = "free")

# replicates within same round and qPCR group
d.at %>%
  group_by(target, exp_round, dpi, inoc, high_N, high_P, q_group) %>%
  summarise(reps = length(unique(sample))) %>%
  filter(reps >1) # 19

# data by virus
d.at.p <- d.at %>%
  filter(target == "PAV" & inoc != "RPV")

d.at.r <- d.at %>%
  filter(target == "RPV" & inoc != "PAV")


#### log-transformed models, uninformative priors ####

# PAV model
m.lu.p <- brm(data = d.at.p, family = gaussian,
              log_conc ~ co * high_N * high_P,
              autocor = cor_ar(~time),
              prior <- c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 1), class = b)),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)

# inspect model
summary(m.lu.p) # coinfection increases concentration
plot(m.lu.p) # convergence among chains
plot(marginal_effects(m.lu.p), points = T)
pp_check(m.lu.p, nsamples = 100) # model distributions slightly higher

# save model
save(m.lu.p, file = "./output/exp-1-analysis-log-uninformative-pav.rda")

# RPV model 
m.lu.r <- update(m.lu.p, newdata = d.at.r)

# inspect model
summary(m.lu.r) # no strong effects
plot(m.lu.r) # convergence among chains
plot(marginal_effects(m.lu.r), points = T)
pp_check(m.lu.r, nsamples = 100) # pretty close

# save model
save(m.lu.r, file = "./output/exp-1-analysis-log-uninformative-rpv.rda")


#### count models, uninformative priors ####

# check mean and variance
mean(d.at.p$conc_rd) #460
var(d.at.p$conc_rd) # much larger variance
mean(d.at.r$conc_rd) # 10683
var(d.at.r$conc_rd) # much larger variance

# PAV model (you can't currently implement autoregressive models, time included as a random effect)
m.cu.p <- brm(data = d.at.p, family = negbinomial(),
              quant_rd ~ offset(mass_ext.mg) + co * high_N * high_P + (1|time),
              prior <- c(prior(normal(460, 100), class = Intercept),
                         prior(normal(0, 1), class = b)),
              iter = 6000, warmup = 1000, chains = 1, cores = 1,
              control = list(adapt_delta = 0.99))

# inspect model
summary(m.cu.p)
yrep.cu.p <- posterior_predict(m.cu.p, draws = 100)
ppc_dens_overlay(log10(d.at.p$quant_rd), log10(yrep.cu.p)) # not a very good fit

# RPV model
m.cu.r <- update(m.cu.p, newdata = d.at.r)

# inspect model
summary(m.cu.r) # positive effect of P
yrep.cu.r <- posterior_predict(m.cu.r, draws = 100)
ppc_dens_overlay(log10(d.at.r$quant_rd), log10(yrep.cu.r)) # not a very good fit


#### log-transformed models, specific priors ####

# use the estimates for coinfection and N addition from CL's models

# PAV model
summary(m.c.p) # mean = 0.23, se = 0.41
summary(m.n.p) # mean = -0.21, se = 0.21
m.li.p <- update(m.lu.p,
                 prior = c(prior(normal(0.23, 0.41), class = b, coef = co),
                            prior(normal(-0.21, 0.21), class = b, coef = high_N),
                            prior(normal(0, 1), class = b, coef = high_P),
                            prior(normal(0, 1), class = b, coef = co:high_N),
                            prior(normal(0, 1), class = b, coef = co:high_P),
                            prior(normal(0, 1), class = b, coef = high_N:high_P),
                            prior(normal(0, 1), class = b, coef = co:high_N:high_P),
                            prior(normal(0, 10), class = Intercept)))

# save model
save(m.li.p, file = "./output/exp-1-analysis-log-informative-pav.rda")

# inspect model
summary(m.li.p) # didn't change estimates too much
plot(m.li.p) # convergence among chains
plot(marginal_effects(m.li.p), points = T)
pp_check(m.li.p, nsamples = 100) # model distributions slightly higher, but closer than model with uninformative priors

# RPV model
summary(m.c.r)
summary(m.n.r)
m.li.r <- update(m.lu.r, 
                  prior = c(prior(normal(-0.15, 0.3), class = b, coef = co),
                            prior(normal(-0.09, 0.26), class = b, coef = high_N),
                            prior(normal(0.36, 0.32), class = b, coef = high_P),
                            prior(normal(-0.24, 0.39), class = b, coef = co:high_N),
                            prior(normal(0, 1), class = b, coef = co:high_P),
                            prior(normal(-0.13, 0.39), class = b, coef = high_N:high_P),
                            prior(normal(0, 1), class = b, coef = co:high_N:high_P),
                            prior(normal(0, 10), class = Intercept)))

# save model
save(m.li.r, file = "./output/exp-1-analysis-log-informative-rpv.rda")

# inspect model
summary(m.li.r) # similar estimates to uninformative priors, tighter confidence intervals
plot(m.li.r) # convergence among chains
plot(marginal_effects(m.li.r), points = T)
pp_check(m.li.r, nsamples = 100) # similar to uninformative priors


#### save data for plotting ####

# merge PAV and RPV
d.out = full_join(d.at.p, d.at.r)

# save file
write_csv(d.out, "./output/exp-1-analysis-data.csv")