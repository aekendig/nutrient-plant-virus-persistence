#### set up ####

# import data
source("./code/exp-1-qPCR-raw-data-processing.R")
# loads tidyverse
# sets working directory to data folder

# load libraries
library(rethinking)
library(tidybayes.rethinking)

# clear all except dataset
rm(list = setdiff(ls(), c("dat")))


#### edit data ####

dpi <- tibble(
  time = 1:8,
  dpi = c(5, 8, 12, 16, 19, 22, 26, 29)
)

dat <- dat %>%
  filter(remove == 0 & material == "shoot") %>%
  mutate(quant_adj = case_when(target == "PAV" & quant_adj < PAVmin ~ 0,
                               target == "RPV" & quant_adj < RPVmin ~ 0,
                               is.na(quant_adj) ~ 0,
                               TRUE ~ quant_adj),
         quant_zero = case_when(quant_adj == 0 ~ 1,
                                TRUE ~ 0)) %>%
  full_join(dpi)


#### examine values ####

# histogram of values
dat %>%
  ggplot() +
  geom_histogram(aes(x = quant_adj, fill = quant_zero)) 

# one value is very high
filter(dat, quant_adj > 1e10) %>% data.frame() # a PAV value much higher than max

# look at all values above max
filter(dat, (target == "RPV" & quant_adj > RPVmax) | (target == "PAV" & quant_adj > PAVmax)) # only 3, remove these too

# remove high values
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
  geom_histogram(aes(x = log(quant_adj)))

# look at mean values
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1), color = "blue", alpha = 0.5) + 
  facet_grid(target ~ inoc, scales = "free") # PAV growth is delayed by coinfection, RPV is enhanced (but more variable)

# look at mean values by nutrient
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # the peak in coinfection is driven by low nutrients, but it is driven by N when PAV is alone

dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "RPV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # it looks like there are multiple peaks in the temporal dynamics and they occur at different times depending on the nutrient and inoculation treatment, highest peaks with P addition

# look at PAV in coinfection
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_cl_boot") +
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_cl_boot", geom = "line") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV" & inoc == "coinfection"), color = "blue", alpha = 0.5) + 
  facet_wrap(~nutrient, nrow = 2, scales = "free") # peaks in the middle with high nutrients (like in single infection), especially P, peaks later with lower


#### overall average titer for successful infections ####

# edit data

d.at <- dat %>%
  filter(quant_zero == 0 &
           inoc %in% c("PAV", "coinfection", "RPV")) %>%
  mutate(co = ifelse(inoc == "coinfection", 1, 0),
         log_quant = log(quant_adj),
         exp_round = round)

## PAV model ##

# data
d.at.p <- d.at %>%
  filter(target == "PAV" &
           inoc != "RPV") %>%
  select(log_quant, q_group, well, exp_round, dpi, high_P, high_N, co) %>%
  data.frame() %>%
  mutate(q_group_n = as.factor(q_group) %>% as.numeric(),
         well_n = as.factor(well) %>% as.numeric(),
         exp_round_n = as.factor(exp_round) %>% as.numeric(),
         dpi_n = as.factor(dpi) %>% as.numeric())

# model
m.at.p <- map2stan(
  alist(
    log_quant ~ dnorm(mu, sigma),
    mu <- a + 
      a_q_group[q_group_n] + 
      a_well[well_n] + 
      a_exp_round[exp_round_n] + 
      a_dpi[dpi_n] + 
      bp * high_P +
      bn * high_N +
      bc * co + 
      bpn * high_P * high_N +
      bpc * high_P * co +
      bnc * high_N * co +
      bnpc * high_P * high_N * co,
    a ~ dnorm(0, 100),
    a_q_group[q_group_n] ~ dnorm(0, sigma_q_group),
    a_well[well_n] ~ dnorm(0, sigma_well),
    a_exp_round[exp_round_n] ~ dnorm(0, sigma_exp_round),
    a_dpi[dpi_n] ~ dnorm(0, sigma_dpi),
    c(bp, bn, bc, bpn, bpc, bnc, bnpc) ~ dnorm(0, 10),
    c(sigma, sigma_q_group, sigma_well, sigma_exp_round, sigma_dpi) ~ dcauchy(0, 1) 
  ),
  data = d.at.p, warmup = 1000, iter = 6000, chains = 3, cores = 2)

# examine model
plot(m.at.p)
precis(m.at.p)
 ## start here - make sure strcutre makes sense and that well should be kept in


#### day and height of highest peak ####

dat %>%
  filter(quant_zero == 0 & (inoc == "coinfection" | inoc == target)) %>%
  group_by(target, inoc, nutrient, dpi) %>%
  summarise(mean_quant = mean(quant_adj),
            se_quant = sd(quant_adj)/sqrt(length(quant_adj))) %>%
  filter(mean_quant == max(mean_quant)) %>%
  select(target, inoc, nutrient, mean_quant, se_quant, dpi)

# need to re-do this, taking into account variation among wells and then among replicates

# look at low coinfection for PAV on day 26
dat %>%
  filter(inoc == "coinfection" & target == "PAV" & dpi == 26) %>%
  select(round, replicate, nutrient, q_group, well, quant_adj) %>%
  arrange(nutrient, round, replicate, q_group) %>%
  data.frame()





#### cumulative titer ####

#### proportion successfully infected ####