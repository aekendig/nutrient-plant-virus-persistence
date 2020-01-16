## Goal: ordinary differential equation model of virus prevalence based on model estimated transmission values

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1
library(cowplot) # version used: 0.9.4
library(brms)  # version used: 2.7.0
library(deSolve) # version used: 1.21

# import data
ptran <- read_csv("./output/pav_transmission_values_same_nuts.csv")
rtran <- read_csv("./output/rpv_transmission_values_same_nuts.csv")

#### edit data ####

# calculate parameters
ptran2 <- ptran %>%
  select(inoculation, nutrient, prev, mechanisms) %>%
  spread(key = inoculation, value = prev) %>%
  mutate(qP = co / single) %>%
  select(-co) %>%
  rename(BP = single)

rtran2 <- rtran %>%
  select(inoculation, nutrient, prev, mechanisms) %>%
  spread(key = inoculation, value = prev) %>%
  mutate(qR = co / single) %>%
  select(-co) %>%
  rename(BR = single)

# combine data
tran_co <- full_join(ptran2, rtran2) %>%
  mutate(interactions = 1)

# create data without interactions
tran_s <- tran_co %>%
  mutate(qP = 1,
         qR = 1) %>%
  mutate(interactions = 0)

# combine data
tran <- full_join(tran_co, tran_s)


#### model ####

simmod = function(parm_dat){
  
  parms = list(BP = parm_dat$BP, qP = parm_dat$qP, BR = parm_dat$BR, qR = parm_dat$qR, Sinit = Sinit, Pinit = Pinit, Rinit = Rinit, Cinit = Cinit, N = Ninit, simtime = simtime)
  
  mymodel = with(as.list(parms), function(t, x, parms){
    
    S = x["S"]
    P = x["P"]
    R = x["R"]
    C = x["C"]
    
    Sdot = -(BP*P + BR*R + qP*BP*(1 - qR*BR)*C + qR*BR*(1 - qP*BP)*C + qP*BP*qR*BR*C)*S/N 
    Pdot = BP*(P + qP*(1 - qR*BR)*C)*S/N - BR*(R + qR*C)*P/N
    Rdot = BR*(R + qR*(1 - qP*BP)*C)*S/N - BP*(P + qP*C)*R/N
    Cdot = BP*(P + qP*C)*R/N + BR*(R + qR*C)*P/N + qP*BP*qR*BR*C*S/N
    
    list(c(Sdot, Pdot, Rdot, Cdot))
  })
  
  xstart = c(S = Sinit, P = Pinit, R = Rinit, C = Cinit)
  
  times = seq(0, simtime, length = simtime)
  
  out = as.data.frame(lsoda(xstart, times, mymodel, parms, hmax=20))
  
  return(out)
}


#### constants ####

Sinit = 4000
Pinit = 1
Rinit = 1
Cinit = 0
Ninit = Sinit + Pinit + Rinit + Cinit
days = 40
simtime = days / 4
thresh_C = 0.5 * Ninit


#### test model ####

# test <- simmod(tran[1, ])
# test %>%
#   mutate(N = S + P + R + C)


#### run model ####

# create columns
tran2 <- tran %>%
  mutate(S = NA,
         P = NA,
         R = NA,
         C = NA,
         time_C = NA)

# add final numbers to each parameter set
for(i in 1:nrow(tran2)){
  
  mod <- simmod(tran2[i, ])
  
  tran2[i, c('S','P', 'R', 'C')] <- mod[simtime, c('S','P', 'R', 'C')]
  # tran2[i, 'time_C'] <- filter(mod, C > thresh_C) %>%
  #   select(time) %>%
  #   min() * 4
  
}
tran2


#### format data ####

prev <- tran2 %>%
  mutate(PAV = P + C,
         RPV = R + C) %>%
  select(mechanisms, nutrient, interactions, PAV, RPV, C) %>%
  gather(key = "state", value = "count", -c(mechanisms, nutrient, interactions)) %>%
  mutate(prev = count / Ninit,
         state = recode(state, C = "Coinfected")) %>%
  mutate(mechanisms = fct_relevel(mechanisms, "density", "non-density"),
         nutrient = fct_relevel(nutrient, "low", "N", "P"),
         state = fct_relevel(state, "PAV", "RPV"))


#### figure settings ####

# palettes
col_pal = c("black", "darkgoldenrod2", "dodgerblue1", "palegreen4")
line_pal = c("solid", "dashed")
shape_pal = c(21, 22, 24)

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


#### figure ####

ggplot(prev, aes(x = nutrient, y = prev)) +
  geom_point(aes(fill = nutrient, shape = state), size = 3) +
  facet_grid(interactions ~ mechanisms) +
  base_theme +
  scale_fill_manual(values = col_pal, name = "Nutrient") +
  scale_shape_manual(values = shape_pal, name = "Infection") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

# coinfection is high with N addition because both viruses do well with N addition. This is driven by non-density-dependent mechanisms
# among-virus interactions have minimal effects on infection prevalence
# compare values to lit