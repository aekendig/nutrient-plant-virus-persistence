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


#### initial conditions ####

Sinit = 1000
Pinit = 1
Rinit = 1
Cinit = 0
Ninit = Sinit + Pinit + Rinit + Cinit
simtime = 25


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
tran <- full_join(ptran2, rtran2)


#### test model ####

test <- simmod(tran[1, ])
test %>%
  mutate(N = S + P + R + C)


#### run model ####

# create columns
tran2 <- tran %>%
  mutate(S = NA,
         P = NA,
         R = NA,
         C = NA)

# add final numbers to each parameter set
for(i in 1:nrow(tran2)){
  
  mod <- simmod(tran2[i, ])
  
  tran2[i, c('S','P', 'R', 'C')] <- mod[simtime, c('S','P', 'R', 'C')]
  
}
# all are coinfected