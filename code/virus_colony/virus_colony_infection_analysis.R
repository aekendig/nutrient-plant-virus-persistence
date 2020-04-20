## Goal: evaluate infection status of inoculum leaves
# % with expected status
# % no infection
# % unexpected infection

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse) # version used: 1.2.1

# import data
dat <- read_csv("./data/virus_colony_infection_data.csv")


#### edit data ####

# examine data
head(dat)
unique(dat$infection_type)
unique(dat$infection_status)
unique(dat$notes)

# look closer at notes
dat %>%
  filter(!is.na(notes)) %>%
  select(infection_status, notes) %>%
  data.frame()

# light band notes
light_band_PAV = c("lighter band, but probably still infects", "light band", "light PAV band", "really light RPV and PAV bands", "really slight PAV band")
light_band_RPV = c("really light RPV and PAV bands", "slight RPV band, contamination?", "slight RPV band")

# add band information
dat2 <- dat %>%
  mutate(PAV_infection = case_when(infection_status %in% c("PAV", "coinfected") & !(notes %in% light_band_PAV) ~ 1,
                                   TRUE ~ 0),
         potential_PAV_infection = case_when(infection_status %in% c("light PAV") | notes %in% light_band_PAV ~ 1,
                                             TRUE ~ 0),
         RPV_infection = case_when(infection_status %in% c("RPV", "coinfected") & !(notes %in% light_band_RPV) ~ 1,
                                   TRUE ~ 0),
         potential_RPV_infection = case_when(infection_status %in% c("light RPV") | notes %in% light_band_RPV ~ 1,
                                             TRUE ~ 0),
         unknown = case_when(is.na(infection_status) ~ 1,
                             TRUE ~ 0),
         contamination = case_when(infection_type == "Healthy" & (PAV_infection == 1 | RPV_infection == 1) ~ 1,
                                             infection_type == "PAV" & RPV_infection == 1 ~ 1,
                                             infection_type == "RPV" & PAV_infection == 1 ~ 1,
                                             TRUE ~ 0),
         potential_contamination = case_when(infection_type == "Healthy" & (potential_PAV_infection == 1 | potential_RPV_infection == 1) ~ 1,
                                             infection_type == "PAV" & potential_RPV_infection == 1 ~ 1,
                                             infection_type == "RPV" & potential_PAV_infection == 1 ~ 1,
                                             TRUE ~ 0))


#### summary tables ####

# across all rounds
dat2 %>%
  group_by(infection_type) %>%
  summarise(n = n(),
            PAV_infection = sum(PAV_infection)/n,
            potential_PAV_infection = sum(potential_PAV_infection)/n,
            RPV_infection = sum(RPV_infection)/n,
            potential_RPV_infection = sum(potential_RPV_infection)/n,
            uknown = sum(unknown)/n)

# proportion potential
0.08/(0.853+0.08)
0.106/(0.773+0.106)

# for each round
dat2 %>%
  group_by(infection_type, round) %>%
  summarise(n = n(),
            PAV_infection = sum(PAV_infection)/n,
            potential_PAV_infection = sum(potential_PAV_infection)/n,
            RPV_infection = sum(RPV_infection)/n,
            potential_RPV_infection = sum(potential_RPV_infection)/n,
            uknown = sum(unknown)/n)

# summary for text
dat2 %>%
  group_by(infection_type) %>%
  summarise(n = n(),
            PAV_infection = (sum(PAV_infection) + sum(potential_PAV_infection))/n,
            RPV_infection = (sum(RPV_infection) + sum(potential_RPV_infection))/n)

# percentage with contamination
sum(dat2$contamination)/nrow(dat2)
sum(dat2$potential_contamination)/nrow(dat2)
