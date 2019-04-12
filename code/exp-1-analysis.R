#### set up ####

# import data
source("./code/exp-1-qPCR-raw-data-processing.R")
# loads tidyverse
# sets working directory to data folder

# clear all except dataset
rm(list = setdiff(ls(), c("dat")))


#### edit data ####

dpi <- tibble(
  time = 1:8,
  dpi = c(5, 8, 12, 16, 19, 22, 26, 29)
)

dat <- dat %>%
  filter(remove == 0 ) %>%
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

# look at mean values
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0), fun.data = "mean_cl_boot") +
  geom_point(data = filter(dat, quant_zero == 1), color = "blue", alpha = 0.5) + 
  facet_grid(target ~ inoc, scales = "free") # PAV growth is delayed by coinfection, RPV is enhanced (but more variable)

# look at mean values by nutrient
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV"), fun.data = "mean_cl_boot") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # the peak in coinfection is driven by low nutrients, but it is driven by N when PAV is alone

dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "RPV"), fun.data = "mean_cl_boot") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "RPV"), color = "blue", alpha = 0.5) + 
  facet_grid(inoc ~ nutrient, scales = "free") # it looks like there are multiple peaks in the temporal dynamics and they occur at different times depending on the nutrient and inoculation treatment, highest peaks with P addition

# look at PAV in coinfection
dat %>%
  ggplot(aes(x = dpi, y = quant_adj)) + 
  stat_summary(data = filter(dat, quant_zero == 0 & target == "PAV" & inoc == "coinfection"), fun.data = "mean_cl_boot") +
  geom_point(data = filter(dat, quant_zero == 1 & target == "PAV" & inoc == "coinfection"), color = "blue", alpha = 0.5) + 
  facet_wrap(~nutrient, nrow = 2, scales = "free") # peaks in the middle with high nutrients (like in single infection), especially P, peaks later with lower


# day of peak
# height of peak
# overall average titer
# proportion successfully infected