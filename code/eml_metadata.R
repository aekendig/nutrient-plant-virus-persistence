## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)

# Put data and code in a folder together to be grabbed by make_eml
# Generate metadata files by editing current ones or generating them (see Github page for tutorial)
# Edit and run this script


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)


#### fill blank spaces in qPCR raw data and save as single file ####

# import all raw data files
# setwd("./data/qPCR_raw")
# qlist <- 
#   list.files(pattern = "qPCR_data.csv") %>%
#   map(~read.csv(., header = T, stringsAsFactors = F))
# 
# # get list of numbers
# qnums <- list.files(pattern = "qPCR_data.csv") %>% 
#   gsub("[^[:digit:]]", "", .)
# setwd("../..")
# 
# # make cycle numeric
# for(i in 1:length(qlist)){
#   qlist[[i]]$cycle <- as.numeric(qlist[[i]]$cycle)
# }
# 
# # convert to a tibble
# qdat <-
#   bind_rows(qlist, .id = "q_group") %>%
#   as_tibble
# 
# # make q_group match file names
# qgroups = tibble(qnums = qnums, 
#                  q_group = unique(qdat$q_group))
# qdat2 <- qdat %>%
#   full_join(qgroups) %>%
#   select(-q_group) %>%
#   rename(q_group = qnums)

# export raw data files
# for(i in 1:44){
#   qname = paste("./temporary/group_", qnums[i], "_qPCR_data.csv", sep = "")
#   write.csv(qlist[[i]], qname, row.names = F)
# }
# manually moved to data folder after checking output

# export completed data file
# write_csv(qdat2, "./data/source_plant_qPCR_data.csv")

#### import templates ####

# list of data files
dlist <- list.files(path = "./data",
                    pattern = ".csv")

# high level text files
# template_core_metadata(path = "./metadata",
#                        license = "CCBY")
# 
# # one for each data table
# template_table_attributes(path = "./metadata",
#                           data.path = "./data",
#                           data.table = dlist)
# 
# # categorical values
# template_categorical_variables(path = "./metadata",
#                                data.path = "./data")
# 
# # look at units
# view_unit_dictionary()


#### data file descriptors ####

dlist

# description list
ddlist <- NA

ddlist[1] <- "comments about controls and standards for each qPCR group"
ddlist[2] <- "includes experimental treatments, measurements taken during harvesting, and molecular analysis information for source plants"
ddlist[3] <- "chlorophyll measurements for source plants"
ddlist[4] <- "qPCR data for source plants"
ddlist[5] <- "experimental treatments and molecular analysis information for receiving plants in transmission trials"

# name list
dnlist <- c("Comments for qPCR", 
            "Source plant experimental and molecular data", 
            "Source plant chlorophyll data", 
            "Source plant qPCR data", 
            "Receiving plant experimental and molecular data")

# print table
# dtable <- data.frame(data = dlist, description = ddlist)
# kable(dtable)


#### code descriptors ####

# list of code files
clist <- list.files(path = "./code",
                    pattern = ".R")

# remove the eml codes
clist <- clist[-2]

# code descripions
cdlist <- c(
  "code to analyze density of viruses within source plants",
  "code to analyze establishment of viruses within source plants",
  "code to create figure of treatment effects on virus establishment and density",
  "code for mathematical model and figures",
  "code to process raw qPCR data files",
  "code to create supplementary figure of comparison between models with and without priors",
  "code to analyze transmission",
  "code to create figure of treatment effects on transmission"
)

# name list
cnlist <- c("Virus density analysis", 
            "Virus establishment analysis", 
            "Figures for virus density and establishment",
            "Mathematical model of disease spread",
            "Processing of raw qPCR data",
            "Supplementary figure comparing models",
            "Transmission analysis",
            "Figures for transmission")

# print table
# ctable <- data.frame(code = clist, desription = cdlist)
# kable(ctable)


#### make eml ####

make_eml(path = "./metadata",
         data.path = "../nutrients_plant_viruses_edi",
         dataset.title = "Soil nitrogen and phosphorus effects on plant virus density, transmission, and species interactions",
         data.table = dlist,
         data.table.name = dnlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         other.entity = clist,
         other.entity.name = cnlist,
         other.entity.description = cdlist,
         temporal.coverage = c("2014-02-10", "2014-08-01"),
         geographic.description = "St. Paul, MN, USA",
         geographic.coordinates = c(44.98, -93.18, 44.98, -93.18),
         maintenance.description = "completed", 
         user.id = "aekendig",
         user.domain = "EDI",
         package.id = "edi.411.2")


#### check warnings ####

eml <- EML::read_eml("./metadata/edi.411.2.xml")
EML::eml_validate(eml)


