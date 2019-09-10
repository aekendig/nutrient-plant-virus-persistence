## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)

# list of data files
dlist <- list.files(path = "./data",
                    pattern = ".csv")


#### import templates ####

# high level text files
template_core_metadata(path = "./metadata",
                       license = "CCBY")

# one for each data table
template_table_attributes(path = "./metadata",
                          data.path = "./data",
                          data.table = dlist)

# categorical values
template_categorical_variables(path = "./metadata",
                               data.path = "./data")

# look at units
view_unit_dictionary()


#### qPCR data templates ####

# manually edited one template
# put quotes are NA's to maintain them and leave the other cells blank
# don't use # symbol

# import templates
qAttTemp <- read.table("./metadata/attributes_group_xx_qPCR_data.txt", sep = "\t", quote = "\"", header = T, fill = T)
qCatTemp <- read.table("./metadata/catvars_group_xx_qPCR_data.txt", sep = "\t", quote = "\"", header = T, fill = T)

# capture qPCR groups
qGroups <- list.files(path = "./data",
                      pattern = "qPCR_data.csv") %>% 
  gsub("[^[:digit:]]", "", .)

# re-write templates for each group
for(i in 1:length(qGroups)){
  attName = paste("./metadata/attributes_group_", qGroups[i], "_qPCR_data.txt", sep = "")
  catName = paste("./metadata/catvars_group_", qGroups[i], "_qPCR_data.txt", sep = "")
  
  write.table(qAttTemp, attName, row.names = F, quote = F, sep = "\t", na = "")
  write.table(qCatTemp, catName, row.names = F, quote = F, sep = "\t", na = "")
}


#### data file descriptors ####

# write same descriptor for each qPCR file
ddlist <- NA

# loop through qPCR files
for(i in 1:length(qGroups)){
  ddlist[i] <- paste("group", qGroups[i], "qPCR data for source plants")
}

# add individual datasets
ddlist[45] <- "comments about controls and standards for each qPCR group"
ddlist[46] <- "sample experiment molecular data - includes experimental treatments, measurements taken during harvesting, and molecular analysis information for source plants"
ddlist[47] <- "chlorophyll measurements for source plants"
ddlist[48] <- "experimental treatments and molecular analysis information for receiving plants in transmission trials"

# print table
dtable <- data.frame(data = dlist, description = ddlist)
kable(dtable)


#### code descriptors ####

# list of code files
clist <- list.files(path = "./code",
                    pattern = ".R")

# code descripions
cdlist <- c(
  "code to create figure of correlation between PAV and RPV in coinfection",
  "code to analyze density of viruses within source plants",
  "code to create figure of treatment effects on virus density",
  "code to to create figure of the relationship between virus density and transmission",
  "code to create metadata",
  "code to analyze infection status of source plants",
  "code to create figure of treatment effects on infection status",
  "code to derive priors for virus density analysis",
  "code to derive priors for transmission analysis",
  "code to process raw qPCR data files",
  "code to create supplementary figure of comparison between models with and without priors",
  "code to create supplementary figure of comparison between models with full and truncated datasets",
  "code to analyze transmission",
  "code to create figure of treatment effects on transmission"
)

# print table
ctable <- data.frame(code = clist, desription = cdlist)
kable(ctable)


#### make eml ####

make_eml(path = "./metadata",
         data.path = "./data",
         dataset.title = "Soil nitrogen and phosphorus effects on plant virus density, transmission, and species interactions",
         data.table = dlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         temporal.coverage = c("2014-02-10", "2014-08-01"),
         geographic.description = "St. Paul, MN, USA",
         geographic.coordinates = NA,
         maintenance.description = "completed", 
         user.domain = "EDI",
         package.id = "edi.1.1")