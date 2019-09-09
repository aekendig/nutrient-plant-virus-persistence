## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)

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

## left off on qPCR run attributes and categories