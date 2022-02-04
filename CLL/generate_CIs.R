# R script for generating the estimates for patients 3, 6, 9, and 21. 
# make sure to be in the CLL directory of the repo.

library(dplyr)
library(pracma)
library(MASS)
library(Rmpfr)


source("confidence_intervals.R")
#list of patients to generate estimates for
pt_list <- list("9","21","6","3") 
generate_CIs_multiple_pts(pt_list)

