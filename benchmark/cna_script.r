# read randomly generated formulae from file and generate asf + csf with cna

library(readr)
#library(Matrix)
#library(lme4)
#library(robustlmm)
library(cna)
library('data.table')

library(foreach)
library(doParallel)

determine_asf_and_csf <- function(number) {
    formula <- lines[number] # condition to generate the truth table
    complete_table <- allCombs(c(2,2,2,2,2,2)) - 1 # 6 columns for the 6 factors
    names(complete_table) <- c("A","B","C","D","E","F") # the column names
    data_set <- selectCases(formula, complete_table)  # generate the truth table
    outputname <- paste("cna_", toString(formatC(number, width = 5, format = "d", flag = "0")), ".txt", sep="") # generate unique file name using the enumeration of the formulae
    options(width= 1000) # max line length before automatic break
    sink(file = outputname) # start to export R output into file 
        
    print(configTable(data_set))
    
    print(
      cna(data_set, # use the truth table data_set
        maxstep=c(5,5,11), # at most 5 conjuncts, 5 disjuncts and 6 factors (+ possibly 5 negated) per atomic solution formula
        only.minimal.msc = TRUE, 
        details = TRUE, 
        acyclic.only = TRUE) |> csf(n.init = 1000000),
      n = 10000) # "n = 10000" with "|> csf ...", otherwise: nsolutions = "all"
    
    #print(
    #  cna(data_set, # use the truth table data_set
    #    maxstep=c(5,5,11), # at most 5 conjuncts, 5 disjuncts and 6 factors (+ possibly 5 negated) per atomic solution formula
    #    only.minimal.msc = TRUE, 
    #    details = TRUE, 
    #    acyclic.only = TRUE),
    #  nsolutions = "all")
    
    sink(file = NULL) # stop exporting into file
}



filename <- "random_formulae.txt"

con <- file(filename,open="r")
lines <- readLines(con)

numCores <- 12
registerDoParallel(numCores)  # use multicores, change to the number of your CPU-cores

set <- c(32,102,110,276,307,343,370,471,708,716,850,1088,1251,1357,1534,1536,1617,1793,2016,2256,2374,2481,2566,2642,2807,2834,3033,3036,3056,3233,3282,3383,3620,3730,3809,4042,4606,4679,4760,4766,4769,5142,5231,5337,5339,5506,5555,5630,5712,6097,6258,6263,6316,6415,6459,6499,6811,6916,6944,7138,7174,7306,7322,7495,7511,7572,7902,7953,8007,8033,8211,8224,8367,8708,8733,8784,8894,8980,9021,9056,9216,9284,9352,9406,9682,9898,9913)
#foreach (i=1:length(lines)) %dopar% {
for (i in set) {
    determine_asf_and_csf(i)
}
