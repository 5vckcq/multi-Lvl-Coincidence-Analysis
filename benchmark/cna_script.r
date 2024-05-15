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
    sink(file = outputname) # start to export R output into file 
    
    print(cna(data_set, # use the truth table data_set
    maxstep=c(5,5,6), # at most 5 conjuncts, 5 disjuncts and 6 factors per atomic solution formula
    only.minimal.msc = TRUE, 
    details = FALSE,
    acyclic.only = TRUE),
    nsolutions = "all")
    
    sink(file = NULL) # stop exporting into file
}



filename <- "random_formulae.txt"

con <- file(filename,open="r")
lines <- readLines(con)

numCores <- 12
registerDoParallel(numCores)  # use multicores, change to the number of your CPU-cores
foreach (i=1:length(lines)) %dopar% {
    determine_asf_and_csf(i)
}
