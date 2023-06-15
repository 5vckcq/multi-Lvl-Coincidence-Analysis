# test -- standard example for preemption
library(cna)

rnames <- c(1,2,3,4) # row names
cnames <- c("A","B","C","D","E") # column names = names of causal factors

data_set <- array(0, dim=c(4,5),dimnames=list(rnames,cnames)) # a table of 16 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row all factors are false

# second row A,C, D,E are true
true <- c("A","C","D","E")
data_set[2,true] <- 1

# third row C, D,E are true
true <- c("C","D","E")
data_set[3,true] <- 1

# fourth row A, B, E are true
true <- c("A","B","E")
data_set[4,true] <- 1

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  ordering = list(c("A","B","C","D"),c("E"))
  ), # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
