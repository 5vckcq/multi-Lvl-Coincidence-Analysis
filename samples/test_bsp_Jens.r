# second test -- a simple causal chain on two levels
library(cna)

rnames <- c(1,2,3,4,5,6,7,8,9) # row names
cnames <- c("A","B","C","D","E","F","G","H","J","K","L","M","N","O","P","Q","R","S","T","U") # column names = names of causal factors

data_set <- array(0, dim=c(9,20),dimnames=list(rnames,cnames)) # a table of 9 rows and 20 columns each referenced by the entries in 
# rnames and cnames

# first row 
true <- c("C","D","F","J","K","L","O","P","Q")
data_set[1,true] <- 1

# second row
true <- c("L","O","T","U")
data_set[2,true] <- 1

# third row
true <- c("C","E","F","H","J","L","N","O","P","Q","S")
data_set[3,true] <- 1

# fourth row
true <- c("B","C","E","F","G","H","J","K","L","O","P","Q","U")
data_set[4,true] <- 1

# fifth row
true <- c("A","C","F","G","J","L","M","N","O","P","Q","R","S","U")
data_set[5,true] <- 1

# sixth row
true <- c("B","C","D","E","F","G","J","L","P")
data_set[6,true] <- 1

# seventh row
true <- c("C","E","F","G","J","L","M","N","P","Q","U")
data_set[7,true] <- 1

# eighth row
true <- c("A","B","C","F","J","K","L","M","N","O","P","Q","R","T")
data_set[8,true] <- 1

# ninth row
true <- c("A","B","C","F","H","J","L","N","O","P","R","S","T","U")
data_set[9,true] <- 1

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,19), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE,
  #ordering = list(c("D","E","F","G","H","I","J"),c("A", "B","C"))
  ), # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
