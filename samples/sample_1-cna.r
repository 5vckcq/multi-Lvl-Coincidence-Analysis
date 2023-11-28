# a simple causal structure on one level
# using cna
library(cna)
rnames <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) # row names
cnames <- c("A","B","C","D","E","F","G") # column names = names of causal factors

data_set <- array(0, dim=c(16,7),dimnames=list(rnames,cnames)) # a table of 16 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row
true <- c("A","B","C","D","E","F","G")
data_set[1,true] <- 1

# second row
true <- c("A","C","D","F")
data_set[2,true] <- 1

# third row
true <- c("B","C","D","F")
data_set[3,true] <- 1

# fourth row
true <- c("A","B","D","E")
data_set[4,true] <- 1

# fifth row
true <- c("A","B","C","E")
data_set[5,true] <- 1

# sixth row
true <- c("A","B","E")
data_set[6,true] <- 1

# seventh row
true <- c("A","C")
data_set[7,true] <- 1

# eighth row
true <- c("A","D")
data_set[8,true] <- 1

# ninth row
true <- c("B","C")
data_set[9,true] <- 1

# tenth row
true <- c("B","D")
data_set[10,true] <- 1

# eleventh row
true <- c("C","D","F")
data_set[11,true] <- 1

# twelfth
true <- c("A")
data_set[12,true] <- 1

# thirdteenth
true <- c("B")
data_set[13,true] <- 1

# fourteenth
true <- c("C")
data_set[14,true] <- 1

# fifteenth
true <- c("D")
data_set[15,true] <- 1

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 


print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per atomic solution formula
  details = FALSE),
  nsolutions = "all") # returns all solutions
sink(file = NULL) # stop exporting into file
