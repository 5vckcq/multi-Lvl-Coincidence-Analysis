##################################
# Taking anti-pregnancy pills, as well as being pregnant carries an increased risk of thrombosis.
# There are the four variables:
# is female (F)
# takes anti-pregnancy pills (A)
# is pregnant (P)
# has high risk of thrombosis (T)
##################################

library(cna)

rnames <- c(1,2,3,4,5) # row names
cnames <- c("F","A","P","T") # column names = names of causal factors

data_set <- array(0, dim=c(5,4),dimnames=list(rnames,cnames)) # a table of 4 rows and 5 columns each referenced by the entries in 
# rnames and cnames

# first row: is female, takes anti-pregnancy pills and has high risk of getting thrombosis
true <- c("F","A","T")
data_set[1,true] <- 1

# second row: is female, pregnant and has high risk of getting thrombosis
true <- c("F","P","T")
data_set[2,true] <- 1

# third row: : is male, takes anti-pregnancy pills and has high risk of getting thrombosis
true <- c("A","T")
data_set[3,true] <- 1

# fourth row: is female, neither takes anti-pregnancy pills, nor pregnant, has no high risk of getting thrombosis
true <- c("F")
data_set[4,true] <- 1

# fifth row: is male, neither takes anti-pregnancy pills, nor pregnant, has no high risk of getting thrombosis


setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(5,5,10), # at most 5 conjuncts, 5 disjuncts and 10 factors per formula
  details = FALSE),
  #ordering = list(c("A","B"),c("D"),c("E"))), # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  nsolutions = "all") # return all solutions

sink(file = NULL) # stop exporting into file
