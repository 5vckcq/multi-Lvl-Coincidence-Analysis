library(cna)

rnames <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) # row names
cnames <- c("em","sm","eml","sml","af","ih") # column names = names of causal factors

data_set <- array(0, dim=c(17,6),dimnames=list(rnames,cnames)) # a table of 4 rows and 5 columns each referenced by the entries in 
# rnames and cnames

# first row: is female, takes anti-pregnancy pills and has high risk of getting thrombosis
true <- c("em","sm","af","ih")
data_set[1,true] <- 1

# second row: is female, pregnant and has high risk of getting thrombosis
true <- c("em","sm")
data_set[2,true] <- 1

# third row: : is male, takes anti-pregnancy pills and has high risk of getting thrombosis
true <- c("sm","af")
data_set[3,true] <- 1

# fourth row: is female, neither takes anti-pregnancy pills, nor pregnant, has no high risk of getting thrombosis
true <- c("sm")
data_set[4,true] <- 1

# fifth row: is male, neither takes anti-pregnancy pills, nor pregnant, has no high risk of getting thrombosis
true <- c("sm","eml","af")
data_set[5,true] <- 1

true <- c("sm","eml")
data_set[6,true] <- 1

true <- c("eml","sml","af","ih")
data_set[7,true] <- 1

true <- c("eml","sml")
data_set[8,true] <- 1

true <- c("sml","af")
data_set[9,true] <- 1

true <- c("sml")
data_set[10,true] <- 1

true <- c("em","sml","af")
data_set[11,true] <- 1

true <- c("em","sml")
data_set[12,true] <- 1

true <- c("em","sm","eml")
data_set[13,true] <- 1

true <- c("em","sm","eml","af","ih")
data_set[14,true] <- 1

true <- c("em","eml","sml","af","ih")
data_set[15,true] <- 1

true <- c("af")
data_set[16,true] <- 1

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
