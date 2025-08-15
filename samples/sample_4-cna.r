# example of a complex mechanism with three levels
# many ambiguities due to several co-extensive factors
# using CNA
library(cna)
library('data.table')

# syntax:
# negation: minuscle of factor (not(A)=a)
# conjunction: "*"
# disjunction: "+"
# implication: "->"
# equivalence: "<->"
formula <- "(P1 <-> P2)*(P1 <-> P3)*(P2 + P3 <-> P4)*(T1 <-> T2)*(T2 <-> T3)*(P1 <-> X1)*(P4 <-> X1)*(T1 <-> X4)*(T3 <-> X4)*(X1*X3 <-> X2)*(X1*X2 <-> X3 )*(X2 + X3 <-> X4)*(X4 <-> S)"  # condition to generate the truth table

complete_table <- allCombs(c(2,2,2,2,2,2,2,2,2,2,2,2,2)) - 1 # 12 columns for the 12 factors

names(complete_table) <- c("X1","X2","X3","X4","T1","T2","T3","X1","X2","X3","X4","S") # the column names

data_set <- selectCases(formula, complete_table)  # generate the truth table

print(data_set)  # print table into console

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(3,2,6), # at most 3 conjuncts, 2 disjuncts and 6 factors per atomic solution formula
  what = "a",  # return only atomic solution formulae (since the structure is too complex to generate csf)
  details = FALSE,
  # this is how the constitution levels get separated
  ), nsolutions = "all") # returns 100 solutions, might be changed to "all" (in quotation marks)
  
print(formula) # since we do not get complex solution formulae from CNA, append our known solution manually
sink(file = NULL) # stop exporting into file
