# ninth test - example from the introductory part of our paper

# syntax:
# negation: minuscle of factor (not(A)=a)
# conjunction: "*"
# disjunction: "+"
# implication: "->"
# equivalence: "<->"

library(cna)
library('data.table')

formula <- "(E1 <-> E2)*(E1 <-> E3)*(E2*E3 <-> E4)*(E5 <-> E6)*(E6 <-> E7)*(F1*F2 <-> F4)*(F3 <-> F5)*(F3 <-> F6)*(F4*F5*F6 <-> F7)*(E1 <-> F1)*(E4 <-> F1)*(E5 <-> F7)*(E7 <-> F7)*(F1*F2*F3 <-> G)*(F7 <-> G)"  # condition to generate the truth table

complete_table <- allCombs(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)) - 1 # 14 columns for the 14 factors

names(complete_table) <- c("E1","E2","E3","E4","E5","E6","E7","F1","F2","F3","F4","F5","F6","F7","G") # the column names

data_set <- selectCases(formula, complete_table)  # generate the truth table

print(data_set)  # print table into console

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print(cna(data_set, # use the truth table data_set
  rm.dup.factors=FALSE, # do not discard logically equivalent factors
  maxstep=c(3,2,6), # at most 3 conjuncts, 2 disjuncts and 6 factors per atomic solution formula
  what = "a",  # return only atomic solution formulae (since the structure is too complex to generate csf)
  details = FALSE,
  ordering = list(c("E1","E2","E3","E4","E5","E6","E7"),c("F1","F2","F3","F4","F5","F6","F7"),c("G")) # ordering(e_1,e_2,...) places e_2 causally downstream to e_1
  # this is how the constitution levels get separated
  ), nsolutions = 100) # returns 100 solutions, might be changed to "all" (in quotation marks)
  
print(formula) # since we do not get complex solution formulae from CNA, we have to append our known solution manually
sink(file = NULL) # stop exporting into file
