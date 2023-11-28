# example of a complex mechanism with three levels
# many ambiguities due to several co-extensive factors
# using QCA
library(QCA)
rnames <- c(1,2,3,4,5,6,7,8) # row names
cnames <- c("E1","E2","E3","E4","E5","E6","E7","F1","F2","F3","F4","F5","F6","F7","G") # column names = names of causal factors

data_set <- array(0, dim=c(8,15),dimnames=list(rnames,cnames)) # a table of 16 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row everything but D is true
true <- c("E1","E2","E3","E4","F1")
data_set[1,true] <- 1

# second row A,E,F,G are true
true <- c("F2")
data_set[2,true] <- 1

# third row A,E,F,H are true
true <- c("E1","E2","E3","E4","F1","F2","F4")
data_set[3,true] <- 1

# fourth row A,B,C,F,G,H,I,J are true
true <- c("F3","F5","F6")
data_set[4,true] <- 1

# fifth row everything is true
true <- c("E1","E2","E3","E4","F1","F3","F5","F6")
data_set[5,true] <- 1

# sixth row A,E,F are true
true <- c("F2","F3","F5","F6")
data_set[6,true] <- 1

# seventh row A,F,G are true
true <- c("E1","E2","E3","E4","E5","E6","E7","F1","F2","F3","F4","F5","F6","F7","G")
data_set[7,true] <- 1

# eighth row all false

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

# this line has to be added manually to pass this information to mLCA
print("ordering = E1, E2, E3, E4, E5, E6, E7 < F1, F2, F3, F4, F5, F6, F7 < G") 
print(causalChain(data_set, ordering = "E1, E2, E3, E4, E5, E6, E7 < F1, F2, F3, F4, F5, F6, F7 < G"), strict = FALSE)

sink(file = NULL) # stop exporting into file
