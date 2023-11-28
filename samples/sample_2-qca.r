# a simple causal-mechanistic structure on two levels featuring negation, disjunction and conjunction
# using QCA, manually adding level information to output
library(QCA)

rnames <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) # row names
cnames <- c("A","B","C","D","E","F","G","H","I","J") # column names = names of causal factors

data_set <- array(0, dim=c(16,10),dimnames=list(rnames,cnames)) # a table of 16 rows and 10 columns each referenced by the entries in 
# rnames and cnames

# first row everything but D is true
true <- c("A","B","C","E","F","G","H","I","J")
data_set[1,true] <- 1

# second row A,E,F,G are true
true <- c("A","E","F","G")
data_set[2,true] <- 1

# third row A,E,F,H are true
true <- c("A","E","F","H")
data_set[3,true] <- 1

# fourth row A,B,C,F,G,H,I,J are true
true <- c("A","B","C","F","G","H","I","J")
data_set[4,true] <- 1

# fifth row everything is true
true <- c("A","B","C","D","E","F","G","H","I","J")
data_set[5,true] <- 1

# sixth row A,E,F are true
true <- c("A","E","F")
data_set[6,true] <- 1

# seventh row A,F,G are true
true <- c("A","F","G")
data_set[7,true] <- 1

# eighth row A,F,H are true
true <- c("A","F","H")
data_set[8,true] <- 1

# ninth row A,D,E,F,G are true
true <- c("A","D","E","F","G")
data_set[9,true] <- 1

# tenth row A,D,E,F,H are true
true <- c("A","D","E","F","H")
data_set[10,true] <- 1

# eleventh row B,D,H,G,I are true
true <- c("B","D","H","G","I")
data_set[11,true] <- 1

# twelfth row A,F are true
true <- c("A","F")
data_set[12,true] <- 1

# thirdteenth row A,D,E,F are true
true <- c("A","D","E","F")
data_set[13,true] <- 1

# fourteenth row D,G are true
true <- c("D","G")
data_set[14,true] <- 1

# fifteenth row D,H are true
true <- c("D","H")
data_set[15,true] <- 1

# sixteenth row only D is true
true <- c("D")
data_set[16,true] <- 1

setwd("..") # export output file to parent folder

sink(file = "r_output.txt") # start to export R output into file 

print("ordering = D, E, F, G, H, I, J < A, B, C") # this line has to be added manually to pass this information to mLCA
print(causalChain(data_set), ordering = "D, E, F, G, H, I, J < A, B, C")

sink(file = NULL) # stop exporting into file
