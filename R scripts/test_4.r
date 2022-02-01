# vierter Test: mit Negationen
library(cna)

rnames <- c(1,2,3,4,5,6,7,8)
cnames <- c("A","B","C","I")

data_set <- array(0, dim=c(8,4),dimnames=list(rnames,cnames))

# notwendige Bedingung erfuellt
wahr <- c("A","I")
data_set[1,wahr] <- 1

# notwendige Bedingung erfuellt, nur ein Konjunkt des Vetos erfuellt
wahr <- c("A","B","I")
data_set[2,wahr] <- 1

# notwendige Bedingung erfuellt, nur ein Konjunkt des Vetos erfuellt
wahr <- c("A","C","I")
data_set[3,wahr] <- 1

# notwendige Bedingung erfuellt, aber Veto ebenso
wahr <- c("A","B","C")
data_set[4,wahr] <- 1

# nur Veto erfuellt
wahr <- c("B","C")
data_set[5,wahr] <- 1

# nur eine Bedingung des Vetos erfuellt
wahr <- c("B")
data_set[6,wahr] <- 1

# nur eine Bedingung des Vetos erfuellt
wahr <- c("C")
data_set[7,wahr] <- 1

# nichts erfüllt -> alles null

setwd("..") # Ausgabe soll in uebergeordneten Ordner erzeugt werden

sink(file = "r_output.txt") # Ausgabe werden ab hier in Datei gespeichert

print(cna(data_set, # unserer Datensatz
  rm.dup.factors=FALSE, # verwerfe Spalten mit identischen Eintraegen nicht
  maxstep=c(5,5,11), # maximal 5 Konjunkte, 5 Disjunkte und 3 Faktoren
  details = FALSE),
  nsolutions = "all") # gib alle Loesungen an

sink(file = NULL) # Ausgabe in Datei endet hier
