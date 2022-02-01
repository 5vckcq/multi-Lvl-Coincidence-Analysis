library(cna)

rnames <- c(1,2,3,4,5,6,7,8)
cnames <- c("F1","F2","F3","F4","F5","F6","F7","E","A","B","C")

data_set <- array(0, dim=c(8,11),dimnames=list(rnames,cnames))

# erste Zeile alles auf 1
data_set[1,] <- 1

# zweite Zeile F1,F3,F5,F6,B auf 1
wahr <- c("F1","F3","F5","F6","B")
data_set[2,wahr] <- 1

# dritte Zeile F1,F2,F4,A auf 1
wahr <- c("F1","F2","F4","A")
data_set[3,wahr] <- 1

# vierte Zeile F1 auf 1
data_set[4,"F1"] <- 1

# fuenfte Zeile F2,F3,F5,F6,B auf 1
wahr <- c("F2","F3","F5","F6","B")
data_set[5,wahr] <- 1

# sechste Zeile F3,F5,F6,B auf 1
wahr <- c("F3","F5","F6","B")
data_set[6,wahr] <- 1

# siebente Zeile F2 auf 1
data_set[7,"F2"] <- 1

# achte Zeile alles auf 0

setwd("..") # Ausgabe soll in uebergeordneten Ordner erzeugt werden

sink(file = "r_output.txt") # Ausgabe werden ab hier in Datei gespeichert

print(cna(data_set, # unserer Datensatz
  rm.dup.factors=FALSE, # verwerfe Spalten mit identischen Eintraegen nicht
  maxstep=c(5,5,11), # maximal 5 Konjunkte, 5 Disjunkte und 11 Faktoren
  what = "a", # zeige nur atomare Loesungsformeln
  details = FALSE,
  ordering = list(c("F1","F2","F3","F4","F5","F6","F7"),c("E","A","B","C"))), # ordering(e_1,e_2,...) setzt e_2 downstream bezueglich e_1, hier werden die Ebenen separiert
  nsolutions = "all") # gib alle Loesungen an

sink(file = NULL) # Ausgabe in Datei endet hier
