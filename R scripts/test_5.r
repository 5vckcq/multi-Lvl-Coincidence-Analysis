# fuenfter Test: Baumgartners Beispiele
library(cna)

setwd("..") # Ausgabe soll in uebergeordneten Ordner erzeugt werden

sink(file = "r_output.txt") # Ausgabe werden ab hier in Datei gespeichert

print(cna(d.women, # unserer Datensatz
  rm.dup.factors=FALSE, # verwerfe Spalten mit identischen Eintraegen nicht
  maxstep=c(5,5,10), # maximal 5 Konjunkte, 5 Disjunkte und 10 Faktoren
  what = "a", # zeige nur atomare Loesungsformeln
  details = FALSE),
  nsolutions = "all") # gib alle Loesungen an

sink(file = NULL) # Ausgabe in Datei endet hier
