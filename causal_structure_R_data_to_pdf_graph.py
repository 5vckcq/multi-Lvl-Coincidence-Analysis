#!/usr/bin/env python3

# Datei: causal_structure_R_data_to_pdf_graph.py

import re                                        # Regex fuer Suche in Strings

def find_causal_factors(st):
    # sucht in String st nach durch ", ", " < " oder "*" separierten Bezeichnern fuer Kausalfaktoren
    # gibt die Liste an Kausalfaktoren zurueck
    
    # loesche "Factors: " aus Textzeile (sofern auftretend)
    st.replace("Factors: ","")
    # loesche Zeilenendensymbol und ggf. Leerzeichen am Zeilenende
    st = re.sub("\r?\n","",st).rstrip()
    # gib Liste der durch ", " oder "< " separierten Eintraege zurueck
    return re.split(",\s*|\s*<\s*|\*", st)

def get_equiv_formula(st):
    # gibt die Teilformeln links und rechts von "<->" in Eingabestring st als Paar von Strings (a,b) zurueck
    a = re.split(" <-> ",st)[0].strip()          # strip() loescht fuehrende Leerzeichen
    b = re.split(" <-> ",st)[1].strip()
    # in cna-Ausgabe haengen rechts noch weitere Informationen, diese von b abtrennen:
    b = re.split("[ \t]",b)[0]
    return (a,b)

def get_components_from_formula(st, factor_list):
    # gibt eine Liste der im String st vorkommenden Elemente aus factor_list zurueck
    component_list = []                          # Deklaration der Ausgabeliste
    for element in factor_list:
        if st.find(element) > -1:
            component_list.append(element)
    
    return component_list

def formula_level(st, level_factor_list):
    # falls alle Faktoren in st vom gleichen Level sind, wird dieses Level zurueck gegeben
    # anderenfalls ist der Rueckgabewert -1
    
    inequal = False
    factors = find_causal_factors(st)
    level = 0
    
    # suche das Level von factors[0]
    for i in range(len(level_factor_list)):
        if factors[0] in level_factor_list[i]:
            level = i
            
    # pruefe ob alle anderen factors vom gleichen Level sind
    for k in range(1,len(factors)):
        if not(factors[k] in level_factor_list[level]):
            inequal = True
        
    if inequal:
        level = -1
        
    return level

def convert_formula(st):
    # wandelt Formelsyntax von String st in Latex-Code um und
    # gibt diesen aus
    
    
    return st

def read_R_file(file_name):
    file_lines = []                              # Deklaration der Liste fuer die Textzeilen
    with open (file_name, 'rt') as text_file:    # Oeffne file_name um Textdatei auszulesen
        for next_line in text_file:              # fuer jede Zeile in file_name
            file_lines.append(next_line)         # fuege sie zur Liste file_lines hinzu
            
    ############################################## 
    # Schritt 1: bestimme Menge der Kausalfaktoren
    ##############################################
    factor_list = []                             # Deklaration der Liste der Kausalfaktoren
    # offenbar steht in Zeile 3 der R-Ausgabe immer entweder "Causal ordering:" oder "Factors:"
    # das sollte spaeter robuster gestaltet werden
    if file_lines[2].find("Causal ordering:") > -1:
        # im Falle von "Causal ordering:" sind Faktoren in Zeile 4 gelistet
        factor_list = find_causal_factors(file_lines[3]) 
        multi_level = True                       # multi_level, wenn dies in R angegeben worden ist
        
    elif file_lines[2].find("Factors:") > -1:
        # im Falle von "Factors:" sind Faktoren in Zeile 3 gelistet
        factor_list = find_causal_factors(file_lines[2])
        multi_level = False
    else:
        print("Programmabbruch keine Kausalfaktoren gefunden")
    
    # weiter falls Kausalfaktoren gefunden worden sind (Liste factor_list ist nicht leer)
    if factor_list:
        ##############################################
        # Schritt 2: suche nach Formeln
        ##############################################
        
        equiv_list = []                          # Deklaration der Liste fuer Aequivalenzformeln
            
        for line in file_lines:
            if line.find("<->") > 0:
                # falls "<->" gefunden wurde, lies Teilformeln links und rechts aus
                # und fuege sie zu equiv_list hinzu
                equiv_list.append(get_equiv_formula(line))
                    
        if not(equiv_list):                      # falls equiv_list leer Programmabbruch
            print("Programmabbruch keine Formeln gefunden")
        else:                                    # sonst weiter mit
            # pruefe ob jeder Kausalfaktor in mindestens einer Formel vorkommt, ansonsten loesche ihn
            # aus factor_list
            for k in range(len(factor_list)-1,-1,-1):
                i = 0
                found = False
                while not(found) and i < len(equiv_list):
                    found = (equiv_list[i][0].find(factor_list[k]) > -1 or equiv_list[i][1] == factor_list[k])
                    i = i + 1
                    
                if not(found):  # das Loeschen in der laufenden for-Schleife sollte aufgrund des Rueckwaertslaufens
                    # keine Probleme machen
                    print("Faktor " + factor_list[k] + " gelöscht, da in keiner Formel auftretend")
                    factor_list.remove(factor_list[k])  
                    print(equiv_list[len(equiv_list) - 1])
                    
            
            # Separierung der Ebenen in den Formeln
            # Vorsicht die folgende Zeile haengt kritisch von der Struktur der Eingabedatei ab!!!
            # Kausalfaktoren muessen in 4. Zeile stehen
            level_count = file_lines[3].count("<")
            # level_count = Anzahl der Levels -1
            
            level_factor_list = []               # Deklaration der separater Faktorlisten
            level_equiv_list = []
            constitution_relation_list = []
            for i in range(level_count + 1):
                if multi_level:
                    st = re.split(" < ",file_lines[3])[i].strip()
                    level_factor_list.append(find_causal_factors(st))
                    level_equiv_list.append([])
                else :
                    level_factor_list.append(find_causal_factors(file_lines[2]))
                    level_equiv_list.append([])
            
            # Entflechtung der Formeln
            for formula in equiv_list : 
                # fuege alle Formeln, die allein aus Elementen der Ebene i bestehen, zu level_equiv_list[i] hinzu
                if formula_level(formula[1], level_factor_list) == formula_level(formula[0], level_factor_list) :
                    print(formula)
                    print(formula_level(formula[0], level_factor_list))
                    level_equiv_list[formula_level(formula[1], level_factor_list)].append(formula)
                    
                elif formula_level(formula[0], level_factor_list) > -1 :
                    # alle Formeln, bei denen rechts ein Element aus einer anderen Ebene mit Elementen links (alle aus gleicher Ebene)
                    # verbunden wird, werde zu constitution_relation_list hinzugefuegt
                    constitution_relation_list.append(formula)
                           
                # alle Formeln mit Ebenenvermischung auf der linken Seite werden verworfen
            
            # ab hier folgend alle Schritte nach Ebenen separiert
            for m in range(level_count + 1):

                ##############################################
                # Schritt 3: suche nach Eingangsfaktoren
                ##############################################
                
                incoming_factor_list = []            # Deklaration der Liste an Eingangsfaktoren
               
                for element in level_factor_list[m]:
                    # Eingangsfaktoren stehen nie rechts in Formeln
                    i = 0
                    found = False
                    while not(found) and i < len(level_equiv_list[m]):
                        found = (level_equiv_list[m][i][1] == element)
                        i = i + 1
                    
                    if not(found):
                        # fuege element zu incoming_factor_list
                        incoming_factor_list.append(element)
                    
            
              
                if not(incoming_factor_list):
                    print("Programmabbruch keine eingehenden Kausalfaktoren gefunden")
                else:
                
                    ##############################################
                    # Schritt 4: entwickle Kausalstruktur iterativ
                    ##############################################
                    
                    # dazu zunaechst Zerlegung von level_factor_list in Faktoren, die nur links in den verbleibenden
                    # Formeln stehen -> upstream_factor_list
                    # bzw. auch rechts von "<->" stehen -> downstream_factor_list
                    upstream_factor_list = []        # Deklaration der Listen
                    downstream_factor_list = []
                
                    # Anfangs besteht upstream_factor_list aus den Elementen von incoming_factor_list,
                    # downstream_factor_list aus allen anderen Elementen von level_factor_list[m]
                    for element in level_factor_list[m]:
                        downstream_factor_list.append(element)
                        
                    for element in incoming_factor_list:
                        upstream_factor_list.append(element)
                        downstream_factor_list.remove(element)
                    
                    # jetzt wird upstream_factor_list sukzessiv ergaenzt um Faktoren, die 
                    # nur von Faktoren abhaengen, die bereits in upstream_factor_list stehen
                    # dazu wird downstream_factor_list von hinten durchgegangen
                    for j in range(len(downstream_factor_list)-1,-1,-1):
                        stop = False
                        k = 0
                        candidate = False
                        # k verweist auf den aktuellen Kandidaten aus downstream_factor_list
                        while not(candidate) and k < len(downstream_factor_list):
                            i = 0
                            candidate = True
                            # es werden alle Formeln dieser Ebene durchsucht, wobei
                            # i auf aktuelle Formel verweist
                            while candidate and i < len(level_equiv_list[m]):
                                # steht der aktuelle Kandidat rechts in der Formel?
                                if downstream_factor_list[k] == level_equiv_list[m][i][1]:
                                    # dann wird geprueft ob alle Faktoren links Elemente von 
                                    # upstream_factor_list sind, dabei:
                                    # level_equiv_list[m][i][0] -> m = Level der Formel
                                    # i = Nummer des Eintrags, 0 = linke Seite der Formel
                                    for element in get_components_from_formula(level_equiv_list[m][i][0], factor_list):
                                        if not(element in upstream_factor_list):
                                            candidate = False
                                
                                i = i + 1
                                
                            k = k + 1
                        
                        # falls in allen Formeln des Levels, alle Faktoren auf der linken Seite
                        # upstream_factors sind, wird downstream_factor_list[k-1] in einen solchen
                        # umgewandelt (k-1, weil zwischenzeitlich k um eins erhoeht wurde)    
                        if candidate :
                            upstream_factor_list.append(downstream_factor_list[k-1])
                                                        
                            # Formeln fuer Output vorbereiten  
                            #print(convert_formula(alle Formeln in denen downstream_factor_list[k-1] rechts steht))
                            
                            # unter restlichen Formeln suchen, ob Redundanzen auftreten (transitiv ueber [k-1] hinweg geschlossen wird)
                            # dann diese Formeln loeschen
                            
                            
                            downstream_factor_list.remove(downstream_factor_list[k-1])
                        else:
                            stop = True
                    

                    print(upstream_factor_list)
                    print(downstream_factor_list)
               
    

if __name__ == '__main__':
    read_R_file("r_output.txt")
