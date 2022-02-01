#!/usr/bin/env python3

from tree import Tree
import itertools
import re
import multiprocessing
import functools


# Umformungsregeln
def commutation_conj(formula) :
    # Kommutativitaet von Konjunkten
    
    q = 0
    start_formula = formula
    
    if (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
        # falls sich die erste Konjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        
        while (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            bracket_term = ""
            for i in range(left_bracket_pos + 1, right_bracket_pos) :
                bracket_term = bracket_term + formula[i]
        
            formula = bracket_term
    
    elif formula.count("(") > 0 :
        # andere Klammerausdruecke werden zur Bearbeitung durch Pseudo-Literale substituiert
        
        while formula.find("(") > -1 :
            # trage Klammerausdruck in ein dictionary ein -> dictionary Bezeichner als Ersatz in der Formel
            # 1. Schritt "( ... )" auslesen und in dictionary speichern
        
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            bracket_term = ""
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i] for i in range(left_bracket_pos, right_bracket_pos + 1)]
            bracket_term = "".join(slist)
            
        
            formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
            
            if q == 0 : # lege dictionary an
                dictionary = {str(q) : bracket_term}
            else :     
                dictionary[str(q)] = bracket_term

            q = q + 1
    
    if formula.find("*") > -1 :
        # Vorgehen: spalte Formel zunaechst nach Disjunkten auf
        # bilde fuer jede Disjunkte separat die moeglichen Permutationen der Konjunkten
        # kombiniere jede Kombination jeder Disjunkten miteinander
        # fuehre dies fuer jede Formel aus eq_list aus
        
        permutation_list = []
        
        disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
            
        partial_permutation_list = [] # fuer jede Formel wird eine partial_permutation_list erstellt
        for disj in disj_list :
            if disj.find("*") > -1 :
                # fuer jede Disjunkte wird eine Unterliste in partial_permutation_list angelegt mit allen
                # Permutationen der Konjunkten in der Disjunkten disj
                # list(itertools.permutations(...) -> Liste der Permuationen der Konjunkten der Disjunkten disj
                partial_permutation_list.append(list(itertools.permutations(re.split("\*", disj))))
                
            else :    
                # falls Disjunkte keinen Konjunktor enthaelt (andere Disjunkte, aber schon), dann
                # trage disj in die Unterliste ein
                partial_permutation_list.append(list(disj))
            
            
        auxiliary_list = [] # fuehrt zusammengehoerige Konjunkte wieder in eine Formel zusammen (z.B. ("A","~B") -> "A*~B")
        for i in range(len(partial_permutation_list)) :                
            auxiliary_list.append([])
            for term in partial_permutation_list[i] :
                
                st = ""
                
                for konj in term :
                    if st == "" :
                        st = st + konj
                    else :
                        st = st  + "*" + konj

                if i < len(partial_permutation_list) - 1 :
                    st = st + " + "
                auxiliary_list[i].append(st)


        # bilde nun die Cartesischen Produkte aus den Varianten der Disjunkten
        for item in itertools.product(*auxiliary_list):
            # .join() ist schneller als Schleife ueber Buchstaben
            #st = ""
            slist = [t for t in item]
            st = "".join(slist)
            
            
            if q > 0 :     # aus Ursprungsformel wurden Klammern substituiert,
                # also Platzhalter "bracket_term1", ... wieder gemaess dictionary ersetzen
                for i in range(q) :
                    bracket_term = dictionary[str(i)]
                    st = st.replace("bracket_term" + str(i), bracket_term, 1) 
            
            elif formula != start_formula :
                # falls nur eine Teilformel ausgewertet wurde, muss st auf start_formula uebertragen werden
                st = start_formula.replace(formula,st,1)
            
            permutation_list.append(st)

        auxiliary_list.clear()
        partial_permutation_list.clear()
            
        return permutation_list
    else : return ""

def commutation_disj(formula) :
    # Kommutativitaet von Disjunkten
    q = 0
    start_formula = formula
    
    if (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
        # falls sich die erste Disjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
            
            
    elif formula.count("(") > 0 :
        # andere Klammerausdruecke werden zur Bearbeitung durch Pseudo-Literale substituiert
        
        while formula.find("(") > -1 :
            # trage Klammerausdruck in ein dictionary ein -> dictionary Bezeichner als Ersatz in der Formel
            # 1. Schritt "( ... )" auslesen und in dictionary speichern 
        
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos, right_bracket_pos + 1)]
            bracket_term = "".join(slist)
            
            formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
            
            if q == 0 : # lege dictionary an
                dictionary = {str(q) : bracket_term}
            else :     
                dictionary[str(q)] = bracket_term

            q = q + 1
    
    
    
    if formula.find(" + ") > -1 :
        # Vorgehen: bestimme alle Disjunkten, bilde alle mgl. Permutationen dieser Disjunkten
        
        # Liste aller mgl. Permutationen von Disjunkten aus Formeln von eq_list
        permutation_list = list(itertools.permutations(re.split("\s*\+\s*", formula)))
        
        return_list = []
        
        # zusammensetzen der Permutationen
        for permut in permutation_list :
            st = ""
            for disj in permut :
                if st == "" :
                    st = st + disj
                else :
                    st = st  + " + " + disj
            
            
            if q > 0 :     # aus Ursprungsformel wurden Klammern substituiert,
                # also Platzhalter "bracket_term1", ... wieder gemaess dictionary ersetzen
                for i in range(q) :
                    bracket_term = dictionary[str(i)]
                    st = st.replace("bracket_term" + str(i), bracket_term, 1) 
            
            elif formula != start_formula :
                # falls nur eine Teilformel ausgewertet wurde, muss st auf start_formula uebertragen werden
                st = start_formula.replace(formula,st,1)    
                
            
            # eintragen in die Ausgabeliste
            return_list.append(st)
        
        permutation_list.clear()
        
        return return_list
    else : return ""

def distribution(formula) :
    # Distributivitaet
     
    q = 0
    while formula.find("(") > -1 :
        # innerhalb von Klammern koennen keine neuen Formeln durch Distribution gefunden werden
        # daher koennen Klammerausdruecke einstweilen durch Pseudo-Literale ersetzt werden
        # trage Klammerausdruck in ein dictionary ein -> dictionary Bezeichner als Ersatz in der Formel
        # 1. Schritt "( ... )" auslesen und in dictionary speichern 
        
        left_bracket_pos = formula.find("(")
        right_bracket_pos = formula.find(")")
        
        # .join() ist schneller als Schleife ueber Buchstaben
        slist = [formula[i]  for i in range(left_bracket_pos, right_bracket_pos + 1)]
        bracket_term = "".join(slist)
        
        
        formula = formula.replace(bracket_term, "bracket_term" + str(q), 1)
        
        if q == 0 :  # lege dictionary an
            dictionary = {str(q) : bracket_term}
        else :     
            dictionary[str(q)] = bracket_term

        q = q + 1
    
    # Formel muss mindestens einen Disjunktor, zwei Konjunktoren enthalten
    if (formula.find(" + ") > -1) and (formula.count("*") > 1) :
                
        # Vorgehen: prueft ob erste n Konjunkte der ersten m Disjunkten identisch sind, wobei n < Anzahl der Konjunkte jeder dieser
        # Disjunkten gelten muss
        # falls ja: Klammer die n Konjunkten aus (A*B*C + A*B*D*E + B ... -> A*B*(C + D*E) + B ...)
        # andere Distributivitaeten koennen ignoriert werden, da sie durch Kommutationen gefunden werden
        
        return_list = []
        
        disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
        counter = 0
        max_disj = len(disj_list)
        conj_list = []
        while counter < len(disj_list) and (max_disj == len(disj_list)) :
            if disj_list[counter].find("*") > -1 :
                conj_list.append(re.split("\*", disj_list[counter]))
            else :
                max_disj = counter
                conj_list.append([])
            counter = counter + 1
        
        # conj_list ist nun eine Liste der Form [[c11, c12, ...], [c21, c22, ... ],  ...] wobei cij die jte Konjunkte der iten 
        # Disjunkten ist
        
        
        # nun: ruecklaufend was ist das groesste n (=Anzahl an gemeinsamen, ersten Konjunkten) der groesstmoeglichen Anzahl m
        # an ersten Disjunkten
        conj_counter = len(conj_list[0]) - 1
        while conj_counter > 0 :
            disj_counter = max_disj
            while disj_counter > 1 :
                condition = True
                for i in range(disj_counter) :
                    if conj_counter < len(conj_list[i]) :  # besitzt die Disjunktion ueberhaupt genuegend viele Konjunkte?
                        for j in range(conj_counter) :
                            if conj_list[0][j] != conj_list[i][j] :
                                condition = False
                    else : condition = False            
                        
                if condition :
                    # stimmen die ersten n Konjunkten der ersten m Disjunkten ueberein -> fuege Distributionsaequivalente zur 
                    # Ausgabeliste hinzu
                    
                    # konjunktiver Vorfaktor
                    
                    # .join() ist schneller als Schleife ueber Buchstaben
                    slist = [conj_list[0][i] + "*(" if i == conj_counter -1 else conj_list[0][i] + "*" for i in range(conj_counter)]
                    st = "".join(slist)
                    
                              
                    # disjunktiver Klammerausdruck
                    for i in range(disj_counter) :
                        for j in range(conj_counter,len(conj_list[i])) :
                            st = st + conj_list[i][j]
                            if j < len(conj_list[i]) - 1:
                                st = st + "*"
                            
                        if i < disj_counter - 1 :
                            st = st + " + "    
                        else :
                            st = st + ")"
                            
                    # weiteren Disjunkte
                    for i in range(disj_counter,len(disj_list)) :
                        st = st + " + " + disj_list[i]
                    
                                
                    if q > 0 : # Ursprungsformel enthielt keine Klammern, daher 
                        for i in range(q) : # Platzhalter "bracket_term1", ... wieder gemaess dictionary ersetzen
                            bracket_term = dictionary[str(i)]
                            st = st.replace("bracket_term" + str(i), bracket_term, 1) 
                        
                    
                    return_list.append(st)
                
                # versuche es mit einer Disjunkten weniger              
                disj_counter = disj_counter - 1
            # versuche es mit einer Konjunkten weniger    
            conj_counter = conj_counter - 1    
        if return_list :
            return return_list
        else : return ""
    else : return ""
    
def substitution(formula, equivalences) :
    # Substitutionen
    
    return_list = []
    
    for equiv in equivalences :
        if formula.find(equiv[0]) > -1 :   # falls die linke Seite einer Aequivalenz in der Formel auftritt
            aux_formula = formula.replace(equiv[0], equiv[1], 1) # erstmaliges Vorkommen durch die rechte Seite ersetzen
            # weitere Vorkommen koennen durch weitere Aufrufe ersetzt werden
            
            if aux_formula.find("(" + equiv[1] + ")") > -1 : # falls nun die Form "(A)" gewonnen wurde, loesche die Klammern um "A"
                aux_formula = aux_formula.replace("(" + equiv[1] + ")", equiv[1])
            
            return_list.append(aux_formula)  
    
    if return_list : return return_list
    else : return ""


def idempotence_disj(formula) :
    # formt den ersten Term der Form A*B*C + A*B*C + ... -> A*B*C + ... um
    
    start_formula = formula
    
    # Umgang mit Klammern
    if (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
        # falls sich die erste Disjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("+")) and (formula.find("+") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
    
    if formula.find("+") > -1 :
        # Idempotenz ist nur moeglich, wenn es mindestens einen Disjunktor gibt
        disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
        
        if disj_list[0] == disj_list[1] :
            # erste und zweite Disjunkte sind identisch -> entferne erste Disjunkte und setze Formel wieder aus disj_list zusammen
            
            slist = [disj_list[i] if i == len(disj_list) - 1 else disj_list[i] + " + " for i in range(1, len(disj_list))]
            st = "".join(slist)
            
            # falls nur erster Klammerausdruck umgeformt wurde, setzen ihn wieder in Gesamtformel ein
            if formula != start_formula :
                st = start_formula.replace(formula,st,1)  
            
            return_list = []
            return_list.append(st)
            return return_list # Rueckgabe in Liste
            
        else : return "" # keine Idempotenz zwischen ersten beiden Disjunkten
        
    else :  return "" # keine Disjunktoren in formula, also gib "" zurueck
    

def idempotence_conj(formula) :
    # formt den ersten Term der Form A*A*B + ... -> A*B + ... um
    
    start_formula = formula
    
    # Umgang mit Klammern
    if (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
        # falls sich die erste Konjunktion in Klammern befindet, wird nur der Klammerterm umgeformt
        while (formula.find("(") < formula.find("*")) and (formula.find("*") < formula.find(")")) :
            # Bestimmung des Strings innerhalb der Klammern
            left_bracket_pos = formula.find("(")
            right_bracket_pos = formula.find(")")
            
            # .join() ist schneller als Schleife ueber Buchstaben
            slist = [formula[i]  for i in range(left_bracket_pos + 1, right_bracket_pos)]
            formula = "".join(slist)
        
    
    if formula.find("*") > -1 :
        # Idempotenz ist nur moeglich, wenn es mindestens einen Konjunktor gibt
        
        valid = True
        if formula.find("+") > -1 :
            # falls die Formel eine Disjunktion ist, spalte in Disjunkte auf
            disj_list = re.split("\s*\+\s*", formula) # Liste der Disjunkten der Formel
            
            if disj_list[0].find("*") > -1 : # erste Disjunkte muss einen Konjuktor enthalten
                conj_list_one = re.split("\*", disj_list[0]) # nur erste Disjunkte ist von Relevanz
            else : valid = False     # sonst kann die Formel nicht in der richtigen Weise idempotent sein
            
        else : # Fall 2: formula ist eine Konjunktion, keine Disjunktion
            conj_list_one = re.split("\*", formula) 
        
        
        if valid :
            # formula enthaelt an der relevanten Stele einen Konjunktor
            
            if conj_list_one[0] == conj_list_one[1] :
                # erste und zweite Konjunkte sind identisch -> entferne erste Konjunkte und setze Formel wieder aus conj_list_one
                # und disj_list zusammen
                
                # Zusammensetzung des Restes der ersten Disjunkten
                slist = [conj_list_one[i] if i == len(conj_list_one) - 1 else conj_list_one[i] + "*" for i in range(1, len(conj_list_one))]
                st = "".join(slist)
                
                if formula.find("+") > -1 :
                    # Hinzufuegen der anderen Disjunkten, falls Ausgangsformel eine Disjunktion ist
                    slist = [" + " + disj_list[i] for i in range(1, len(disj_list))]
                    st = st + "".join(slist)
                
                # falls nur erster Klammerausdruck umgeformt wurde, setzen ihn wieder in Gesamtformel ein
                if formula != start_formula :
                    st = start_formula.replace(formula,st,1)  
            
                return_list = []
                return_list.append(st)
                return return_list # Rueckgabe in Liste
            
            
            else : return "" # keine Idempotenz zwischen ersten beiden Konjunkten

        else : return "" # Konjunktoren an falschen Stellen
        
    else :  return "" # keine Konjunktoren in formula, also gib "" zurueck

# Hilfsfunktion fuer parallele Ablaeufe
def smap(f):
    return f()

def search_formula(formula,formula_list) :
    # ermittelt alle aequivalenten Ausdruecke die mittels sechs Umformungsoperationen aus formula entstehen koennen und 
    # prueft, ob daraus eine Formel aus formula_list generiert werden kann
    # gibt entsprechenden Wahrheitswert zurueck
    #
    # Umformungsoperationen sind:
    # 1) Substitution eines Ausdrucks mittels der Aequivalenzrelationen aus formula_list
    # 2) Distributivumformung: A*B*C + A*B*~C*D + A*B*G + B ... -> A*B*(C + ~C*D + G) + B ...
    # 3) Kommutativitaet der Disjunktion
    # 4) Kommutativitaet der Konjunktion
    # 5) Idempotenz der Disjunktion
    # 6) Idempotenz der Konjunktion
    
    # Baum wird angelegt
    root = Tree(formula[0]) 
    
    # Auspaltung von formula_list
    formula_list_eq = []
    formula_trans_list = []
    for f in formula_list :
        # a) Formeln mit gleichem Term auf der rechten Seite
        if f[1] == formula[1] :
            formula_list_eq.append(f)
        
        # b) Formeln mit Termen niedrigerer Ordnung auf der rechten Seite koennen fuer Umformungen genutzt werden 
        else :  # (in Hauptfkt. wurde bereits geprueft, dass dies alle Formeln sind, die die erste Bedingung nicht erfuellen)
            formula_trans_list.append(f)
    
    formula_list.remove(formula)
    
    # Abgleichsliste, sodass keine doppelten Elemente im Baum vorkommen
    full_set = {formula[0]}
    
    # erlaubte Regeln (vertices)
    rules = {substitution, distribution, commutation_disj, commutation_conj, idempotence_disj, idempotence_conj}

    
    print("Baum "+ str(root) + " läuft durch")
    
    # Abbruch der Schleifen, sobald ein Ergebnis gefunden
    stop = not(formula_list_eq)   # sollte formula_list_eq leer sein, brauchen keine Umformungen untersucht zu werden,
    # es kann nie eine aequivalente gefunden werden


    for tree in root.dfs() :
        
        for formula_t in formula_list :
            # vergleich der Aequivalenten von formula mit den anderen Formeln der Formelliste
            if formula_t[0] == str(tree) : # falls Gleichheit gefunden wurde, sind wir fertig
                 print(formula_t[0] + " gefunden.")
                 stop = True
                 break  # Break aus innerer for-Schleife
            
        if stop : break # Break aus aeusserer for-Schleife
       

                
        
        func1 = functools.partial(substitution, str(tree),formula_trans_list)
        func2 = functools.partial(distribution, str(tree))
        func3 = functools.partial(commutation_disj, str(tree))
        func4 = functools.partial(commutation_conj, str(tree))
        func5 = functools.partial(idempotence_disj, str(tree))
        func6 = functools.partial(idempotence_conj, str(tree))
        
        pool = multiprocessing.Pool(processes=6)
        res = pool.map(smap, [func1, func2, func3, func4, func5, func6]) # parallele Anwendung der vier Umformungsoperationen
        pool.close()
        pool.join()
                      
        
        
        
        for sub_list in res :
            if not(sub_list == "") : 
                for j in sub_list :  # Regeln geben nur dann "" zurueck, wenn sie nicht auf die untersuchte Formel anwendbar sind
                    if not(j in full_set) : # sonst pruefe, ob umgeformte Formel nicht bereits im Baum enthalten ist
                        full_set.add(j)     # falls nicht ergaenze sie, ebenfalls zur Abgleichsliste
                        tree.add_child(Tree(j))
        
                                
    return stop

def main() :
    # nur zur Demonstration der Hauptfunktion search_formula
    
    formula = "E*G*~D*G*H + E*G*H + E*A" # die wird manipuliert
    
    # ein paar Vergleichsformeln - u.a. fuer Substitutionen
    st = "H*I"
    st1 = "G"
    st2= "E*A"
    st3 = "D"
    
    # da jede Formel eine Aequivalenz ist, bestehen sie aus zwei Teilen (linke Seite, rechte Seite)
    formula_1 =(st,st1)
    formula_2 = (st2,st3)
    
    # diese Formel soll rauskommen
    st4 = "D + E*G*~D*H + E*G*H"
    st5 = "J"
    formula_3 = (st4,st5) 

    formula_list = []
    formula_list.append(formula_1)
    formula_list.append(formula_2)
    formula_list.append(formula_3)
    print(search_formula(formula, formula_list))
    
if __name__ == '__main__':
    main()
