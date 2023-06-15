#!/usr/bin/env python3
#
# file: obtain_equivalence_formulae.py

# DESCRIPTION

import re                 # regex for complex search patterns in strings
import pandas as pd       # for reading csv files that contain truth tables
import itertools          # itertools provides functions to obtain all permutations of a string and Cartesian products of lists
import multiprocessing    # multiprocessing and functools for multicore usage
from auxiliary_functions import get_components_from_formula, get_factor_level, get_factor_order, get_equiv_formula

def list_to_string(in_list):
    # converts a nested list of the form out_list[DISJUNCT][CONJUNCT] or 
    # a list simple list of the form out_list[CONJUNCT] into a string
    # which connects disjuncts by " + " and conjuncts by "*"
    
    if in_list: # list is not empty
        if type(in_list[0]) == list: # list is nested
            return ' + '.join(['*'.join(disj) for disj in in_list])
        else: # list of strings
            return "*".join(in_list)
    else:
        return ""
    

def string_to_list(st):
    # converts a string of disjunctive normal form into a nested list
    # out_list[DISJUNCT][CONJUNCT]
    out_list = []
    disj_list = re.split("\s*\+\s*", st)
    for disj in disj_list:
        out_list.extend([re.split("\*", disj)])

    return out_list

def negate_formula(in_formula):
    # negates each entry in in_formula

    neg_formula = {}
    # convert the input formula in a list of its disjuncts
    if in_formula:
        for i in in_formula:
            neg_formula[i] = not(in_formula[i])  
    
    return neg_formula
    
def find_effects(formula, factor_list):
    # returns the list of effects in a formula generated from a configuration table
    #
    # causal factors are effects only if they do not satisfy either of three conditions
    # 1) the causal factor is one in every line of the configuration table (in this case they are irrelevant)
    # 2) it is zero in every line of the configuration table (same as 1)
    # 3) two lines in the table differ only by the value of the factor (this means they might only be first causes)
    
    # start with the full list of causal factors and reduce it accordingly to 1)-3) until only effects remain
    effect_list = [x for x in factor_list]
    for i in range(len(effect_list)-1,-1,-1):
        cond = True
        # first test: appears effect_list[i] in every formula (and its negation nowhere)?
        for term in formula:
            if not(effect_list[i] in term):
                cond = False
                break

        if not(cond):
            cond = True
            # second test: appears the negation of effect_list[i] in every formula?
            st = "~" + effect_list[i]
            for term in formula:
                if not(st in term):
                    cond = False
                    break
            
            if not(cond):
                for term in formula:
                    for sec_term in formula:
                        if sec_term == term:
                            sec_cond = False
                        else:
                            sec_cond = True
                            for fac in term:
                                # check whether all conditions are not satisfied for fac from term:
                                # (i) fac is effect_list[i] and sec_term its negation
                                # (ii) sec_term contains effect_list[i] and fac is its negation
                                # (iii) fac is contained in sec_term
                                if not((fac == effect_list[i] and st in sec_term) or (fac == st and effect_list[i] in sec_term) or fac in sec_term):
                                    sec_cond = False
                                    break
                            
                        if sec_cond:
                            break        
                    if sec_cond:
                        cond = sec_cond
                        break
        if cond:
            # delete the causal factor if either of the three exclusion criteria is true
            print(effect_list[i] + " discarded. It has no causal relevance for any other causal factor.")
            del effect_list[i]
            
    return effect_list
    
def get_instance_formula_to_factor(in_formula, factor, level_factor_list_order):
    # the instance function of a factor is obtained by removing that factor or its negation from each of the conjunct
    # and setting its value to True or False respectively
    output = {}
    
    full_formula = [x for x in in_formula]
    
    factor_level = get_factor_level(factor, level_factor_list_order)
    factor_order = get_factor_order(factor, level_factor_list_order)

    
    for lvl in range(len(level_factor_list_order)):
        for order in range(len(level_factor_list_order[lvl])):
            # removing all factors of different level from the instance function
            # reduces the further calculations significantly
            # however, an additional test will be necessary to check if the truncated instance formula
            # are correct (possible error: full instance function would be A + B <-> C, truncated function with A,C < B
            # becomes A <->, C which is wrong)
            if not(factor_level == lvl or factor_level == lvl + 1):
                for fac in level_factor_list_order[lvl][order]:
                    neg_fac = "~" + fac
                    for term in range(len(full_formula)-1,-1,-1):
                        if fac in full_formula[term]:
                            full_formula[term].remove(fac)
                        elif neg_fac in full_formula[term]:
                            full_formula[term].remove(neg_fac)
                            #full_formula[term] = full_formula[term].replace(factor,"")

            
    neg_fac = "~" + factor
    for term in full_formula:
        if factor in term:
            aux = [x for x in term]
            aux.remove(factor)
            output[list_to_string(aux)] = True
        elif neg_fac in term:
            aux = [x for x in term]
            aux.remove(neg_fac)
            output[list_to_string(aux)] = False
    
    # test the instance formula for wrong entries due to truncation
    delete_list = [] # list of wrong entries in output
    for disj_ins in output:
        if output[disj_ins]:
            neg_fac = "~" + factor
        else:
            neg_fac = factor

        for disj in full_formula:    
            found_all = True
            for term in string_to_list(disj_ins)[0]:
                if not(term in disj):
                    found_all = False
                    break
            if found_all:
                # check whether factor appears in disj as part of full_formula as it appears in the preliminary instance function
                # if not: add the entry to delete_list and remove it from the instance function
                if neg_fac in disj:
                    delete_list.append(disj_ins)
                    
    if delete_list:
        # delete_list is non empty
        # remove the list entries from the dictionary output
        delete_list = list(set(delete_list)) # get rid of potential dublicates
        for entry in delete_list:
            del output[entry]            
    
    return output

def reduce_term_by(term, literal):
    #
    if term == literal:
        # if term is atomic 
        # return empty string because no reducted term exists
        return ""
    else:
        # term is a complex        
        st = literal + "*"
        if term.find(st) > 0 or term.find(st) == -1:
            # literal is not first conjunct
            # then cut string of conjunctions to the left
            st = "*" + literal
        return term.replace(st,"")
        
def contains_term(original_term, comparison_term):
    # in: original_term - a string
    #     comparison_term - another string
    # out: truth value whether original_term is contained in comparison_term
    #      "contains" can imply that e.g. "A*C" is contained in "A*B*C"
    
    if original_term == "":
        return False
    else:
        aux_list = string_to_list(original_term)
        all_found = True
        for fac in aux_list[0]:
            if not(fac in string_to_list(comparison_term)[0]):
                all_found = False
                break
        
        return all_found

def absorb_terms(arg): 
    # in: tuple of arguments: 
    # arg[0] - is a nested list arg[0][LIST OF CONJUNCTS][CONJUNCTS]
    # arg[1] - is a nested list, which should include arg[0]; arg[1][LIST OF DISJUNCTS][LIST OF CONJUNCTS][CONJUNCTS]
    # out: the nested listed reduced by the absorption rule
    # absorption rule: in the list of disjunctions, discard disjuncts that can be formed by conjuncting terms to other disjuncts of the list
    # so: if all elemts of a list of conjuncts are completely contain in another, discard the larger list
    return [[disj for disj in arg[1] if all([x in conj_2 for x in conj] for conj in disj)] for conj_2 in arg[0]]


def distribution(formula):
    # applies the distribution rule on a logical formula given as a string as often as possible
    # then applies the absorption rule to simplify the resulting formula as often as possible
    # returns the simplified formula as a string
    
    # as long as formula has a conjunctor right before or after a bracket
    if (formula.find(")*") > -1 or formula.find("*(") > -1):
                
        formula = formula[:-1] # get rid of trailing ")"
        conj_list = re.split("\)\*", formula) # list of conjuncts of formula
        conj_list = [conj[1:] for conj in conj_list] # get rid of leading "("
        
        disj_list = [] # list of disjuncts per conjunct
        for conj in conj_list:
            disj_list.append(re.split("\s*\+\s*", conj))        
        # disj_list is a list [[d11, d12, ...], [d21, d22, ... ],  ...] with dij being the j-th disjunct in conjunct i
        
        # rebuild formula
        formula = list_to_string([[*x] for x in itertools.product(*disj_list)])
        
        disj_list.clear()
        conj_list.clear()
        
        
        disj_list = re.split("\s*\+\s*", formula) # list of disjuncts of new formula
        
        #conj_list = [list(set(re.split("\*", disj))).sort() for disj in disj_list] # doesn't work
        
        conj_list = []
        set_disjuncts = set() # set in place of a list automatically discards duplicates
        for disj in disj_list:
            if not(disj in set_disjuncts):
                a = list(set(re.split("\*", disj)))
                a.sort()
                conj_list.append(a)
                set_disjuncts.add(disj)
        
        new_list = []
        for disj in conj_list:
            if not(disj in new_list):
                new_list.append(disj)   
        
        # simplify by absorption (a*b*c*d*e + a*c*d <-> a*c*d)        
        
        # start working on all CPUs
        arguments = ([f, new_list] for f in new_list)
        with multiprocessing.Pool() as pool:
	        # call the function for each item in parallel
            absorbed_terms = pool.map(absorb_terms, arguments)
        
        pool.close()
        pool.join()

        new_list = [i for i in new_list if i not in absorbed_terms]        

        new_list.sort()
        # translate formula encoded in the nested list of disjuncts of conjuncts into a string
        formula = list_to_string(new_list)

    return formula        
    
def get_prime_implicants(instance_formula, factor, level_factor_list):
    # prime implicants are obtained by comparing the min-terms of positive instance function
    # with those of the negative instance function 
    # - if a section of one positive min-term is not part of any negative min term,
    #   that positive min-term can be reduced to this section
    # - if every section is a part of at least one negative min term, the considered term is a prime factor
    #   of the positive instance function
    
    
    # get number of True terms in instance_formula
    num_entries = 1
    for disj in instance_formula:
        if num_entries < len(get_components_from_formula(disj, level_factor_list)):
            num_entries = len(get_components_from_formula(disj, level_factor_list))
    
    
    reduced_term_list = []
    prime_imp_list = []
    # prepare sub-lists of reduced_term_list
    for k in range(num_entries,-1,-1):
        reduced_term_list.append([])

    
    for term in instance_formula:
        if instance_formula[term]:
            reduced_term_list[num_entries].append(term)

    for k in range(num_entries,-1,-1):
        for term in reduced_term_list[k]:
            # only check term if not already known as prime implicator
            if not(term in prime_imp_list):
                any_lit_found = True
                term_list = string_to_list(term)[0]
                # reduce term by one of its literals and check whether the reduced fragments is not contained in any negative term
                for lit in term_list: 
                    
                    # check only fragments that are neither already known to be not contained in any negative term
                    # (that are not listed in reduced_term_list[k-1]), nor empty strings
                    if not(reduce_term_by(term, lit) in reduced_term_list[k-1]) and reduce_term_by(term, lit) != "":
                        contained = False
                        for neg_term in instance_formula:
                            # loop over all negative terms of the instance formula
                            if not(instance_formula[neg_term]): # only continue with False neg_terms
                                if contains_term(reduce_term_by(term, lit), neg_term):
                                    contained = True
                                    break # break from for-loop over negative terms after the fragment has been found in one negative term

                        any_lit_found =  any_lit_found and contained
                        if not(contained): # if a fragment resulting from deleting a literal from term is not contained in any
                            # negative instance term, add the fragment to the list of terms that is gradually shortened and checked to
                            # for being a prime implicator
                            reduced_term_list[k-1].append(reduce_term_by(term, lit))

                    
                if any_lit_found: # if all fragments resulting from deleting a literal from term are contained in negative instance terms,
                    # then term is a prime implicator
                    prime_imp_list.append(term)
                    #print('\x1b[0;36;40m' + str(term) + " added to PI list.\x1b[0m") #dummy-print
    
    # add atomic prime implicants
    if reduced_term_list[1]:
        for at_term in reduced_term_list[1]:
            if not(at_term in prime_imp_list):
                prime_imp_list.append(at_term)
   
    for n_pi in range(len(prime_imp_list)-1,-1,-1):
        for pi in prime_imp_list:
            if contains_term(pi, prime_imp_list[n_pi]) and not(pi == prime_imp_list[n_pi]):
                del prime_imp_list[n_pi]
                break
    
    return prime_imp_list

def get_rdnf(pi_list, formula, factor_list):
    # obtains reduced disjunctive normal form by Petrick's algorithm
    # in: pi_list - list of strings = list of all prime implicants
    #     formula - dictionary of strings with all values either True or False
    #               = keys are the min terms of the formula, True or False their
    #                 truth values
    # out: solution_list - list of strings - contains all solutions each item is the string of a reduced disjunctive normal form of formula
    
    solutions_list = []
    out_formula = ""
    
    # define list of essential prime implicants
    e_pi_list = []
    # check for essential prime implicants
    uncovered_terms = {}
    for term in formula:
        if formula[term]:
            aux_list = []
            for lit in pi_list:
                if contains_term(lit,term):
                    aux_list.append(lit)
            
            if not(aux_list):
                # no pi found that makes this min-term true, should never happen
                print("Problem: No prime implicant has been found for " + str(term))
            else:
                if len(aux_list) == 1:
                    # term is only covered by one prime implicant -> this is an essential prime implicant
                    if not(aux_list[0] in e_pi_list):
                        e_pi_list.append(aux_list[0])
                else:
                    # several prime implicants cover term
                    uncovered_terms[term] = aux_list
    

    # simplest case would be that the disjunction of essential prime implicants covers al min terms
    # check if this the case 
    all_covered = True
    for term in uncovered_terms:
        covered = False
        for epi in e_pi_list:
            if epi in term:
                covered = True
                break
        all_covered = all_covered and covered
        if not(all_covered):
            break
            
                
    # the corresponding formula:
    for epi in e_pi_list:
        out_formula = out_formula + " + " + epi
    out_formula = out_formula[3:] # remove the leading " + "
    
    if all_covered:
        # first case: disjunction of all essential prime implicants covers all min terms
        # there is only one solution
        solutions_list.append(out_formula)
    else:
        # disjunction of essential prime implicants does not cover all min terms
        # proceed in four steps:
        # A) for each min term uncovered by the essential prime implicants: form the disjunction of the PI that cover it
        # B) form the product of all such disjunctions
        # C) transform it into disjunctive normal form
        # D) simplify it by applying idempotence and absorption law
        # -> each disjunct of the remaining formula is one solution + to be added to the disjunction of essential PIs
        
        # create a dictionary for PIs
        dict_pi = {}
        for i in range(len(pi_list)):
            abbr = "PI" + str(i)
            dict_pi[pi_list[i]] = abbr
            
        aux_formula = "("
        
        for u_term in uncovered_terms:
            # check whether one essential PI covers u_term
            e_covered = False
            for epi in e_pi_list:
                if epi in uncovered_terms[u_term]:
                    e_covered = True
                    break
            if not(e_covered):
                for PI in uncovered_terms[u_term]:
                    if aux_formula[-1] == "(":
                        # first term of a product
                        aux_formula = aux_formula + dict_pi[PI]
                    else:
                        # all further terms
                        aux_formula = aux_formula + " + " + dict_pi[PI]
                # end of disjunctions
                aux_formula = aux_formula + ")*("
        # remove trailing "*("
        aux_formula = aux_formula[:-2]

        aux_formula = distribution(aux_formula)

        # every disjunct of aux_formula constitutes one solution
        sol_list = re.split("\s*\+\s*",aux_formula)
        
        for sol in sol_list:
            # in each solution "*" are to be changed into " + " (part of Petrick's algorithm)
            st = sol.replace("*"," + ")
            for i in range(len(pi_list)-1,-1,-1):
                # replace the PI-placeholders by their values
                st = st.replace("PI"+str(i),pi_list[i])
            if out_formula != "":
                st = out_formula + " + " + st
            elif st[0] == "(":
                # if there is no contribution by ePI, remove the leading and trailing bracket
                st = st[1:]
                if len(sol_list) == 1:
                    # remove 
                    st = st[:-1]
            elif len(st) > 0 and st[-1] == ")":
                st = st[:-1]
            if len(st) > 1 and st[-2] == "+":
                # remove trailing " + "
                st = st[:-3]

            solutions_list.append(st)        
    
    # remove solutions that are disjunctions of other solutions (if A <-> B then A + C <-> B should not be in the list of solutions)
    for n_sol in range(len(solutions_list)-1,-1,-1):
        list1 = get_components_from_formula(solutions_list[n_sol], factor_list)
        for sol in solutions_list:
            list2 = get_components_from_formula(sol, factor_list)
            if len(list2) < len(list1):
                cond = True
                for fac in list2:
                    if not(fac in list1):
                        cond = False
                        break
                
                if cond:
                    del solutions_list[n_sol]
                    break
                            
    return solutions_list

def get_truth_table_from_file(file_path):
    # reads a csv file into a data frame
    # in: file_path - string - path to csv-file
    # out: factor_list - list of strings - list of column heads in csv-file = list of causal factors
    #      formula - string - logical formula derived from csv-truth table
    
    # different strings that will be interpreted as True
    true_values = [1, "1", "T", "t", "w", "W", "true", "True", True]
    # different strings that will be interpreted as False    
    false_values = [0, "0", "F", "f", "false", "False", False]
    
    # read csv file into df, assume that the column separators are one character from sep
    df = pd.read_csv(file_path, sep='[:,;|_]', engine='python')    

    factor_list = df.columns.tolist()
    level_factor_order_list = [] # declare nested lists of causal factors by constitution level and causal ordering:
    # level_factor_order_list[LEVEL][ORDER][FACTOR]
    level_factor_order_list.append([])   # declare lists for zeroth level
    level_factor_order_list[0].append([])  # zeroth order
    
    order_information = "" # variable to mark the row that contains the information on the causal and constitution ordering
    
    formula = ""
    # go through the data frame and get the respective logical formula
    for index, row in df.iterrows():
        formula_row = ""
        for col in factor_list:
            if row[col] in true_values:
                formula_row = formula_row + "*" + col
            elif row[col] in false_values:
                formula_row = formula_row + "*~" + col
            elif row[col] == "<" or row[col] == "<<":
                 # this row contains information on the causal and constitutional separation of the causal factors
                 order_information = row
                
        if formula_row != "":
            # add the disjunctor for the next row
            formula = formula + formula_row + " + "
    
    if not(type(order_information) == str):
        # categorise the causal factors
        level = 0 # start with level zero
        order = 0 # and order zero
        for col in factor_list:
            if str(order_information[col]).find("<<") > -1:
                # add new level
                level_factor_order_list.append([])
                level = level + 1
                level_factor_order_list[level].append([]) # add zeroth order to the new level
                order = 0
            elif str(order_information[col]).find("<") > -1:
                # add new causal phase for the current level
                level_factor_order_list[level].append([])
                order = order + 1
            
            level_factor_order_list[level][order].append(col)
    else:
        # no order information given -> all factors are of zeroth order and zeroth level
        level_factor_order_list[0][0].extend(factor_list)
    
    # remove "*" at the beginning of each disjunct
    formula = formula.replace(" *"," ")
    # remove the leading "*" and the trailing " + "
    formula = formula[1:-3]
    return level_factor_order_list, factor_list, formula

         
def read_data_from_csv(file_path):
    # (1) reads the csv file under file_path
    #     converts the truth table into a logical formula via get_truth_table_from_file
    #     determines the causal factors and any given information on constitution levels and causal order
    # (2) derives instance formulae for all causal factors that might be effects from the obtained total formula
    # (3) obtains the prime implicants for each instance formula
    # (4) transforms the instance formulae into disjunctive normal forms by means of the prime implicants
    # output:
    # abort - Boolean value, False in case that any error occurred (invalid input), True otherwise
    # level_factor_list - a nested list of causal factors subdivided by constitution levels
    # list_equiv_tuple - list of pairs, each pair corresponds to an obtained equivalence formulae
    #                    in form of: element[0] <-> element[1]
    # order_input_information - nested list of pairs, for each constitution level, the order information are translated into 
    #                           binary relations between factors (fac[0], fac[1]) means: fac[0] < fac[1]
    
    list_equiv_formula = []
    list_equiv_tuple = []
    level_factor_list = []
    order_input_information = []
    
    # determine the list of causal factors and the min-term formula from the truth table
    level_factor_order_list, factor_list, formula_st = get_truth_table_from_file(file_path)  


    # translate the information on the causal ordering within each level into the nested list order_input_information
    for lvl in range(len(level_factor_order_list)):
        order_input_information.append([]) # append a new sublist per constitution level
        for id_order in range(len(level_factor_order_list[lvl])):
            for fac in level_factor_order_list[lvl][id_order]:
                for id_order_2 in range(id_order + 1,len(level_factor_order_list[lvl])):
                    # run over all higher orders
                    for fac_2 in level_factor_order_list[lvl][id_order_2]:
                        new_pair = (fac, fac_2)
                        order_input_information[lvl].append(new_pair)
    
    
    if factor_list and formula_st != "":
        # transform the formula from string into a nested list 
        # elements are the disjuncts of the min-term formula as lists of the conjuncts each disjunct
        # formula[DISJUNCT][CONJUNCT]
        formula = string_to_list(formula_st)
        
        # determine which causal factors might be effects       
        effects_list = find_effects(string_to_list(formula_st),factor_list)
        
        # check for co-extensive factors - only for one factor of each set of co-extensive factors,
        # the prime implicants have to be determined
        list_of_coextensives = [] # this becomes a nested list: every sublist contains factors that are mutually coextensive
        for i in range(len(effects_list)-1):
            for j in range(i+1,len(effects_list)):
                co_ext = True
                for disj in formula:
                    neg_i = "~" + effects_list[i]
                    neg_j = "~" + effects_list[j]
                    if (((effects_list[i] in disj) and not(effects_list[j] in disj)) or
                    ((effects_list[j] in disj) and not(effects_list[i] in disj)) or
                    ((neg_i in disj) and not(neg_j in disj)) or ((neg_j in disj) and not(neg_i in disj))):
                        co_ext = False # two factors are not coextensive if one or its negation appears in one disjunct but the other
                        break          # does not
                
                if co_ext:
                    if not(list_of_coextensives): # if list of coextensives is still empty
                        list_of_coextensives.append([]) # append an empty sublist
                        list_of_coextensives[0].append(effects_list[i])
                        list_of_coextensives[0].append(effects_list[j])
                    else: # list is non-empty search for a sublist that contains effects_list[i]
                        new_list = True
                        for sublist in list_of_coextensives:
                            if ((effects_list[i] in sublist) or (effects_list[j] in sublist)): 
                                # at least one of both factors is already contained in one sublist of coextensive factors
                                new_list = False
                                if not(effects_list[i] in sublist): # and effects_list[i] is not,
                                    sublist.append(effects_list[i]) # then add effects_list[i] to sublist
                                elif not(effects_list[j] in sublist): # other case effects_list[j] is not contained,
                                    sublist.append(effects_list[j]) # then add it to sublist
                                break
                        if new_list: # neither factor is already contained in any sublist
                            list_of_coextensives.append([]) # create a new sublist
                            list_of_coextensives[-1].append(effects_list[i]) # append effects_list[i]
                            list_of_coextensives[-1].append(effects_list[j]) # and effects_list[j]

        # remove all but one factor of each set of coextensive factors
        if list_of_coextensives: # list is not empty
            for sublist in list_of_coextensives:
                for i in range(1,len(sublist)):
                    print('Factor ' + str(sublist[i]) +  " is coextensive with " + str(sublist[0]))
                    effects_list.remove(sublist[i])
        
        for fac in effects_list:
            # for each of these factors: 
            # determine its instance formula (the min-term equivalence formula to fac)
            i_formula = get_instance_formula_to_factor(string_to_list(formula_st), fac, level_factor_order_list)

            
            # determine the prime implicants of this instance formula
            pi_list = get_prime_implicants(i_formula, fac, level_factor_order_list)

            if pi_list:
                # list of prime implicants is non-empty
                # obtain all possible transformations of the instance formula into the reduced disjunctive normal form
                for sol in get_rdnf(pi_list, i_formula, level_factor_order_list):
                    # since sol is equivalent to fac, append the equivalence-operator and the second equivalent (fac)
                    sol = sol + " <-> " + fac
                    list_equiv_formula.append(sol)
        
        # add equivalence relations for coextensive factors
        if list_of_coextensives:
            for sublist in list_of_coextensives:
                for i in range(1,len(sublist)):
                    # replace sublist[0] by sublist[i] and vice versa in the equivalence formulae for sublist[0]
                    for formula in list_equiv_formula:
                        if formula[-1] == sublist[0]:
                            st = formula[:-1].replace(sublist[i],sublist[0]) # the new formula is the old one without the last
                            st = st + sublist[i] # expression and with all occurences of sublist[i] replaced by sublist[0]
                            # then append sublist[i] as second equivalent of the equivalence formula
                            list_equiv_formula.append(st) # add the new formula to the formulae list
        
        if list_equiv_formula:
            abort = False
            # transform the list of strings into a list of pairs, such that element[0] <-> element[1]
            list_equiv_tuple = [get_equiv_formula(formula) for formula in list_equiv_formula]
        else:
            # obtained list of equivalence formulae is empty
            abort = True
        
        
        # create level_factor_list
        for lvl in range(len(level_factor_order_list)):
            level_factor_list.append([])
            for order in level_factor_order_list[lvl]:
                level_factor_list[lvl].extend(order)
    else:
        # factor_list is empty
        abort = True
    
    return abort, level_factor_list, list_equiv_tuple, order_input_information
                   
def main():
    # main function
    file_path = "samples/Beispieldatensatz_Data_Science_1.csv"
    _, fa_l_list, list_equiv_formula, _ = read_data_from_csv(file_path)
    
    
    

    
                     
if __name__ == '__main__':
    # start main() when executing this script file
    main()
