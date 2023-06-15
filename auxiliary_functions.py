#!/usr/bin/env python3

# file: auxiliary_functions.py

import re                          # regex for complex search patterns in strings
import itertools                   # itertools provides functions to obtain all permutations of a string and Cartesian products of lists

# auxiliary functions:
def powerset(in_set):
    # returns the powerset of the given input set in_set
    aux_list = list(in_set)
    sec_aux_list = []
    for i in range(2**len(aux_list)):
        sub_list = []
        for j in range(len(aux_list)):
            if i & 2**j:
                sub_list.append(aux_list[j])
                
        sec_aux_list.append(sub_list)
    return sec_aux_list

def sort_by_second(in_tuple):
    # returns the second element of tuple (used for sort function by second tuple-element)
    return in_tuple[1]

def list_comparison(list1, list2):
    # compares two nested lists
    # returns True iff the sorted lists with sorted sublists are equal
    for subl in list1:
        subl.sort()
    for subl in list2:
        subl.sort()
    list1.sort()
    list2.sort()
    return list1 == list2
    
def find_causal_factors(st) :
    # looks in string st for denominators of causal factors, which are separated by ", " or " < "
    # returns a (possibly empty) list of causal factors
    
    # deletes "Factors: " from line (if it occurs)
    st = st.replace("Factors: ","")
    
    # deletes end-of-line-symbol ("\n") and spaces at the end of line if necessary
    st = re.sub("\r?\n","",st).rstrip()
    
    # returns the list of components of st that were separated by ", " or " < "
    return re.split(",\s*|\s*<\s*", st)
    

def get_causal_prefactors(factor, formula_list, factor_list) :
    # returns a list of direct and indirect causal prefactors of a given factor relative to a list of causal relations
    
    return_list = []
    
    
    for formula in formula_list :
        if formula[1] == factor :
            # if the factor we are interested in is the target factor of the relation formula,            
            # then add all factors that appear on the left side of formula to the return list
            return_list.extend(get_components_from_formula(formula[0], factor_list))
            
            for pfac in get_components_from_formula(formula[0], factor_list) :
                # get the indirect prefactors recursively
                return_list.extend(get_causal_prefactors(pfac, formula_list, factor_list))
    
    return_list = list(set(return_list))  # get rid of duplicates
    return return_list

    
def get_equiv_formula(st):
    # returns the leftside and rightside partial formulae of the "alpha <-> beta" from the input string st
    # as pairs of strings (a,b)
    a = re.split(" <-> ",st)[0].strip()          # strip() removes leading spaces
    b = re.split(" <-> ",st)[1].strip()
    
    # conversion of the negation syntax (in cna by minuscle) such that "a" -> "~A"
    # 1. step: add "~" before each minuscle, which is either
    # a) at the beginning of a formula
    # b) follows a conjunctor
    # c) follows a disjunctor
    a = re.sub(r'^([a-z])',  r'~\1', a)
    # explanation:  "sub" replaces each instance of a minuscle (expressed by "[a-z]")
    # by itself plus the prefix "~",
    # if it has been found at the first position of the string (implicated by "^")

    # b) if following a "*", the letter will be placed behind "*~"
    a = re.sub(r'\*([a-z])',  r'*~\1', a)
    # The regex expression "\*" picks the star symbol "*" from the string.

    # c) if following " + ", the letter will be placed behind "+ ~"
    a = re.sub(r'\s\+\s([a-z]+)',  r' + ~\1', a)
    # in regex "\s" corresponds to spaces, "\+" to "+"

    # 2. step replacement of the minuscle that follow to "~" by majuscle
    a = re.sub(r'(~[a-z]+)', lambda pat: pat.group(1).upper(), a)
    
    
    # The lines of the cna output contain further stuff, we can get rid off it:
    b = re.split("[ \t]",b)[0]
    return (a,b)

def get_components_from_formula(st, factor_list):
    # returns a list of the elements of factor_list that appear in the input string st
    # the returned list is empty if no factor from factor_list appears in st or factor_list is empty, no list of strings or no list at all
    
    component_list = []                          # declaration of the list
    
    # since we use several nested lists of causal factors, we have to treat all cases separately
    # - the factors are elements of the list as in factor_list from main -> first case
    # - the factors are elements of the elements of the list as in level_factor_list from main -> second case
    # - the factors are elements of elements of the elements of the list as in level_factor_list_order from main -> third case

    if isinstance(factor_list[0], str) :
        # first case: check the elements from factor_list
        for element in factor_list:
            if st.find(element) > -1:
                component_list.append(element)
    
    elif isinstance(factor_list[0], list) :
        if isinstance(factor_list[0][0], str) :
            # second case: traverse the sublists of factor_list and check for occurrences in st 
            for m in range(len(factor_list)) :
                for element in factor_list[m] :
                    if st.find(element) > -1 :
                        component_list.append(element)
                        
        elif isinstance(factor_list[0][0], list):
            if isinstance(factor_list[0][0][0], str) :
                # third case: traverse the subsublists of factor_list and check for occurrences in st
                for m in range(len(factor_list)) :
                    for o in range(len(factor_list[m])) :
                        for element in factor_list[m][o] :
                            if st.find(element) > -1 :
                                component_list.append(element) 
    
    return component_list

def get_factor_level(factor, level_factor_list):
    # returns the level of factor in the nested factor list that has either the form level_factor_list[LEVEL][FACTOR] 
    # or level_factor_list[LEVEL][ORDER][FACTOR]
    # if list is empty the factor has not been found return -1
    
    level = -1
    
    if not(level_factor_list):
        # level factor list is empty
        return level
    else:
        # check whether level_factor_list is further separated by causal orders
        multi_order = isinstance(level_factor_list[0], list)
        
        if multi_order:
            for m in range(len(level_factor_list)) :
                lvl_found = False
                for i in range(len(level_factor_list[m])) :
                    if factor in level_factor_list[m][i] :
                        level = m
                        lvl_found = True
                        break
                if lvl_found:
                    break        
        
        else:
            for i in range(len(level_factor_list)) :
                if factor in level_factor_list[i] :
                    level = i
                    break
        
        return level        

def get_formula_level(st, level_factor_list):
    # if all factors in st are of the same level, get_formula_level returns this level,
    # otherwise it returns -1
    
    if level_factor_list:
        # level_factor_list is non-empty
        inequal = False
        factors = get_components_from_formula(st, level_factor_list) # list of factors that occur in st
        level = -1
        
        if factors:
            # list of factors is non-empty
            # Since we use different forms of nested lists, we have to distinguish the different cases:
            # case 1: multi_order = True -> level_factor_list is subdivided into levels and causal orders
            # case 2: multi_order = False -> level_factor_list is only subdivided into levels
            multi_order = isinstance(level_factor_list[0], list)
    
    
            # first case
            if multi_order :
                # get the level of the first factor
                level = get_factor_level(factors[0], level_factor_list)
                if level > -1:
                    # check if all remaining factors have the same level as factors[0]
                    for k in range(1,len(factors)):
                        inequal = True
                        # for-loop over the orders
                        for i in range(len(level_factor_list[level])) :
                            if factors[k] in level_factor_list[level][i] :
                                inequal = False
                                break  # after we have found the factor in the right level list, we can stop the loop over the orders
        
                    if inequal:
                        level = -1
    
            # second case
            else :
                # get the level of the first factor
                level = get_factor_level(factors[0], level_factor_list)
                if level > -1:    
                    # check if all remaining factors have the same level as factors[0]
                    for k in range(1,len(factors)):
                        if not(factors[k] in level_factor_list[level]) :
                            inequal = True
                            break  # we can stop after we have found the first factor that is not of the same level as factors[0]
        
                    if inequal:
                        level = -1
    
    else:
        # level_factor_list is empty
        level = -1
        
    return level

def get_factor_order(factor, factor_list):
    # determines the order of a causal factor relative to a factor list of the form
    # list[LEVEL][ORDER][FACTORS], if factor is contained in FACTORS, return the corresponding value of ORDER in the respective level
    # otherwise return -1
    order = -1
    for m in range(len(factor_list)) :
        for o in range(len(factor_list[m])) :
            if factor in factor_list[m][o] :
                order = o
                break       # we can stop after we have 
        if order != -1 : break  # determined the order
        
    return order
    
def get_formula_order(formula, factor_list):
    # determines the order of a formula = the highest order of the factors occurring in it
    # The order of the factors is determined by factor_list using the function get_factor_order.
    # if all factors have a determinable order, the maximum value is returned,
    # otherwise -1
    
    order = -1
       
    
    for fac in get_components_from_formula(formula, factor_list):
        if order < get_factor_order(fac, factor_list) :
            order = get_factor_order(fac, factor_list)
                
    return order
    
    
def get_ordered_dnf(formula):
    # input: a string
    # applies an alphabetic order to a disjunctive normal form formula by applying the commutation rules, such that
    # (a) all conjuncts in every conjunction are sorted alphabetically
    # (b) all disjuncts in every disjunction are sorted alphabetically 
    # e.g. A < A + B + D < A + C + D < D
    # returns the transformed formula as string
    
    if not(type(formula) == str):
        # if formula is not of type string, treat it like an empty string
        formula = ""
                                         
    # split formula into disjuncts
    disj_list = re.split("\s\+\s", formula)
                                        
    # split disjuncts into conjuncts
    conj_list = [re.split("\*", disj) for disj in disj_list] # nested list conj_list[DISJUNCT][CONJUNCT IN DISJUNCT]
    for disj in conj_list:
        # sort each conjunct
        disj.sort()
    # sort each disjunct
    conj_list.sort()                                        
    # reconstruct the formula
    new_formula = ""
    for disj in conj_list:
        for conj in disj:
            new_formula = new_formula + conj + "*"
        # remove trailing "*"
        new_formula = new_formula[:-1]
        new_formula = new_formula + " + "
                                        
    # remove trailing " + "
    new_formula = new_formula[:-3]
    return new_formula

def get_clusters(formula_list, factor_list):
    # determines clusters of causally connected factors
    # two factors are causally connected, iff they either appear in the same formula of formula_list or they are connected via intermediate factors
    # returns the list of clusters, each cluster is a list of the contained factors
    
    # list for output 
    new_list_of_connected = []
    # each cluster is a sublist
    new_list_of_connected.append([])
    
    for fac in factor_list:
        if new_list_of_connected == [[]]:
            # always add first factor to zeroth cluster
            new_list_of_connected[0].append(fac)
            # add also all factors that are connected with fac to the new cluster
            for formula in formula_list:
                
                if (formula[1] == fac) or (fac in get_components_from_formula(formula[0], factor_list)):                    
                    # if fac is part of this formula
                    # add the right side
                    if not(formula[1] in new_list_of_connected[0]):
                            new_list_of_connected[0].append(formula[1])
                    # and remaining factors from the left side
                    for sec_fac in get_components_from_formula(formula[0], factor_list):
                        if not(sec_fac in new_list_of_connected[0]):
                            new_list_of_connected[0].append(sec_fac)

        else:
            already_contained = False
            for cluster in new_list_of_connected:
                # is fac already an element of a cluster?
                if fac in cluster:
                    # then stop
                    already_contained = True
                    break
            
            if not(already_contained):
                cluster_found = False
                for formula in formula_list:
                    if formula[1] == fac:
                        # fac is the right-side term
                        # look for clusters that contain any factor from the left side
                        for sec_fac in get_components_from_formula(formula[0], factor_list):
                            for cluster in new_list_of_connected:
                                if sec_fac in cluster:
                                    cluster_found = True
                                    if not(fac in cluster):
                                        cluster.append(fac)

                    elif fac in get_components_from_formula(formula[0], factor_list):
                        # fac is one part of left side-term
                        # then consider the right-side factor
                        for cluster in new_list_of_connected:
                            if formula[1] in cluster:
                                cluster_found = True
                                if not(fac in cluster):
                                    cluster.append(fac)
                        # and the other factors from the left side
                        for sec_fac in get_components_from_formula(formula[0], factor_list):
                            for cluster in new_list_of_connected:
                                if sec_fac in cluster:
                                    cluster_found = True
                                    if not(fac in cluster):
                                        cluster.append(fac)

                if not(cluster_found):
                    # create a new cluster
                    new_list_of_connected.append([])
                    # add fac to the new cluster
                    new_list_of_connected[len(new_list_of_connected)-1].append(fac)
                    # add also all factors that are connected with fac to the new cluster
                    for formula in formula_list:
                        if formula[1] == fac:
                            # fac is the term on the right side
                            for sec_fac in get_components_from_formula(formula[0], factor_list):
                                if not(sec_fac in new_list_of_connected[len(new_list_of_connected)-1]):
                                    new_list_of_connected[len(new_list_of_connected)-1].append(sec_fac)
                        
                        elif fac in get_components_from_formula(formula[0], factor_list):    
                            # fac is part of the term on the left side
                            # then add the right-side term
                            if not(formula[1] in new_list_of_connected[len(new_list_of_connected)-1]):
                                new_list_of_connected[len(new_list_of_connected)-1].append(formula[1])
                            # and the other left-side factors
                            for sec_fac in get_components_from_formula(formula[0], factor_list):
                                if not(sec_fac in new_list_of_connected[len(new_list_of_connected)-1]):
                                    new_list_of_connected[len(new_list_of_connected)-1].append(sec_fac)
                                    
            else:
                # fac is already element of one cluster
                # add all factors that are connected to fac to this cluster                        
                for formula in formula_list:
                    if formula[1] == fac:
                        # fac is the right side term
                        # so add all factors from the left side that are not already elements of this cluster
                        for sec_fac in get_components_from_formula(formula[0], factor_list):
                            if not(sec_fac in cluster):
                                cluster.append(sec_fac)
                    elif fac in get_components_from_formula(formula[0], factor_list):
                        # check the right side factor
                        if not(formula[1] in cluster):
                            cluster.append(formula[1])
                            
                        # and consider the remaining factors on the left side
                        for sec_fac in get_components_from_formula(formula[0], factor_list):
                            if not(sec_fac in cluster):
                                cluster.append(sec_fac)
    
    # merge clusters in case that they have elements in common
    if len(new_list_of_connected) > 1:
        for num_cluster in range(len(new_list_of_connected)-1,-1,-1):
            for sec_num_cluster in range(len(new_list_of_connected)):
                if num_cluster != sec_num_cluster:
                    elements_in_common = False
                    for fac in new_list_of_connected[num_cluster]:
                        if fac in new_list_of_connected[sec_num_cluster]:
                            elements_in_common = True
                            break
                    
                    if elements_in_common:
                        # add all factors from entry num_cluster to sec_num_cluster
                        for fac in new_list_of_connected[num_cluster]:
                            if not(fac in new_list_of_connected[sec_num_cluster]):
                                new_list_of_connected[sec_num_cluster].append(fac)
                                
                        # delete the cluster num_cluster
                        del new_list_of_connected[num_cluster]
                        break # leave for-loop of sec_num_cluster                          

    return new_list_of_connected        
