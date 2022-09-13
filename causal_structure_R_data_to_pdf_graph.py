#!/usr/bin/env python3

# file: causal_structure_R_data_to_pdf_graph.py

# proceeds in eight steps


import os                       # operating system interfaces is required to find the files of the own path
import sys                      # system-specific parameters and functions, needed to get optional script arguments
import codecs                   # for en- and decoding of strings (esp. to get tex-files in utf-8)
import re                       # regex for complex search patterns in strings
import jinja2                   # Latex interface
from jinja2 import Template
import itertools                # itertools provides functions to obtain all permutations of a string and Cartesian products of lists
import time                     # for measuring the running time

# syntax definitions for Latex expressions
latex_jinja_env = jinja2.Environment(
	block_start_string = '\BLOCK{',
	block_end_string = '}',
	variable_start_string = '\VAR{',
	variable_end_string = '}',
	comment_start_string = '\#{',
	comment_end_string = '}',
	line_statement_prefix = '%%',
	line_comment_prefix = '%#',
	trim_blocks = True,
	autoescape = False,
	loader = jinja2.FileSystemLoader(os.path.abspath('.')))


# auxiliary function:
def powerset(in_set):
    # returns the powerset of a given input set in_set
    aux_list = list(in_set)
    sec_aux_list = []
    for i in range(2**len(aux_list)):
        sub_list = []
        for j in range(len(aux_list)):
            if i & 2**j:
                sub_list.append(aux_list[j])
                
        sec_aux_list.append(sub_list)
    return sec_aux_list

def list_comparison(list1, list2):
    # compares two nested lists of three levels (lists of lists of tuples)
    # returns True iff the sorted lists with sorted sublists are equal
    for subl in list1:
        subl.sort()
    for subl in list2:
        subl.sort()
    list1.sort()
    list2.sort()
    return list1 == list2

def is_transitive(formula_list, factor_list):
    order_factor_list = []
    # add (empty) zeroth order
    order_factor_list.append([])
    for fac_num in range(len(factor_list)-1,-1,-1):
        first_order = True
        for formula in formula_list:
            if formula[1] == factor_list[fac_num]:
                first_order = False
                break
        
        if first_order:
            order_factor_list[0].append(factor_list[fac_num])
    
    if not(order_factor_list[0]):
        # if no factor of order zero has been found
        return False, order_factor_list  
    elif not(factor_list):
        # if the incomming list is empty
        return True, order_factor_list  
    else:
        # proceed to obtain the order of the non-zeroth order factors by using the formulae in formula_list
        ###########################################
        # step 6: define causal order iteratively #
        ###########################################
                    
        # create a total ordering of all causal factors of one level into their causal order
        # this is done iteratively starting with the yet obtained order zero by categorising all further factors
            
        # create a list of all yet non-categorised factors (= those that are causally down_stream to the considered factors)
        downstream_factor_list = []
                
        # initially downstream_factor_list consists of those elements of level_factor_list[m], that are not of order 0
        for element in factor_list:
            if not (element in order_factor_list[0]):
                downstream_factor_list.append(element)
                        
           
        # successively add to the list level_factor_list_order[m] those factors which appear on the right side of causal relation
        # whose left side factors are all already contained in level_factor_list_order[m]
        # this is done by four nested lists:
        # 1) a while loop that runs over the indexes of the elements of downstream_factor_list
        # it may have to pass the same element multiple times since it might be necessary that other factors are to
        # be classified first
        # 2) a for loop over all causal equivalence formulae
        # checks whether in all formulae where the considered element stands on the right side (as a singular term)
        # all terms on the left side have an order assigned, therefore a third loop has to been run over
        # 3) a for loop over the causal factors that appear on the left side of the current formula
        # 4) a for loop over all orders in level_factor_list_order to check whether the elements from 3) are already contained
        # in one sublist

        # traverse downstream_factor_list regressively (loop 1)
        j = len(downstream_factor_list) - 1
        while j > -1 :
            classifiable = False    # Boolean value whether a causal order can be assigned to downstream_factor_list[j] given
            # the current level_factor_list_order
            order = 0               # the order that will be given to downstream_factor_list[j]
                
            # loop 2 - for loop over all causal formulae
            for formula in formula_list:
                if formula[1] == downstream_factor_list[j] :
                    # go through all formulae where element j is the right side term
                    # then check whether in all of these formulae every factor on the left side has already an order assigned
                    # if so, element j gets the max order + 1
                    # if not, continue with the next element of downstream_factor_list
                        
                    order = 0 # order has to be reset to zero
                        
                    # loop 3 - for loop over all factors on the left side of formula
                    for fac in get_components_from_formula(formula[0], factor_list):
                        fac_is_listed = False      # is there already an order assigned to the currently considered factor fac?
                            
                        # loop 4 - for loop over all orders
                        for o in range(len(order_factor_list)) :
                            if fac in order_factor_list[o] :
                                fac_is_listed = True
                                if order < o + 1 :
                                    order = o + 1  # the order of downstream_factor_list[j] is at least one higher than that of fac
                                      
                                break              # break from loop over orders, if fac has already been found
                        # end of loop 4 over orders        
                            
                        if not(fac_is_listed) :
                            classifiable = False   # if one source factor has no assigned order, the target factor is not (yet)
                            # classifiable
                            break                  # break from loop over factors, since one non-categorised factor suffices
                                
                        else :
                            classifiable = True    # an order can be assigned to target factor j (given the current information)
                                
                    # end of loop 3 over factors of formula[0]
                        
                    if not(classifiable) :
                        break                      # break from loop over formulae after one has been found that turns out that
                        # factor j is unclassifiable by now
                            
            # end of loop 2 over all causal relations
                
            if classifiable :
                if len(order_factor_list) < order :
                    # the determined order of j is by at least 2 higher than that of the highest elements of the current list
                    # this should never happen - and possibly can't
                    print("Error in determining the order of the causal factors")
                    abort = True
                    break                          # break from while loop, something went wrong with level_factor_list_order
                else :
                    if len(order_factor_list) == order :
                        # factor j is the first one of its order
                        # a new sublist is to be created
                        order_factor_list.append([])
                       
                    # j to the ordered factor list    
                    order_factor_list[order].append(downstream_factor_list[j])     
                        
                    # and delete it from the unordered one
                    del downstream_factor_list[j]
                        
                    # proceed with the last element of downstream_factor_list
                    j = len(downstream_factor_list) - 1
                            
            else :
                j = j - 1                          # proceed with the next factor in downstream_factor_list
                        
        # end of loop 1 over the elements of downstream_factor_list
            
            
                       
        if downstream_factor_list:
            # If downstream_factor_list is not empty, there are some causal factors that cannot be categorised
            # into level_factor_list_order.
                
            # test-wise set the causal order of all remaining factors to max order + 1
            circular = True   
            order_factor_list.append([])
            
            for fac in downstream_factor_list:
                order_factor_list[len(order_factor_list) - 1].append(fac)
                    
        else :
            circular = False
    
    return not(circular), order_factor_list          
                
        

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
        for element in factor_list :
            if st.find(element) > -1 :
                component_list.append(element)
    
    elif isinstance(factor_list[0], list) :
        if isinstance(factor_list[0][0], str) :
            # second case: traverse the sublists of factor_list and check for occurrences in st 
            for m in range(len(factor_list)) :
                for element in factor_list[m] :
                    if st.find(element) > -1 :
                        component_list.append(element)
                        
        elif isinstance(factor_list[0][0], list) :
            if isinstance(factor_list[0][0][0], str) :
                # third case: traverse the subsublists of factor_list and check for occurrences in st
                for m in range(len(factor_list)) :
                    for o in range(len(factor_list[m])) :
                        for element in factor_list[m][o] :
                            if st.find(element) > -1 :
                                component_list.append(element) 
    
    return component_list

def get_formula_level(st, level_factor_list):
    # if all factors in st are of the same level, get_formula_level returns this level,
    # otherwise it returns -1
    
    inequal = False
    factors = get_components_from_formula(st, level_factor_list) # list of factors that occur in st
    level = -1
    
    # Since we use different forms of nested lists, we have to distinguish the different cases:
    # case 1: multi_order = True -> level_factor_list is subdivided into levels and causal orders
    # case 2: multi_order = False -> level_factor_list is only subdivided into levels
    multi_order = isinstance(level_factor_list[0], list)
    
    
    # first case
    if multi_order :
        # get the level of factors[0]
        for m in range(len(level_factor_list)) :
            for i in range(len(level_factor_list[m])) :
                if factors[0] in level_factor_list[m][i] :
                    level = m
            
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
        # get the level of factors[0]
        for i in range(len(level_factor_list)) :
            if factors[0] in level_factor_list[i] :
                level = i
            
        # check if all remaining factors have the same level as factors[0]
        for k in range(1,len(factors)):
            if not(factors[k] in level_factor_list[level]) :
                inequal = True
                break  # we can stop after we have found the first factor that is not of the same level as factors[0]
        
        if inequal:
            level = -1
        
    return level

def get_factor_order(factor, factor_list) :
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


def minimise_constitution_relations(level_factor_list_order, level, level_equiv_list, constitution_relation_list, color_map, mode, color_index):
    # DESCRIPTION
    
    level_count = len(level_factor_list_order)
    new_constitution_list = []
    return_list = [] 
    m = level
    
    for o in range(len(level_factor_list_order[m])) :
        # loop over all causal orders of level m

        for fac in level_factor_list_order[m][o] :
            # loop over all factors of order o and level m
                
            # create an auxiliary list of constitution relations with respect to fac
            auxiliary_list = []
            for c in constitution_relation_list :
                # loop over all constitution relations
                if c[1] == fac :
                    for l_fac in get_components_from_formula(c[0], level_factor_list_order) :
                        if not(l_fac in auxiliary_list) :
                            # add the factors on the left side of c to the auxiliary list if fac stands on its right side
                            auxiliary_list.append(l_fac) 
                    
            
            if len(auxiliary_list) > 2 :
                # further steps are only necessary if auxiliary_list contains more than two entries
                # 1) find the highest order -> factors of this order a kept for sure
                # 2) determine the lowest order that should be kept
                # a) if fac is of order 0 (on level m), factors of downto order 0 (on level m-1) should be kept
                # b) if fac is of a higher order, its lowest factors must not be in constitution relation
                # with fac's causal pre-factors
                    
                max_order = 0
                    
                # determine the value of max_order
                for l_fac in auxiliary_list :
                    if get_factor_order(l_fac, level_factor_list_order) > max_order :
                        max_order = get_factor_order(l_fac, level_factor_list_order)

                            
                if o > 0 :
                    # if the considered factor fac on the higher level is not an incomming factor, check whether its alleged
                    # constitution factors already constitute a causally upstream factor on fac's level 
                    # -> if it has been found remove it from the prospective list for fac
                        
                    # a list of all direct and indirect causal prefactors of fac
                    prefactor_list = get_causal_prefactors(fac, level_equiv_list[m], level_factor_list_order[m])
                        
                    for l in range(len(auxiliary_list) - 1, -1, -1) :
                        for pfac in prefactor_list :
                            pair = (auxiliary_list[l], pfac)
                            if pair in new_constitution_list :
                                del auxiliary_list[l]
                                break   # break from loop over prefactors, since we are done with this element fron auxiliary list
                                                   
                # otherwise if the considered factor on the higher level is an incomming factor, do not restrict the lowest order
                # of factors from auxiliary_list                         

            # carry over the entries from auxiliary_list into new_constitution_list    
            for lfac in auxiliary_list :
                entry = (lfac, fac)
                if not(entry in new_constitution_list) :
                    new_constitution_list.append(entry) 
                        
                    if "color" in mode:
                        # record constitution relation in color_map
                        color_map["draw"][lfac] = "color" + str(color_index)
                        color_map["text"][fac] = "color" + str(color_index)
                
                
            if auxiliary_list :
                # if the color has been used for fac, the next factor gets a new color index        
                color_index = color_index + 1
                if color_index > 11 : color_index = 0 # after 11 colors, use the first one again
                        
                        
                
            # clear the auxiliary list    
            auxiliary_list.clear()                    
                    
                    
    # now discard constitution relations to terms that are middle terms of causal chains whose
    # upstream and downstream factors are also in a constitution relation with the considered higher level factor
        
    # completely new loops are necessary in order to make the first one run with correct results
    # (some relations that will be discarded now, are necessary before in order to test whether some constitution relations
    # overlap with those of causal pre-factors
        
    for o in range(len(level_factor_list_order[m])) :
        # loop over all causal orders of level m

        for fac in level_factor_list_order[m][o] :
            # loop over all factors of order o and level m                
            # create a new auxiliary list of constitution relations with respect to fac
            auxiliary_list = []
                                    
            for entry in new_constitution_list :
                # obtain the results from the loop above
                if entry[1] == fac :
                    auxiliary_list.append(entry[0])
                    
            if len(auxiliary_list) > 2 :
                for l in range(len(auxiliary_list) - 1, -1, -1) :
                        
                    # get the level of auxiliary_list[l]
                    level = 0
                    for lvl in range(len(level_factor_list_order)) :
                        if any(auxiliary_list[l] in o for o in level_factor_list_order[lvl]) :
                            # in this if-clause we check whether auxiliary_list[l] is element in any sublist of
                            # level_factor_list_order[lvl] (the sublists are the order lists, which contain the causal factors)
                                
                            level = lvl
                            break          # break from the lvl loop, after it has been obtained

                    if get_causal_prefactors(auxiliary_list[l], level_equiv_list[lvl], auxiliary_list) :
                        # if factor l has causal prefactors within auxiliary_list
                          
                        # and if it is also a causal prefactor of another element in auxiliary list,
                        # it should be removed
                        for lfac in auxiliary_list :
                            if auxiliary_list[l] in get_causal_prefactors(lfac, level_equiv_list[lvl], auxiliary_list) :
                                del auxiliary_list[l]
                                break      # break from loop over the elements of auxiliary_list
                            

            # carry over the entries from auxiliary_list into return_list    
            for lfac in auxiliary_list :
                entry = (lfac, fac)
                if not(entry in return_list) :
                    return_list.append(entry) 
                
            # clear the auxiliary list    
            auxiliary_list.clear()
    return return_list, color_map, color_index

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
            # always add first factor to first cluster
            new_list_of_connected[0].append(fac)
            # add also all factors that are connected with fac to the new cluster
            for formula in formula_list:
                if (formula[1] == fac) or (fac in get_components_from_formula(formula[0], factor_list)):
                    for sec_fac in get_components_from_formula(formula[0], factor_list):
                        if not(sec_fac in new_list_of_connected[len(new_list_of_connected)-1]):
                            new_list_of_connected[len(new_list_of_connected)-1].append(sec_fac)
                    
                    # in case that fac is part of the left-side term of formula
                    # add the right side
                    if not(formula[1] in new_list_of_connected[len(new_list_of_connected)-1]):
                            new_list_of_connected[len(new_list_of_connected)-1].append(formula[1])
                    # and remaining factors from the left side
                    for sec_fac in get_components_from_formula(formula[0], factor_list):
                        if not(sec_fac in new_list_of_connected[len(new_list_of_connected)-1]):
                            new_list_of_connected[len(new_list_of_connected)-1].append(sec_fac)

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
                        # then add all terms from the left side
                        for sec_fac in get_components_from_formula(formula[0], factor_list):
                            for cluster in new_list_of_connected:
                                if sec_fac in cluster:
                                    cluster_found = True
                                    if not(fac in cluster):
                                        cluster.append(fac)

                    elif fac in get_components_from_formula(formula[0], factor_list):
                        # fac is one part of left side-term
                        # then add the right-side term
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


def rearrange_level_factor_list(in_level_factor_list_order, level_equiv_list, constitution_relation_list, mode):
    # prepares the factor list for plotting such that the nodes are arranged to minimise crossings of vertices
    # 1) factors of same level and order zero are grouped when belonging to the same constitution relation
    # 2) factors of same level and subsequent orders are arranged such that arrow crossing in the causal graphs becomes minimised
    # (in a very rudimentary way)
    # 3) discards constitution relations for middle terms
    
    
    # adapt level_factor_list_order to partial solution
    # first: create level_factor_list
    level_factor_list = []
    counter = 0
    for lvl in in_level_factor_list_order:
        level_factor_list.append([])
        for order in lvl:
            for fac in order:
                level_factor_list[counter].append(fac)        
        counter = counter + 1
    
    _,_, _, level_factor_list_order = determine_factor_order(level_factor_list, level_equiv_list)
    
    new_level_factor_list_order = []
    
    # performing 1)

    dictionary = {} # create a dictionary, for each lower level factor, the upper level factor it is a part of will be added
    
    for m in range(len(level_factor_list_order) - 1) :
        # highest level has not to be considered
        
        new_level_factor_list_order.append([]) # append a sublist for each level
        new_level_factor_list_order[m].append([]) # in the sublist for level m append a subsublist for the first order factors
        
        if constitution_relation_list :        
            # constitution_relation_list is not empty
            # use the constitution relations to group the factors that belong to the same upper factor
            
            for fac in level_factor_list_order[m][0] :
                # for each factor of first order
                constitute = ""  # the upper level factor fac is a part of
                for formula in constitution_relation_list :
                    if fac in get_components_from_formula(formula[0], level_factor_list_order) :
                        constitute = formula[1]
                        break # break from formula loop (one factor might appear in several constitution relations, but only one can be 
                        # considered for constructing an ordering)
            
                dictionary[fac] = constitute # add a new dictionary entry
    
            # fill first order of new_level_factor_list_order according to dictionary            
            
            # get first key of dictionary and add it as first element into new_level_factor_list_order[m][0]

            c_fac = ""   # dictionary[c_fac] will be used as comparison value to find the other factors with the same target
            for fac in dictionary :
                if get_formula_level(fac, level_factor_list_order) == m :
                    if dictionary[fac] != "" :
                        # start value of c_fac should be some factor with dictionary[c_fac] != "" if possible
                        c_fac = fac 
                        break   # end for loop after first dictionary entry with an assigned upper level factor has been found
            
            if c_fac == "" :
                # if no factor of first order is linked to any factor of an upper level, start with the first one from the original list
                c_fac = level_factor_list_order[m][0][0] 
            
            new_level_factor_list_order[m][0].append(c_fac)
            
            counter = 0
            while len(level_factor_list_order[m][0]) > len(new_level_factor_list_order[m][0]) :
                # as long as the dictionary contains more entries than the new factor list
                
                if not(level_factor_list_order[m][0][counter] in new_level_factor_list_order[m][0]) :
                    # if the factor with index counter is not yet in the new list
                    # add, if
                    if c_fac == "" :
                        # there is no current upper level target assigned
                        new_level_factor_list_order[m][0].append(level_factor_list_order[m][0][counter])
                        # set the upper level target to that of the indexed factor
                        c_fac = level_factor_list_order[m][0][counter]
                        # reset counter to be sure that we start at the beginning of level_factor_list_order[m][0]
                        counter = -1
                        
                    elif dictionary[level_factor_list_order[m][0][counter]] == dictionary[c_fac] :
                        # the target of the indexed factor is the same as that of c_fac
                        new_level_factor_list_order[m][0].append(level_factor_list_order[m][0][counter])
                
                # adjust the counter for the next run through the loop    
                if counter == len(level_factor_list_order[m][0]) - 1 :
                    counter = 0
                    c_fac = ""
                else :
                    counter = counter + 1
            
        # end constitution_relation_list is not empty                    
        
        else :
            # constitution_relation_list is empty
            # keep the order of factors from level_factor_list_order[m][0]
            
            for fac in level_factor_list_order[m][0] :
                new_level_factor_list_order[m][0].append(fac)

    # end for loop over the levels        
    
    max_level = len(level_factor_list_order) - 1    
    new_level_factor_list_order.append([]) # append sublist for the highest level    
    new_level_factor_list_order[max_level].append([]) # append subsublist for the zeroth order on the highest level    
    
    for fac in level_factor_list_order[max_level][0] :
        new_level_factor_list_order[max_level][0].append(fac)
        
    # performing 2)
    for m in range(len(level_factor_list_order)) :
        for o in range(len(level_factor_list_order[m]) - 1) : # go through all orders but the last
        
            new_level_factor_list_order[m].append([]) # append a subsublist for order o+1 on level m

            for fac in new_level_factor_list_order[m][o] :
                for formula in level_equiv_list[m] : 
                    if fac in get_components_from_formula(formula[0], level_factor_list_order) and (get_factor_order(formula[1], level_factor_list_order) == o+1) and not(formula[1] in new_level_factor_list_order[m][o+1]) :
                        # if the considered factor appears on the left side of formula
                        # AND the factor on formula's right side is of the subsequent order
                        # AND that factor is not in new_level_factor list yet
                        
                        new_level_factor_list_order[m][o+1].append(formula[1]) # then add this factor

            # end of loop over new_level_factor_list_order[m][o]
            
        # end of loop over orders

    # end of loop over levels
    
    # start with 3) preparing the final constitution_relation_list
    # preparing the color_map for colored graphs
    color_index = 0 # set index of first color
    # in case that the plot mode is color, a color map is created, with different specifications for text color and node color
    # this is done here, since some constitution relations to be discarded are needed for determining the colors of all nodes
    color_map = { "draw" : {}, "text" : {}}
    
    # standard color for all factors is black
    for m in range(len(new_level_factor_list_order)) :
        for o in range(len(new_level_factor_list_order[m])) :
            for e in range(len(new_level_factor_list_order[m][o])) :
                color_map["draw"][new_level_factor_list_order[m][o][e]] = "black"
                color_map["text"][new_level_factor_list_order[m][o][e]] = "black"
    
    new_constitution_relation_list = []
    for m in range(1,len(new_level_factor_list_order)) :
        partial_const_list, color_map, color_index = minimise_constitution_relations(new_level_factor_list_order, m, level_equiv_list, constitution_relation_list, color_map, mode, color_index)
        new_constitution_relation_list.extend(partial_const_list)
    
    
    if "color" in mode:
        for lvl in new_level_factor_list_order:
            for order in new_level_factor_list_order[m]:
                for fac in order:
                    colored = False
                    for formula in constitution_relation_list :
                        if formula[1] == fac:
                            # record constitution relation in color_map
                            for lfac in get_components_from_formula(formula[0], new_level_factor_list_order):
                                color_map["draw"][lfac] = "color" + str(color_index)
                            color_map["text"][formula[1]] = "color" + str(color_index)
                            colored = True  
                            
                    if colored:        
                        # if the color has been used for fac, the next factor gets a new color index        
                        color_index = color_index + 1
                        if color_index > 11 : color_index = 0 # after 11 colors, use the first one again    
    
    return new_level_factor_list_order, new_constitution_relation_list, color_map

def convert_causal_relation(formula, level_factor_list_order, tex_code, color, color_map) :
    # translates a formula of causal relations into TikZ-Latex code
    # returns the code as string
    
    st = ""  # output code
    
    ###########################################
    # determine the syntactic type of formula #
    ###########################################
    # assumption: the cna output contains only the following types of formula:
    # 1) equivalence between causal factors (A <-> B) 
    # 2) non-equivalence between causal factors (~A <-> B)
    # 3) equivalence of one causal factor with conjunctions of factors ( A*B* ... <-> E)
    # 4) equivalence of one causal factor with conjunctions of factors and negated factors ( A*~B* ... <-> E)
    # 5) equivalence of one causal factor with disjunctions of factors ( A + B + ... <-> E)
    # 6) equivalence of one causal factor with disjunctions of factors and negated factors ( A + ~B + ... <-> E)
    # 7) equivalence of one causal factor with disjunctions of at least one conjunct of factors and possibly negated factors
    #    ( A + ~B*C + ... <-> E)
    
    
    level = get_formula_level(formula[0], level_factor_list_order)
    
    #######################################
    # equivalences with (negated) factors #
    #######################################
    
    for o in range(len(level_factor_list_order[level])) :
        for fac in level_factor_list_order[level][o] :
            if fac == formula[0] :  # the left side of formula equals one factor from the factor list
                st = "\draw[->, " + color + "] (" + fac + ".east) -- (" + formula[1] + ".west);"
                break    # stop after the factor has been found
            elif formula[0] == "~" + fac : # the left side of formula equals the negation of one factor from the factor list
                color_neg = color_map["draw"][fac]
                st = "\\node[neg, " + color_neg + "] (" + fac + "neg) at ([xshift=\LNeg]" + fac + ".south east) {};\n\draw[->, " + color + "] (" + fac + "neg) -- (" + formula[1] + ");"
                break # stop after the factor has been found
                
                
    if st == "" :
        if formula[0].find("+") > -1 :
        
            ################
            # disjunctions #
            ################
            
            # create a list of the possibly complex disjuncts of formula
            disjunctor_list = re.split("\s*\+\s*", formula[0])  
            
            
            for disj in disjunctor_list : 
                # running through this list
                # each disjunct is tested whether it is
                # A) a simple causal factor
                # B) a conjunction
                # C) a negated causal factor
                # each case is treated separately

                
                if disj in get_components_from_formula(formula[0], level_factor_list_order) :
                    # case A: the discunct is one causal factor
                    
                    # a straight arrow is drawn from source factor to target factor
                    st = st + "% simple disjunction\n"
                    st = st + "\draw[->, " + color + "] (" + disj + ".east) to (" + formula[1] + ".west);\n"
                
                elif formula[0].find("*") > -1 :
                    # case B: the disjunct is a conjunction
                    st = st + "% complex disjunction\n"
                    
                    # first, the conjuncts are connected by curved lines meeting in one junction point
                    # placed one right (with a slight upward shift) to factor of the highest causal order
                    # second, this junction point is connected with the target factor by a straight arrow like
                    # in case A)
                    conjunctor_list = re.split("\*", disj)
                    
                    # set the junction of the conjuncts
                    # place it beside the (first) conjunct of the highest causal order (the factor that is most to the right in the graph)
                    # find the corresponding node of that conjunct -- f_fac
                    cross_point = ""                    # name of node of the junction of the conjunctions
                    for fac in get_components_from_formula(disj, level_factor_list_order) :
                        cross_point = cross_point + fac
                        if fac == get_components_from_formula(disj, level_factor_list_order)[0] :
                            f_fac = fac
                        else :
                            if get_factor_order(fac, level_factor_list_order) > get_factor_order(f_fac, level_factor_list_order) :
                                f_fac = fac
                    
                    cross_point = cross_point + formula[1]
                    
                    if get_factor_order(f_fac, level_factor_list_order) < get_factor_order(formula[1], level_factor_list_order) :
                        # this is the normal non-circular case
                        position = "at ([xshift=\hDisjConj, yshift=\\vDisjConj]" + f_fac + ".east)"
                        circular = False
                    else :
                        # circular case, f_fac has the same horizontal coordinate as the target factor,
                        # so the junction should not be placed to the right of f_fac but to the left
                        
                        position = "at ([xshift=\hcDisjConj, yshift=\\vDisjConj]" + f_fac + ".west)"
                        circular = True
                    
                    # Attention: It might happen that several disjuncts of conjuncts meet at the same factor f_fac,
                    # therefore we have to check whether the position of the junction node has to be shifted.
                    q = 1
                    while (tex_code.find(position) > -1) or (st.find(position) > -1) :  
                        # this position has already been specified in earlier vertices (tex_code) or this one
                        # -> shift it above by \tDisjConj
                        if circular :
                            position = "at ([xshift=\hcDisjConj, yshift={\\vDisjConj + " + str(q) + "*\\tDisjConj}]" + f_fac + ".west)"
                        else :
                            position = "at ([xshift=\hDisjConj, yshift={\\vDisjConj + " + str(q) + "*\\tDisjConj}]" + f_fac + ".east)"
                        
                        q = q + 1
                        
                    st = st  + "% junction of the conjuncts\n\\node[aux, " + color + "] (" + cross_point + "aux) " + position + " {};\n% partial arrows from the conjuncts to the junction\n"
                        
                    
                    
                        
                    

                    
                    for conj in conjunctor_list :
                    
                        # now connect the conjuncts with the junction
                        if conj[0] == "~" and conj[1:] in get_components_from_formula(disj, level_factor_list_order) :
                            # case B i) the conjunct is a negated factor
                            color_neg = color_map["draw"][conj[1:]]
                            st = st  + "\\node[neg, " + color_neg + "] (" + conj[1:] + "neg) at ([xshift=\LNeg]" + conj[1:] + ".south east) {};\n"
                            st = st + "\draw[conjunctonsegment, " + color + "] (" + conj[1:] + "neg) to (" + cross_point + "aux);\n"
                            
                        elif conj in get_components_from_formula(disj, level_factor_list_order) :
                            # case B ii) the conjunct is a mere factor
                            st = st + "\draw[conjunctonsegment, " + color + "] (" + conj + ".east) to (" + cross_point + "aux);\n"
                            
                        else :
                            # Is there anything else that might happen??
                            print("The disjunction " + formula[0] + "  could not be plotted because its substructure was not recognized(1).")

                    # connect the junction with the target factor
                    
                    # experimental label above the connection for very convoluted graphs
                    st_conj = "$"
                    for conj in conjunctor_list :
                        if conj[0] == "~" :
                            st_conj = st_conj + "\\neg " + conj[1:] + "\cdot "                    
                        else :
                            st_conj = st_conj + conj + "\cdot "
                        
                    st_conj = st_conj[:-6] + "$"
                    st = st + "% arrow from junction to target factor\n\draw[->, " + color + "] (" + cross_point + "aux) -- (" + formula[1] + ".west) node[draw=none,text=black,fill=none,font=\\tiny,pos=0,sloped,above=\LabelDist] {\scalebox{.3}{" + st_conj + "}};\n"
                
                elif (disj[0] == "~") and (disj[1:] in get_components_from_formula(formula[0], level_factor_list_order)) :
                    # case C: the disjunction is a negated factor
                    # = the first character is "~" and the further characters correspond to one element from level_factor_list_order
                    
                    # assumption: one disjunctive chain can contain either a mere or its negation (otherwise above "elif" has to been
                    # changed to "if")
                    
                    st = st + "% negated disjunct\n"
                    color_neg = color_map["draw"][disj[1:]]
                    st = st  + "\\node[neg, " + color_neg + "] (" + disj[1:] + "neg) at ([xshift=\LNeg]" + disj[1:] + ".south east) {};\n"
                    st = st + "\draw[->, " + color + "] (" + disj[1:] + "neg) to (" + formula[1] + ".west);\n"
                
                
                else :
                    # disjunction is of a different form (is there any further possible??)
                    
                    print("The disjunction " + formula[0] + "  could not be plotted because its substructure was not recognized.")

            
            
        elif formula[0].find("*") > -1 :
            
            ################
            # conjunctions #
            ################
            
            # procedure is the same as conjuncts in disjunctions
            # the only difference is that we here do not deal with a subformula of formula[0], but the whole
            
            # place junction of the conjuncts
            st = "% junction of the conjuncts\n\\node[aux, " + color + "] (" + formula[1] + "aux) at ([xshift=\LConj]" + formula[1] + ".west) {};\n% partial arrows from the conjuncts to the junction\n"
            
            # plot the arrows from the conjuncts to the junction
            for conj in get_components_from_formula(formula[0], level_factor_list_order) :

                if formula[0].find("~" + conj) > -1 :
                    # case A: the factor conj appears negated n formula
                    
                    # assumption: a conjunction chain can only contain a factor or its negation
                    # (otherwise the subsequent "else" has to be changed into a new "if")
                    color_neg = color_map["draw"][conj]
                    st = st  + "\\node[neg, " + color_neg + "] (" + conj + "neg) at ([xshift=\LNeg]" + conj + ".south east) {};\n"
                    st = st + "\draw[conjunctonsegment, " + color + "] (" + conj + "neg) to (" + formula[1] + "aux);\n"
                else :
                    # case B: factor conj occurs non-negated
                    
                    st = st + "\draw[conjunctonsegment, " + color + "] (" + conj + ".east) to (" + formula[1] + "aux);\n"
            
            # draw arrow from junction to target factor
            # experimental with tiny label above vertex
            st = st + "% arrow from junction to target factor\n\draw[->, " + color + "] (" + formula[1] + "aux) -- (" + formula[1] + ") node[draw=none, text=black, fill=none, font=\\tiny, above=\LabelDist, pos=0, sloped] {\scalebox{.3}{$" + formula[0].replace("*", "\cdot ").replace("~", "\\neg ") + "$}};\n"
            
        else :
            # this should never happen
                        
            ##############################
            # structure not recognisable #
            ##############################
            
            print(formula[0] + " -> " + formula[1] + "  could not be drawn in because its structure was not recognized.")
            
    return st
    
def convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list, color) :
    # converts a formula of constitution relations into TikZ-Latex code
    # returns the code as string
    
    st = ""
    
    # constitution relations are drawn differently depending on whether they are to the left or to the right of the upper level factor
    c_left = False
    c_right = False
    
    # check whether it is a left- or rightside relation
    for f in constitution_relation_list :
        if (formula[1] == f[1]) and (formula[0] != f[0]) :
            # is there a further constitution relation to the same causal factor which includes factors of higher causal order
            # than those from formula? -> if true it is a leftside relation
            # if there is no further constitution relation it is neighter left- nor rightside
            # if there further relations but of lower order -> rightside relation
            if get_formula_order(formula[0], level_factor_list_order) < get_formula_order(f[0], level_factor_list_order) :
                c_left = True
            elif get_formula_order(formula[0], level_factor_list_order) > get_formula_order(f[0], level_factor_list_order) :
                c_right = True
    
    # draw one connecting line toward formula[1] for each causal factor in formula[0]
    # (usually there should only be one factor in formula[0])
    for fac in get_components_from_formula(formula[0], level_factor_list_order) :           
        if c_left and not(c_right) :
            # case 1: leftside relation
            st = st + "\draw[crelationleft, " + color + "] (" + fac + ".north west) to (" + formula[1] + ".south);\n"
    
        elif not(c_left) and c_right :
            # case 2: rightside relation
            st = st + "\draw[crelationright, " + color + "] (" + fac + ".north east) to (" + formula[1] + ".south);\n"
    
        else: 
            # case 3: otherwise
            st = st + "\draw[crelationstraight, " + color + "] (" + fac + ".north) to (" + formula[1] + ".south);\n"
    
    return st

# main functions start here
def read_R_file(file_name):
    # gets the relevant data of the cna output and returns in some lists
    # level_factor_list -- list of causal factors separated into one sublist for each constitution level,
    # level_equiv_list -- list of causal relations separated into one sublist for each constitution level,
    # constitution_relation_list -- list of constitution relations
    
    
    file_lines = []                              # declaration of a list to include the lines of text of the R output file
    with open (file_name, 'rt') as text_file:    # open file_name
        for next_line in text_file:              # for each line in that file do
            file_lines.append(next_line)         # append it to the list file_lines
            
    abort = False                                # since corrupted/unexpected input cannot be used, we check at several stages
                                                 # whether we may continue or not - this what the variable abort is for
           
    ######################################## 
    # step 1: determine the causal factors #
    ########################################
    factor_list = []                             # declaration of factor list
    
    
    # Attention the following might change if the formatting of the cna output changes
    for i in range(len(file_lines)) : # search for the list of causal factors in the R output
        if file_lines[i].find("Causal ordering:") > -1:
            # case 1: if the factors are divided into several levels, cna prints "Causal ordering:"
            # then the factors are listed in the subsequent line
            factor_list = find_causal_factors(file_lines[i+1]) 
            file_line_factors = i + 1            # Later it will be helpful to know the line where the factors are listed. 
            break                                # leave for-loop after the line has been found
        
        elif file_lines[i].find("Factors:") > -1:
            # case 2: the R input does not include a separation of causal factors into different levels,
            # then the output contains "Factors:" followed by the causal factors in the same line
            factor_list = find_causal_factors(file_lines[i])
            file_line_factors = i
            break
    
    level_count = file_lines[file_line_factors].count("<") + 1 # level_count = number of constitution levels
    
    if not(factor_list) :  # factor_list is empty
        print("Abort no causal factors have been found in " + file_name)
        abort = True
        return abort, "", "", ""
    
    # continue if factor_list is non-empty
    else :
        #########################
        # step 2: find formulae #
        #########################
        
        equiv_list = []                          # declaration of the list for atomic solution formulae
            
        for line in file_lines :
            if line.count("<->") == 1 :          # cna's atomic solution formulae 
                # exactly one "<->" has been found in the line
                # read the partial formulae on its left and right side and add them to equiv_list
                equiv_list.append(get_equiv_formula(line))
                


                    
        if not(equiv_list) :  # if equiv_list is empty 
            print("Abort no formula has been found in " + file_name)
            abort = True
            return abort, "", "", ""
        else:                                    # otherwise continue
            # check whether each causal factor appears in some atomic formula, otherwise remove it from factor_list

            for k in range(len(factor_list)-1,-1,-1):
                i = 0
                found = False
                while not(found) and i < len(equiv_list):
                    found = (equiv_list[i][0].find(factor_list[k]) > -1 or equiv_list[i][1] == factor_list[k])
                    i = i + 1
                    
                if not(found):  # since the for-loop is regressive, it should be no problem to remove the elements from the list
                    # within the loop
                    print("Factor " + factor_list[k] + " has been discarded, since it does not occur in any formula.")
                    factor_list.remove(factor_list[k])  
                    
            
            
            #####################################
            # step 3: categorising the formulae #
            #####################################
            
            # separating into the constitutive levels
            
            
            
            level_factor_list = []               # declaration of new lists
            level_equiv_list = []
            constitution_relation_list = []
            
            for i in range(level_count):
                if level_count > 1: # multi-level case
                    st = re.split(" < ",file_lines[file_line_factors])[i].strip()
                else :   # single-level case
                    st = file_lines[file_line_factors]
                    
                level_factor_list.append(find_causal_factors(st))
                level_equiv_list.append([])
                
            for formula in equiv_list : 
                # add all formulae that contain only factors from level i to level_equiv_list[i]
                if get_formula_level(formula[1], level_factor_list) == get_formula_level(formula[0], level_factor_list) :
                    level_equiv_list[get_formula_level(formula[1], level_factor_list)].append(formula)
                    
                elif get_formula_level(formula[0], level_factor_list) > -1 :
                    # alle Formeln, bei denen rechts ein Element aus einer anderen Ebene mit Elementen links (alle aus gleicher Ebene)
                    # verbunden wird, werde zu constitution_relation_list hinzugefuegt
                    # assumption: only constitution relations with level difference of one are maintained
                    
                    if get_formula_level(formula[0], level_factor_list) == get_formula_level(formula[1], level_factor_list) - 1 :
                        constitution_relation_list.append(formula)
                        
                           
                # all further formulae will not be considered any longer
    return abort, level_factor_list, level_equiv_list, constitution_relation_list
    
    
    
def determine_factor_order(level_factor_list, level_equiv_list):    
    # uses the list of causal relations in level_equiv_list to determine a total causal ordering of the causal factors
    # in level_factor_list for each level separately
    # returns a Boolean value whether the process has been aborted due to some error, another Boolean whether this
    # ordering permits a unique order and the new factor list level_factor_list_order that is nested twice by level and causal order
    
    level_count = len(level_factor_list)
    abort = False
    unique = True 
    
            
    level_factor_list_order = []  # declaration of a new list for causal factors with one sublist for each level that contains one
    # sublist for each causal order in this level, wherein we find the causal factors
    # order = 0 -> incoming factors, 
    # order = i -> this is a target factor, all of its source factors are of order < i and at least one is of order i - 1
            
    # starting from here, all steps will be executed separately for each level
    for m in range(level_count):
                
        ##########################################
        # step 5: determine the incoming factors #
        ##########################################
                
                
        level_factor_list_order.append([]) 
        level_factor_list_order[m].append([])  # add the empty list for factors of order 0
                

        for element in level_factor_list[m]:
            # incoming factors do never appear on the right side of equivalence formulae
            i = 0
            found = False
            while not(found) and i < len(level_equiv_list[m]):
                found = (level_equiv_list[m][i][1] == element)
                i = i + 1
                    
            if not(found):
                # add element to incoming_factor_list
                level_factor_list_order[m][0].append(element)
        
        if not(level_factor_list_order[m][0]):
            # no incoming factors have been found
            # Exist circular relations? (A<->B; B<->A)
            # yes -> structure is non-unique
            # no -> cannot continue with input data
            transitive = True
            
            for element in level_factor_list[m] :
                e_circular = True
                for formula in level_equiv_list[m] :
                    if formula[1] == element :
                        if not((formula[1],formula[0]) in level_equiv_list[m]) :
                            e_circular = False
                  
                if e_circular :
                    transitive = False
            
            unique = transitive
            if unique :
                # no incoming factors found but no circularity determined -> abort
                print("Abort, no incoming causal factors have been found.")
                abort = True
                return abort, unique, not(transitive), level_factor_list_order    
            else:
                # create provisional zeroth order of all factors
                for element in level_factor_list[m]:
                    level_factor_list_order[m][0].append(element)
                
        else:
            # level_factor_list_order[m][0] is not empty
                        
                
            transitive, level_factor_list_order[m] = is_transitive(level_equiv_list[m], level_factor_list[m]) 
                
               
    return abort, unique, not(transitive), level_factor_list_order
    
    
    
    
    
    
def find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list) :
    # bestimmt Kausal- und Konstitutionsstruktur aus den eingegebenen Listen von Kausalfaktoren, Kausalrelationen und 
    # Konstitutionsbeziehungen
    # Ausgabe erfolgt in Form von (moeglichst) minimalen Kausal- und Konstitutionsrelationen
    
    
    level_count = len(level_factor_list_order)
    solutions_list = []
    
    for m in range(level_count) :
        #########################################################
        # Schritt 6: falsche Konstitutionsbeziehungen verwerfen #
        #########################################################
        # vor den Kausalrelationen durchfuehren um volle Zahl (auch redundante) an Aequivalenzformeln
        # fuer Umformungen nutzen zu koennen       
                                
                
                    
        #################################################
        # Schritt 7: redundante Kausalformeln streichen #
        #################################################
        # Kausalformeln, die durch Einsetzen von anderen Formeln in einander gebildet werden,
        # sollen verworfen werden  
        solutions_list.append([])       

        if len(level_equiv_list[m]) > 1:
            # check whether there are multiple causal paths to the same factor -> all but one should be discarded
            list_fac = []
            list_redundant_fac = []
            for formula in level_equiv_list[m]:
                if formula[1] in list_fac and not(formula[1] in list_redundant_fac):
                    list_redundant_fac.append(formula[1])
                elif not(formula[1] in list_fac):
                    list_fac.append(formula[1])
            
            if len(list_redundant_fac) > 0:
                list_redundant_equiv = []
                list_unique_equiv = []
                # add possibly redundant formulae to list_redundant_equiv
                counter = 0
                for fac in list_redundant_fac:
                    list_redundant_equiv.append([])
                    for formula in level_equiv_list[m]:
                        if formula[1] == fac:
                            list_redundant_equiv[counter].append(formula)
                    counter = counter + 1
                
                for formula in level_equiv_list[m]:
                    redundant = False
                    for i in range(len(list_redundant_equiv)):
                        if formula in list_redundant_equiv[i]:
                            redundant = True
                            break
                    
                    if not(redundant):
                        list_unique_equiv.append(formula)
                
  
                    
                unique_equiv = []
                unique_equiv.append(list_unique_equiv)            
                
                
                # move circular formulae from list_redundant_equiv to circular_list
                circular_list = []
                for fac in list_redundant_equiv:
                    for formula in fac:
                        if len(get_components_from_formula(formula[0], level_factor_list_order)) == 1:
                            #candidates are formulae of the form A<->B, now we have to look for B<->A
                            for sec_fac in list_redundant_equiv:
                                for sec_formula in sec_fac:
                                    if formula[1] == sec_formula[0] and sec_formula[1] == formula[0]:
                                        # pair of circular formulae found
                                        if not(formula in circular_list):
                                            circular_list.append(formula)
                                        if not(sec_formula in circular_list):
                                            circular_list.append(sec_formula)
                
                # since it might happen that list_redundant_equiv becomes empty for some factor, delete all empty entries
                for fac in list_redundant_equiv:
                    for num in range(len(fac)-1,-1,-1):
                        if fac[num] in circular_list:
                            del fac[num]
                

                
                for num in range(len(list_redundant_equiv)-1,-1,-1):
                    if list_redundant_equiv[num] == []:
                        del list_redundant_equiv[num]
                        
                
                list_circular_factors = []
                for formula in circular_list:
                    for fac in get_components_from_formula(formula[0], level_factor_list_order):
                        if not(fac in list_circular_factors):
                            list_circular_factors.append(fac)
                    if not(formula[1] in list_circular_factors):
                            list_circular_factors.append(formula[1])
                
                 
                # form the powerset of the circular formulae                
                new_circular_list = powerset(set(circular_list))                
                
                for num in range(len(new_circular_list)-1,-1,-1):
                    # accept only transitive solutions (non-circular)
                    transitive, _ = is_transitive(new_circular_list[num],list_circular_factors)
                    if not(transitive):
                        del new_circular_list[num]
                    
                    else:
                        # discard all solutions with multiple causal paths towards one factor (every path is considered to be complete, so there can't be more than one)
                        test_list = []
                        single_path = True
                        for formula in new_circular_list[num]:
                            if formula[1] in test_list:
                                single_path = False
                                break
                            else:
                                test_list.append(formula[1])
                    
                        if not(single_path):
                            del new_circular_list[num]
                        
                        else:    
                            # discard all solutions where some factor from list_circular_factors is missing
                            complete = True
                            for fac in list_circular_factors:
                                f_complete = False
                                for formula in new_circular_list[num]:
                                    if (fac in get_components_from_formula(formula[0], list_circular_factors)) or fac == formula[1]:
                                        f_complete = True
                                        break
                                if not(f_complete):
                                    complete = False
                                    break
                    
                            if not(complete):
                                del new_circular_list[num]
                            
                            else:
                                # discard solutions where initially connected factors become separated
                                
                                # list of clusters = causally connected factors
                                list_of_connected = get_clusters(circular_list, list_circular_factors)
                                
                                # second step do the same again but this time with the particular solution new_circular_list[num]
                                # in case of a valid solution the lists of clusters should be equal
                                # new list of clusters = causally connected factors
                                new_list_of_connected = get_clusters(new_circular_list[num], list_circular_factors)
                                                    
                                
                                # check whether new_list_of_connected is equal to list_of_connected
                                # this is done by directly comparing the sorted lists with sorted sublists      
                                
                                if not(list_comparison(list_of_connected, new_list_of_connected)) :  
                                    del new_circular_list[num]
              
                aux_list = list(itertools.product(*list_redundant_equiv,new_circular_list,unique_equiv))
                
                
                sec_aux_list = []
                counter = 0
                for sol in aux_list:
                    sec_aux_list.append([])
                    for term in sol:
                        if isinstance(term, list):
                            for el in term:
                                sec_aux_list[counter].append(el)
                        else:
                            sec_aux_list[counter].append(term)
                    counter = counter + 1

                solutions_list[m].extend(sec_aux_list)
            else:
                # there is no factor with several causal paths
                # search for circular formulae
                circular_list = []
                for formula in level_equiv_list[m]:
                    if len(get_components_from_formula(formula[0], level_factor_list_order)) == 1:
                        #candidates are formulae of the form A<->B, now we have to look for B<->A
                        for sec_formula in level_equiv_list[m]:
                            if formula[1] == sec_formula[0] and sec_formula[1] == formula[0]:
                                # pair of circular formulae found
                                if not(formula in circular_list):
                                    circular_list.append(formula)
                                if not(sec_formula in circular_list):
                                    circular_list.append(sec_formula)
                
                if len(circular_list) > 0:
                    print("Circularity of two factors -- code not yet written")
                    # circularity of maximal 2 factors (if there were more, the list of formulae would also contain several causal paths for one factor)
                else:
                    solutions_list[m].append(level_equiv_list[m]) 
        else:
            # there is only one formula in level_equiv_list[m] -> it will surely not be redundant
            solutions_list[m].append(level_equiv_list[m])             

    # Schleife ueber Konstitutionsebenen endet hier
    
    
    thr_aux_list = [[*row] for row in itertools.product(*solutions_list)] # transform tuples back into a list

    # removing duplicates
    frt_aux_list = []
    for it in thr_aux_list:
        if not(it in frt_aux_list):
            frt_aux_list.append(it)  
    
    
    # removing solutions that do only differ in the order of the formulae
    final_list = []
    # always consider the first solution as unique
    final_list.append(frt_aux_list[0])

    for sol in frt_aux_list:
        unique = True
        for solu in final_list:
            count = 0
            sol_found = True
            for lvl in sol:
                term_found = True
                for term in lvl:
                    if not(term in solu[count]):
                        term_found = False
                        break
                count = count + 1
            
                if not(term_found): # if one term has not been found, it is a candidate for a new unique solution
                    sol_found = False
                    break # leave the for-loop over 
                
            if sol_found:
                unique = False
                break            
        # add new unique formula to final_list
        if unique:
            final_list.append(sol)
                       
    
    return final_list
# Ende find_structure                   


def print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, constitution_relation_list, color_map) :
    # Prepares the TikZ code for plotting one solution
    # a) places the causal factors as nodes separated by causal order (horizontally) and constitution level (vertically)
    # b) adds the causal and constitution relations as vertices between nodes
    
    ######################################
    # step 8: preparing the output files #
    ######################################
    
    tex_code = "% placement of the nodes\n"
    
    ######################################################
    # step 8a) - placement of the nodes = causal factors #
    ######################################################
    
    # [possible improvement]
    # - start with the level thas has the highest number of causal orders (highest horizontal length)
    # here: start with level 1
    
    placement = ""

    for m in range(len(level_factor_list_order)) :
        # add some tex-comments in order to increase the readability of the tex-code
        tex_code = tex_code  + "% factors of level " + str(m) + ":\n"
        
        max_num_factors_order = len(level_factor_list_order[m][0])
        for o in range(len(level_factor_list_order[m])) :
            
            # placement of the factors in their causal order
            # factors of the same order are placed on top of each other
            
            tex_code = tex_code  + "% causal order " + str(o) + ":\n"
            for e in level_factor_list_order[m][o] :
                
                # add a line to tex_code in which the node is placed, its name is the same as the one of the causal factor
                # and it is displayed on a label
                tex_code = tex_code + "\\node[draw=" + color_map["draw"][e] + ", text=" + color_map["text"][e] + "] " + placement + " (" + e + ") {$" + e +"$};\n"
                
                if o == 0 :
                    # highlight incoming factors
                    tex_code = tex_code + "\hilightsource{" + e + "};\n"
                
                
                if not(any(e in get_components_from_formula(formula[0], level_factor_list_order) for formula in level_equiv_list[m])) :
                    # outgoing factors = those that have no outgoing arrows
                    # They do not appear on the complex side of a causal relation.
                    
                    # highlight outgoing factors
                    tex_code = tex_code + "\hilighttarget{" + e + "};\n"
                
                # prepare the variable placement for the next factor
                placement = "[above= \LvDist of " + e + "]"
                # the next factor of the same level and causal order will be positioned above by \LvDist
                
            # factors of the subsequent order -> next factor will be placed to the right of the bottom factor of the current order
            

            if level_factor_list_order[m][o] != []:
                placement = "[right= \LhDist of " + level_factor_list_order[m][o][0] + "]"
            
            if len(level_factor_list_order[m][o]) > max_num_factors_order :
                max_num_factors_order = len(level_factor_list_order[m][o])
                
        # factors of the subsequent level -> next factor will be shifted upwards by \iLvDist
        # more precisely: \iLvDist  + height of the current level (= max_num_factors_order * \LvDist) above the first factor of this level
        placement = "[above= {" + str(max_num_factors_order) + "*\LvDist + " + str(max_num_factors_order - 2) + "* \HeightNode  + \iLvDist} of " + level_factor_list_order[m][0][0] + "]"
        
    ##########################################################
    # step 8 b) - plot the causal and constitution relations #
    ##########################################################
    tex_code = tex_code  + "\n% causal relations\n"
    for m in range(len(level_equiv_list)) : 
        tex_code = tex_code  + "% of level "  + str(m) + "\n"
        for formula in level_equiv_list[m] :
            tex_code = tex_code  + "% formula: "  + formula[0] + " <-> " + formula[1] + "\n"
            
            # standard color is black 
            color = "black"
            if color_map["draw"][formula[1]] == color_map["draw"][get_components_from_formula(formula[0], level_factor_list_order)[0]] :
                # if the color map entry of target node and the first source node (any other would do it likewise) are identical
                # use their color for plotting the causal relation
                color = color_map["draw"][formula[1]]
                
            tex_code = tex_code + convert_causal_relation(formula, level_factor_list_order, tex_code, color, color_map) + "\n\n"
    
    tex_code = tex_code  + "\n% constitution relations\n"        
    for formula in constitution_relation_list :
        tex_code = tex_code  + "% formula: "  + formula[0] + " <-> " + formula[1] + "\n"
        
        # standard color is gray
        color = "gray"
        if color_map["text"][formula[1]] != "black" :
            color = color_map["text"][formula[1]]
        
        tex_code = tex_code + convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list, color) + "\n"  
    
    
    return tex_code   
      
def create_pdf(tex_code_table, latex_template_file, total_solutions):
   # creates a pdf of tex_code_table
   # Latex_Template assumes that tex_code_table is a table of pairs of a number (used to reference the index of the element) and
   # a string that contains valid TikZ code.
   
    output_file = "output_graph.tex"
    # defining the template (already prepared file)
    template = latex_jinja_env.get_template(latex_template_file)
    if total_solutions != len(tex_code_table):
        render = template.render(data = tex_code_table, maxnumber = str(len(tex_code_table)) + " (of " + str(total_solutions) + " in total)")
    else:
        render = template.render(data = tex_code_table, maxnumber = len(tex_code_table))
    
    # save the generated string as tex file
    f = open(output_file, 'wb')
    f.write(render.encode('utf-8'))
    f.close()
    
    
    # compile this tex file with pdflatex -- requires the Latex compiler to be installed on the executing system
    with codecs.open(str(output_file), "w","utf-8") as letter:
        letter.write(render);
        letter.close();
        os.system("pdflatex -interaction=batchmode " + str(output_file))
        print('File ' + str(output_file[:-4]) + '.pdf created.')
    
    # remove automatically generated log files
    os.remove(str(output_file[:-4]+".log"))
    os.remove(str(output_file[:-4]+".aux"))

def create_separtate_formula_list(formula_list):
   # if the optional parameter "-fl" or "--fulllist" has been set,
   # create also a text file listing all possible causal structures as tex-formulae
   
    output_file = "output_formula_list.txt"
    
    # write the formulae to output_file
    f = open(output_file, 'w')
    for formula in formula_list:
     f.write(formula)
     f.write('\n')  # line break after each formula
    f.close()      
      
                   
def main() :
    # main function
    input_file = "r_output.txt"
    latex_template_file = "Latex_Template.tex"
    start_time = time.time()
    
    level_factor_list = []               # declaration of the lists
    level_factor_order_list = []
    level_equiv_list = []
    constitution_relation_list = []
    
    # steps 1 - 3 start function read_R_file -> converts cna output into lists that are sorted by constitution level
    # of causal factors (level_factor_list), causal relations (level_equiv_list) and one list of constitution relations
    # (constitution_relation_list), if the cna output is not as expeceted stop the procedure with abort = True
    if os.path.exists(input_file) :
        abort, level_factor_list, level_equiv_list, constitution_relation_list = read_R_file(input_file)
    else :
        abort = True
        print("Error: Expected input data file " + input_file + " has not been found in path folder.")
    
    if not(abort) :
        # continue with steps 4 and 5 that determine the causal order of the factors
        abort, unique, circular, level_factor_list_order = determine_factor_order(level_factor_list,level_equiv_list)

        if not(abort) :
                # if the causal order of the factors is uniquely determinable:
                # steps 6 and 7 minisation of the causal and constitution relations
                
                mode = []  # list of special output modes depending on optional parameters (see below)
            
                # there are two plot modes possible:
                # "bw" - black/white
                # "color" - in color
                # The plot mode can be specified when running the script by adding "-c" or "-bw" respectively.
                mode.append("bw") # Standard mode is black/white.
            
                if any(arg == "-c" or arg == "--color" for arg in sys.argv) : 
                    # sys.argv is the list of arguments given when executing the script.
                    mode.append("color")                          # e.g. python script.py -c (The script ifself is one element of sys.argv.)
                    mode.remove("bw")
                elif any(arg == "-bw" or arg == "--blackwhite" for arg in sys.argv) : 
                    # Extra case for black/white, just in case that the standard mode will be changed.
                    mode.append("bw")
                
            
                # further option: Exports a second pdf-file containing the full list of possible causal structures as formulae            
                if any(arg == "-fl" or arg == "--fulllist" for arg in sys.argv) : 
                    # sys.argv is the list of arguments given when executing the script.
                    mode.append("fulllist")                       # e.g. python script.py -c -fl (The script ifself is one element of sys.argv.)
                    separate_formula_list = []                    # list of formulae in tex-code
                
                
                
                tex_table = []
                counter = 0
                
                # create a local solution_term_list for the reduced solution_list
                solution_term_list = find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list)
                
                total_solutions = len(solution_term_list)
                print("Number of solutions " + str(total_solutions))
                print("Time needed to find all solutions " + str(time.time() - start_time))                
                
                if len(solution_term_list) > 1000 :   # limiting the output to 1000 solutions
                    # remove this if-clause if more solutions are required
                    print("More than 100 solutions obtained, plotting only the first 1000.")
                    for i in range(len(solution_term_list) - 1, 999, -1) :
                        del solution_term_list[i]
                
                for sol in solution_term_list:
                    # level_equiv_list has to be adapted to sol
                    new_level_equiv_list = []
                    for i in range(len(level_factor_list)) :
                        new_level_equiv_list.append([])
        
                    for eq_lvl in sol:
                        for formula in eq_lvl :
                            new_level_equiv_list[get_formula_level(formula[0], level_factor_list)].append(formula)

                    # continue with steps 5 and 6 that determine the causal order of the factors
                    #abort, unique, circular, level_factor_list_order = determine_factor_order(level_factor_list, level_equiv_list)
                    if not(abort) :
                        # step 7 minimalisation of the constitution relations
                        # in the same step the colors of nodes and node texts are defined                                
                        #new_constitution_relation_list, color_map = find_structure(level_factor_list_order, level_equiv_list, constitution_relation_list, mode)
                          
                        # rearrange the factors in the factor list for improved placement in the plot
                        # also get rid of unnecessary constitution graphs (only the outer left and outer right part of the structure constituting a higher level factor should be drawn)
                        new_level_factor_list_order, new_constitution_relation_list, color_map = rearrange_level_factor_list(level_factor_list_order, new_level_equiv_list, constitution_relation_list, mode)
                        
                        # step 8: graphical output as a graph in pdf                        
                      
                        # generating the tex-code 
                        try:
                            st = print_structure_in_tikz_plot(new_level_factor_list_order, new_level_equiv_list, new_constitution_relation_list, color_map)
                        except:
                            st = ""
                            
                        # subtitle of the graph will be the formula in tex-math syntax
                        subtitle = "$"
                        for lvl in sol:
                            for term in lvl:
                                subtitle = subtitle + "(" + term[0].replace("*", " \cdot ").replace("~", "\\neg ") + "\leftrightarrow " + term[1] + ")\cdot"
                        
                        subtitle = subtitle[:-5] + "$"  # remove the "\cdot" at the end of the last term
                        
                        if "fulllist" in mode:          # when in fulllist mode, append the formula to the list to be exported
                            separate_formula_list.append("$" + subtitle + "$")
                          
                        subtitle = "\\tiny " + subtitle # formulae might be quite long, so subtitle should be written in tiny
                            
                        circular = False
                        if circular :
                            subtitle = "circular causal structure\\\\[4mm]" + subtitle
                        
                        # counter enumerates the solutions
                        counter = counter + 1
                        entry = (counter, st, subtitle)
                        tex_table.append(entry)
                     
                        
             
            
                # after one entry for each solution has been generated in tex_table compile the pdf
                if tex_table :
                    # if tex_table is non-empty
                    if os.path.exists(latex_template_file) :
                        create_pdf(tex_table, latex_template_file, total_solutions)
                    else :
                        abort = True
                        print("Error: Cannot plot the graphs since the expected template file " + latex_template_file + " does not exist in path folder.")
                    if "fulllist" in mode:
                        # create file
                        create_separtate_formula_list(separate_formula_list)
        
            
                print("Total execution time " + str(time.time() - start_time))
        else :
            # solution_list is empty
            # no solution survived selection of valid solutions
            print("No valid complex solution formula has been found in " + input_file + ".") 
          
                    
if __name__ == '__main__':
    main()
