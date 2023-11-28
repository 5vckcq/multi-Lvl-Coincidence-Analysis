#!/usr/bin/env python3

# file: mLCA.py

# This script generates multi-level causal-mechanistic models from Boolean data tables that are provided with
# tiered list of causal factors (assignment to constitutive levels).
# It generates causal graphs for each unique solution. These are finally exported
# into a graph in Latex TikZ-code. Therefore latex and the tikz library have to be installed on the system running this script. 

# The data tables can be processed in three different ways: A) as csv-files with Boolean data (for details on formatting requirements
# see obtain_equivalence_formulae.py), B) reading the atomic solutions from the output of the R-package cna,
# C) or from the output of the QCA package. If any of the two latter options is used, it is assumed that causal factors pertaining to
# different levels are separated by the causal ordering relation "<".
# Relations between factors of different levels are not considered as causal but constitution relations.

# It proceeds in eight main steps:
# step 1: obtains lists of the causal factors with level assignment, identifies the equivalence formulae
#         and categorises them into constitution relations and causal relations of different levels,
#         already discards all formulae that do not fit in any of these categories (function read_input)
# step 2: obtains the list of all causal structures that are compatible with the causal relations (function find_structures)
# step 3: determines a non-strict total order of the causal factors of each constitutive level
# step 4: prepares the lists for the graphical output (grouping of related causal factors, discarding of some constitution relations)
# step 5: translating the obtained structures into a graphical output via Latex


import os                          # operating system interfaces is required to find the files in the local path
import sys                         # system-specific parameters and functions, needed to get optional script arguments
import codecs                      # for en- and decoding of strings (esp. to write tex-files in utf-8)
import re                          # regex for complex search patterns in strings
import itertools                   # itertools provides functions to obtain all permutations of a string and Cartesian products of lists
import time                        # for measuring the run time

from auxiliary_functions import *
# further package that contains functions for deriving equivalence formulae from a given truth table
from obtain_equivalence_formulae import read_data_from_csv 
from plot_graph import *

def is_transitive(formula_list, factor_list):
    # checks whether the list of formula constitutes a transitive relation of the involved causal factors
    # e.g. A->B, B->C, but not A->B, B->C,C->A
    # returns the truth value
    # and order_factor_list = a list of causal phases/orders
    # if factor_list is empty, returns True and the empty list

    # a new list of the form order_factor_list[ORDER][FACTOR]
    order_factor_list = []
    # add (empty) zeroth order
    order_factor_list.append([])
    
    # run through factor_list and check for each factor
    # whether it is of zeroth order -> add it to order_factor_list[0]
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
        # stop the function return False and the empty list
        return False, order_factor_list  
        
        
    elif not(factor_list):
        # if the incomming list is empty
        # return True and the empty list
        return True, order_factor_list  
        
    else:
        # the usual case: the input list factor_list is non-empty and some of its elements have been
        # identified as of zeroth order
        # proceed to obtain the order of the non-zeroth order factors by using the formulae in formula_list
        ###################################
        # define causal order iteratively #
        ###################################
                    
        # create a total ordering of all causal factors of one level in their causal order
        # this is done iteratively starting with the yet obtained order zero by categorising all further factors
            
        # create a list of all still non-categorised factors (= those that are causally down stream to the considered factors)
        downstream_factor_list = []
                
        # initially downstream_factor_list consists of those elements of level_factor_list[m], that are not of order 0
        for element in factor_list:
            if not (element in order_factor_list[0]):
                downstream_factor_list.append(element)
                        
           
        # successively add to the list level_factor_list_order[m] those factors which appear on the right side of causal relation
        # whose left side factors are all already contained in level_factor_list_order[m]
        # this is done through four nested loops:
        # 1) a while loop that runs over the indexes of the elements of downstream_factor_list
        # it may have to pass the same element multiple times since it might be necessary to classify other factors first
        # 2) a for loop over all causal equivalence formulae
        # searches for formulae in which the considered factor is the right-side (atomic) term
        # check whether all factors appearing in the left side term have an order assigned, done using a third loop:
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
                
            circular = True                       
            # set the causal order of all remaining factors to max order + 1
            order_factor_list.append([])
            
            for fac in downstream_factor_list:
                order_factor_list[len(order_factor_list) - 1].append(fac)
   
        else :
            circular = False
    
    return not(circular), order_factor_list          
                
def reduce_structural_redundancy(factor_list, formula_list):
    # input: list of causal factors
    #        list of equivalence relations between causal factors
    #
    # reduces the list of equivalence relations to dissolve structural redundancies
    # -> all factors that are causally connected through formula_list remain causally connected,
    #    while a linear ordering of the causal factors can be defined (which is assumed not to be possible for formula_list)
    #
    # returns: truth-value of success
    #          nested list of causal factors (sub-divided by causal order) if True; initial factor_list if False
    #          reduced formula list
    
    # list of original causal connections for comparison    
    original_list_of_causally_connected = get_clusters(formula_list, factor_list)
    for cluster in original_list_of_causally_connected:
        cluster.sort()
    original_list_of_causally_connected.sort()
    
    non_circ_solution = False
    
    for redundant_formula in formula_list:
        local_formula_list = [formula for formula in formula_list if formula != redundant_formula]
        # remove one formula from formula_list
        
        local_cluster = get_clusters(local_formula_list, factor_list)
        for cluster in local_cluster:
            cluster.sort()
        local_cluster.sort()
        
        if local_cluster == original_list_of_causally_connected:
            # after removing the formula, all factors remain connected as before
            non_circular, factor_list = is_transitive(local_formula_list, factor_list)
            if non_circular:
                # and circularity has been resolved
                # then stop -> reduced solution has been found
                non_circ_solution = True
                break
            
            else:
                # still circular -> recursion with reduced formula list
                factor_list = [order for level in factor_list for order in level] # flatten factor_list (running is_transitive above
                                                                                  # created a nested factor_list[ORDER)[FACTOR]
                                                                                  # for the recursion step a simple factor_list is needed
                non_circ_solution, factor_list, local_formula_list = reduce_structural_redundancy(factor_list, local_formula_list)
    
    return non_circ_solution, factor_list, local_formula_list       


def minimise_constitution_relations(level_factor_list_order, level, level_equiv_list, constitution_relation_list, color_map, mode, color_index):
    # transforms constitution relations and discards inaccurate constitution relations
    # constitution relations are inaccurate if the lower order factors are in fact constituents of an upstream higher order factor
    # relative to the higher order factor that appears in the constitution relation
    # steps:
    # step 1 - constitution relations are rewritten, complex formulae in the constituents are split up into separate constitution
    #          relations (e.g. ('A*B','C') becomes ('A','C'),('B','C')
    #          in this step, an individual list for each higher order factor for which constitution relations exist is created
    # step 2 - A) find the highest order of constituents -> factors of this order are kept for sure
    #          B) determine the lowest order that should be kept
    #             i) if the higher level factor is of order 0 (on level m), factors of downto order 0 (on level m-1) should be kept
    #             ii) if the higher level factor is of a higher order, its lowest factors must not be in constitution relation
    #                 with the factor's causal pre-factors
    # step 3 - discard constitution relations to terms that are middle terms of causal chains whose
    #          upstream and downstream factors are also in a constitution relation with the considered higher level factor 
    
    
    level_count = len(level_factor_list_order)
    new_constitution_list = []
    return_list = [] 
    m = level
    
    for o in range(len(level_factor_list_order[m])) :
        # loop over all causal orders of level m
        
        # step 1: create a list of all constituents (= factors that appear on the left side of a constitution relation)
        # for each factor which appears as right side factor in constitution relations
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
                # step 2: only necessary if auxiliary_list contains more than two entries
                # 1) find the highest order -> factors of this order are kept for sure
                # 2) determine the lowest order that should be kept
                # a) if fac is of order 0 (on level m), factors of downto order 0 (on level m-1) should be kept
                # b) if fac is of a higher order, its lowest factors must not be in constitution relation
                # with fac's causal pre-factors
                    
                max_order = 0
                    
                # determine the value of max_order
                for l_fac in auxiliary_list:
                    if get_factor_order(l_fac, level_factor_list_order) > max_order :
                        max_order = get_factor_order(l_fac, level_factor_list_order)

                            
                if o > 0 :
                    # if the considered factor fac on the higher level is not an incomming factor, check whether its alleged
                    # constitution factors already constitute a causally upstream factor on fac's level 
                    # -> if it has been found remove it from the prospective list for fac
                        
                    # a list of all direct and indirect causal prefactors of fac
                    prefactor_list = get_causal_prefactors(fac, level_equiv_list[m], level_factor_list_order[m])
                        
                    for l in range(len(auxiliary_list) - 1, -1, -1) :
                        for pfac in prefactor_list:
                            pair = (auxiliary_list[l], pfac)
                            if pair in new_constitution_list:
                                del auxiliary_list[l]
                                break   # break from loop over prefactors, since we are done with this element fron auxiliary list
                                                   
                # otherwise if the considered factor on the higher level is an incomming factor, do not restrict the lowest order
                # of factors from auxiliary_list                         

            # carry over the entries from auxiliary_list into new_constitution_list    
            for lfac in auxiliary_list :
                entry = (lfac, fac)
                if not(entry in new_constitution_list) :
                    new_constitution_list.append(entry) 
                    
                    # record constitution relation in color_map
                    if "color" in mode:
                        color_map["draw"][lfac] = "color" + str(color_index)
                        color_map["text"][fac] = "color" + str(color_index)
                
                
            if auxiliary_list :
                # if the color has been used for fac, the next factor gets a new color index        
                color_index = color_index + 1
                if color_index > 11 : color_index = 0 # after 11 colors, use the first one again
                
            # clear the auxiliary list    
            auxiliary_list.clear()                        
                    
    # step 3: discard constitution relations to terms that are middle terms of causal chains whose
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


def arrange_factors(level_factor_list, level_equiv_list, constitution_relation_list, mode):
    # step 3:
    # determines the (non-strict) transitive order of the causal factors
    # step 4:
    # prepares the factor list for plotting such that the nodes are arranged to minimise crossings of vertices
    # A: factors of same level and order zero are grouped when belonging to the same constitution relation
    # B: factors of same level and subsequent orders are arranged such that arrow crossing in the causal graphs becomes minimised
    # (in a very rudimentary way)
    # C: discards constitution relations for middle terms
    
    ##########################################
    # step 3: ordering of the causal factors #
    ##########################################
        
    # adapt the causal order to particular solution    
    level_factor_list_order, level_equiv_list = determine_factor_order(level_factor_list, level_equiv_list)
    
    ################################################################################
    # step 4: rearranging the factors for optimised placement in the output graphs #
    ################################################################################
        
    new_level_factor_list_order = []
    
    ########################################################################################################
    # step 4A: group factors of order zero and same level if they belong to the same constitution relation #
    ########################################################################################################
    
    dictionary = {} # create a dictionary, for each lower level factor, the upper level factor it is a constituent of will be added
    
    for m in range(len(level_factor_list_order) - 1) :
        # highest level has not to be considered as constituents
        
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
        
    ##############################################################################################################################
    # step 4B: factors of same level and subsequent orders are arranged such that arrow crossing in the causal graphs is reduced #
    ##############################################################################################################################
        
    for m in range(len(level_factor_list_order)) :
        fac_counter = 0
        for o in range(len(level_factor_list_order[m]) - 1) : # go through all orders but the last
        
            new_level_factor_list_order[m].append([]) # append a subsublist for order o+1 on level m

            for fac in new_level_factor_list_order[m][o] :
                for formula in level_equiv_list[m] : 
                    if fac in get_components_from_formula(formula[0], level_factor_list_order) and (get_factor_order(formula[1], level_factor_list_order) == o+1) and not(formula[1] in new_level_factor_list_order[m][o+1]) :
                        # if the considered factor appears on the left side of formula
                        # AND the factor on formula's right side is of the subsequent order
                        # AND that factor is not in new_level_factor list yet
                        
                        new_level_factor_list_order[m][o+1].append(formula[1]) # then add this factor
                        fac_counter = fac_counter + 1

            # end of loop over new_level_factor_list_order[m][o]
            
        # end of loop over orders
        
        # add all factors that have not been categorised in any order into max order + 1
        total_level_factor_list = []
        for order in level_factor_list_order[m]:
            for fac in order:
                total_level_factor_list.append(fac)
                
        if fac_counter < len(total_level_factor_list):
            for order in new_level_factor_list_order[m]:
                for fac in order:
                    if fac in total_level_factor_list:
                        total_level_factor_list.remove(fac)
            
            if total_level_factor_list:
                new_level_factor_list_order[m].append([])
                new_level_factor_list_order[m][len(new_level_factor_list_order[m])-1].extend(total_level_factor_list)     
                       
    # end of loop over levels
    
    ##################################################
    # step 4C: reduce the constitution_relation_list #
    ##################################################
    # preparing the color_map for colored graphs
    color_index = 0 # set index of first color
    # in case that the plot mode is color, a color map is created, with different specifications for text color and node color
    # this is done here, since some constitution relations to be discarded are needed for determining the colors of all nodes
    color_map = { "draw" : {}, "text" : {}}
    
    # standard color for all factors is black
    for m in range(len(new_level_factor_list_order)) :
        for o in range(len(new_level_factor_list_order[m])) :
            for e in range(len(new_level_factor_list_order[m][o])):
                color_map["draw"][new_level_factor_list_order[m][o][e]] = "black"
                color_map["text"][new_level_factor_list_order[m][o][e]] = "black"
    
    new_constitution_relation_list = []
    for m in range(1,len(new_level_factor_list_order)) :
        partial_const_list, color_map, color_index = minimise_constitution_relations(new_level_factor_list_order, m, level_equiv_list, constitution_relation_list, color_map, mode, color_index)
        new_constitution_relation_list.extend(partial_const_list)  
    
    return new_level_factor_list_order, new_constitution_relation_list, color_map, level_equiv_list



# main functions start here
def read_input(file_name):
    # obtains the relevant data from the R output of either QCA or CNA and returns two lists:
    # level_factor_list -- list of causal factors separated into one sublist for each constitutive level,
    # equiv_list -- list of equivalence formulae
        
    # functionality depends on the output syntax of both R-packages and is adapted to CNA version 3.5.4;
    # QCA version 3.21 
    # this function has to be adjusted to any changes in the output formatting of these packages    
    
    abort = False                                    # since corrupted/unexpected input cannot be used, we check at several stages
                                                     # whether we may continue or not - this what the variable abort is for
    
    file_lines = []                                  # declaration of a list to include the lines of text from the R output file
    try:
        with open (file_name, 'rt') as text_file:    # open file_name
            for next_line in text_file:              # for each line in that file do
                file_lines.append(next_line)         # append it to the list file_lines
                
    except OSError:                                  # OSError - OS-error while opening the file
        abort = True
        print("Error: Expected input data file " + file_name + " has not been found in path folder.")
        return abort, "", ""
    except UnicodeError:                             # UnicodeError - file cannot be read as text
        abort = True
        print("Error: Input data file " + file_name + " is invalid. It contains non-unicode symbols.")
        return abort, "", ""
    except EOFError:                                 # EOFError - reached end-of-file before finding any data
        abort = True
        print("Error: Input data file " + file_name + " appears to be empty.")
        return abort, "", ""
    
    if len(file_lines) < 1:
        abort = True
        print("Error: Input data file " + file_name + " appears to be empty.")
        return abort, "", ""
    
    
    ###################################################################
    # step 1 A: check whether input conforms CNA or QCA output syntax #
    ###################################################################
    input_from_cna = False
    input_from_qca = False
    
    if file_lines[0].find("--- Coincidence Analysis (CNA) ---") > -1:
        # CNA output starts with line "--- Coincidence Analysis (CNA) ---"
        input_from_cna = True
    
    else:
        for line in file_lines:
            if line[0] == "M" and line.count("<->") == 1:
                # lines with equivalence formulae in QCA start with "M"
                input_from_qca = True
                break
    
    abort = (input_from_cna == False) and (input_from_qca == False)
    if abort:
        print("Input file does not conform to syntax of CNA, nor to QCA. Please use ")
        return abort, "", ""
           
   ##########################################
   # step 1 B: determine the causal factors #
   ##########################################
    factor_list = []                             # declaration of factor list
    st = ""
    
    if input_from_cna:  
        # Attention the following might change if the formatting of the cna output changes
        for i in range(len(file_lines)):  # search for the list of causal factors in the R output
            if file_lines[i].find("Causal ordering:") > -1:
                # case 1: if the factors are divided into several levels, cna prints "Causal ordering:"
                # then the factors are listed in the subsequent line
                st = file_lines[i+1].replace("Factors: ","") # deletes "Factors: " from line (if it occurs)
                file_line_factors = i + 1                     # It will be helpful to know the line where the factors are listed. 
                break                                         # leave for-loop after the line has been found
        
            elif file_lines[i].find("Factors:") > -1:
                # case 2: the R input does not include a separation of causal factors into different levels,
                # then the output contains "Factors:" followed by the causal factors in the same line
                st = file_lines[i].replace("Factors: ","") # deletes "Factors: " from line (if it occurs)
                file_line_factors = i
                break
                
        level_count = file_lines[file_line_factors].count("<") + 1 # level_count = number of constitutive levels
        
    elif input_from_qca:
        # by default QCA does not return a list of causal factors
        # check whether it has been added manually
        file_line_factors = -1
        for i in range(len(file_lines)):
            if file_lines[0].find("ordering") > -1:
                file_line_factors = i
                break
        
        if file_line_factors > -1:
            # case A: it has been manually added with the ordering information (using "print("ordering = ..."" in the R-script)
            st = re.split("=",file_lines[file_line_factors])[1]        # keep only text after " ... ordering ="
            st = st.strip().rstrip()                                   # remove leading and trailing spaces
            st = st[:-1]                                               # remove trailing '"'
            level_count = file_lines[file_line_factors].count("<") + 1 # level_count = number of constitutive levels
        else:
            # case B the causal factors are to be read from the returned formulae
            level_count = 0                                            # This implies that no information on level assignment is included.
            aux_fac_list = []
            for line in file_lines:
                if line.count("<->") == 1: 
                    # "<->" symbolises equivalence operator
                    aux_str = re.split(":",line)[1]  # QCA output lines start with "Mxx:", we have to get rid of this enumeration
                    formula = get_equiv_formula(aux_str)
                    # read all factors from formula, add them to aux_fac_list if they aren't already elements
                    # 1) right-side term (is always atomic)
                    if not(formula[1] in aux_fac_list):
                        aux_fac_list.append(formula[1])
                    # 2) factors from left-side term (always disjunctive normal forms)
                    # decompose formula by first obtaining list of all disjuncts
                    disj_list = re.split("\s\+\s", formula[0]) # create list of all disjuncts of left term from formula
                    for disj in disj_list:
                        # split every disjunct into its conjuncts (which are atomic or negations of atomic terms)
                        conj_list = re.split("\*", disj)
                        for element in conj_list:
                            if element[0] == "~":
                                element = element[1:]          # remove negators
                            
                            if not(element in aux_fac_list):
                                aux_fac_list.append(element)   # add factor if not already included to aux_fac_list
                    
                
            st = ", ".join(aux_fac_list)                       # in order to align aux_fac_list to the other cases, combine it
                                                               # string, use ", " as separator between factors
    
    factor_list = find_causal_factors(st)
    
    if not(factor_list) :  # factor_list is empty
        print("Abort no causal factors have been found in " + file_name)
        abort = True
        return abort, "", ""
    
    # continue if factor_list is non-empty
    else :
        ##########################
        # step 1C: find formulae #
        ##########################
        
        equiv_list = []                          # declaration of the list for atomic solution formulae
            
        for line in file_lines :
            if line.count("<->") == 1:          # CNA's/QCA's atomic solution formulae 
                # exactly one "<->" has been found in the line
                # read the partial formulae on its left and right side and add them to equiv_list
                if input_from_qca:
                    line = re.split(":",line)[1]  # QCA output lines start with "Mxx:", we have to get rid of this enumeration
                equiv_list.append(get_equiv_formula(line)) 
     
        if not(equiv_list) :  # if equiv_list is empty 
            print("Abort no formula has been found in " + file_name)
            abort = True
            return abort, "", ""
        else:                                    # otherwise continue
            # check whether each causal factor appears in at least one atomic formula, otherwise remove it from factor_list

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
            
                 
            # separate factor list by constitutive level
            level_factor_list = []   # declaration of new, nested list of the form level_factor_list[LEVEL][FACTOR]
            if level_count > 1: # multi-level case
                for i in range(level_count):
                    if input_from_cna:
                        st = re.split(" < ",file_lines[file_line_factors])[i].strip()
                    elif input_from_qca:
                        st = re.split("=",file_lines[file_line_factors])[1] # remove leading " ... ordering =" from line
                        st = st.strip().rstrip()                            # remove leading and trailing spaces
                        st = st[:-1]                                        # remove trailing '"'
                        st = re.split(" < ",st)[i].strip()
                    level_factor_list.append(find_causal_factors(st))
            else :   # single-level case
                level_factor_list.append(factor_list)                       # just use factor_list
            
            # remove all factors that have already been removed from factor_list
            for lvl in level_factor_list:
                for id_fac in range(len(lvl)-1,-1,-1):
                    if not(lvl[id_fac] in factor_list):
                        del lvl[id_fac]
            # check whether any level has an empty factor list -> in this case remove the level
            for i in range(len(level_factor_list)-1,-1,-1):
                if not(level_factor_list[i]):
                    del level_factor_list[i]
    
    return abort, level_factor_list, equiv_list    
            
            
def categorise_formulae(equiv_list, level_factor_list):
    # separates the equivalence relations in causal relations (subdivided by their constitutive level) and constitution relations
    # all further equivalence relations from equiv_list are discarded
    # input:
    # equiv_list a list of 2-tuples of strings (element [0] corresponds to the possibly complex left side of "<->",
    # element [1] to the atomic right side term
    # output:
    # level_factor_list -- list of causal factors separated into one sublist for each constitutive level,
    # equiv_list -- list of equivalence formulae
                
    level_equiv_list = [] # nested: level_equiv_list[LEVEL][FORMULA] (preparation for step 2)
            
    # create level-sublists
    for lvl in level_factor_list:
        level_equiv_list.append([])
            
    #####################################
    # step 2: categorising the formulae #
    #####################################
    # divide equiv_list into constitution and causal relations, with the latter separated the constitutive levels
    # declaration of new lists
    constitution_relation_list = []
            
                
    for formula in equiv_list : 
        # add all formulae that contain only factors from level i to level_equiv_list[i]
        if get_formula_level(formula[1], level_factor_list) == get_formula_level(formula[0], level_factor_list):
            level_equiv_list[get_formula_level(formula[1], level_factor_list)].append(formula)
                    
        elif get_formula_level(formula[0], level_factor_list) > -1 :
            # all formulae, which relate the element on the right side with factors of one different level on the left side
            # are added to constitution_relation_list
            # assumption: only constitution relations with level difference of one are maintained
                    
            if get_formula_level(formula[0], level_factor_list) == get_formula_level(formula[1], level_factor_list) - 1 :
                constitution_relation_list.append(formula)
                        
                           
        # all further formulae will not be considered any longer

    return level_equiv_list, constitution_relation_list
    
    
    
def determine_factor_order(level_factor_list, level_equiv_list):    
    # uses the list of causal relations in level_equiv_list to determine a total causal ordering of the causal factors
    # in level_factor_list for each level separately
    # returns the new factor list level_factor_list_order that is nested twice by level and causal order
    # level_factor_list_order[LEVEL][ORDER][FACTOR]
    
    level_count = len(level_factor_list)
    
            
    level_factor_list_order = []  # declaration of a new list for causal factors with one sublist for each level that contains one
    # sublist for each causal order in this level, wherein we find the causal factors
    # order = 0 -> incoming factors, 
    # order = i -> this is a target factor, all of its source factors are of order < i and at least one is of order i - 1
            
    # starting from here, all steps will be executed separately for each level
    for m in range(level_count):
                
        #############################################
        # step 4: determine the factor causal order #
        #############################################
               
        level_factor_list_order.append([]) 
        level_factor_list_order[m].append([])  # add the empty list for factors of order 0
                
                       
        non_circular, level_factor_list_order[m] = is_transitive(level_equiv_list[m], level_factor_list[m])
       
        
        if not(non_circular):
            non_circular, level_factor_list_order[m], level_equiv_list[m] = reduce_structural_redundancy(level_factor_list[m], level_equiv_list[m])
       
       
            if not(non_circular):
                print("Cannot determine causal ordering for " + str(level_factor_list_order[m]))

    return level_factor_list_order, level_equiv_list
    
    
    
    
def find_structures(in_level_factor_list, in_level_equiv_list, mode):
    # This functions combines the causal relations to solutions for the underlying multi-level structure.
    # Each solution consists of a minimal set of causal relations to causally connect every causal factor.
    # find_structures returns a list with all valid solutions in form of a list of solutions, which are lists of constitutive levels,
    # that are lists of causal relations
    # solutions_list[SOLUTION][LEVEL][CAUSAL RELATION]
    # It is the third step of the overall procedure.
    #
    # Since combining all causal relations that can logically be generated from the truth table
    # can lead to various causal relations for the same effect, as well as formulae where an effect
    # might appear as cause of its own cause, some formulae have to be ignored when creating
    # the causal structure of the mechanism. This happens with the following steps:
    # for each consitution level separate steps:
    # A: find circular sub-structures
    # B: discard as many causal relations as necessary to get rid of all simple circularities
    # C: find causal relations that have the same effect
    # D: compose the list of all valid solutions (it consists of all combinations of valid solution
    #    each consisting of one formula per effect of step A combined with one possible resolution of
    #    circularities and those formulae that are the common core of all valid structures)
    # combining the partial solutions for each constitutive level
    # E: Cartesian product of the partial solutions - also discard dublicates of solutions and solutions
    #    that only differ in the order of formulae 
    
    ##########################################################################
    # step 3: single out valid and unique solutions for the causal structure #
    ##########################################################################
    
    # test whether the causal factors of each level are completely connected
    # if not create additional virtual levels of the connected clusters
    level_factor_list = []
    virtual_level_dict = {} # dictionary that assigns the ordinal number of the corresponding virtual levels to each real level (e.g. virtual_level_dict[2] = [4, 5, 6])
    vl_counter = 0
    for lvl in range(len(in_level_factor_list)):
        level_factor_list.extend(get_clusters(in_level_equiv_list[lvl], in_level_factor_list[lvl]))
        num_clust_lvl = len(get_clusters(in_level_equiv_list[lvl], in_level_factor_list[lvl])) # number of causally unconnected clusters of causal factors in constitutive level lvl
        if num_clust_lvl > 1:
            virtual_level_dict[lvl] = []
            for i in range(0, num_clust_lvl):
                virtual_level_dict[lvl].append(vl_counter + i)
            vl_counter = vl_counter + num_clust_lvl
            
        else:
            v_lvl = lvl + vl_counter - 1
            virtual_level_dict[lvl] = [v_lvl]
    
    level_count = len(level_factor_list)
    
    if len(in_level_factor_list) == level_count:
        # no unconnected levels found
        level_equiv_list = in_level_equiv_list
    else:
        # additional levels created
        level_equiv_list = []
        for orig_lvl in range(len(virtual_level_dict)):
            if len(virtual_level_dict[orig_lvl]) == 1:
                # this level has no virtual levels
                level_equiv_list.append(in_level_equiv_list[orig_lvl]) # create new sublist for level lvl
                
            else:
                # this level has virtual levels -- subdivide in_level_equiv_list according to the factors
                for v_lvl in virtual_level_dict[orig_lvl]:
                    level_equiv_list.append([])
                    for formula in in_level_equiv_list[orig_lvl]:
                        if formula[1] in level_factor_list[v_lvl]:
                            level_equiv_list[-1].append(formula)
                    
    
    solutions_list = []
    
    # all steps for analysing the causal structures can and will be done
    # for each constitutive level separately
    for m in range(level_count):
        
        solutions_list.append([])       
        if len(level_equiv_list[m]) > 1:
            
            ###############################################
            # step 3A: find circular causal sub-structure #
            ###############################################
            # since the structures will be developped along the direction of causation from first causes
            # towards indermediate effects until the final effect(s)
            # a causal ordering in each constitutive level will be erected
            # in order to do so, circularities in the causal relations have to be resolved
            # this will not be possible in all cases - only resolvable circularities will be dealt with here
            # this is done with the next step
                                
            # find circular formulae keep them in circular_list
            circular_list = []
            for formula in level_equiv_list[m]:
                if len(get_components_from_formula(formula[0], level_factor_list)) == 1:
                    #candidates are formulae of the form A<->B, now we have to look for B<->A
                    # B<->A might not exist, even though A<->B exists, in case that a causal downstream has been defined and A < B
                    for sec_formula in level_equiv_list[m]:
                        if formula[1] == sec_formula[0] and sec_formula[1] == formula[0]:
                            # pair of circular formulae found
                            if not(formula in circular_list):
                                circular_list.append(formula)
                            if not(sec_formula in circular_list):
                                circular_list.append(sec_formula)
                           
                                
                                 
            ###############################################################################
            # step 3B: find partial sets that exhibit a transitive order and are complete #
            ###############################################################################
            # in order to resolve simple circularities, only subsets of these formulae can be retained
            # therefore, the powerset of all formulae that account for circular sub-structures
            # is searched for elements that allow for erecting a transitive causal order and
            # are complete
            # elements of the powerset are retained iff they
            # i) allow for a transitive ordering
            # ii) contain at most one causal relation per effect
            # iii) contain all factors that are listed as part of circular substructures (list_circular_factors)
            # iv) do not break causal connections between factors that are causally connected in the full list
            #     of causal relations
                
            # create a list of all causal factors of circular sub-structures
            list_circular_factors = []
            for formula in circular_list:
                for fac in get_components_from_formula(formula[0], level_factor_list):
                    if not(fac in list_circular_factors):
                        list_circular_factors.append(fac)
                if not(formula[1] in list_circular_factors):
                        list_circular_factors.append(formula[1])
           
                
            # form the powerset of the circular formulae                 
            new_circular_list = powerset(set(circular_list))
                                
            for num in range(len(new_circular_list)-1,-1,-1):
                # i) accept only transitive solutions (non-circular)
                transitive, _ = is_transitive(new_circular_list[num],list_circular_factors)
                if not(transitive):
                    del new_circular_list[num]
                    
                else:
                    # ii) discard all solutions with multiple causal paths towards one factor
                    # (every path is considered to be complete, so there can't be more than one)
                    effect_list = []
                    cause_list = {}   # dictionary cause_list[EFFECT] = [LIST OF CAUSES] for all effects with multiple causal paths
                    single_path = True
                    for formula in new_circular_list[num]:
                        if formula[1] in effect_list:
                            single_path = False
                            cause_list[formula[1]].append(formula[0])
                        else:
                            effect_list.append(formula[1])
                            cause_list[formula[1]] = [formula[0]]
                  
                    # remove all keys in cause_list which have only one-elementary lists
                    for effect in list(cause_list.keys()):
                        if len(cause_list[effect]) == 1:
                            del cause_list[effect]
                    
                    keep_solution_step_ii = True

                    if not(single_path):
                        
                        if "simple" in mode:
                            # in simple mode discard these solutions
                            keep_solution_step_ii = False
                            del new_circular_list[num]
                            
                        else:
                            # complex structures between coextensive factors should be generated
                            for id_formula in range(len(new_circular_list[num])-1,-1,-1):
                                if new_circular_list[num][id_formula][1] in cause_list:
                                    # if the effect of id_formula is listed in the cause-effect-dictionary
                                    # delete the respective term, such that only the one-directional terms of new_circular_list[num] survive
                                    del new_circular_list[num][id_formula]
                        
                            # apply special rules: 
                            # combine co-extensive factors to possible conjunction or disjunction formulae according to the rules
                            # (1) for all "atomic causes" of one effect exist also causal chains
                            # (1) A <-> C, B <-> C => A*B <-> C
                        
                            # (2) between the "atomic causes" of one effect exist no direct causal relations
                            # (2a) A <-> C, B <-> C, D <-> C => A + B + D <-> C
                            # (2b) A <-> C, B <-> C, D <-> C => A*B*D <-> C
                            # (2c) A <-> C, B <-> C, D <-> C => A*B + D <-> C (same with A or B in place of D)
                            # (2d) A <-> C, B <-> C, D <-> C => A*D + B*D <-> C (same with A,B,D interchanged)
                            for effect in cause_list:
                            
                            
                                # (1) A <-> C, ..., B <-> C => A*...*B <-> C
                                formula_left = ''
                                for cause in cause_list[effect]:
                                    formula_left = formula_left + cause + "*"
                                formula_left = formula_left[:-1] # remove the trailing "*"
                                # formula_left will be ordered to facilitate later comparison for duplicates
                                formula_left = get_ordered_dnf(formula_left)
                                formula = (formula_left, effect)
                                new_circular_list[num].append(formula)
                            
                        
                    if keep_solution_step_ii:    
                        # iii) discard all solutions where some factor from list_circular_factors is missing
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
                            keep_solution_step_ii = False
                            
                        else:
                            # iv) discard solutions where initially connected factors become separated
                                
                            # list of clusters = causally connected factors
                            list_of_connected = get_clusters(circular_list, list_circular_factors)
                                
                            # second step do the same again but this time with the particular solution new_circular_list[num]
                            # in case of a valid solution the lists of clusters should be equal
                            # new list of clusters = causally connected factors
                            new_list_of_connected = get_clusters(new_circular_list[num], list_circular_factors)
                                                    
                                
                            # check whether new_list_of_connected is equal to list_of_connected
                            # this is done by directly comparing the sorted lists with sorted sublists      
                                
                            if not(list_comparison(list_of_connected, new_list_of_connected)):
                                del new_circular_list[num]
                                keep_solution_step_ii =False
                    
                                    
            
            
            # in case that additional rule (2a) has been applied, add also the corresponding solution for (2b)
            # and check for (2c) and (2d)
            # That means: Look for all formulae in new_circular_list that contain a disjunctor.
            # Every disjunctor is the result of an application of rule (2a).
            # Applying (2b) is to replace all disjunctors by conjunctors.
            # Applying (2c) is to construct every combination of disjunctors and conjunctors (so replacing some disjunctors by conjunctors).
            # Applying (2d) is to append further conjuncts to each disjunct.
            
            # a list for additional formulae to new_circular_list - it will be added later
            additional_circular = []
            for sol in new_circular_list:
                loc_term_list = []
                for id_term in range(len(sol)):
                    if sol[id_term][0].find("*") > -1:
                        new_formula = (sol[id_term][0].replace("*"," + "), sol[id_term][1])
                        loc_term_list.append(new_formula)
                    else:
                        loc_term_list.append(sol[id_term])
                additional_circular.append(loc_term_list)
            
            for sol in additional_circular:
                dis_terms = [] # list of all terms of sol that contain disjunctors
                for term in sol:
                    if term[0].find("+") > -1:
                        dis_terms.append(term)
                        
                
                if len(dis_terms) > 0:
                    # proceed with sol only if it contains at least one disjunctive formula
                    
                    local_sol = [] # nested list for each formula local_sol[TERM][FORMULA]
                    # first level: number of the term of the original formula that es being modified
                    # second level: is a list of all possible terms it can be changed into
                    
                    for i in range(len(dis_terms)):
                        # for-loop over the number of terms that contain disjunctors
                        list_of_factors = get_components_from_formula(dis_terms[i][0], level_factor_list)
                        
                        local_sol.append([]) # add a new sublist for the i-th term of sol      
                        # now fill this sublist with all combinations of disjunctors and conjunctors in the respective term:
                        for j in range(2**dis_terms[i][0].count("+")):
                            # create 2^#disjunctors entries
                            # formulae will be numbered by a binary scheme
                            # junctors no.
                            # ++...+++  0
                            # ++...++*  1
                            # ++...+*+  2
                            # ++...+**  3
                            # ...
                            # **...*** 2^#disjunctors - 1
                            
                            # start all with first factor - junctors and further factors will be added incrementally
                            local_sol_f = (list_of_factors[0], dis_terms[i][1])
                            local_sol[i].append(local_sol_f)
                        
                        for id_for in range(len(local_sol[i])):
                            # loop over all formulae
                            counter_for = id_for
                            for id_fac in range(1,len(list_of_factors)):
                                # for loop over the factors in the formula, except the first one, which has been already
                                # added to every new_formula
                                counter_fac = len(list_of_factors) - id_fac - 1 # reverse the numbering of factors (second becomes
                                # last etc.)
                                if counter_for >= 2**counter_fac:
                                    # decide whether the factor id_fac is to be joined via conjunction or disjunction
                                    local_sol[i][id_for] = (local_sol[i][id_for][0] + "*" + list_of_factors[id_fac], dis_terms[i][1])
                                    counter_for = counter_for - 2**counter_fac
                                else:
                                    local_sol[i][id_for] = (local_sol[i][id_for][0] + " + " + list_of_factors[id_fac], dis_terms[i][1])
                        
                        # By now (2b) and (2c) are done. All necessary formulae became elements of local_sol[i].
                        # Continuing with (2d):
                        
                        # Form every possible disjunct between atomic and maximal conjunct (A, A*B, A*B*C, ...)
                        # But take care that no disjunct becomes a subset of another disjunct (NOT A*B + B + ...)
                        # Also take substitutions by further formulae into account.
                        for formula in local_sol[i]:
                            # run over all formula in local_sol[i] to check for eligible disjunctive formulae
                            if formula[0].find("+") > -1:
                                # only proceed with formula if it contains at least one disjunctor
                                
                                f_disj_list = re.split("\s\+\s", formula[0]) # create list of all disjuncts of formula
                                
                                # nested list f_conj_list[DISJUNCT][CONJUNCT IN DISJUNCT]
                                f_conj_list = conj_list = [re.split("\*", disj) for disj in f_disj_list] 
                                new_disj_list_2d = [] # list of additional terms due to rule (2d)
                                # this list will be nested new_disj_list_2d[DISJUNCT][CONJUNCT]
                                
                                for id_disj in range(len(f_disj_list)):
                                    new_disj_list_2d.append([]) # create empty entry for next disjunct
                                    
                                    # the list of factors that can be added as further conjuncts
                                    fac_to_be_added = [fac for fac in list_of_factors if not(fac in f_conj_list[id_disj])]
                                    
                                    # the list of all possible forms between atomic and the maximal conjunct is determined
                                    # by using the Cartesian product of the present factor and the powerset of fac_to_be_added
                                    aux_list_2d = [[sublist_1, sublist_2] for (sublist_1, sublist_2) in itertools.product([f_conj_list[id_disj]], powerset(fac_to_be_added))]
                                    # structure of this set: e.g. for A + B + C -> [[[['A'], []], [['A'], ['B]], [['A'], ['C']], [['A'], ['B', 'C']]], [[['B'], []] ...] ... ]
                                    
                                    # flatten the inner lists (e.g. [['A'], []] -> ['A'] and [['A'], ['B','C']] -> ['A','B','C']
                                    sec_aux_list_2d = []
                                    for id_element in range(len(aux_list_2d)):
                                        sec_aux_list_2d.append([])
                                        for ij in range(2):
                                            for subelement in aux_list_2d[id_element][ij]:
                                                sec_aux_list_2d[id_element].append(subelement)
                                        sec_aux_list_2d[id_element].sort()
                                    
                                    # add these newly obtained disjuncts to the list of complete terms
                                    new_disj_list_2d[id_disj].extend(sec_aux_list_2d)
                                    
                                # the totality of new DNF formulae is the Cartesian product of all variants for each disjunct
                                sec_new_disj_list_2d = [list(x) for x in list(itertools.product(*new_disj_list_2d))]
                                    
                                # now discard all invalid formulae = one disjunct is a subset of another disjunct
                                # e.g.: 
                                # A*B + A*B*C
                                for e_counter in range(len(sec_new_disj_list_2d)-1,-1,-1):
                                    # run through the list and check for every term individually
                                    invalid = False
                                    for d_counter_1 in range(len(sec_new_disj_list_2d[e_counter])):
                                        # for every disjunct in the term
                                        for d_counter_2 in range(len(sec_new_disj_list_2d[e_counter])):
                                            # compare with each other disjunct
                                            if d_counter_1 != d_counter_2:
                                                # whether it is a subset of the other
                                                if set(sec_new_disj_list_2d[e_counter][d_counter_1]).issubset(sec_new_disj_list_2d[e_counter][d_counter_2]):
                                                    # if True, then discard this term
                                                    invalid = True
                                                    break
                                                    
                                                
                                        if invalid:
                                            break
                                    if invalid:
                                        del sec_new_disj_list_2d[e_counter]
                                # add formulae in sec_new_disj_list to local_sol[i]
                                for new_term in sec_new_disj_list_2d:
                                    # convert each new_term into a string of logical formula
                                    str_formula = ""
                                    for disj in new_term:
                                        if str_formula != "":
                                            str_formula = str_formula + " + "
                                        
                                        for conj in disj:
                                            str_formula = str_formula + conj + "*"
                                            
                                        str_formula = str_formula[:-1]
                                    
                                    
                                    
                                    # add the newly obtained term to local_sol[i] if it is not already contained
                                    compl_formula = (str_formula, dis_terms[i][1])
                                    if not(compl_formula in local_sol[i]):
                                        local_sol[i].append(compl_formula)
                        
                        
                        
                    # form Cartesian of all sublists of local_sol and the unchanged terms
                    sol_base = [x for x in sol if x[0].find("+") == -1] # list of all terms of sol that do not contain disjunctors

                    # and are the common core of every result for sol
                    if len(sol_base) > 0:
                        # add the Cartesian product of sol_base and all sublists of local_sol to
                        # additional_circular if they are not already contained therein
                        additional_circular.extend([[*line] for line in itertools.product(*local_sol,sol_base) if not([*line] in additional_circular)])
                    else:
                        # in case that the common core is empty, add the Cartesian product of all sublists of local_sol to
                        # additional_circular if they are not already contained therein
                        additional_circular.extend([[*line] for line in itertools.product(*local_sol) if not([*line] in additional_circular)])

            # add additional terms without duplicates to new_circular_list          
            # reconstruct new
            for sol in new_circular_list:
                sol.sort(key=sort_by_second)
                
            for sol in additional_circular:
                sol.sort(key=sort_by_second)
                if not(sol in new_circular_list):
                    new_circular_list.append(sol)
                                
            ########################################################################################
            # step 3C: find mutually exclusive causal relations (same effect but different causes) #
            ########################################################################################
            # each causal relation is supposed to be a complete causal path towards the effect (= the right side term of a causal relation)
            # hence, there cannot be two or more causal relations with the same effect
            # each solution can only contain one of this set of causal relations
            
            # check whether there are multiple causal paths to the same factor -> all but one should be discarded
            list_fac = []
            list_redundant_equiv = []
            list_unique_equiv = []
            
            # list_redundant_fac will contain the effects for which several causal relations exist
            list_redundant_fac = []
            for formula in level_equiv_list[m]:
                if not(formula in circular_list) and (formula[1] in list_fac and not(formula[1] in list_redundant_fac)):
                    list_redundant_fac.append(formula[1])
                elif not(formula in circular_list) and not(formula[1] in list_fac):
                    list_fac.append(formula[1])
            
            # the complement set will be unique_equiv
            unique_equiv = []
            
            # find causal relations of which has to be chosen
            if len(list_redundant_fac) > 0:
                # there are causal factors for which causal relations have to be discarded
                
                # add possibly redundant formulae to list_redundant_equiv
                counter = 0
                for fac in list_redundant_fac:
                    list_redundant_equiv.append([])
                    for formula in level_equiv_list[m]:
                        if not(formula in circular_list) and formula[1] == fac:
                            list_redundant_equiv[counter].append(formula)
                    counter = counter + 1
                
                for formula in level_equiv_list[m]:
                    redundant = False
                    for i in range(len(list_redundant_equiv)):
                        if formula in list_redundant_equiv[i]:
                            redundant = True
                            break
                    
                    if not(redundant) and not(formula in circular_list):
                        list_unique_equiv.append(formula)
                
  
                # unique_equiv contains all causal relations that appear in every solution
                # = the only causal relation for the respective effect
                if list_unique_equiv:
                    unique_equiv.append(list_unique_equiv) 
            
            else:
                # there are no redundant formulae
                # all formulae from level_equiv_list[m] that are not elements of circular_list should go into unique_equiv
                unique_equiv = [[formula for formula in level_equiv_list[m] if not(formula in circular_list)]]
            
            #####################################################################
            # step 3D: compose a list of solutions for the constitutive level m #
            #####################################################################                
            # this will be the Cartesian product of list_redundant_equiv, new_circular_equiv and unique_equiv
            # in case that one or two of these lists are empty, different cases have to be distinguished:
            
            if new_circular_list:
                if list_redundant_equiv:
                    if unique_equiv:
                        aux_list = list(itertools.product(*list_redundant_equiv,new_circular_list,unique_equiv))
                    else:
                        if list_redundant_equiv:
                            aux_list = itertools.product(*list_redundant_equiv,new_circular_list)

                        else:
                            aux_list = []
                            aux_list.extend(new_circular_list)
                else:
                    if unique_equiv:
                        aux_list = list(itertools.product(new_circular_list,unique_equiv))
                    else:
                        aux_list = []
                        aux_list.extend(new_circular_list)
            else:
                if list_redundant_equiv:
                    if unique_equiv:
                        aux_list = list(itertools.product(*list_redundant_equiv,unique_equiv))
                    else:
                        aux_list = list(*list_redundant_equiv)
                else:
                    aux_list = []
                    aux_list.extend(unique_equiv)
                
            # avoid problems with Cartesian products of lists and tuples
            # make sure that every component is a list
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
            
            # equations from redundant_equiv might have the same effect as formulae from the circular group
            # in this case the formula from redundant_equiv has to be removed from the solution
            for sol in solutions_list[m]:
                delete_list = []
                for f_formula in sol:
                    for s_formula in sol:
                        if f_formula != s_formula and f_formula[1] == s_formula[1]:
                            if len(f_formula[0]) > len(s_formula[0]):
                                # retain the shorter formula (the one that stems from the circular group)
                                delete_list.append(f_formula)
                            else:
                                delete_list.append(s_formula)
                
                delete_list = list(set(delete_list)) # get rid of duplicates in delete_list

                for formula in delete_list:
                    sol.remove(formula)

            # check solutions for completeness every causal factor that apppears in level_equiv_list
            # has to appear in at least one causal formula
            for i in range(len(solutions_list[m])-1,-1,-1):
                all_found = True
                
                for fac in level_factor_list[m]:
                    # test whether fac appears in level_equiv_list[m]
                    fac_required = False
                    if level_equiv_list[m]:
                        for formula in level_equiv_list[m]:
                            if formula[1] == fac or fac in get_components_from_formula(formula[0], level_factor_list):
                                fac_required = True
                                break
                    
                    if fac_required:
                        fac_found = False
                        for formula in solutions_list[m][i]:
                            if formula[1] == fac or fac in get_components_from_formula(formula[0], level_factor_list):
                                fac_found = True
                                break
                        if not(fac_found):
                            all_found = False
                            break
                        
                if not(all_found):
                    del solutions_list[m][i]
        else:
            # there exists only one formula in level_equiv_list[m] or no formula at all -> it will surely not be redundant or circular
            solutions_list[m].append(level_equiv_list[m])             

    # for-loop over constitutive levels ends here
       
    #####################################################################
    # step 3E: combine the partial solutions of the constitutive levels #
    #####################################################################    

    # merge virtual levels into their original form
    if len(in_level_factor_list) != len(level_factor_list):
        # virtual levels have been created
        orig_solutions_list = []
        for orig_lvl in range(len(in_level_factor_list)):
            if len(virtual_level_dict[orig_lvl]) == 1:
                # this constitutive level has no virtual levels
                orig_solutions_list.append(solutions_list[virtual_level_dict[orig_lvl][0]])
            else:
                # this constitutive level has virtual levels to be merged
                new_sublist = []
                for v_lvl in virtual_level_dict[orig_lvl]:
                    new_sublist.append(solutions_list[v_lvl])
                orig_solutions_list.append([])
                orig_solutions_list[-1] = [list(x) for x in itertools.product(*new_sublist)]
                
                # flatten the second level of orig_solutions_list (by introducing new variable loc_list, which become the set-sum of orig_solutions_list[j] subsets)
                for j in range(len(orig_solutions_list[-1])):
                    loc_list = []
                    for i in range(len(orig_solutions_list[-1][j])):
                        loc_list = loc_list  + orig_solutions_list[-1][j][i]
                
                    orig_solutions_list[-1][j] = loc_list
               
        final_list = [[*row] for row in itertools.product(*orig_solutions_list)] # transform tuples back into a list   
    
    else:
        # no virtual levels created
        # construct the Cartesian product of the partial solutions
        final_list = [[*row] for row in itertools.product(*solutions_list)] # transform tuples back into a list   
    
    # due to previous sorting, no dublicates are in the list
    # high computational effort for searching for duplicates (with taking into account the possibility of commutations)
    # is unnecessary
    dublicate_list = []
    for i in range(len(final_list)-1,-1,-1):
        if final_list[i] in dublicate_list:
            del final_list[i]
        else:
            dublicate_list.append(final_list[i])
            
    final_list.sort()
    return final_list
# End of find_structures                  


   
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
    latex_template_file = "Latex_Template.tex"
    start_time = time.time()
    
    level_factor_list = []               # declaration of the factor list    
    order_input_information = []
    abort = False
    
    # sys.argv is the list of arguments given when executing the script.
    # e.g. "python script.py --csv" contains two arguments "script.py" (the script name ifself is one element of sys.argv) and "--csv"
    for i in range(len(sys.argv)):
        if sys.argv[i] == "--r-import":
            # "--r-import" is the command to use the text output of a r-script (cna or QCA) 
            if i < len(sys.argv):
                input_file = sys.argv[i+1] # assume that the next argument is path to the r-output file
    
                # step 1 function read_input -> converts cna/QCA output into lists of causal factors and equivalence relations
                # if the output is not as expected, stop the procedure with abort = True
                if os.path.exists(input_file):
                    abort, level_factor_list, equiv_list = read_input(input_file)
                    if not(abort):
                        # classify the equivalence relations into causal relations (level_equiv_list) subdivided by constitutive level
                        # or into constitution relations
                        level_equiv_list, constitution_relation_list = categorise_formulae(equiv_list, level_factor_list)
                        
                else:
                    # input_file does not exist
                    abort = True
                    print("Error: Expected input data file " + input_file + " has not been found in path folder.")
            
            else:
                # no input_file specified
                abort = True
                print("Error: Missing expected file path after \"--r-import\" command.")
                
            
            break # break from loop over script parameters after "--r-import" has been found
        
        
        elif sys.argv[i] == "--csv": 
            # Read truth table from csv-file and generate equivalence formulae
            if i < len(sys.argv):
                input_file = sys.argv[i+1] # assume that the next argument is path to the csv table
                if os.path.exists(input_file):
                    # execute function from obtain_equivalence_formulae.py
                    abort, level_factor_list, equiv_list, order_input_information = read_data_from_csv(input_file)
                        
                    if not(abort):
                        # classify the equivalence relations into causal relations (level_equiv_list) subdivided by constitutive level
                        # or into constitution relations
                        level_equiv_list, constitution_relation_list = categorise_formulae(equiv_list, level_factor_list)

                    else:
                        print("Error while reading the input file " + input_file)
                
                else:
                    # input_file does not exist
                    abort = True
                    print("Error: Expected input data file " + input_file + " has not been found in path folder.")
                
            
            else:
                # no input_file specified
                abort = True
                print("Error: Missing expected file path after \"--csv\" command.")
            
            break # break from loop over script parameters after "--csv" has been found
    
    if not(level_factor_list) and not(abort):
        # level_factor_list is still empty -> probably neither "--csv", nor "--r-import" has been found as argument
        # in case that neither "--csv", nor "--r-import" parameter has been specified raise error, due to
        # missing input data
        abort = True
        print("Error: Expected input data not found\n Use either \"--csv FILEPATH_TO_CSV_TABLE\" or \"--r-import FILEPATH_TO_R_OUTPUT\" to specify the path to the input data.")    
    
    if not(abort) :                
        # check whether special options like color mode or output of the full list of formulae has been activated
        mode = []  # list of special output modes depending on optional parameters (see below)
            
        # prepare these special modes
        # there are two plot modes possible:
        # "bw" - black/white
        # "color" - in color
        # The plot mode can be specified when running the script by adding "-c" or "-bw" respectively.
        mode.append("color") # Standard mode is color.
            
        if any(arg == "-c" or arg == "--color" for arg in sys.argv) : 
            # sys.argv is the list of arguments given when executing the script.
            if not("color" in mode):
                mode.append("color")                          # e.g. python script.py -c (The script ifself is one element of sys.argv.)
        elif any(arg == "-bw" or arg == "--blackwhite" for arg in sys.argv) : 
            # further case for black/white
            if not("bw" in mode):
                mode.append("bw")
                # remove the standard mode
                mode.remove("color")
                
            
        # further option: Exports a second pdf-file containing the full list of possible causal structures as formulae            
        if any(arg == "-fl" or arg == "--fulllist" for arg in sys.argv) : 
            # sys.argv is the list of arguments given when executing the script.
            mode.append("fulllist")                       # e.g. python script.py -c -fl (The script ifself is one element of sys.argv.)
            separate_formula_list = []                    # list of formulae in tex-code
                
        # simple mode does not derive complex structures between co-extensive factors
        # results of simple mode should be the same as those of cna
        if any(arg == "-c" or arg == "--complex" for arg in sys.argv) : 
            pass                                          # e.g. python script.py -c (The script ifself is one element of sys.argv.)
        else:                                             # default behaviour: simple mode without complex structures between
            mode.append("simple")                         # co-extensive factors

        # prepare the list for the tex-code        
        tex_table = []
                
        # step 3 generate all possible solutions for the underlying causal structure
        solution_term_list = find_structures(level_factor_list, level_equiv_list, mode)
                
        total_solutions = len(solution_term_list)
        print("Time needed to find all " + str(total_solutions) + " solutions " + str(round(time.time() - start_time,2)) + " seconds.")                
         
        # prepare each solution individually for graphical output
        sol_counter = 0
  
        for sol in solution_term_list:            
            # prepare a local version of level_equiv_list
            new_level_equiv_list = []
            for i in range(len(level_factor_list)) :
                new_level_equiv_list.append([])
       
            for eq_lvl in sol:
                for formula in eq_lvl :
                    new_level_equiv_list[get_formula_level(formula[0], level_factor_list)].append(formula)

            # step 4: determine the causal order and
            # step 5: arrange the factors for improved placement in the plot
            # also get rid of unnecessary constitution graphs
            # (only the outer left and outer right part of the structure constituting a higher level factor should be drawn)
            new_level_factor_list_order, new_constitution_relation_list, color_map, new_level_equiv_list = arrange_factors(level_factor_list, new_level_equiv_list, constitution_relation_list, mode)
            
            # check for uncategorised causal factors
            aux_factor_list = [] # list of factors that appear in new_level_equiv
            for lvl in new_level_equiv_list:
                for formula in lvl:                    
                    aux_factor_list.extend(get_components_from_formula(formula[0], level_factor_list))
                    aux_factor_list.append(formula[1])
            
            aux_factor_list = list(set(aux_factor_list))
            
            # are all factors of aux_factor_list in some level and order of new_level_factor_list_order?
            found_all = True
            for factor in aux_factor_list:
                found = False
                for lvl in new_level_factor_list_order:
                    for order in lvl:
                        for fac in order:
                            if factor == fac:
                                found = True
                                break
                        if found: break
                    if found: break
                    
                if not(found):
                    found_all = False
                    total_solutions = total_solutions - 1 # do not count this solution
                    break
                    
            # continue only with solutions where all factors are categorised
            if found_all:
                # check whether the intra-level causal ordering conforms the input on the causal ordering
                
                # translate the order information from new_level_factor_list_order into a nested list (by levels) of
                # binary order relations (fac[i], fac[j]) iff fac[i] is of same level but lower order than fac[j]
                order_relations = []
                for lvl in range(len(new_level_factor_list_order)):
                    order_relations.append([]) # append a new sublist per constitutive level
                    for id_order in range(len(new_level_factor_list_order[lvl])):
                        for fac in new_level_factor_list_order[lvl][id_order]:
                            for id_order_2 in range(id_order + 1,len(new_level_factor_list_order[lvl])):
                                # run over all higher orders
                                for fac_2 in new_level_factor_list_order[lvl][id_order_2]:
                                    new_pair = (fac, fac_2)
                                    order_relations[lvl].append(new_pair)
                
                order_preserved = True
                
                for lvl in range(len(order_relations)):
                    # check whether there are any re-specified causal order information
                    if lvl <= len(order_input_information) and order_input_information:
                        # the input information is consistent with the obtained causal order of the particular solution
                        # iff the set of the input order relations is a subset of the obtained order relations (per level)
                        if not(set(order_input_information[lvl]).issubset(order_relations[lvl])):
                            order_preserved = False
                            total_solutions = total_solutions - 1 # do not count this solution
                            break
                
                
                # limiting the output to 1000 solutions
                # remove the respective condition in the if-clause if more solutions are required
                if order_preserved and sol_counter < 1000:
                    #############################################
                    # step 6: graphical output as a graph in pdf                         
                    # generating the tex-code 
                    sol_counter = sol_counter + 1
                    try:
                        st = print_structure_in_tikz_plot(new_level_factor_list_order, new_level_equiv_list, new_constitution_relation_list, color_map)
                    except:
                        st = ""
                            
                    # subtitle of the graph will be the formula in tex-math syntax
                    subtitle = convert_formula_to_tex_code(new_level_equiv_list)
                    subtitle = "\\tiny " + subtitle # formulae might be quite long, so subtitle should be written in tiny
                        
                    # counter enumerates the solutions
                    entry = (sol_counter, st, subtitle)
                    tex_table.append(entry)
            
        # after one entry for each solution has been generated in tex_table compile the pdf
        if tex_table:
            # if tex_table is non-empty
            if os.path.exists(latex_template_file) :
                create_pdf(tex_table, latex_template_file, total_solutions)
                print("Number of solutions " + str(total_solutions))        
                if total_solutions > 1000:
                    print("More than 1000 solutions obtained, plotting only the first 1000.")
            else :
                abort = True
                print("Error: Cannot plot the graphs since the expected template file " + latex_template_file + " does not exist in path folder.")
            if "fulllist" in mode:
                # when in fulllist mode, add the formulae to the list to be exported
                for sol in solution_term_list:
                    separate_formula_list.append("$" + convert_formula_to_tex_code(sol) + "$")       
                # create file
                create_separtate_formula_list(separate_formula_list)
             
            # print the duration of the total run time into the console
            print("Total run time " + str(round(time.time() - start_time,2)) + " seconds.")  
            
        else:
            # solution_list is empty
            # no solution survived selection of valid solutions
            print("No valid complex solution formula has been found in " + input_file + ".") 
                    
if __name__ == '__main__':
    # start main() function when the py-file is executed
    main()
