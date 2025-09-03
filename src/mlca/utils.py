#!/usr/bin/env python3

# file: utils.py

"""utility functions, mostly string and list comparisons
"""

import re                          # regex for complex search patterns in strings
import itertools                   # itertools provides functions to obtain all permutations of a string and Cartesian products of lists

def powerset(in_set: set) -> set:
    """Returns the powerset of the input in_set.

    Parameters
    __________
    in_set : set

    Returns
    _______
    set
        powerset of in_set = the set of all sets that can be formed of in_set and its elements
    """
    aux_list = list(in_set)
    sec_aux_list = []
    for i in range(2**len(aux_list)):
        sub_list = []
        for j in range(len(aux_list)):
            if i & 2**j:
                sub_list.append(aux_list[j])
                
        sec_aux_list.append(sub_list)
    return sec_aux_list

def list_comparison(list1: list, list2: list) -> bool:
    """Compares two nested lists. Returns True iff the sorted lists with sorted sublists are equal.

    Parameters
    __________
    list1 : list of lists
    list2 : list of lists

    Returns
    _______
    bool
        True if both lists are equivalent, else false

    """
    for subl in list1:
        subl.sort()
    for subl in list2:
        subl.sort()
    list1.sort()
    list2.sort()
    return list1 == list2

def contains_term(original_term: str, comparison_term: str) -> bool:
    """Checks whether comparison_term contains all substrings of original_term,
    possibly not in one piece, but spaced out over original_term.

    Parameters
    __________
    original_term: str
        string which is searched for
    comparison_term: str
        string in which is searched for comparison_term

    Returns
    _______
    bool
        truth value of whether original_term contains comparison_term
    """

    #Contains may mean that e.g. "A*C" is contained in "A*B*C".
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

def flatten_nested_list(in_list: list) -> list:
    """Flattens an homogenous list up to two times in case that it is a nested list.

    Parameters
    __________
    in_list : list
        list that will be transformed
        It is assumed that in_list is homogenous insofar that if the first element
        is a list, all further elements are lists, too, and if the first element
        is not a list, the further elements are neither.

    Returns
    _______
    list
        flattened list, up to two nested levels are removed and their elements
        are direct elements of the returned list
        empty list in case that in_list is not a list
    """

    if type(in_list) == list:
        if not in_list:
            return in_list
        else:
            if type(in_list[0]) == list:
                flat_list = [x for sub_list in in_list for x in sub_list]
                if type(flat_list[0]) == list:
                    flat_list = [x for sub_list in flat_list for x in sub_list]
                return flat_list
            else:
                return [x for x in in_list]
    else:
        # in_list is not a list
        return []
    
def find_causal_factors(st: str) -> list:
    """Returns all substrings of st that are separated by ", " or " < ".

    Parameters
    __________
    st : str

    Returns
    _______
    list
        (possibly empty) list of substrings of st separated by ", " or " < "
    """
    
    # deletes end-of-line-symbol ("\n") and spaces at the end of line if necessary
    st = re.sub("\r?\n","",st).rstrip()
    
    # returns the list of components of st that were separated by ", " or " < "
    return re.split(r',\s*|\s*<\s*', st)
    

def get_causal_prefactors(factor: str, formula_list: list, factor_list: list) -> list :
    """Returns a list of prefactors to a given factor.

    Prefactors are defined recursively. A prefactor of first order is a factor that
    appears in the first element of a tuple whose second element is equal to the
    given factor. Factors that appear in the first element of tuples whose second
    element is a prefactor of order n are prefactors of order n+1.

    Parameters
    __________
    factor : str
        factor for which the list of its prefactors is created
    formula_list : list of 2-tuples of str
        tuples represent causal relations: first element represents the cause, second element the effect
    factor_list: list of str
        list considered factors

    Returns
    _______
    list of str
        list of prefactors to the given factor
    """
    
    return_list = []
    
    
    for formula in formula_list:
        if formula[1] == factor:
            # if the factor we are interested in is the target factor of the relation formula,            
            # then add all factors that appear on the left side of formula to the return list
            return_list.extend(get_components_from_formula(formula[0], factor_list))
            
            for pfac in get_components_from_formula(formula[0], factor_list) :
                # get the indirect prefactors recursively
                return_list.extend(get_causal_prefactors(pfac, formula_list, factor_list))
    
    return_list = list(set(return_list))  # get rid of duplicates
    return return_list

def find_effects(formula: list, factor_list: list) -> list:
    """Determines which elements from factor_list are dependent variables (effects)
    given the dependencies expressed in formula.

    Causal factors are effects only if they do not satisfy either of three conditions:

    1) the causal factor is one in every line of the configuration table (in this case they are irrelevant)
    2) it is zero in every line of the configuration table (same as 1)
    3) two lines in the table differ only by the value of the factor (this means they might only be first causes)

    These conditions follow M. Baumgartner (2009) "Uncovering deterministic causal structures: a Boolean
    approach", p. 83.

    Returns the list of effects.

    Parameters
    __________
    formula: list of str
        list of string, it is assumed that each element is a conjunctive formula of
        the factors appearing from factor_list
    factor_list: list of str
        list of the potential variables of formula

    Returns
    _______
    list of str
        list of elements from factor_list that do not satisfy any of the conditions 1)-3)
    """

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
                        cond = True
                        break
        if cond:
            # delete the causal factor if either of the three exclusion criteria is true
            #print(effect_list[i] + " discarded. It has no causal relevance for any other causal factor.")
            del effect_list[i]

    return effect_list

def list_to_string(in_list: list) -> str:
    """Converts a nested list of the form in_list[DISJUNCT][CONJUNCT] or
    a list simple list of the form in_list[CONJUNCT] into a string
    which connects disjuncts by " + " and conjuncts by "*"

    Parameters
    __________
    in_list: list of lists of str or list of str

    Returns
    _______
    str
        String elements of list of lowest order are connected by '*',
        sublists by ' + ', if in_list is empty, return ''
    """

    if in_list: # list is not empty
        if type(in_list[0]) == list: # list is nested
            return ' + '.join(['*'.join(disj) for disj in in_list])
        else: # list of strings
            return "*".join(in_list)
    else:
        return ""


def string_to_list(st: str) -> list:
    """Converts a string into nested list:
    ' + ' separates sublists, '*' elements of the sublists

    Parameters
    __________
    st: str
        String, expected to express a DNF with disjunctor ' + ' and
        conjunctor '*'

    Returns
    _______
    list of lists of str
        nested list of form out_list[DISJUNCT][CONJUNCT]
    """
    out_list = []
    disj_list = re.split(r'\s*\+\s*', st)
    for disj in disj_list:
        out_list.extend([re.split(r'\*', disj)])

    return out_list
    
def get_equiv_formula(st: str) -> tuple:
    """Transforms a string into a tuple of strings (a,b) with the following characteristics:

    * The substring " <-> " within the input st separates a from b. If st does not contain " <-> ", b=''. If
      st contains multiple " <-> ", b includes everything between the the first two instances of the separator.
    * Any white spaces or substrings followed by white spaces in a are removed.
    * If st contains any lower letter, they are replaced by the '~' appended by the corresponding upper letter.
    * Any trailing white spaces in b or substrings following white spaces in b are removed.

    This procedure transforms output from cna into tuples (a, b), such that a + " <-> " + b
    is a well-formed equivalence relation with a being a disjunctive normal form and b an atomic expression.

    Parameters
    __________
    st : str

    Returns
    _______
    tuple of str
    """

    a = re.split(" <-> ",st)[0].strip()          # strip() removes leading spaces
    # in case that the line starts with some unnecessary stuff, followed by spaces, capture only content
    # behind white space
    if bool(re.search(r'\s{2,}', a)):
        a = re.split(r'\s{2,}', a)[1]        
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
    b = re.split(r'[ \t]',b)[0]
    return (a,b)

def get_components_from_formula(st: str, factor_list: list) -> list:
    """Returns a list of the elements of factor_list that appear in the input string st.
    The returned list is empty if no factor from factor_list appears in st or factor_list is empty,
    not a list of strings or not a list at all.

    Parameters
    __________
    st : str
        string which is searched for elements from factor_list
    factor_list : list of lists of lists of str or list of lists of str or list of str
        list whose elements are searched for in st

    Returns
    _______
    list of str
        elements of factor_list that have been found in st or empty list
    """
    
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

def get_factor_level(factor: str, level_factor_list: list) -> int:
    """Returns the index of the sublist of the nested list level_factor_list which contains factor.
    If factor is element in several sublists, the index of the first sublist is returned.
    If the list is empty or the factor has not been found in any sublist return -1.

    Parameters
    __________
    factor : str
        string which is searched for in the sublists of level_factor_list
    level_factor_list : list of lists of str or list of lists of lists of str
        nested list whose lowest level sublists are searched for factor

    Returns
    _______
    int
        index of sublist which contains factor, -1 if factor is not contained in
        any sublist or list is empty
    """
    
    # level_factor_list has either the form level_factor_list[LEVEL][FACTOR]
    # or level_factor_list[LEVEL][ORDER][FACTOR]
    level: int = -1
    
    if not(level_factor_list):
        # level factor list is empty
        return level
    else:
        # check whether level_factor_list is further separated by causal orders
        multi_order: bool = isinstance(level_factor_list[0], list)
        
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

def get_formula_level(st: str, level_factor_list: list) -> int:
    """Searches in string st for elements of sublists of level_factor_list.
    If all elements found are elements of the same sublist of level_factor_list,
    return the index of this sublist, elsewise, return -1.


    Parameters
    __________
    st : str
        string which is tested for containing only elements of the same sublist
    level_factor_list : list of lists of str or list of lists of lists of str
        nested list whose sublists are went through looking for substrings of str

    Returns
    _______
    int
        index of the sublist to which all elements from
    """

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

def get_factor_order(factor: str, factor_list: list) -> int:
    """Searches in sublists of second order of the nested list factor_list
    for factor. Returns the index of the first order sublist in whose sublist
    factor has been found. If factor is contained in several sublists,
    the lowest index is returned. If factor is not an element in any sublist, return -1.


    Parameters
    __________
    factor : str
        string which is searched for
    factor_list : list of lists of lists of str
        nested list whose sublists are went through looking for substrings of str

    Returns
    _______
    int
        index of the sublist of factor_list whose sublist contains factor
    """
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
    
def get_formula_order(formula: str, factor_list: list) -> int:
    """Applies get_components_from_formula on formula and searches for all obtained
    substrings in sublists of factor_list. If all substrings have been found, returns
    the highest index of a sublist that contains a substring, elsewise, returns -1.

    Parameters
    __________
    formula : str
        string whose substrings are searched for in sublists of factor_list
    factor_list : list of lists of lists of str
        nested list whose sublists are went through looking for substrings of str

    Returns
    _______
    int
        highest index of sublist that contains substring of formula
    """
    # specific use in mLCA
    # determines the order of a formula = the highest order of the factors occurring in it
    # The order of the factors is determined by factor_list using the function get_factor_order.
    # if all factors have a determinable order, the maximum value is returned,
    # otherwise -1
    
    order = -1
       
    
    for fac in get_components_from_formula(formula, factor_list):
        if order < get_factor_order(fac, factor_list) :
            order = get_factor_order(fac, factor_list)
                
    return order
    
    
def get_ordered_dnf_string(formula: str) -> str:
    """Alphabetically orders all disjuncts and conjuncts in a string, assuming
    ' + ' as disjunctor and '*' as conjunctor. Returns the string with the
    ordered formula. If formula is not a string, return ''.

    Parameters
    __________
    formula : str
        string expected to contain a disjunctive normal form with disjunctor ' + '
        or conjunctor '*'

    Returns
    _______
    str
        string with alphabetically ordered formula
    """
    # applies an alphabetic order to a disjunctive normal form formula by applying the commutation rules, such that
    # (a) all conjuncts in every conjunction are sorted alphabetically
    # (b) all disjuncts in every disjunction are sorted alphabetically 
    # e.g. A < A + B + D < A + C + D < D
    # returns the transformed formula as string
    
    if not(type(formula) == str):
        # if formula is not of type string, treat it like an empty string
        formula = ""
                                         
    # split formula into disjuncts
    disj_list = re.split(r'\s\+\s', formula)
                                        
    # split disjuncts into conjuncts
    conj_list = [re.split(r'\*', disj) for disj in disj_list] # nested list conj_list[DISJUNCT][CONJUNCT IN DISJUNCT]
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


def get_ordered_dnf_list(formula: list) -> list:
    """Alphabetically orders all disjuncts and conjuncts in a list, assuming
    a nested list whose elements on the first level are disjuncts and elements
    on the second level are conjuncts.

    Parameters
    __________
    formula : list of list of str
         nested list of string

    Returns
    _______
    list of list of str
        alphabetically ordered list with alphabetically ordered sublists
    """

    if not(type(formula) == list):
        # if formula is not of a list, return the empty list
        return []
    elif len(formula) < 1: # formula is empty
        return []
    elif any(type(element) != list for element in formula):
        formula.sort()
        return formula
    else:
        for sublist in formula:
            # sort each disjunct
            sublist.sort()
        # sort the whole disjunction
        formula.sort()
        return formula

def get_clusters(formula_list: list, factor_list: list) -> list:
    """Determines clusters of elements from factor_list that are connected via formulae
    from formula_list. Two elements are connected if they appear in the same formula,
    or if they both appear in formula within chain of formulae in which every formula
    shares at least one element with any other formula of that chain.
    Returns the list of clusters as a nested list.

    Parameters
    __________
    formula_list : list of 2-tuples of str
        list of 2-tuples of strings that a checked for connecting elements from
        factor_list
    factor_list : list of str or list of lists of str or list of lists of lists of str
        possibly nested lists of strings

    Returns
    _______
    list of lists of str
        list of clusters of factors, each cluster contains all factors that are mutually
        connected, returns a list containing the empty list if formula_list is empty or
        not a list
    """
    # determines clusters of causally connected factors
    # two factors are causally connected, iff they either appear in the same formula of formula_list or they are connected via intermediate factors
    # returns the list of clusters, each cluster is a list of the contained factors
    

    if not factor_list:
        return [[]] # result if factor_list is empty
    elif not(type(factor_list) == list):
        return [[]] # result if factor_list is no list
    else:
        # flatten factor_list in case that the list is nested
        factor_list = flatten_nested_list(factor_list)
        if type(factor_list[0]) == list:
                factor_list = [x for sub_list in factor_list for x in sub_list]
                if type(factor_list[0]) == list:
                    factor_list = [x for sub_list in factor_list for x in sub_list]
    
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

def get_truthvalue(formula: list, assignment: dict) -> bool:
    """Evaluates the truth value of a logical formula in DNF based on
    a specified assignment of truth values to its variables.
    If the formula is syntactically ill-formed or assignment does not
    assign a truth value to all occuring variables, the function returns
    False, otherwise it returns the truth value of the function.

    Parameters
    __________
    formula : list of str or list of lists of str
        Possibly nested list of string. Each string is expected to
        contain one factor or a negated factor. Elements of sublists
        are connected via conjunctors, the sublists via disjunctors.
    assignment : dictionary of str -> bool
        Dictionary that assigns Boolean values to string entries.

    Returns
    _______
    bool
        True if formula is well-formed, a truth value is assigned to
        all variables and formula is True under this assignment.
        False if (i) formula is ill-formed, (ii) assignment does not
        contain truth value assignments to all variables of formula, or
        (iii) formula is False under the given assignment.
    """
    # recursive function to return the truth value of formula under the given assignment of truth values to variables
    # returns False in case of syntax error in formula or if not all variables are assigned to a truth value
    # assignment is expected to be a dictionary whose entries encompass the variables occuring in formula with Boolean values
    # assumptions for formula corresponding to mLCA-style: 
    # (1) if formula is a list with more than one element, each element is interpreted as a conjunct of a chain of conjunctions
    # (2) if formula is a 1-elemental list whose element is a 2-tuple or formula is a 2-tuple, both tuple elements are interpreted as connected by an equivalence
    # (2a) left side term is a disjunctive normal form string variable
    # (2b) right side term is atomic
    # (3) if formula is a string, it is a disjunctive normal form with disjunctor ' + ',
    #     conjunctor '*' and negator '~'
    if type(assignment) == dict:
        if type(formula) == list and len(formula) > 1:
            # case (1) -> conjunction
            truthvalue = get_truthvalue(formula[0], assignment) and get_truthvalue(formula[1:], assignment)
        elif type(formula) == list and len(formula) == 1:
            # case (2) -- list with only one element, assumed to be a tuple
            if type(formula[0]) == tuple and len(formula[0]) == 2:
                # -> logical equivalence
                truthvalue = get_truthvalue(formula[0][0], assignment) == get_truthvalue(formula[0][1], assignment)
            else:
                # type error
                print('syntax/type error')
                truthvalue = False
        elif type(formula) == tuple and len(formula) == 2:
            # case (2) -- 2-tuple -> logical equivalence
            truthvalue = get_truthvalue(formula[0], assignment) == get_truthvalue(formula[1], assignment)
        elif type(formula) == str and formula in assignment:
            # case (2a) - atomic term
            if type(assignment[formula]) == bool:
                truthvalue = assignment[formula]
            else:
                print('type error in dictionary')
                truthvalue = False

        elif type(formula) == str and len(formula) > 1:
            # case (2b) - complex string
            if formula.find(' + ') > -1:
                # split disjunction
                terms = formula.split(' + ', 1)
                truthvalue = get_truthvalue(terms[0], assignment) or get_truthvalue(terms[1], assignment)
            elif formula.find('*') > -1:
                # split conjunction
                terms = formula.split('*', 1)
                truthvalue = get_truthvalue(terms[0], assignment) and get_truthvalue(terms[1], assignment)
            elif formula.find('~') == 0:
                # negator (as main operator) can only be in the first position
                truthvalue = not(get_truthvalue(formula[1:], assignment))
            else:
                # syntax error
                truthvalue = False
        else:
            # syntax/type error
            truthvalue = False
        
    else:
        print('type error: no dict')
        truthvalue = False
    return truthvalue

def create_assignments(variable_list: list) -> list:
    """Creates a list of all possible truth value assignments for variable_list.
    Truth value assignments are stored in dictionaries whose keys are the
    variable names and whose values their truth values.

    Parameters
    __________
    variable_list : list of str
        list of variables

    Returns
    _______
    list of dict
        list of value assignments in form of dictionaries whose keys are the
        variable names and whose values their truth values
    """
    # generate all cases for Boolean variables
    assignment_list = []
    for bool_value in itertools.product([True, False], repeat = len(variable_list)):
        assignment_list.append(dict(zip(variable_list, bool_value)))
    return assignment_list
    
def count_true(function_list: list, factor_list: list) -> int:
    """Counts the number of truth value assignments to the variables in factor_list that make the logical
    formula true. The formula is constructed by connecting both elements of each tuple in function_list
    with an equivalence operator and combining all such equivalence expressions using conjunctors.

    Parameters
    __________
    function_list : list of 2-tuples of str
        List of 2-tuples of which the first element is expected to be a DNF using the syntactical rules:
        * '*' represents a conjunctor
        * '~' represents a negation
        * ' + ' represents a disjunctor

        The second element of each tuple is expected to be an atomic term.

    factor_list : list of lists of lists of str or list of lists of str or list of str
        (nested) list of factors that are possibly variables of the logical function from function_list

    Returns
    _______
    int
        Amount of True values in the truth table for the given logical function function_list.
        Returns 0, in case that function_list contains syntax errors, or if function_list
        or factor_list is not a list.
    """
    # input: logical formula given as list of tuples
    # counts for how many assignments of the variables occuring in function_list the function becomes true
    # returns: zero if syntax/type errors in function
    #          otherwise the number of possible assignments such that function_list becomes true
        
    counter = 0 # counts in how many cases function_list is True
    
    # determine variables occuring in formula
    st = ''
    for term in function_list:
        st = st + '(' + term[0] + '<->' + term[1] + ')*'
    st = st[:-1]
    
    variable_list = get_components_from_formula(st, factor_list)
    
    # generate all cases for Boolean variables
    assignment_list = create_assignments(variable_list)
    
    for assignment in assignment_list:
        if get_truthvalue(function_list, assignment):
            counter = counter + 1
    return counter

def get_coextensive_factors(factor_list: list, formula: list, respect_levels: bool = False) -> list:
    """Determines clusters of coextensive factors by invoking the function determine_coextensive_clusters.
    If respect_levels is set to True, it will be executed for each level separately, and factors from
    different level will not be put into the same cluster of coextensive factors.

    Parameters
    __________
    factor_list: list of str or list of lists of lists of str
        possibly nested list of factors, it is assumed that the list is nested if
        respect_levels is set to True
    formula: list of list of str
        nested list representing a DNF, first level list contains disjuncts,
        second level the conjuncts of each disjunct, it is assumed that
        the string elements are either the factors from factor_list or their negations
    respect_levels: bool, optional
        determines whether coextensive factors of different levels are still clustered
        together (respect_levels=False) or not (respect_levels=True)

    Returns
    _______
    list of lists of str
        a nested list with one sublist per cluster of coextensive factors, the sublists
        contain all factors that are coextensive with each other
    """

    if respect_levels and factor_list and isinstance(factor_list[0], list) and \
       isinstance(factor_list[0][0], list):
        level_list = []
        for index, lvl in enumerate(factor_list):
            level_list.append([])
            level_list[index] = determine_coextensive_clusters(flatten_nested_list(lvl), formula)
        return level_list
    elif factor_list:
        return determine_coextensive_clusters(flatten_nested_list(factor_list), formula)
    else:
        return []

def determine_coextensive_clusters(factor_list: list, formula: list) -> list:
    """Determines clusters of coextensive factors in factor_list from a DNF-formula.
    Returns a nested list of clusters of coextensive factors.

    Parameters
    __________
    factor_list: list of str
        list of factors
    formula: list of list of str
        nested list representing a DNF, first level list contains disjuncts,
        second level the conjuncts of each disjunct, it is assumed that
        the string elements are either the factors from factor_list or their negations

    Returns
    _______
    list of lists of str
        a nested list with one sublist per cluster of coextensive factors, the sublists
        contain all factors that are coextensive with each other
    """
    list_of_coextensives = [] # this becomes a nested list: every sublist contains factors that are mutually coextensive

    for i in range(len(factor_list)-1):
        for j in range(i+1,len(factor_list)):
            co_ext = True
            for disj in formula:
                neg_i = "~" + factor_list[i]
                neg_j = "~" + factor_list[j]
                if (((factor_list[i] in disj) and not(factor_list[j] in disj)) or \
                   ((factor_list[j] in disj) and not(factor_list[i] in disj)) or \
                   ((neg_i in disj) and not(neg_j in disj)) or ((neg_j in disj) and not(neg_i in disj))):
                    co_ext = False # two factors are not coextensive if one or its negation appears in one disjunct but the other
                    break          # does not

            if co_ext:
                if not(list_of_coextensives): # if list of coextensives is still empty
                    list_of_coextensives.append([]) # append an empty sublist
                    list_of_coextensives[0].append(factor_list[i])
                    list_of_coextensives[0].append(factor_list[j])
                else: # list is non-empty search for a sublist that contains factor_list[i]
                    new_list = True
                    for sublist in list_of_coextensives:
                        if ((factor_list[i] in sublist) or (factor_list[j] in sublist)):
                            # at least one of both factors is already contained in one sublist of coextensive factors
                            new_list = False
                            if not(factor_list[i] in sublist): # and factor_list[i] is not,
                                sublist.append(factor_list[i]) # then add factor_list[i] to sublist
                            elif not(factor_list[j] in sublist): # other case factor_list[j] is not contained,
                                sublist.append(factor_list[j]) # then add it to sublist
                            break
                    if new_list: # neither factor is already contained in any sublist
                        list_of_coextensives.append([]) # create a new sublist
                        list_of_coextensives[-1].append(factor_list[i]) # append factor_list[i]
                        list_of_coextensives[-1].append(factor_list[j]) # and factor_list[j]
    return list_of_coextensives
