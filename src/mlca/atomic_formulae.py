#!/usr/bin/env python3
#
# file: atomic_formulae.py
"""
Functions to derive atomic solution formulae from Boolean coincidence data tables
"""
import copy               # for deep-copy of lists
import re                 # regex for complex search patterns in strings
import pandas as pd       # for reading csv files that contain truth tables
import itertools          # itertools provides functions to obtain all permutations of a string and Cartesian products of lists
import multiprocessing    # multiprocessing and functools for multicore usage
import suspension_search as ss
from utils import get_components_from_formula, get_factor_level, get_factor_order, get_equiv_formula, list_to_string, \
                  string_to_list, contains_term, flatten_nested_list, find_effects, get_coextensive_factors

def get_instance_formula_to_factor(in_formula: list, factor: str, level_factor_list_order: list) -> dict:
    """Derives the instance function for factor from the formula in_formula.

    The instance function of a factor is obtained by removing that factor or its negation from each of the
    conjuncts in in_formula and setting its value to True or False respectively.
    For computational effiency, formulae that mix different levels and do not satisfy the conditions
    on constitution relations are reduced such that factors from other levels than the level of factor
    are discarded.

    Returns the instance formula as a dictionary in which the formulae (in form of lists) are the keys,
    and the respective truth values the corresponding dicitonary values.

    Parameters
    __________
    in_formula: list of lists of str
        nested list of strings, expected to have the form [DISJUNCTS][CONJUNCTS]
        each conjunct is either a factor from level_factor_list_order or a negated
        factor (='~' + factor)
    factor: str
        factor whose instance function is to be generated
    level_factor_list_order: list of lists of lists of str
        nested list of factors, form [LEVEL][CAUSAL_ORDER][FACTOR]
        expected to contain all

    Returns
    _______
    dict
        dictionary representing the instance function for factor, keys are formulae in form of lists,
        the values are Boolean
    """

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

def reduce_term_by(term: str, literal: str) -> str:
    """Removes the string literal from the string term. If term starts
    with literal, the string literal + '*' is removed from term, otherwise
    the string '*' + literal.

    It is used to remove literal from a chain of conjuncts in the string term.

    Parameters
    __________
    term: str
        string to be reduced, expected to express a conjunctive
        formula with conjunctor '*'
    literal: str
        string which should be removed from term

    Returns
    _______
    str
        term reduced by literal + '*' or '*' + literal
    """

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

def absorb_terms(arg: tuple) -> list:
    """Applies the absorption rule a*b*c*d*e + a*c*d <-> a*c*d to simplify a DNF which is
    passed to this function in form of a nested list.

    Returns the simplified DNF as nested list.

    Parameters
    __________
    arg: tuple
        arg[0] - is a nested list arg[0][LIST OF CONJUNCTS][CONJUNCTS]
        arg[1] - is a nested list, arg[1][LIST OF DISJUNCTS][LIST OF CONJUNCTS][CONJUNCTS],
                 it should include arg[0]

    Returns
    _______
    list of lists of str
        nested list of the form [DISJUNCTS][CONJUNCTS] corresponding to
        arg[1] after applying the absorption rule
    """

    # absorption rule: in the list of disjunctions, discard disjuncts that can be formed by conjuncting terms to other disjuncts of the list
    # so: if all elemts of a list of conjuncts are completely contained in another, discard the larger list
    return [[disj for disj in arg[1] if all([x in conj_2 for x in conj] for conj in disj)] for conj_2 in arg[0]]


def distribution(formula: str) -> str:
    """Applies the distribution rule on a logical formula given as a string as often as possible,
    then applies the absorption rule to simplify the resulting formula as often as possible.

    Returns the simplified formula as a string.

    Parameters
    __________
    formula: str
        expected to express a DNF with disjunctor ' + ' and conjunctor '*'

    Returns
    _______
    str
        simplified formula after applying distribution and absorption rules
    """

    # as long as formula has a conjunctor right before or after a bracket
    if (formula.find(")*") > -1 or formula.find("*(") > -1):

        formula = formula[:-1] # get rid of trailing ")"
        conj_list = re.split(r'\)\*', formula) # list of conjuncts of formula
        conj_list = [conj[1:] for conj in conj_list] # get rid of leading "("

        disj_list = [] # list of disjuncts per conjunct
        for conj in conj_list:
            disj_list.append(re.split(r'\s*\+\s*', conj))
        # disj_list is a list [[d11, d12, ...], [d21, d22, ... ],  ...] with dij being the j-th disjunct in conjunct i

        # rebuild formula
        formula = list_to_string([[*x] for x in itertools.product(*disj_list)])

        disj_list.clear()
        conj_list.clear()


        disj_list = re.split(r'\s*\+\s*', formula) # list of disjuncts of new formula

        #conj_list = [list(set(re.split(r'\*', disj))).sort() for disj in disj_list] # doesn't work

        conj_list = []
        set_disjuncts = set() # set in place of a list automatically discards duplicates
        for disj in disj_list:
            if not(disj in set_disjuncts):
                a = list(set(re.split(r'\*', disj)))
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

def get_prime_implicants(instance_formula: dict, factor: str, level_factor_list: list) -> list:
    """Applies the Quine-McCluskey algorithm to obtain prime implicants by comparing the
    min-terms of positive instance function with those of the negative instance function:

    - if a section of one positive min-term is not part of any negative min term,
      that positive min-term can be reduced to this section
    - if every section is a part of at least one negative min term, the considered term is a prime factor
      of the positive instance function

    Parameters
    __________
    instance_formula: dict
        dictionary with keys: conjunctive formulae in form of lists,
        values: truth value
    factor: str
        factor for whose instance function the prime implicants are derived
    level_factor_list: list of lists of str
        nested list of factors, expected to contain all variables in
        instance_function and factor

    Returns
    _______
    list of str
        list of prime implicants for the instance formula for factor
    """


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

def get_rdnf(pi_list: list, formula: dict, factor_list: list) -> list:
    """Obtains reduced disjunctive normal form using Petrick's algorithm.

    Parameters
    __________
    pi_list: list of str
        list of prime implicants
    formula: dict
        dictionary of strings with all values either True or False
        = keys are the min terms of the formula, True or False their
        truth values
    factor_list: list of str
        expected to express a DNF with disjunctor ' + ' and conjunctor '*'

    Returns
    _______
    list of str
        contains all solutions each item is the string of a reduced disjunctive normal form of formula
    """

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
        sol_list = re.split(r'\s*\+\s*',aux_formula)

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

def get_truth_table_from_file(file_path: str) -> tuple:
    """Reads a csv file into a pandas data frame with Boolean entries.

    Returns the list factors (=column heads) and a string expressing
    the logical formula that corresponds to the truth table, when
    all columns of a row are connected by conjunctors and the rows
    by disjunctors.

    Parameters
    __________
    file_path: str
        path to csv-file

    Returns
    _______
    list of str
        list of column heads in csv-file = list of factors
    str
        logical formula derived from csv-truth table
    """

    # different strings that will be interpreted as True
    true_values = [1, "1", "T", "t", "w", "W", "true", "True", True]
    # different strings that will be interpreted as False
    false_values = [0, "0", "F", "f", "false", "False", False]

    # read csv file into df, assume that the column separators are one character from sep
    df = pd.read_csv(file_path, sep='[:,;|_]', engine='python', index_col = False)

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

def create_factor_ordering(level_factor_order_list: list) -> list:
    # create level_factor_list
    output = []
    for lvl in range(len(level_factor_order_list)):
        output.append([])
        for order in level_factor_order_list[lvl]:
            output[lvl].extend(order)
    return output

def suspension_search_asf(level_factor_order_list: list, formula_st: str) -> list:
    """Determines the list of atomic solution formulae using a breadth first
    suspension tree search.

    Parameters
    __________
    level_factor_order_list: list of lists of lists of str
        nested list of factors by levels and by causal orders,
        form: level_factor_order_list[LEVEL][ORDER][FACTOR]
    formula_st: str
        logical formula from which the atomic solution formulae
        are to be derived,

        notation: "*" - conjunctor, "~" - negator,
        " + " - disjunctor

    Returns
    _______
    list of tuple (str, str)
        list of equivalence relations in form of 2-tuples
        element[0] - DNF formula; element[1] - atomic
    """

    data_table = []
    for index, line in enumerate(string_to_list(formula_st)):
        data_table.append({})
        for factor in flatten_nested_list(level_factor_order_list):
            if "~" + factor in line:
                data_table[index][factor] = False
            elif factor in line:
                data_table[index][factor] = True

    coextensive_factor_list = get_coextensive_factors(level_factor_order_list, string_to_list(formula_st), respect_levels=True)

    coextensive_list_ignore_levels = get_coextensive_factors(level_factor_order_list, string_to_list(formula_st), respect_levels=False)
    constitution_coextensive_list = [(x, y) for cluster in coextensive_list_ignore_levels for x in cluster for y in cluster \
        if get_factor_level(x, level_factor_order_list)+1==get_factor_level(y, level_factor_order_list)]
    # list of pairs (x,y) with x being an coextensive factor to y and level(x)+1==level(y)

    effects_list = find_effects(string_to_list(formula_st),flatten_nested_list(level_factor_order_list))

    nested_effects_list = [] # create a nested list
    for lvl_index, lvl in enumerate(level_factor_order_list):
        nested_effects_list.append([])
        nested_effects_list[lvl_index].append([])
        for fac in flatten_nested_list(lvl):
            if fac in effects_list:
                nested_effects_list[lvl_index][0].append(fac)

    causes_list = [] # list of first causes
    for lvl_index, lvl in enumerate(level_factor_order_list):
        if lvl_index > 0:
            remove_list = [x for level in level_factor_order_list for order in level for x in order if lvl != level ]
            reduced_data_table = ss.reduce_data_table(data_table, remove_list)
            new_formula = []
            for line_index, line in enumerate(string_to_list(formula_st)):
                new_formula.append([])
                for lit in line:
                    if any(fac == lit or "~" + fac == lit for fac in flatten_nested_list(level_factor_order_list[lvl_index])):
                        new_formula[line_index].append(lit)

            nested_effects_list[lvl_index][0] = find_effects(new_formula, flatten_nested_list(level_factor_order_list[lvl_index]))

        else:
            reduced_data_table = copy.deepcopy(data_table) # deepcopy makes also copies of the elements which are dictionaries
        stop = False
        counter = 0
        causes_list.append([])

        if len(flatten_nested_list(lvl)) == 1:
            nested_effects_list[lvl_index][0].append(flatten_nested_list(lvl)[0])
            stop = True
            causes_list[lvl_index].append([])

        while not(stop):
            causes_list[lvl_index].append([])

            for fac in flatten_nested_list(lvl):
                if not(fac in nested_effects_list[lvl_index][counter]):
                    causes_list[lvl_index][counter].append(fac)


            if len(causes_list[lvl_index][counter]) == 0 or (counter > 0 and ([x for x in causes_list[lvl_index][counter] if x not in causes_list[lvl_index][counter-1]] == [])):
                # if there are no first causes
                if counter == 0:
                    causes_list[lvl_index] = [flatten_nested_list(level_factor_order_list[lvl_index])]
                else:
                    del causes_list[lvl_index][-1] # remove last element
                    del causes_list[lvl_index][-1] # remove second-last element
                    del nested_effects_list[lvl_index][-1]
                    del nested_effects_list[lvl_index][-1]
                stop = True
            else:
                remove_list = [x for x in causes_list[lvl_index][counter]]
                for i in range(len(level_factor_order_list)):
                    if i != lvl_index:
                        remove_list.extend(nested_effects_list[i][0])
                reduced_data_table = ss.reduce_data_table(data_table, remove_list)
                nested_effects_list[lvl_index].append([])
                new_formula = []
                for line_index, line in enumerate(string_to_list(formula_st)):
                    new_formula.append([])
                    for lit in line:
                        if any(fac == lit or "~" + fac == lit for fac in nested_effects_list[lvl_index][counter]):
                            new_formula[line_index].append(lit)
                nested_effects_list[lvl_index][counter+1] = find_effects(new_formula, nested_effects_list[lvl_index][counter])

            counter += 1

    # remove all coextensive factors but one
    reduced_factor_list = copy.deepcopy(level_factor_order_list) # deepcopy makes also copies of the elements which are lists
    for lvl_index, lvl in enumerate(coextensive_factor_list):
        for cluster in lvl:
            for ind in range(1,len(cluster)):
                for i in range(len(nested_effects_list[lvl_index])):
                    if cluster[ind] in nested_effects_list[lvl_index][i]:
                        nested_effects_list[lvl_index][i].remove(cluster[ind])
                    for level in reduced_factor_list:
                        for order in level:
                            if cluster[ind] in order:
                                order.remove(cluster[ind])

    max_conj = len(flatten_nested_list(causes_list))
    max_disj = len(data_table)
    suspension_acc = 0.2

    equiv_relations = {}

    for lvl in range(len(nested_effects_list)):
        for i in range(len(nested_effects_list[lvl])):

            for target_factor in nested_effects_list[lvl][i]:
                # Create a root node to start with
                root_value = []
                root = ss.Node(root_value, level=-1)

                # adapt max_disj and max_conj for target_factor
                # max_disj should not be larger than the number of cases in which target_factor
                # is True
                true_cases = len([x for x in data_table if x[target_factor] == True])
                if true_cases < max_disj:
                    local_max_disj = true_cases
                else:
                    local_max_disj = max_disj

                #number_factors
                #if
                local_max_conj = max_conj

                #"""
                successful, found_nodes, last_find = ss.suspension_bfs(root, target_factor, reduced_factor_list, data_table, \
                    target_factor_level=get_factor_level(target_factor, level_factor_order_list), active=causes_list[lvl][i], \
                        max_disj=local_max_disj, max_conj=local_max_conj, max_depth=500, suspension_acc=suspension_acc)

                solutions_list = []
                if successful:
                    for node in found_nodes:
                        solutions_list.append(node.value)
                elif found_nodes:
                    for node in found_nodes:
                        solutions_list.append(node.value)

                else:
                    pass

                if not target_factor in equiv_relations:
                    # if it is the first run for target_factor
                    equiv_relations[target_factor] = solutions_list
                else:
                    # if there already exists an entry for target_factor
                    equiv_relations[target_factor].extend(solutions_list)



    # remove empty keys whose value is [] from equiv_relations
    remove_list = [key for key in equiv_relations if equiv_relations[key]==[]]
    for key in remove_list:
        del equiv_relations[key]

    # introduce relations for coextensive factors
    # 1) add DNF <-> B for all DNF <-> A if level(B)==level(A)
    # 2) add DNF <-> B for all DNF <-> A if ((level(B)+1==level(A) and level(DNF)==level(B)) OR level(B)==level(A)+1 and level(DNF)==level(A))
    # 3) go through every DNF and replace any combination of instances of "A" by "B" if level(A) == level(B)
    # 4) add B <-> A if level(A)==level(B)
    # 5) add B <-> A if level(A)==level(B)+1

    for lvl in coextensive_factor_list:
        for cluster in lvl:
            for fac_a in cluster:
                # case 1)
                if fac_a in equiv_relations:
                    for fac_b in cluster:
                        if get_factor_level(fac_a, level_factor_order_list) == get_factor_level(fac_b, level_factor_order_list):
                            equiv_relations[fac_b] = equiv_relations[fac_a]
                # case 4)
                for fac_b in cluster:
                    if get_factor_level(fac_a, level_factor_order_list) == get_factor_level(fac_b, level_factor_order_list) and \
                       fac_a != fac_b:
                        if not(fac_a in equiv_relations):
                            equiv_relations[fac_a] = []
                        if not [[fac_b]] in equiv_relations[fac_a]:
                            equiv_relations[fac_a].append([[fac_b]])


    # case 2)
    for pair in constitution_coextensive_list:
        if pair[0] in equiv_relations:
            for formula in equiv_relations[pair[0]]:
                if get_components_from_formula(list_to_string(formula), level_factor_order_list):
                    if get_factor_level(get_components_from_formula(list_to_string(formula), level_factor_order_list)[0], level_factor_order_list) == \
                       get_factor_level(pair[0], level_factor_order_list):
                        if not pair[1] in equiv_relations:
                            equiv_relations[pair[1]] = []
                        if not formula in equiv_relations[pair[1]] and not(list_to_string(formula) == pair[1]):
                            equiv_relations[pair[1]].append(formula)

        if pair[1] in equiv_relations:
            for formula in equiv_relations[pair[1]]:
                if get_components_from_formula(list_to_string(formula), level_factor_order_list):
                    if get_factor_level(get_components_from_formula(list_to_string(formula), level_factor_order_list)[0], level_factor_order_list) == \
                       get_factor_level(pair[0], level_factor_order_list):
                        if not pair[0] in equiv_relations:
                            equiv_relations[pair[0]] = []
                        if not formula in equiv_relations[pair[0]] and not(list_to_string(formula) == pair[0]):
                            equiv_relations[pair[0]].append(formula)
        # case 5)
        if not(pair[1] in equiv_relations):
            equiv_relations[pair[1]] = []
        if not([[pair[0]]]) in equiv_relations[pair[1]]:
            equiv_relations[pair[1]].append([[pair[0]]])


    # case 3)
    add_list = []
    for lvl in coextensive_factor_list:
        for cluster in lvl:
            for fac_a in cluster:
                for fac_b in cluster:
                    if fac_a != fac_b:
                        for key in equiv_relations:
                            if key != fac_a and key != fac_b:
                                for formula in equiv_relations[key]:
                                    # construct every combination of instances of fac_a replaced by fac_b in formula
                                    new_formula_list = ss.replace_instances_all_combs(list_to_string(formula), fac_a, fac_b)
                                    for new_formula in new_formula_list:
                                        if not(string_to_list(new_formula) in equiv_relations[key]):
                                            pair = (new_formula, key)
                                            if not pair in add_list:
                                                # add new combinations to add_list
                                                add_list.append(pair)

    output_list = ss.convert_dict_to_pair_list(equiv_relations)
    output_list.extend(add_list)

    return output_list


def read_data_from_csv(file_path: str, mode: str = "td") -> tuple:
    """Main function of atomic_formulae.py.

    (1) reads the csv file under file_path
        converts the truth table into a logical formula via get_truth_table_from_file
        determines the causal factors and any given information on constitution levels and causal order
    (2) derives instance formulae for all causal factors that might be effects from the obtained total formula
    (3) obtains the prime implicants for each instance formula
    (4) transforms the instance formulae into disjunctive normal forms by means of the prime implicants

    Parameters
    __________
    file_path: str
        path to csv-file
    mode: str, optional
        specifies the method of deriving the atomic solution formulae,
        If mode="bu", the bottom-up approach using the suspension tree search
        algorithm is used, otherwise the traditional Petrick algorithm (top down).

    Returns
    _______
    bool
        Boolean value for error tracking, False in case that any error occurred (invalid input), True otherwise
    list of lists of str
        a nested list of causal factors subdivided by constitution levels
    list of tuples of str
        list of pairs, each pair corresponds to an obtained equivalence formulae
        in form of: element[0] <-> element[1]
    list of lists of tuples
        nested list of pairs, for each constitution level, the order information are translated into
        binary relations between factors (fac[0], fac[1]) means: fac[0] < fac[1]
    """
    list_equiv_formula = []
    list_equiv_tuple = []
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

    if mode == "td":
        # top down approach using Petrick's algorithm

        if factor_list and formula_st != "":
            # transform the formula from string into a nested list
            # elements are the disjuncts of the min-term formula as lists of the conjuncts each disjunct
            # formula[DISJUNCT][CONJUNCT]
            formula = string_to_list(formula_st)

            # determine which causal factors might be effects
            effects_list = find_effects(string_to_list(formula_st),factor_list)

            # check for co-extensive factors - only for one factor of each set of co-extensive factors,
            # the prime implicants have to be determined
            list_of_coextensives = get_coextensive_factors(effects_list, formula)

            # keep only one factor per level of each set of coextensive factors
            # only keep the last element (will be factor of highest order within that level)
            if list_of_coextensives: # list is not empty
                for sublist in list_of_coextensives:
                    if len(sublist) > 1:
                        for level in level_factor_order_list:
                            found_one = False
                            for order in reversed(level):
                                for factor in order:
                                    if factor in sublist:
                                        if not(found_one):
                                            found_one = True
                                        else:
                                            effects_list.remove(factor)

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
                for effect in effects_list:
                    for sublist in list_of_coextensives:
                        if effect in sublist:
                            for fac in sublist:
                                if not(fac in effects_list) and not (fac == effect) and (get_factor_level(fac, level_factor_order_list) == get_factor_level(effect, level_factor_order_list)):
                                    # replace effect by fac and vice versa in the equivalence formulae for effect
                                    for formula in list_equiv_formula:
                                        lgth = len(effect)
                                        if formula.endswith(effect):
                                            # skip formula if cause-term contains factor of higher causal order than fac
                                            skip = False
                                            # case 1: fac is cause and effect of higher order than fac
                                            if (formula.find(fac) > -1) and (get_factor_order(fac, level_factor_order_list) < get_factor_order(effect, level_factor_order_list)):
                                                skip = True
                                            # case 2: another factor of higher order than fac is among causes
                                            if not(skip):
                                                for index, order in enumerate(level_factor_order_list[get_factor_level(fac, level_factor_order_list)]):
                                                    if index > get_factor_order(fac, level_factor_order_list):
                                                        for f in order:
                                                            if (f != effect) and (formula.find(f) > -1):
                                                                skip = True
                                                                break # break from for-loop over factors
                                                    if skip:
                                                        break # break from for-loop over orders

                                            if not(skip):
                                                st = formula[:-lgth].replace(fac,effect) # the new formula is the old one without the last
                                                st = st + fac # expression and with all occurences of fac replaced by effect
                                                # then append fac as second equivalent of the equivalence formula
                                                list_equiv_formula.append(st) # add the new formula to the formulae list
                            break # break from for-loop over sublists

            if list_equiv_formula:
                abort = False
                # transform the list of strings into a list of pairs, such that element[0] <-> element[1]
                list_equiv_tuple = [get_equiv_formula(formula) for formula in list_equiv_formula]
            else:
                # obtained list of equivalence formulae is empty
                abort = True

            level_factor_list = create_factor_ordering(level_factor_order_list)

        else:
            # factor_list is empty
            abort = True

        return abort, level_factor_list, list_equiv_tuple, order_input_information

    elif mode=="bu":
        # alternative mode: suspension tree search

        list_equiv_tuple = suspension_search_asf(level_factor_order_list, formula_st)
        abort = len(list_equiv_tuple) == 0
        level_factor_list = create_factor_ordering(level_factor_order_list)

        return abort, level_factor_list, list_equiv_tuple, order_input_information

if __name__ == '__main__':
    # when executing this script file
    print('Either run the command line application cli.py or the graphical user interface gui.py.')
