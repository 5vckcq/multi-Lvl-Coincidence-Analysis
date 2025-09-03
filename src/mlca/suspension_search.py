#!/usr/bin/env python3

# file: suspension_search.py

"""Alternative search algorithm for atomic solution formulae using
a suspension tree search.
"""

from collections import deque
import itertools                            # provides function for Cartesian product
import multiprocessing                      # multiprocessing and functools for multicore usage
from operator import itemgetter, attrgetter # provide efficient sorting functions

from utils import create_assignments, get_truthvalue, list_to_string, string_to_list, get_ordered_dnf_list, \
                  get_factor_level, flatten_nested_list, contains_term, get_components_from_formula, find_effects, \
                  get_coextensive_factors

class Node(object):
    """Node of a tree-structure

    Parameters
    __________
    value: list of lists of str
        logical formula in DNF represented as nested list
    name: str
        name of the node
    level: int
        constitutive level the node belongs to
    accuracy: float
        accuracy of how well value performs as equivalent to the target
    recall: float
        recall of how well value performs as equivalent to the target
    specificity: float
        specificity of how well value performs as equivalent to the target
    suspended: bool
        suspended nodes are postponed when transpassing the tree
    children: list of Node
        list of the node's children
    """

    def __init__(self, value: list, name: str = "root", level: int = 0, accuracy: float = 0, recall: float = 0, \
                 specificity: float = 0, suspended: bool = False, children: list = None) -> None:
        """Initialises a new tree.

        Parameters
        __________
        value: list of lists of str
            logical formula in DNF represented as nested list
        name: str, optional
            name of the node
        level: int, optional
            constitutive level the node belongs to
        accuracy: float, optional
            accuracy of how well value performs as equivalent to the target
        recall: float, optional
            recall of how well value performs as equivalent to the target
        specificity: float, optional
            specificity of how well value performs as equivalent to the target
        suspended: bool, optional
            suspended nodes are postponed when transpassing the tree
        children: list of Node, optional
            list of the node's children
        """
        self.value = value
        self.suspended = suspended
        self.name = name
        self.level = level
        self.accuracy = accuracy
        self.recall = recall
        self.specificity = specificity
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)

    def __str__(self) -> str:
        """String conversion of Node object returns the name parameter of the node.

        Returns
        _______
        str
            self.name
        """
        return self.name

    def add_child(self, child_node: 'Node') -> None:
        """Adds child_node as child to the currenct node.

        Parameters
        __________
        child_node: Node
            node that is to be added as child
        """
        assert isinstance(child_node, Node)
        self.children.append(child_node)

    def get_all_nodes(self) -> list:
        """Returns the list of all descendant nodes.

        Returns
        _______
        list of Node
            list of all descendant nodes
        """
        list_of_nodes = []
        list_of_nodes.extend(self.children)
        for child in self.children:
            list_of_nodes.extend(child.get_all_nodes())
        return list_of_nodes

    def create_new_nodes(self, active_nodes: list, data_table: dict, target_factor: str, created_nodes: list, \
                         suspended: bool = False, target_factor_level: int = 0, max_disj: int = 0, max_conj: int = 0) -> list:
        """Creates new nodes as children of the current node.

        Constructs all DNF that can be formed of the current node combined with any of the nodes from
        active_nodes as children of the current node.
        Suspends new children if the accuracy of new formulae is not better than those of the parent nodes or
        if the parameter suspended is set to True, or if the new node can be simplified by absorbtion law
        to a node with accuracy=1.
        Nodes to formulae which are no minimal DNF are not created (conjunctions must not contain the
        same literal more than once, disjunctions must not contain the same disjunct more than once).
        The children are sorted by descending accuracy and increasing character length to continue
        with the most promising nodes first.

        Parameters
        __________
        actice_nodes: list of Node
            list of all non-suspended nodes
        data_table: list of dict (str, bool)
            truth table in form of a list of dictionaries, each row corresponds to
            one list element, each element is dictionary with the same keys (the factors)
            and Boolean values
        target_factor: string
            name of factor for whose equivalent DNF is searched for
        created_nodes: list of Node
            list of all previously created nodes
        suspended: bool, optional
            if suspended is set to True, all newly created nodes will be suspended
        target_factor_level: int, optional
            constitutive level of target_factor
        max_disj: int, optional
            maximum number of disjunctions in each DNF-formula
            if max_disj=0, there is no upper limit
        max_conj: int, optional
            maximum number of conjunctions per disjunct in each DNF-formula
            if max_conj=0, there is no upper limit

        Returns
        _______
        list of Node
            list of all newly created nodes
        """
        out_list = []
        if (not isinstance(self, Node) and len(active_nodes) < 1) or get_accuracy(self.value, data_table, target_factor) == 1.0:
            # unexpected behaviour, at least an empty root element should exist
            # and new nodes should be addable from active_nodes
            # also skip nodes that are already equivalent to target_factor
            return []
        elif len(active_nodes) > 0 and self.value == []:
            # first level: filling first level of nodes below root
            literals = [y for x in flatten_nested_list(active_nodes) for y in [x, "~"+x] if x != target_factor] # adding negated factors to list
            list_to_add = [] # list of triples of values and their accuracy for target_factor and the factor's level'
            for lit in literals:
                if lit[0] == "~":
                    level = get_factor_level(lit[1:], active_nodes)
                else:
                    level = get_factor_level(lit, active_nodes)

                # only add factors of (1) same level as target_factor (=> causal relations)
                # or (2) one level below level of target_factor (=> constitution relations)
                if target_factor_level == level or (target_factor_level > 0 and target_factor_level == level + 1):
                    list_to_add.append(([[lit]], get_accuracy([[lit]], data_table, target_factor), level))

            # sort elements by accuracy for target_factor in descending order to start with the most promising literals
            list_to_add.sort(key=itemgetter(1), reverse=True)
            for value, acc, level in list_to_add:
                new_node = Node(value, name=str(value)[3:-3], level=level, accuracy=acc, recall=get_recall(value, data_table, target_factor),\
                                specificity=get_specificity(value, data_table, target_factor))
                if suspended:
                    new_node.suspended = True
                out_list.append(new_node)
                self.add_child(new_node)
            return out_list
        else:
            # higher level, connecting nodes by conjunctors or disjunctors
            new_nodes = [] # store triples of (new_node, accuracy, second_parent) (first parent is the current node)
            for anc_node in active_nodes:
                if self.level == anc_node.level:
                    # only complex terms of factors from the same level are meaningful
                    if not anc_node.value == self.value and len(anc_node.value) == 1 and anc_node.accuracy < 1.0:
                        #add possible conjunctions between anc_node and anc_node2
                        # if anc_node != self and anc_node is no disjunction
                        for disj in self.value:
                            if not(any(conj2 in disj for conj2 in anc_node.value[0]) or any("~" + conj2 in disj for conj2 in anc_node.value[0]) or \
                                   any("~" + conj in anc_node.value[0] for conj in disj)):
                                # no conjunct from anc_node2 is already part of the conjunction, neither any negation of factor that is in the conjunction
                                if max_conj == 0 or len(disj) + len(anc_node.value[0]) < max_conj + 1:
                                    # maximal conjunction length not surpassed

                                    # replace old disj-term by new disj + "*" + anc_node[0] and reorder disjuncts and conjuncts alphabetically
                                    self_value_string = list_to_string(self.value)
                                    disj_string = list_to_string(disj)
                                    new_string = list_to_string([[*disj, *anc_node.value[0]]])
                                    if self_value_string.find(disj_string) == 0:
                                        # disj is first disjunct
                                        new_value = get_ordered_dnf_list(string_to_list(new_string + self_value_string[len(disj_string):]))
                                    elif len(self_value_string.replace(disj_string, "")) == self_value_string.find(disj_string):
                                        # disj is last term
                                        new_value = get_ordered_dnf_list(string_to_list(self_value_string[:-len(disj_string)] + new_string))
                                    else:
                                        # disj is middle term
                                        new_value = get_ordered_dnf_list(string_to_list(self_value_string.replace("+ " + disj_string + " +", "+ " + new_string + " +")))
                                    new_nodes.append((new_value, get_accuracy(new_value, data_table, target_factor), anc_node))


                    # next step add further disjunctive terms:
                    if self.recall < 1.0 and anc_node.recall < 1.0 and \
                        (max_disj == 0 or len(self.value) + len(anc_node.value) < max_disj + 1):
                        if (len(anc_node.value) == 1 and not(anc_node.value[0] in self.value)):
                            # only add terms of one disjunct
                            # avoid appending disjuncts that are already part of current node
                            # (if one factor appears in several disjuncts, introduce first the other conjuncts)
                            new_value = get_ordered_dnf_list([*self.value, *anc_node.value])
                            new_nodes.append((new_value, get_accuracy(new_value, data_table, target_factor), anc_node))

            # sort new_nodes by descending accuracy to continue with the most promising elements first
            new_nodes.sort(key = lambda y: (y[1], -len(list_to_string(y[0]))), reverse=True) # second key guarantees that for same accuracy,
            # shorter expressions come first

            # add new nodes to out_list
            for value, acc, sec_parent in new_nodes:
                if not(any(node.value == value for node in created_nodes)):
                    new_node = Node(value, name=list_to_string(value), level=self.level, accuracy=acc,\
                                    recall=get_recall(value, data_table, target_factor), specificity=get_specificity(value, data_table, target_factor))
                    to_be_created = True
                    if suspended:
                        new_node.suspended = True
                    elif not((new_node.recall > self.recall and new_node.recall > sec_parent.recall) or \
                         (new_node.specificity > self.specificity and new_node.specificity > sec_parent.specificity)) or \
                         not(new_node.accuracy > self.accuracy or new_node.accuracy > sec_parent.accuracy):
                        # suspend if new node is not better than both current node and second parent
                        new_node.suspended = True
                    elif len(new_node.value) > 1 and any(any(len(node.value) == len(new_node.value) and node.accuracy == 1.0 \
                         and all(contains_term(list_to_string(disj2), list_to_string(disj1)) or any(disj2 == disj for disj in new_node.value) \
                         for disj2 in node.value) for node in created_nodes) for disj1 in new_node.value):
                        # remove disjunctions for which all disjuncts are either equal to or contain all disjuncts of an already found equivalent
                        to_be_created = False
                    elif any(node.accuracy == 1. and len(new_node.value) > len(node.value) and all(disj in new_node.value for disj in node.value) for node in created_nodes):
                        # enforce minimal necessity cf. Baumgartner (2009) "Uncovering Deterministic Causal Structures: A Boolean Approach" p. 4
                        to_be_created = False
                    elif any(node.accuracy == 1. and len(new_node.value) > 1 and len(node.value) == 1 and any(new_node.value == \
                             get_ordered_dnf_list([[*node.value[0], conj], [*node.value[0], "~" + conj]]) for disj in new_node.value for conj in disj) for node in created_nodes):
                        # new node has the form X*A + X*~A with X.accuracy=1
                        new_node.suspended = True
                    elif any(any(len(node2.value) == 1 and node2.accuracy == 1. and contains_term(list_to_string(node2.value), list_to_string(disj)) \
                         for node2 in created_nodes) for disj in value):
                        # enforce minimal sufficiency cf. Baumgartner (2009) "Uncovering Deterministic Causal Structures: A Boolean Approach" p. 4
                        to_be_created = False

                    if to_be_created:
                        out_list.append(new_node)
                        self.add_child(new_node)
                        created_nodes.append(new_node)

            return out_list

def get_accuracy(formula: list, data_table: list, target_factor: str) -> float:
    """Returns the accuracy of formula as equivalent to target_factor.
    Accuracy is defined as accuracy = TP + TN / (P + N) with the abbreviations
    number of true positives (TP), number of true negatives (TN), positives (P) and
    negatives (N).
    If formula is not a list or an empty list or P+N=0, the function returns 0.

    Parameters
    __________
    formula: list of lists of str
        nested list representing a DNF-formula
    data_table: list of dict (str, bool)
        truth table in form of a list of dictionaries, each row corresponds to
        one list element, each element is dictionary with the same keys (the factors)
        and Boolean values
    target_factor: str
        name of the factor whose values determine the true values (TP and TN)

    Returns
    _______
    float
        ratio TP + TN / (P + N), or zero if formula is not a list or empty,
        or P+N=0
    """
    if not(isinstance(formula, list)) or formula == []:
        return 0
    else:
        count_correct = 0
        count_wrong = 0
        for assignment in data_table:
            if get_truthvalue((list_to_string(formula), target_factor), assignment):
                count_correct += 1
            else:
                count_wrong += 1
        if count_correct + count_wrong == 0:
            count_wrong = 1 # avoid division by zero in accuracy formula
        return (count_correct/(count_correct + count_wrong))

def get_recall(formula, data_table, target_factor):
    """Returns the recall of how well formula functions as equivalent to target_factor.
    Recall is defined as recall = TP/P with the abbreviations
    number of true positives (TP), and number of positive values (P).
    If formula is not a list or an empty list or P=0, the function returns 0.

    Parameters
    __________
    formula: list of lists of str
        nested list representing a DNF-formula
    data_table: list of dict (str, bool)
        truth table in form of a list of dictionaries, each row corresponds to
        one list element, each element is dictionary with the same keys (the factors)
        and Boolean values
    target_factor: str
        name of the factor whose values determine the true values (TP)

    Returns
    _______
    float
        ratio TP / P, or 0 if formula is not a list or empty,
        or 1 if P=0
    """
    if not(isinstance(formula, list)) or formula == []:
        return 0
    else:
        P = 0
        TP = 0

        for assignment in data_table:
            if get_truthvalue(target_factor, assignment):
                P += 1
                if get_truthvalue(list_to_string(formula), assignment):
                    TP += 1
        if P == 0:
            # avoid division by zero
            # if there are no cases for target_factor being True,
            # recall = 1 if TP=P=0
            TP = 1
            P = 1
        return (TP/P)

def get_specificity(formula, data_table, target_factor):
    """Returns the specificity of how well formula as equivalent to target_factor.
    Specificity is defined as recall = TN/N with the abbreviations
    number of true negatives (TN), and number of negative values (N).
    If formula is not a list or an empty list or N=0, the function returns 0.

    Parameters
    __________
    formula: list of lists of str
        nested list representing a DNF-formula
    data_table: list of dict (str, bool)
        truth table in form of a list of dictionaries, each row corresponds to
        one list element, each element is dictionary with the same keys (the factors)
        and Boolean values
    target_factor: str
        name of the factor whose values determine the true values (TN)

    Returns
    _______
    float
        ratio TN / N, or 0 if formula is not a list or empty,
        or 1 if N=0
    """
    if not(isinstance(formula, list)) or formula == []:
        return 0
    else:
        N = 0
        TN = 0

        for assignment in data_table:
            if not(get_truthvalue(target_factor, assignment)):
                N += 1
                if not(get_truthvalue(list_to_string(formula), assignment)):
                    TN += 1
        if N == 0:
            # avoid division by zero
            TN = 1
            N = 1
        return (TN/N)

def replace_instances_all_combs(original_string: str, target: str, replacement: str) -> list:
    """Returns the list of strings which can be generated by parially replacing instances of the target
    string in original_string by replacement. For n instances of target, the resulting list will
    contain 2**n elements.

    Parameters
    __________
    original_string: str
        string in which instances of target will be partially replaced
    target: str
        string whose instances will be partially replaced
    replacement: str
        string which will replace target

    Returns
    _______
    list of str
        list of all combinations of replacements or not for every instance of target in original_string
    """

    # create a list of combinations where each instance of the target string is replaced with its replacement or kept as it is
    combinations = []
    index = 0

    while index < len(original_string):
        if original_string.startswith(target, index):
            # if target string is found, add its replacement and itself to options
            combinations.append([replacement, target])
            index += len(target)  # Move past the length of the target string
        else:
            # otherwise, add the current character as an option
            combinations.append([original_string[index]])
            index += 1

    # Generate all combinations of these options
    result = [''.join(combination) for combination in itertools.product(*combinations)]

    return result

def convert_dict_to_pair_list(original_dict: dict) -> list:
    """Converts a dictionary whose values are lists of lists into a list of pairs.
    Each pair has the form (content, key) with key being the respective key from
    the dictionary and content the string generated by list_to_string on each
    element of original_dict[key].

    Parameters
    __________
    original_dict: dict
        a dictionary whose values are lists

    Returns
    _______
    list of tuples
        list which contains the pairs (list_to_string(formula), key) for each key in original_dict
        and each formula in original_dict[key]
    """
    result = [(list_to_string(equiv), key) for key in original_dict for equiv in original_dict[key]]
    return result

def suspension_bfs(root: Node, target: str, causes_list: list, data_table: list, target_factor_level: int = 0, active: list = [], \
                   max_depth: int = 0, counter: int = 0, max_disj: int = 0, max_conj: int = 0, suspension_acc: float = 0.1, threshold: float = 0.9) -> tuple:
    """Dscepr

    Parameters
    __________
    root: Node
        root node of the tree
    target: str
        target factor whose equivalents are searched for
    causes_list: list of str
        list of factors that will populate the first level of tree,
        together with their negations
    data_table: list of dict (str, bool)
        truth table in form of a list of dictionaries, each row corresponds to
        one list element, each element is dictionary with the same keys (the factors)
        and Boolean values
    target_factor_level: int, optional
        constitutive level of the target factor
    active: list, optional
        list of non-suspended nodes
    max_depth: int, optional
        maximal number of nodes to visit during the search
    counter: int, optional
        intial value for the counter of how many nodes have been visited
        used for recursive calls
    max_disj: int, optional
        maximum number of disjunct per DNF formula,
        value zero is treated as no upper limit
    max_conj: int, optional
        maximum number of conjuncts per disjunct,
        value zero is treated as no upper limit
    suspension_acc: float, optional
        limit value for accuracy, nodes below this limit get suspended
    threshold: float, optional
        accuracy threshold for accepted formulae in case that no formula
        with accuracy=1 is found

    Returns
    _______
    bool
        True if DNF formulae that are equivalent to target_factor have been
        found, False otherwise
    list of Node
        if first returned value is True, the list of nodes contains all
        nodes with accuracy=1 that have been found
        if first returned value is False, list contains all found nodes
        with accuracy > threshold
    int
        number of visited nodes until the last node that is contained in
        the returned list has been created
    """

    if root is None:
        return False, [], 0

    created_nodes = [root]
    active_nodes = root.get_all_nodes()
    suspended_nodes = []
    ancestors = []
    queue = deque([(root, ancestors, active_nodes, suspended_nodes)])  # Store tuple of (current node, list of ancestors, active nodes, suspended nodes)
    equivalent_list = []  # list to store all formulae that are equivalent to the target factor or close enough
    good_enough_list = [] # list for nodes that are not perfect but accuracy > threshold
    last_element_added = 0 # counts in which circle the last element has been added to equivalent_list

    while queue:
        current_node, ancestors, active_nodes, suspended_nodes = queue.popleft()
        counter += 1

        if max_depth != 0 and counter > max_depth:
            break


        # set suspension conditions
        if current_node != root and current_node.accuracy < suspension_acc or \
           current_node.recall < threshold and len(current_node.value) == max_disj:
               # condition 1: node's accuracy is too low
               # condition 2: recall is below threshold, but maximum number of disjuncts is reached
            current_node.suspended = True
            if current_node in active_nodes:
                active_nodes.remove(current_node)
            if not(current_node in suspended_nodes):
                suspended_nodes.append(current_node)
            continue

        # temporarily suspend this node's processing
        if current_node.suspended:
            continue

        # check if the current node meets the target
        if current_node.accuracy == 1.0 and not(current_node in equivalent_list):
            equivalent_list.append(current_node)  # add to found nodes
            last_element_added = counter
            # suspend this node since it needs no further extension,
            # neither should it be used to extend other nodes
            if current_node in active_nodes:
                active_nodes.remove(current_node)
                suspended_nodes.append(current_node)
        elif current_node.accuracy > threshold:
            good_enough_list.append(current_node)
            if not(equivalent_list):
                last_element_added = counter

        # populating the tree
        if current_node == root and active_nodes == []:
            # in the first run, active_nodes is replaced by the list of factors -> the first children will be the set of literals
            current_node.create_new_nodes(causes_list, data_table, target, created_nodes, suspended=False, target_factor_level=target_factor_level, max_disj=max_disj, max_conj=max_conj)
            for child in current_node.children:
                created_nodes.append(child)
                if active:
                    if not(child.name in active or child.name[1:] in active):
                        child.suspended = True

        elif current_node != root:
            current_node.create_new_nodes(active_nodes, data_table, target, created_nodes, suspended=False, target_factor_level=target_factor_level, max_disj=max_disj, max_conj=max_conj)


        for child in current_node.children:
            if child.suspended or child.accuracy < suspension_acc:
                child.suspended = True
                suspended_nodes.append(child)
                if child in active_nodes:
                    active_nodes.remove(child)
            else:
                if current_node == root:
                    active_nodes.append(child)
                if child.accuracy == 1.:
                    equivalent_list.append(child)
                    last_element_added = counter
            if current_node != root:
                queue.append((child, ancestors + [current_node], active_nodes, suspended_nodes))
            else:
                queue.append((child, ancestors, active_nodes, suspended_nodes))

    if equivalent_list:
        return True, equivalent_list, last_element_added  # return all found nodes
    elif good_enough_list:
        return False, good_enough_list, last_element_added
    elif suspended_nodes and counter < max_depth:
        # After finishing BFS check if no solution was found, reevaluate suspensions.
        # loop through all nodes and reactivate them if suspended
        def reactivate_suspended(node):
            node.suspended = False
            for child in node.children:
                reactivate_suspended(child)

        reactivate_suspended(root)
        return suspension_bfs(root, target, causes_list, data_table, target_factor_level=target_factor_level, max_depth=max_depth, counter=counter, max_disj=max_disj, max_conj=max_conj, suspension_acc=suspension_acc/2., threshold=threshold)
    else:
        return False, good_enough_list, last_element_added


def reduce_data_table(data_table: list, remove_list: list) -> list:
    """Removes elements from remove_list from the keys in the list of dictionaries data_table.
    Returns the new list of dictionaries without the keys listed in remove_list. In case that
    reduction of keys leads to identical dictionaries, only one is kept.

    Parameters
    __________
    data_table: list of dict (str, bool)
        truth table in form of a list of dictionaries, each row corresponds to
        one list element, each element is dictionary with the same keys (the factors)
        and Boolean values
    remove_list: list of str
        list of keys to be removed from the dictionaries

    Returns
    _______
    list of dict
        data_table reduced by the keys listed in remove_list and without duplicate elements
    """
    new_table = [] # create new list of lines for the truth table

    # copy data for all keys that are not in factors_to_remove
    for index, line in enumerate(data_table):
        new_table.append({})
        for key in line:
            if not(key in remove_list):
                new_table[index][key] = line[key]

    # look for identical lines and delete them
    for i in range(len(new_table)-1,-1,-1):
        for j in range(len(new_table)):
            if new_table[i] == new_table[j] and i != j:
                del new_table[i]
                break

    return new_table
