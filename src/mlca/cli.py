#!/usr/bin/env python3

# file: cli.py
"""
This script generates multi-level causal-mechanistic models from Boolean data tables that are provided with
tiered list of causal factors (assignment to constitutive levels).
It generates causal graphs for each unique solution. These are finally exported
into a graph in Latex TikZ-code.

It proceeds in seven main steps:

1. obtains lists of the causal factors with level assignment, identifies the equivalence formulae
2. categorises them into constitution relations and causal relations of different levels,
   already discards all formulae that do not fit in any of these categories (function read_input)
3. obtains the list of all causal structures that are compatible with the causal relations (function find_structures)
4. determines a non-strict total order of the causal factors of each constitutive level
5. derives constitution relations between factors of different levels
6. prepares the lists for the graphical output (grouping of related causal factors, discarding of some constitution relations)
7. translates the obtained structures into a graphical output via Latex
"""

import codecs                      # for en- and decoding of strings (esp. to write tex-files in utf-8)
import os                          # operating system interfaces is required to find the files in the local path
import sys                         # system-specific parameters and functions, needed to get optional script arguments
import time                        # for measuring the run time

__all__ = ("main")

from utils import get_formula_level, get_components_from_formula, get_causal_prefactors, flatten_nested_list
# further file that contains functions for deriving equivalence formulae from a given truth table:
from atomic_formulae import read_data_from_csv
# functions for plotting the results:
import plot_graph
import mlca


def main(input_file="", input_type="", force_mode="") -> str:
    """Function to generate multi-level causal mechanistic models from data tables

    Returns
    _______
    str
       name of created pdf file
    """
    latex_template_file = "Latex_Template.tex"#str(Path("..").resolve().joinpath("config").joinpath("Latex_Template.tex"))
    start_time = time.time()
    
    level_factor_list = []               # declaration of the factor list    
    order_input_information = []
    abort = False

    pdf_file = ""
    
    # sys.argv is the list of arguments given when executing the script.
    # e.g. "python script.py --csv" contains two arguments "script.py" (the script name ifself is one element of sys.argv) and "--csv"
    # check whether special options like color mode or output of the full list of formulae has been activated
    mode = []  # list of special output modes depending on optional parameters (see below)
    # prepare these special modes
    # there are two plot modes possible:
    # "bw" - black/white
    # "color" - in color
    # The plot mode can be specified when running the script by adding "-c" or "-bw" respectively.
    mode.append("color") # Standard mode is color.

    if any(arg == "-c" or arg == "--color" for arg in sys.argv) or "color" in force_mode:
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
    if any(arg == "-fl" or arg == "--fulllist" for arg in sys.argv) or "fulllist" in force_mode:
        # sys.argv is the list of arguments given when executing the script.
        mode.append("fulllist")                       # e.g. python script.py -c -fl (The script ifself is one element of sys.argv.)
        separate_formula_list = []                    # list of formulae in tex-code

    # simple mode does not derive complex structures between co-extensive factors
    # results of simple mode should be the same as those of cna
    if any(arg == "-c" or arg == "--complex" for arg in sys.argv) or "complex" in force_mode:
        pass                                          # e.g. python script.py -c (The script ifself is one element of sys.argv.)
    else:                                             # default behaviour: simple mode without complex structures between
        mode.append("simple")                         # co-extensive factors

    # selecting mode of deriving asf, either "td" (default) or "bu"
    mode.append("td")
    if any(arg == "-td" for arg in sys.argv) or "td" in force_mode:
        if not("td" in mode):
            mode.append("td")
    elif any(arg == "-bu" for arg in sys.argv) or "bu" in force_mode:
        if not("bu" in mode):
            mode.append("bu")
            # remove the standard mode
            mode.remove("td")

    for i in range(len(sys.argv)):
        if sys.argv[i] == "--r-import" or input_type == "R":
            # "--r-import" is the command to use the text output of a r-script (cna or QCA) 
            if i < len(sys.argv) or input_file != "":
                if i < len(sys.argv) and sys.argv[i] == "--r-import":
                    input_file = sys.argv[i+1] # assume that the next argument is path to the r-output file
    
                # step 1 function read_input -> converts cna/QCA output into lists of causal factors and equivalence relations
                # if the output is not as expected, stop the procedure with abort = True
                if os.path.exists(input_file):
                    abort, level_factor_list, equiv_list = mlca.read_input(input_file)
                    if not(abort):
                        # classify the equivalence relations into causal relations (level_equiv_list) subdivided by constitutive level
                        # or into constitution relations
                        level_equiv_list, constitution_relation_list = mlca.categorise_formulae(equiv_list, level_factor_list)
                        
                else:
                    # input_file does not exist
                    abort = True
                    print("Error: Expected input data file " + input_file + " has not been found in path folder.")
            
            else:
                # no input_file specified
                abort = True
                print("Error: Missing expected file path after \"--r-import\" command.")
                
            
            break # break from loop over script parameters after "--r-import" has been found
        
        
        elif sys.argv[i] == "--csv" or input_type == "csv": 
            # Read truth table from csv-file and generate equivalence formulae
            if i < len(sys.argv) or input_file != "":
                if sys.argv[i] == "--csv" and i < len(sys.argv):
                    input_file = sys.argv[i+1] # assume that the next argument is path to the csv table
                if os.path.exists(input_file):
                    # execute function from obtain_equivalence_formulae.py
                    if "bu" in mode:
                        abort, level_factor_list, equiv_list, order_input_information = read_data_from_csv(input_file, mode="bu")
                    else:
                        abort, level_factor_list, equiv_list, order_input_information = read_data_from_csv(input_file)
                        
                    if not(abort):
                        # classify the equivalence relations into causal relations (level_equiv_list) subdivided by constitutive level
                        # or into constitution relations
                        print(str(round(time.time() - start_time,2)) + " seconds needed to find " + str(len(equiv_list)) + " equivalences.")
                        level_equiv_list, constitution_relation_list = mlca.categorise_formulae(equiv_list, level_factor_list)

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
        # prepare the list for the tex-code        
        tex_table = []

        # step 3 generate all possible solutions for the underlying causal structure
        solution_term_list = mlca.find_structures(level_factor_list, level_equiv_list, mode)

        # step 4: determine the causal order
        complete_sol_list = []
        for sol in solution_term_list:
            noncircular, list_of_pairs = mlca.reduce_redundancies(level_factor_list, sol)
            if noncircular:
                for pair in list_of_pairs:  # pair is a 2-tuple of [0]: level factor order list and [1]: list of causal relations        
                    if not(pair in complete_sol_list):
                        complete_sol_list.append(pair) # complete_sol_list is a list of tuples [0]: factor list nested by level and causal ordering,
                                                       # [1]: causal relations list, nested by level
            
        # check whether some solutions are fragment of others -> ignore them
        for i in range(len(complete_sol_list)-1,-1,-1):
            for sol in complete_sol_list:
                if sol != complete_sol_list[i]:
                    fragment = True
                    for lvl in range(len(complete_sol_list[i][1])):            
                        for formula in complete_sol_list[i][1][lvl]:
                            if not(formula in sol[1][lvl]):
                                fragment = False
                                break # break from for-loop over formulae once one has been found that is not included in sol
                        if not(fragment):
                            break # break from for-loop over levels once one formula has been found that is not included in sol
                    if fragment:
                        del complete_sol_list[i] # remove fragment from list of solutions
                        break # break from for-loop over solutions once one has been found that includes all formulae of complete_sol_list[i]
       
         
        # prepare each solution individually for graphical output
        sol_counter = 0
        total_solutions = len(complete_sol_list)

        for sol in complete_sol_list:
            # prepare a local version of level_equiv_list
            new_level_equiv_list = []
            for i in range(len(level_factor_list)) :
                new_level_equiv_list.append([])
       
            for eq_lvl in sol[1]:   # sol[1] is a nested list of equivalence relations that constitutes one causal model
                for formula in eq_lvl :
                    new_level_equiv_list[get_formula_level(formula[0], level_factor_list)].append(formula)

            # step 5: arrange the factors for improved placement in the plot
            # also get rid of unnecessary constitution graphs
            # (only the outer left and outer right part of the structure constituting a higher level factor should be drawn)
            new_level_factor_list_order, new_constitution_relation_list, color_map, new_level_equiv_list = mlca.arrange_factors(sol[0], new_level_equiv_list, constitution_relation_list, mode) # sol[0] is the
            # corresponding level_factor_list_order to sol[1] (the level_equiv_list)
            
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
                    # the input information has to be consistent with the obtained causal order of the particular solution
                    for pre_order_rel in order_input_information[lvl]:
                        if pre_order_rel[1] in get_causal_prefactors(pre_order_rel[0], sol[1][lvl], flatten_nested_list(sol[0][lvl])):
                            # pre_order_rel expresses that element [1] must not be a cause of element [0]
                            # thus, remove sol, if it does not conform to this
                            order_preserved = False
                            total_solutions = total_solutions - 1 # do not count this solution
                            break # break from for-loop over predefined order relations
                    if not(order_preserved):
                        break # break from for-loop over levels


                        """ removed
                        # the input information is consistent with the obtained causal order of the particular solution
                        # iff the set of the input order relations is a subset of the obtained order relations (per level)
                        if not(set(order_input_information[lvl]).issubset(order_relations[lvl])):
                            order_preserved = False
                            total_solutions = total_solutions - 1 # do not count this solution
                            break
                        """
                
                # limiting the output to 1000 solutions
                max_solutions_plot = 1000
                if order_preserved and sol_counter < max_solutions_plot:
                    #############################################
                    # step 6: graphical output as a graph in pdf                         
                    # generating the tex-code 
                    sol_counter = sol_counter + 1
                    try:
                        st = plot_graph.print_structure_in_tikz_plot(new_level_factor_list_order, new_level_equiv_list, new_constitution_relation_list, color_map)
                    except:
                        st = ""
                            
                    # subtitle of the graph will be the formula in tex-math syntax
                    subtitle = plot_graph.convert_formula_to_tex_code(new_level_equiv_list)
                    subtitle = "\\tiny " + subtitle # formulae might be quite long, so subtitle should be written in tiny
                        
                    # counter enumerates the solutions
                    entry = (sol_counter, st, subtitle)
                    tex_table.append(entry)

            
        # after one entry for each solution has been generated in tex_table compile the pdf
        if tex_table:
            # if tex_table is non-empty
            pdf_file = plot_graph.create_pdf(tex_table, latex_template_file, total_solutions)
            print("Number of solutions " + str(total_solutions))
            if total_solutions > 1000:
                print("More than 1000 solutions obtained, plotting only the first 1000.")
            if "fulllist" in mode:
                # when in fulllist mode, add the formulae to the list to be exported
                for sol in solution_term_list:
                    separate_formula_list.append("$" + convert_formula_to_tex_code(sol) + "$")       
                # create file
                mlca.create_separate_formula_list(separate_formula_list)
             
            # print the duration of the total run time into the console
            print("Total run time " + str(round(time.time() - start_time,2)) + " seconds.")  
            
        else:
            # solution_list is empty
            # no solution survived selection of valid solutions
            print("No valid complex solution formula has been found in " + input_file + ".") 

    return pdf_file
if __name__ == '__main__':
    # start main() function when the py-file is executed
    main()
