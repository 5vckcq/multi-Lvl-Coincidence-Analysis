#!/usr/bin/env python3

# generates causal models from cna-asfs and stores them into text files

import multiprocessing    # multiprocessing and functools for multicore usage
import re
from os import listdir, stat
from os.path import isfile, join
import mLCA


def tuple_to_string(tup):
    # converts formulae from mLCA which are treated as tuples into string
    # first element is left side of an equivalence formula, second the right equivalent
    # adds parentheses around the formula
    st = '(' + tup[0] + ' <-> ' + tup[1] + ')'
    return st

def get_models_from_asf(inputfile):
    path_to_mlca_output = './mlca_output/'
    
    print("start with " + inputfile)
    # mLCA generation of causal models for asf from inputfile
    inputfile = './cna_output/' + inputfile
    abort, level_factor_list, equiv_list = mLCA.read_input(inputfile)

    if not abort:
        level_equiv_list, constitution_relation_list = mLCA.categorise_formulae(equiv_list, level_factor_list)
        solution_term_list = mLCA.find_structures(level_factor_list, level_equiv_list)
        
        complete_sol_list = []
        for sol in solution_term_list:            
            non_circ, list_of_pairs = mLCA.reduce_redundancies(level_factor_list, sol)
            
            if non_circ:
                for red_sol in list_of_pairs:
                    if not(red_sol in complete_sol_list):
                        complete_sol_list.append(red_sol) # complete_sol_list is a list of tuples 
                        # [0]: factor list nested by level and causal ordering,
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
        
        # convert solution list into list of strings comparable to cna-output
        list_of_solutions = []
        for sol in complete_sol_list:
            st = ""
            for formula in sol[1][0]:
                st = st + tuple_to_string(formula) + '*'
            list_of_solutions.append(st[:-1])

        # export into a new text file    
        file_name = path_to_mlca_output + "mLCA_" + re.split("_",inputfile)[2]
        with open(file_name, 'w') as file:
            file.write("\n".join(list_of_solutions))
    else:
        # if no solution could be found, create empty file
        file_name = path_to_mlca_output + "mLCA_" + re.split("_",inputfile)[2]
        with open(file_name, 'w') as file:
            file.write('') 
    return


def main():
    path_to_cna_output = './cna_output'
    
    
    
    # list of all files in path_to_cna_output
    files_list = filter( lambda file_name: isfile 
                       (join(path_to_cna_output, file_name)), 
                        listdir(path_to_cna_output) ) 
  
    
    # sort list by file size  
    files_list = sorted( files_list, 
                        key =  lambda size: stat 
                       (join(path_to_cna_output, size)).st_size) 
  
                       
    """                       
    files_list = [path_to_cna_output + "/" + file for file in listdir(path_to_cna_output) if isfile(join(path_to_cna_output, file))]    
    
    
    """
    # start working on 11 CPUs    
    with multiprocessing.Pool(11) as pool:
        # call the function for each item in parallel
        pool.map(get_models_from_asf, files_list)
    pool.close()
    pool.join()
    """
    
    # single core variant   
    for file in files_list:
        get_models_from_asf(file)
    """
    return
    
    

if __name__ == '__main__':
    # start main()
    main()
