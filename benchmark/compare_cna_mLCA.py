#!/usr/bin/env python3

# file: compare_cna_mLCA.py

import re # complex string search patterns
import os.path
from auxiliary_functions import get_equiv_formula, get_components_from_formula, count_true


def read_list_from_cna(file_cna):
    # reads list of causal models from cna output - which includes further lines, such as atomic solution formulae and mere text - and aligns the syntax to mLCA
    # If cna found complex solutions, it list them after the atomic solution formulae. Then all lines before the csf can simply ignored.
    # In some cases cna returns "Same as asf", in others it list the atomic solution formulae as complex solution formulae, both cases have to be dealt with separately.
    
    lines_cna = []
    with open (file_cna, 'rt') as text_file:    # open cna-output
        for next_line in text_file:             # add content line by line to
            lines_cna.append(next_line)         # list lines_cna
    
    asf_list = []
    no_complex_solution = True                  # case csf are emitted but identical to asf
    same_as_asf = False                         # case cna returns "Same as asf"
    for line in lines_cna:
        if line.find('Same as asf') > -1:
            # cna found no complex solution:
            same_as_asf = True
            no_complex_solution = False
            break # break from for loop after "Same as asf" has been found in one line
        elif line.count('<->') > 1:
            no_complex_solution = False
            break # break from for loop after first complex solution has been found
    
    
    
    if same_as_asf or no_complex_solution:
        # find atomic solutions in lines
        for aline in lines_cna:
            if aline.count('<->') == 1:
               tup = get_equiv_formula(aline)
               st = '(' + tup[0] + ' <-> ' + tup[1] + ')'
               asf_list.append(st)
        
        if no_complex_solution:
            asf_list = list(set(asf_list)) # remove potential duplicates
            # sort according to effect (the second to last element of the string - the last one being ')')
            asf_list.sort(key=lambda x:x[-2])
            
        lines_cna = asf_list
    else:
        # else - cna generated a list of csf - delete all lines that do not include csf
        for index in range(len(lines_cna)-1,-1,-1):
            if lines_cna[index].count('<->') < 2:
                del lines_cna[index]
        # format remaining lines properly
        if lines_cna: 
            # if there are csf
            for index in range(len(lines_cna)-1,-1,-1):
                details = re.split(r'\s{2,}' ,lines_cna[index])
                exhaustiveness = 0.0
                faithfulness = 0.0
                for i in range(len(details)-2):
                    if details[i].find("TRUE") > -1 or details[i].find("FALSE") > -1:
                        exhaustiveness = float(details[i+1])
                        faithfulness = float(details[i+2])
                        break
                
                
                if exhaustiveness < 1.0 or faithfulness < 1.0:
                    #print(file_cna + " " + lines_cna[index]) #dummy
                    del lines_cna[index]
                    
                else:
                    # remove everything before first '(' and after last ')'
                    lines_cna[index] = '(' + re.split(r'\((.*)', lines_cna[index])[1]
                    # split the line at every ")"
                    aux_list =  re.split(r'\)', lines_cna[index])
                    # reconnect all elements except the last and add the ")"
                    lines_cna[index] = ')'.join(aux_list[:-1]) + ')'
                
                
                    # replace lower letters by '~' + upper letter
                    # 1. step: add "~" before each minuscle, which is either
                    # a) follows a '('-bracket
                    # b) follows a conjunctor
                    # c) follows a disjunctor
                    # a):
                    lines_cna[index] = re.sub(r'\(([a-z])',  r'(~\1', lines_cna[index])
                    # b) if following a "*", the letter will be placed behind "*~"
                    lines_cna[index] = re.sub(r'\*([a-z])',  r'*~\1', lines_cna[index])
                    # The regex expression "\*" picks the star symbol "*" from the string.
   
                    # c) if following " + ", the letter will be placed behind "+ ~"
                    lines_cna[index] = re.sub(r'\s\+\s([a-z]+)',  r' + ~\1', lines_cna[index])
                    # in regex "\s" corresponds to spaces, "\+" to "+"

                    # 2. step replacement of the minuscle that follow to "~" by majuscle
                    lines_cna[index] = re.sub(r'(~[a-z]+)', lambda pat: pat.group(1).upper(), lines_cna[index])
        
    return lines_cna
    
def read_list_from_mlca(file_mlca):
    # reads lists of causal models from text file
    
    lines_mlca = []   
    with open (file_mlca, 'rt') as text_file:    
        for next_line in text_file:             
            lines_mlca.append(next_line)          
    
    # re-grouping the conjuncts, such that the effects are ordered alphabetically (this is how cna does it)
    for index in range(len(lines_mlca)):
        # remove trailing '\n'
        lines_mlca[index] = re.sub("\r?\n","",lines_mlca[index]).rstrip()
        # remove trailing ')'
        lines_mlca[index] = lines_mlca[index][:-1]
        #print(lines_mlca[index])
        # spilt formulae by ')*'
        conj_list = re.split(r'\)\*', lines_mlca[index])
        # sort after last element
        conj_list.sort(key=lambda x:x[-1])
        #print(conj_list)
        # re-compose the conjunctive formula
        lines_mlca[index] = ')*'.join(conj_list) + ')'

    return lines_mlca
        
    
def compare(file_cna, file_mlca):
    lines_cna = read_list_from_cna(file_cna)
    lines_mlca = read_list_from_mlca(file_mlca)
    
    cna_more = False
    mlca_more = False
    
    if (set(lines_cna) - set(lines_mlca) != set()) and not(len(lines_mlca) == 0):
        diff = list(set(lines_cna) - set(lines_mlca))
        
        if diff:
            # (C) check whether all causal factors that appear in asfs also appear in complex solution
            for index in range(len(diff)-1,-1,-1):
                if get_components_from_formula(diff[index], ['A', 'B', 'C', 'D', 'E', 'F']) != get_components_from_formula(lines_mlca[0], ['A', 'B', 'C', 'D', 'E', 'F']):
                    del diff[index]
        
        if diff:
            cna_more = True
            out_string = str(diff)
            out_string = re.sub(r',\s', r'\n',out_string)
            #print('cna but not mlca:\n' + out_string + '\n set ' + str(file_mlca)) # solutions in cna, missing in mlca
    if (set(lines_mlca) - set(lines_cna) != set()):
        mlca_more = True
        out_string = str(set(lines_mlca) - set(lines_cna))
        out_string = re.sub(r',\s', r'\n',out_string)
        #print('mlca but not cna:\n' + out_string + '\n set ' + str(file_mlca)) # solutions in mlca, missing in cna
        
    
    return cna_more, mlca_more, len(lines_cna), len(lines_mlca)
    
    

if __name__ == '__main__':
    path_to_cna_output = './cna_output_final/'
    path_to_mlca_output = './mlca_output/'
    
    cna_more_counter = 0
    mlca_more_counter = 0
    nums_sol_same = 0
    nums_sol_more_m = 0
    count_empty = 0
    diff = 0
    same_counter = 0
    
    out_string = ""
    
    for i in range(10001):
        filename_mlca = 'mLCA_' + str(i).zfill(5) + '.txt'
        filename_cna = 'cna_' + str(i).zfill(5) + '.txt'
        
        if os.path.isfile(path_to_mlca_output + filename_mlca) and os.path.isfile(path_to_cna_output + filename_cna):
            cna_more, mlca_more, num_solsc, num_solsm = compare(path_to_cna_output + filename_cna, path_to_mlca_output +  filename_mlca)
            if cna_more:
                cna_more_counter = cna_more_counter + 1
            if mlca_more:
                mlca_more_counter = mlca_more_counter + 1
                nums_sol_more_m = nums_sol_more_m + num_solsm
                diff = diff + num_solsm - num_solsc
                out_string = out_string + filename_mlca + ", "
            if not(cna_more) and not(mlca_more):
                same_counter = same_counter + 1
                nums_sol_same = nums_sol_same + num_solsc
            if num_solsc == 0:
                count_empty = count_empty + 1

    print("cna but not mlca: " + str(cna_more_counter))                
    print("mlca but not cna: " + str(mlca_more_counter))
    print(out_string)
    print("total number of different models: " + str(diff))
    print("total number of models in data sets with different results: " + str(nums_sol_more_m))
    print("same: " + str(same_counter))
    print("total number of models in data sets with same results: " + str(nums_sol_same))
    print("number of data sets without causal models: " + str(count_empty))
