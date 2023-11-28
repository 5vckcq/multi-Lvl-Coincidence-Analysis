#!/usr/bin/env python3

# file: plot_graph.py

import os                          # operating system interfaces is required to find the files in the local path
import re                          # regex for complex search patterns in strings
import jinja2                      # Latex interface
import codecs                      # for en- and decoding of strings (esp. to write tex-files in utf-8)

from auxiliary_functions import get_components_from_formula, get_factor_level, get_factor_order, get_formula_level, get_formula_order

# syntax definitions for Latex expressions
latex_jinja_env = jinja2.Environment(
	block_start_string = '\BLOCK{',
	block_end_string = '}',
	variable_start_string = '\VAR{',
	variable_end_string = '}',
	autoescape = False,
	loader = jinja2.FileSystemLoader(os.path.abspath('.')))

def convert_formula_to_tex_code(solution):
    # converts a full solution for a structure into its corresponding logical formula in tex-syntax
    # solution is expected to be a list of levels, which is a list of formulae, each of which is 2-tuple: first element is the possibly
    # complex left-side term of an equivalence formula, the second element is the atomic right-side term
    # solution[LEVEL][FORMULA][(LEFT TERM, RIGHT TERM)]
    tex_code = "$"
    for lvl in solution:
        for term in lvl:
            tex_code = tex_code + "(" + term[0].replace("*", " \cdot ").replace("~", "\\neg ") + "\leftrightarrow " + term[1] + ")\cdot"
                        
    tex_code = tex_code[:-5] + "$"  # remove the "\cdot" at the end of the last term
    return tex_code
    
def convert_causal_relation(formula, level_factor_list_order, tex_code, color, color_map):
    # translates a formula of causal relations into TikZ-Latex code
    # returns the code as string
    
    st = ""  # output code
    ###########################################
    # determine the syntactic type of formula #
    ###########################################
    # assumption: there are only the following types of formula (all DNFs):
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
                           
    if st == "":
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

                if disj in get_components_from_formula(formula[0], level_factor_list_order):
                    # case A: the discunct is one causal factor
                    
                    # in order to prevent overlapping horizontal arrows                    
                    # get positions of disj and formula[1] in their respective level_factor_list_order sublists
                    # first position i of disj
                    found = False
                    for order in level_factor_list_order[level]:
                       for i in range(len(order)):
                           if order[i] == disj:
                               found = True
                               break
                       if found: break
                    
                    # now position j of formula[1]
                    found = False
                    for order in level_factor_list_order[level]:
                       for j in range(len(order)):
                           if order[j] == formula[1]:
                               found = True
                               break
                       if found: break

                    if i == j:
                        # a straight arrow is drawn from source factor.north east to target factor.west
                        st = st + "% simple disjunction with shifted starting point\n"
                        st = st + "\draw[->, " + color + "] (" + disj + ".north east) to (" + formula[1] + ".west);\n"
                    else:
                        # a straight arrow is drawn from source factor to target factor
                        st = st + "% simple disjunction\n"
                        st = st + "\draw[->, " + color + "] (" + disj + ".east) to (" + formula[1] + ".west);\n"
                
                elif disj.find("*") > -1 :
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
                    else:
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
            
            print(formula[0] + " -> " + formula[1] + "  could not be drawn because the causal structure has not been recognized.")
            
    return st
    
    
def convert_constitution_relation(formula, level_factor_list_order, constitution_relation_list, color) :
    # converts a formula of constitution relations into TikZ-Latex code
    # returns the code as string
    
    st = ""
    
    # constitution relations are drawn differently depending on whether they are to the left or to the right of the upper level factor
    c_left = True
    c_right = True
    
    # check whether it is a left- or rightside relation
    for f in constitution_relation_list :
        if (formula[1] == f[1]) and (formula[0] != f[0]) :
            # is there a further constitution relation to the same causal factor which includes factors of higher causal order
            # than those from formula? -> if true it is a leftside relation
            # if there is no further constitution relation it is neighter left- nor rightside
            # if there further relations but of lower order -> rightside relation
            if get_formula_order(formula[0], level_factor_list_order) < get_formula_order(f[0], level_factor_list_order) :
                c_right = False
            elif get_formula_order(formula[0], level_factor_list_order) > get_formula_order(f[0], level_factor_list_order) :
                c_left = False
    
    # draw one connecting line toward formula[1] for each causal factor in formula[0]
    # (usually there should only be one factor in formula[0])
    for fac in get_components_from_formula(formula[0], level_factor_list_order) :           
        if c_left and not(c_right) :
            # case 1: leftside relation
            st = st + "\draw[crelationleft, " + color + "] (" + fac + ".north west) to (" + formula[1] + ".south);\n"
    
        elif not(c_left) and c_right :
            # case 2: rightside relation
            st = st + "\draw[crelationright, " + color + "] (" + fac + ".north east) to (" + formula[1] + ".south);\n"
    
        elif c_left and c_right: 
            # case 3: single component
            st = st + "\draw[crelationstraight, " + color + "] (" + fac + ".north) to (" + formula[1] + ".south);\n"
    
    return st
    
def print_structure_in_tikz_plot(level_factor_list_order, level_equiv_list, constitution_relation_list, color_map) :
    # Prepares the TikZ code for plotting one solution
    # a) places the causal factors as nodes separated by causal order (horizontally) and constitutive level (vertically)
    # b) adds the causal and constitution relations as vertices between nodes
    
    ######################################
    # step 1: preparing the output files #
    ######################################
    
    tex_code = "% placement of the nodes\n"
    
    ###################################################
    # step 2: placement of the nodes = causal factors #
    ###################################################
    
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
        
    ######################################################
    # step 3: plot the causal and constitution relations #
    ######################################################
    tex_code = tex_code  + "\n% causal relations\n"
    for m in range(len(level_equiv_list)): 
        tex_code = tex_code  + "% of level "  + str(m) + "\n"
        for formula in level_equiv_list[m]:

            tex_code = tex_code  + "% formula: "  + formula[0] + " <-> " + formula[1] + "\n"
            
            # standard color is black 
            color = "black"
            if color_map["draw"][formula[1]] == color_map["draw"][get_components_from_formula(formula[0], level_factor_list_order)[0]] :
                # if the color map entry of target node and the first source node (any other would do it likewise) are identical
                # use their color for plotting the causal relation
                color = color_map["draw"][formula[1]]
                
            tex_code = tex_code + convert_causal_relation(formula, level_factor_list_order, tex_code, color, color_map) + "\n\n"
    
    tex_code = tex_code  + "\n% constitution relations\n"        
    for formula in constitution_relation_list:
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
