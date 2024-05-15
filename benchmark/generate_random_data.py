#!/usr/bin/env python3

# file: generate_random_data.py

# generates random conjuncted DNFs

import random # pseudo-random number generator
import string # provides list of letters

def generate_random_dnf(factor_list, N_d, N_c):
    m = random.randrange(1,N_d) # number of disjuncts (between 1 and N_d-1)
    st = ""
    for i in range(m): # loop over disjuncts
        if N_c < len(factor_list):
            l = random.randrange(1,N_c) # number of conjuncts in disjunct i (between 1 and N_c-1)
        else:
            l = random.randrange(1,len(factor_list)) # if N_c > number of literals in factor_list, use length of that list instead
        
        list_of_conjuncts = random.sample(factor_list, k=l) # select l different elements from factor_list
        list_of_conjuncts.sort()

        for ii in range(len(list_of_conjuncts)):
            negated = random.randrange(0,2) # if negated == 0 -> positive literal, if negated == 1 -> negated literal
            if negated == 0:
                st = st + list_of_conjuncts[ii] + "*"
            else:
                st = st + str.lower(list_of_conjuncts[ii]) + "*" # negation syntax in cna is to use lower case letters for negated factors
        st = st[:-1] + " + "                                     # delete trailing "*"

    return st[:-3]   # return formula as string (without trailing " + ")


def generate_random_asfs(N_f, N_d, N_c):
    formula = "" # formula as string
    
    if N_f < len(string.ascii_uppercase):
        factor_list = list(string.ascii_uppercase)[:N_f] # define factor_list as the first N_f letters of the alphabet
    else:
        factor_list = list(string.ascii_uppercase) # N_f > 26 is treated like N_f = 26
        
    random.seed(None) # start pseudo-random number seed with system time
    
    # determine number of equivalences to be combined into one complex (conjunction of equivalence formulae)
    n = random.randrange(1,N_f-1) # generates random number 1 <= n < N_f-1 -> we get at least two atomic terms per complex formula (start counting with 0)
    # and less than number of factors - 1 (at least one factor has to be an incoming factor)
    
    # start loop for all equivalence formulae
    for i in range(n):
        # every equivalence formula has the form DNF <-> literal
        # the first the atomic term of the first equivalence formula is the last element of factor_list, the one to the second formula the second-to-last element etc.
        # this is no restriction of the algorithm since the letters are meaningless
        
        atomic_term = factor_list[N_f-1-i]
        remaining_factor_list = factor_list[:N_f-i-1]
        complex_term = generate_random_dnf(remaining_factor_list, N_d, N_c)
        formula = formula + "(" + complex_term + " <-> " + atomic_term + ")*"

    return formula[:-1] # return formula as string without trailing "*"
    
    

if __name__ == '__main__':
    # start 10000 generate_random_asfs for 6 factors, less than five disjuncts and less than five conjuncts per DNF
    st = ""
    for i in range(10000):
       st = st + generate_random_asfs(6, 5, 5) + "\n"
    with open('random_formulae.txt', 'w') as f:
        f.write(st)
