# CausalExplanationSimplicity

## causal_structure_truth_table_to_pdf_graph.py

Python script that generates the possible causal diagrams either from a truth table in csv-format or the R output of QCA or cna using the latex template Latex_Template.tex. The template file must be in the same folder as the script file. 

**This is the main file (to be executed).**

The input data has to be specified when running the script either by adding the parameter "--csv FILEPATH" or "--r-import FILEPATH".

The script can be executed with the optional argument "-bw" (alternatively "--blackwhite"), which changes the output from colored graphs into black/white. The argument "-c" (resp. "--color") forces the color mode, which is currently set as standard.

If the tex-formated formulae should be exported into a separate file, the optional argument "-fl" (or "--fulllist") should be set.

In general, all causal structures that are logically compatible with the input data are generated. In case of co-extensive causal factors, this results in enormous amounts of solutions. If only the most simple solutions are wanted, use the option "-s" (or "--simple"). Then, the results should be the same as those of cna.

Running the script requires besides a Python3 environment, a latex distribution installed on the system running this script, as well as the jinja2 and itertools Python3 libraries and the TikZ library for latex.


## Latex_Template.tex

This template sets the environment for plotting the obtained mechanisms. It can be modified at will but must not be renamed.


## samples

A folder that contains some samples either truth tables in csv-format or R-scripts that make use of the cna-package. Some of the examples model multi level mechanisms. The content of this folder is not required for running the script.

**test_1**
> causal structure from (Harbecke, 2018) with an added upper level

**test_2**
> a causal chain with two levels, a very simple multi-level structure 

**test_3**
> a variation of test_2 which involves disjunctions, conjunctions and negations

**test_4**
> a simple single level structure

**test_5**
> highly ambigious single level structure, taken from "Critical Tension: Sufficiency and Parsimony in QCA" (Dusa 2019)

**test_6**
> causal structure from the introductory part of "An Algorithm for Boolean Constitutive Inference for Mechanistic Explanation and Potential Measures of Simplicity"
