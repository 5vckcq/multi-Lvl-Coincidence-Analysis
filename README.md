# multi-Lvl-Coincidence-Analysis

multi-Lvl-Coincidence-Analysis (mLCA) is a Python package for causal-mechanistic modeling. The core function of the package is to derive causal-mechanistic models from coincidence data tables. The resulting models are visualized using Ti<em>k</em>Z.


## Installation

Install the required Python packages from environment.yml using conda. Running the main script mlca.py requires a latex distribution installed on the system, as well as the Ti<em>k</em>Z library for latex.

## Usage

We provide an (experimental) graphical interface. To run it type:
```
$ python mlca-gui.py
```
You can also use the command line application. In this case you have to specify the path to your data table by replacing 'PATH' in
```
$ python mLCA.py --csv PATH
```
Alternatively, when using mLCA alongside the R packages cna or QCA, you have to save the output from these packages in a textfile (see the R-scripts in the samples folder for an instruction on how to do this) and indicate the path to this file with 
```
$ python mLCA.py --r-import PATH
```
If any causal ordering has been defined in cna or QCA, it will be reinterpreted as level separation in mLCA.

Further optional arguments can be appended:
* "-bw" (alternatively "--blackwhite") changes the output from colored graphs into black/white. The argument "-c" (resp. "--color") forces the color mode, which is currently set as default.
* If the tex-formated formulae should be exported into a separate file, the optional argument "-fl" (or "--fulllist") should be set.
* By default complex relations between co-extensive variables are excluded, which makes mLCA for single level models comparable to cna. If all possible models should be generated, this can be done by adding the argument "-c" or "--complex".

### Preparing the data tables

The table head of the csv-table should contain the labels for the causal factors.
The following entries in the cells in the table body are interpreted as True: "1", "T", "t", "w", "W", "true", "True". The strings "0", "F", "f", "false", "False" are interpreted as False.
Admissible column separators are ":", ",", ";", "|" and "_".
  
It is possible to include information on the causal order and on constitutive levels. In order to do so, add the following separators after the last row of truth values:
* The entry "<<" is interpreted as level-separator. The respective factor and those of the columns to its right are assigned to a higher level than those factors corresponding to the columns to its left.
* The character "<" marks the causal stream within the same constitutive level. Variables in columns to the left cannot be effects to the variable associated with this column or to its right.

### Customizing the graphical output

The causal-mechanistic hypergraphs are generated using the tex-template file Latex_Template.tex. You can modify it at will to change all further outputs. You can also modify the generated output_graph.tex to customize an individual hypergraph.

## Documentation

cf. comments in source code

## samples

The folder samples contains some examples for truth tables in csv-format or R-scripts that make use of the cna-package. Some of the examples model multi-level mechanisms.

**sample_1**
> single level causal structure

**sample_2**
> a simple two-level structure 

**sample_3**
> a complex three-level structure

## License

The code in this repository, including all code samples is released under the [GPL-3.0 license](LICENSE).
