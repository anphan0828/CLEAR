# CLEAR
Concise List Enrichment Analysis Reducing Redundancy 


## CLEAR.R

`CLEAR.R` contains the main implementation of the CLEAR method. The script defines the `CLEAR()` function, which performs Bayesian gene set enrichment analysis while reducing redundancy among gene sets.

The function takes gene-level statistics (either test statistics or p-values) and a gene set annotation list as input, and returns posterior probabilities indicating whether each gene set is enriched. The model uses an MCMC sampling procedure to estimate enrichment configurations and related parameters.

This file serves as the core implementation of the method and can be sourced directly in R by running:

source("CLEAR.R")

After sourcing the script, users can call the `CLEAR()` function on their own data.


## CLEAR_paper

The `CLEAR_paper/` folder contains the scripts used to reproduce the analyses and figures presented in the CLEAR manuscript. These scripts include the generation of simulated datasets, benchmarking of CLEAR against other enrichment methods, and the code used to generate the plots and summary statistics reported in the paper.

## Example

The `example/` folder contains a small tutorial demonstrating how to run CLEAR with example data. It includes the input data, the CLEAR function script, and a runnable tutorial script that performs the analysis and saves the results.

- **`CLEAR.R`**  
  The main implementation of the CLEAR method. The tutorial script sources this file to run the model.

- **`example_data.rds`**  
  Example gene-level statistics used in the tutorial. This file contains:
  - `genes` – gene identifiers  
  - `p_values` – simulated p-values  
  - `t_stat` – simulated test statistics  

- **`GO_ecoli.rds`**  
  Gene Ontology (GO) gene set annotations for *E. coli*, used as the pathway database in the example.

- **`run_example.R`**  
  A runnable tutorial script that demonstrates how to apply CLEAR using both test statistics and p-values. The script runs several model configurations and saves the results.

- **`results/`**  
  A folder automatically created by the tutorial script to store the output `.rds` result files.
