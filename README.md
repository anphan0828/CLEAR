# CLEAR
Concise List Enrichment Analysis Reducing Redundancy 


## CLEAR.R (under `R` folder)

- **`R/CLEAR.R`** contains the main implementation of the CLEAR method. The script defines the `CLEAR()` function, which performs Bayesian gene set enrichment analysis while reducing redundancy among gene sets. The tutorial script sources this file to run the model.

The function takes gene-level statistics (either test statistics or p-values) and a gene set annotation list as input, and returns posterior probabilities indicating whether each gene set is enriched. The model uses an MCMC sampling procedure to estimate enrichment configurations and related parameters.

This file serves as the core implementation of the method and can be sourced directly in R by running:
```R
source("R/CLEAR.R")
```
After sourcing the script, users can call the `CLEAR()` function on their own data.

- **`R/CLEAR_tutorial.R`** : A runnable tutorial script that demonstrates how to apply CLEAR using both test statistics and p-values. The script runs several model configurations and saves the results.

## CLEAR_paper

The `CLEAR_paper/` folder contains the scripts used to reproduce the analyses and figures presented in the CLEAR manuscript. These scripts include the generation of simulated datasets, benchmarking of CLEAR against other enrichment methods, and the code used to generate the plots and summary statistics reported in the paper. ``CLEAR_paper/00README.md`` file describes in detail how to reproduce and/or extend the analysis of CLEAR and other methods on real datasets.

## Example

The `example/` folder contains example data to run CLEAR. It includes the input data, annotation data of *E. coli*, and example result files.

- **`example_data.rds`**  
  Example gene-level statistics used in the tutorial. This file contains:
  - `genes` – gene identifiers  
  - `p_values` – simulated p-values  
  - `t_stat` – simulated test statistics  

- **`GO_ecoli.rds`**  
  Gene Ontology (GO) gene set annotations for *E. coli*, used as the pathway database in the example.

- **`results/`**  
  A folder automatically created by the tutorial script to store the output `.rds` result files.
