## Setting up environment
### Build a mamba environment 
- Create environment and install packages: ```mamba create -n clear_mamba_env -f environment.yaml```.
- Activate the environment: ```mamba activate clear_mamba_env```
- Run scripts in the **root directory** following instructions below


## Analyzing human RNA-seq and microarray datasets 

### 0. Retrieve datasets from packages
- `0_de_tcga_array.R`: collects 15 TCGA datasets from GSEABenchmarkeR and perform DESeq2
- `0_process_geo_data.R`: collects 24 microarray datasets from KEGGdzPathwaysGEO and perform limma 

### 1. Collect datasets and annotation data
```bash scripts_real_data_analysis/1_generate_dataset_log.sh <real_data_folder> <gaf_file> <obo_file>```
- `real_data_folder` (example: real_data/tcga_data) contains TSV files with 2 columns gene symbols and gene-level statistics. If the TSV file contains more than 2 columns, you must provide an R script for data preprocessing when running the methods.
- `gaf_file` (example: ext_data/03Nov24goa_human.gaf) contains GO annotations in GAF format. This file can be downloaded from GO Data Archive
- `obo_file` (example: ext_data/24go.obo) contains GO structure in OBO format.
- You can also provide annotation in the form of 2-column TSV with gene symbols and GO term. Make sure these GO terms have been propagated (genes annotated to child GO terms are also annotated to parent GO terms).
- `1_datasets_log.ndjson` stores all datasets that we used and the corresponding dataset-specific annotation file.

### 2. Fill in methods
- All methods we used are listed in `2_methods_log.ndjson`. This file was created manually, any parameters being used must be listed here.
- Example entry: ```{"ID_method": "63-1","method":"clear_v5/run_CLEAR_tnormal.R","method_params":{"i": 1000000, "b": 500000, "c":20,"C":500,"s":1}}```. Here the method script path is listed as `clear_v5/run_CLEAR_tnormal.R`, with 1,000,000 iterations (`"i": 1000000`), 500,000 burn-in steps (`"b": 500000`), only use gene sets between size 20 and 500 (`"c":20,"C":500`), and seed 1 (`"s":1`). These parameters follow the arguments defined in the method script.

### 3. Run methods and metrics
- Metric scripts are listed in `3_metrics_log.ndjson`. This file was created manually, any parameters being used must be listed here.

### 4. Construct run list
- To construct a run list: ```python3 4_prerun.py --datasets 1:15 --methods 63-1,64-1,65-1,66-1,66,67,69,70 --metrics 1:3```. If no dataset/method/metric is specified, default run list includes all.

### 5. Run the methods
- To run the dataset-method-metric combinations: ```python 5_run_selected.py 0 1 --data_prepare $data_prepare``` (where data_prepare is an R script that is used to preprocess the gene file if more than 2 columns are in the gene file)
- All result files (containing gene set and p-value or posterior probability) are stored in results/ folder. All metric scores are stored in `main_result.txt` file. This file keeps a record of which dataset/method/metric combination have been run, and running ```python3 4_prerun.py``` will check for existing results as stored in this file and exclude from run list.

### 6. Visualize results
- Data are parsed and plots are created using `6_clear_notebook.ipynb` file.