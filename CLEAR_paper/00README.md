## Setting up environment

### Build a mamba environment (from project root)

- Create environment and install packages:
    `mamba create -n clear_mamba_env -f environment.yaml`
- Activate environment:
    `mamba activate clear_mamba_env`

### Define working directories

From project root (`CLEAR/`), define paths first:

```bash
PAPER_DIR="CLEAR_paper"
REAL_DATA_DIR="CLEAR_paper/real_data"
SCRIPT_DIR="CLEAR_paper/scripts_real_data_analysis"
```

There are two execution contexts in this workflow:

- **Dataset retrieval scripts (Step 0):** these scripts can be skipped if you download pre-processed data (TCGA RNA-seq and GEO microarray) from our [Figshare repository](https://doi.org/10.6084/m9.figshare.31891432).
- **Benchmarking framework scripts (Steps 1-6):** logs/results are in `CLEAR_paper/real_data`, while driver scripts are in `CLEAR_paper/scripts_real_data_analysis`.
All steps listed in the framework should be run from project root (`CLEAR/`), unless otherwise specified.

Important:

- `4_prerun.py` and `5_run_selected.py` now support `--dir` (or `-d`) to point to the real-data directory containing logs/results.
- Keep `1_datasets_log.ndjson`, `2_methods_log.ndjson`, `3_metrics_log.ndjson`, `main_result.txt`, `4_run_list`, `results/`, and `saved_data/` in `CLEAR_paper/real_data`.


## Analyzing human RNA-seq and microarray datasets

### (Optional) 0. Retrieve datasets from package sources

- `scripts_real_data_analysis/0_de_tcga_array.R`: retrieves 15 TCGA datasets from `GSEABenchmarkeR` and runs DESeq2.
- `scripts_real_data_analysis/0_process_geo_data.R`: retrieves 24 microarray datasets from `KEGGdzPathwaysGEO` and runs limma.

Example:

```bash
cd "${PAPER_DIR}" && Rscript scripts_real_data_analysis/0_process_geo_data.R
```

Notes:
- These scripts can be skipped if you download pre-processed data from the FigShare repository
- These scripts write outputs under `CLEAR_paper/real_data/tcga_data/` and `CLEAR_paper/real_data/geo_microarray_data/`.
- `0_de_tcga_array.R` takes positional arguments `(SLURM_ARRAY_TASK_ID, ARRAY_SIZE)` for array execution.
    For a local single-worker run, use:
    `(cd "${PAPER_DIR}" && Rscript scripts_real_data_analysis/0_de_tcga_array.R 0 1)`


### 1. Build dataset log and dataset-specific annotations

Make sure dataset log is created inside `CLEAR_paper/real_data` (where benchmarking logs are stored). Build dataset log and dataset-specific annotation with datasets in `tcga_data` folder:

```bash
cd "${REAL_DATA_DIR}" && bash ${SCRIPT_DIR}/1_generate_dataset_log.sh \
    tcga_data \
    ext_data/03Nov24goa_human.gaf \
    ext_data/24go.obo
```

Arguments:

- `real_data_folder`: folder containing input gene-level TSV files.
    Expected format is two columns TSV (gene symbol, statistic), with headers. If more columns are present, use `--data_prepare` in Step 5.
- `gaf_file`: GO annotation file in GAF format.
- `obo_file`: GO ontology structure in OBO format.

Outputs:

- `1_datasets_log.ndjson`: one JSON line per dataset (stored in `CLEAR_paper/real_data`).
- `real_data/saved_data/save_<ID>_0to100000_annotations.tsv`: dataset-specific propagated GO annotations.

Example line in `1_datasets_log.ndjson`:

```json
{"ID_dataset": 1, "gene_list": "real_data/tcga_data/BLCA.tsv", "annotations": "save_1_0to100000_annotations.tsv"}

```


### 2. Fill in method log (`2_methods_log.ndjson`)

- The framework reads `2_methods_log.ndjson` from the directory passed through `--dir` (default resolves to `CLEAR_paper/real_data`).
- Each method entry includes:
    - `ID_method`: method ID (string; can include repeat suffix, e.g. `63-1`)
    - `method`: path to script
    - `method_params`: argument dictionary (without the leading `-`)

Example CLEAR entry:

```json
{"ID_method":"63-1","method":"methods/run_CLEAR_tnormal.R","method_params":{"i":1000000,"b":500000,"c":20,"C":500,"s":1}}
```

Example custom method entry (new R method):

```json
{"ID_method":"200","method":"methods/run_my_method.R","method_params":{"c":20,"C":500,"alpha":0.1}}
```


### 3. Fill in metric log (`3_metrics_log.ndjson`)

- The framework reads `3_metrics_log.ndjson` from the directory passed through `--dir` (default resolves to `CLEAR_paper/real_data`).
- Each metric entry includes `ID_metric`, `metric` script path, and `metric_params`.

Example existing metric entry:

```json
{"ID_metric":1,"metric":"scripts_real_data_analysis/3_run_phenorelevance_clear.py","metric_params":{"use_gene_file":1,"p":"ext_data/GO_BP.cp"}}
```

Example custom metric entry:

```json
{"ID_metric":10,"metric":"scripts_real_data_analysis/3_run_my_metric.py","metric_params":{"n":50}}
```


### 4. Construct run list

Build combinations of dataset/method/metric that have not yet been recorded in `main_result.txt`. These commands should be run from project root.

```bash
python3 "${SCRIPT_DIR}/4_prerun.py" \
    --dir "${REAL_DATA_DIR}" \
    --datasets 1:15 \
    --methods 63-1,64-1,65-1,66-1,66,67,69,70 \
    --metrics 1:3
```
Arguments:
- `--dir` / `-d`: folder containing `1_datasets_log.ndjson`, `2_methods_log.ndjson`, `3_metrics_log.ndjson`, `main_result.txt`, and `results/`.
- If no `--datasets`, `--methods`, or `--metrics` are provided, all entries are considered.

Examples:

- Run only one dataset/method/metric combination set:
    `python3 ${SCRIPT_DIR}/4_prerun.py -d ${REAL_DATA_DIR} --datasets 1 --methods 67 --metrics 2`
- Use repeated method IDs (each method ID is run as a separate method with different seed):
    `python3 ${SCRIPT_DIR}/4_prerun.py -d ${REAL_DATA_DIR} --methods 63-1:10,64-1:10 --metrics 4`


### 5. Run selected jobs

```bash
python3 "${SCRIPT_DIR}/5_run_selected.py" 0 1 --dir "${REAL_DATA_DIR}"
```

With preprocessing script for >2-column gene files:

```bash
python3 "${SCRIPT_DIR}/5_run_selected.py" 0 1 \
    --dir "${REAL_DATA_DIR}" \
    --data_prepare "${SCRIPT_DIR}/data_prepare_tcga.R"
```

Notes:

- `--dir` / `-d` selects where `4_run_list`, `main_result.txt`, `results/`, and `saved_data/` are read/written.
- Method scripts and metric scripts are resolved using paths in the method/metric logs, so you can launch `5_run_selected.py` from any directory.
- `--data_prepare` is specific to type of input dataset. In this project `data_prepare_tcga.R` is used for TCGA RNA-seq datasets and `data_prepare_geo.R` is used for GEO microarry datasets. 

Outputs:

- `results/<ID_dataset>_<ID_method>_result.tsv`: method result tables.
- `main_result.txt`: metric outputs and run-completion log.


### 6. Visualize and summarize

- Main result parsing and plotting notebook:
    `CLEAR_paper/scripts_real_data_analysis/6_clear_notebook.ipynb`
- Runtime/memory benchmark helper script:
    `CLEAR_paper/scripts_real_data_analysis/7_runtime_analysis.R`


## Brief HPC/SLURM notes

- The run scripts are array-job oriented:
    - `0_de_tcga_array.R` and `5_run_selected.py` expect `(SLURM_ARRAY_TASK_ID, ARRAY_SIZE)` positional arguments.
    - Local serial execution is still possible with `0 1` (meaning all jobs will be run sequentially one at a time)
