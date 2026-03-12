import os
import subprocess
import pandas as pd
import pickle as cp
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PAPER_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(PAPER_DIR)
DEFAULT_REAL_DATA_DIR = os.path.join(PAPER_DIR, "real_data")


def resolve_existing_path(path_value, data_dir):
    """Resolve a path from logs/args against common project roots."""
    if path_value is None:
        return None
    if os.path.isabs(path_value):
        return path_value

    candidates = [
        path_value,
        os.path.join(data_dir, path_value),
        os.path.join(PAPER_DIR, path_value),
        os.path.join(PROJECT_ROOT, path_value),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return os.path.abspath(candidate)
    return os.path.abspath(os.path.join(PROJECT_ROOT, path_value))


# Parsing SLURM_ARRAY_TASK_ID and array size
parser = argparse.ArgumentParser(description='Taking in array ID and figuring out what combination to take')
parser.add_argument('SLURM_ARRAY_TASK_ID', type=int, help='An integer from the array manager')
parser.add_argument('ARRAY_SIZE', type=int, help='Size of the array')
parser.add_argument('--dir', '-d', type=str, default=DEFAULT_REAL_DATA_DIR,
                    help='Directory containing run files (4_run_list, main_result.txt, results/, saved_data/, log files).')
# parser.add_argument('--data_folder', type=str, default="data", help='Path to data folder')
parser.add_argument('--data_prepare', type=str, default=None, help='Path to data preprocessing script')
args = parser.parse_args()
SLURM_ID = args.SLURM_ARRAY_TASK_ID
ARRAY_SIZE = args.ARRAY_SIZE
DATA_DIR = os.path.abspath(args.dir)

if not os.path.isdir(DATA_DIR):
    raise FileNotFoundError(f"Data directory not found: {DATA_DIR}")

run_list_path = os.path.join(DATA_DIR, "4_run_list")
main_result_path = os.path.join(DATA_DIR, "main_result.txt")
results_dir = os.path.join(DATA_DIR, "results")
saved_data_dir = os.path.join(DATA_DIR, "saved_data")

os.makedirs(results_dir, exist_ok=True)
if not os.path.exists(main_result_path):
    with open(main_result_path, 'w') as f:
        f.write('ID_dataset\tID_method\tID_metric\tvalue\n')

# data_folder = args.data_folder
if args.data_prepare is not None:
    data_prepare = resolve_existing_path(args.data_prepare, DATA_DIR)
else:
    data_prepare = None

with open(run_list_path, "rb") as f:
    run_list = cp.load(f)

# Actually assign to job array after constructing run list
print(f"Assigning {len(run_list)} jobs to SLURM_ID {SLURM_ID} with ARRAY_SIZE {ARRAY_SIZE}")
selected_list = [run_list[id] for id in range(len(run_list)) if ((id % ARRAY_SIZE) == SLURM_ID)]
for each in selected_list:
    (dataset, method, metric) = each
    
    ID_dataset = str(dataset["ID_dataset"])
    # gene_list = f'{data_folder}/{dataset["gene_list"]}'
    gene_list = resolve_existing_path(str(dataset["gene_list"]), DATA_DIR)
    annotations = os.path.join(saved_data_dir, dataset["annotations"])
    
    ID_method = str(method["ID_method"])
    method_script = resolve_existing_path(method["method"], DATA_DIR)
    method_params = []
    for param, value in method["method_params"].items():
        method_params.append(f'-{param}')
        method_params.append(value)
    result_filename = os.path.join(results_dir, ID_dataset + "_" + ID_method + "_result.tsv")
    if not os.path.exists(result_filename):
        if method_script.endswith(".py"):
            if data_prepare is None:
                print("Running command: " + " ".join(["python3", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename] + [str(x) for x in method_params]))
                p = subprocess.run(["python3", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename] + [str(x) for x in method_params], cwd=PROJECT_ROOT)
            else:
                print("Running command: " + " ".join(["python3", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename, "-p", data_prepare] + [str(x) for x in method_params]))
                p = subprocess.run(["python3", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename, "-p", data_prepare] + [str(x) for x in method_params], cwd=PROJECT_ROOT)
        elif method_script.endswith(".R"):
            if data_prepare is None:
                p = subprocess.run(["Rscript", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename] + [str(x) for x in method_params], cwd=PROJECT_ROOT)
            else:
                p = subprocess.run(["Rscript", method_script, "-g", gene_list, "-a", annotations, "-o", result_filename, "-p", data_prepare] + [str(x) for x in method_params], cwd=PROJECT_ROOT)
        if p.returncode != 0:
            print(p.stderr)
            print(f'{ID_dataset}_{ID_method} failed')
            continue # DO NOT execute metric if method fails
            # exit(1) # DO NOT EXIT because needs to run later combination
    # Run metric only if method completes successfully
    ID_metric = str(metric["ID_metric"])
    metric_script = resolve_existing_path(metric["metric"], DATA_DIR)
    metric_params = []
    for param, value in metric["metric_params"].items():
        metric_params.append(f'-{param}')
        if isinstance(value, str) and ('/' in value or value.endswith(('.tsv', '.cp', '.txt', '.gaf', '.obo'))):
            metric_params.append(resolve_existing_path(value, DATA_DIR))
        else:
            metric_params.append(value)
    #print(" ".join(["python3", metric_script, "-i", result_filename, "-a", annotations] + [str(x) for x in metric_params]))
    if "-use_gene_file" in metric_params:
        index = metric_params.index("-use_gene_file")
        if index + 1 < len(metric_params):
            metric_params[index + 1] = gene_list # remove the flag for "use_gene_file", replace with path to gene list that is required for certain metrics
            metric_params[index] = "-g" # change name to match flag in metric script
        else:
            print("use_gene_file flag is present but no value is provided")
        print("Metric running command: " + " ".join(["python3", metric_script, "-i", result_filename, "-a", annotations] + [str(x) for x in metric_params]))
    p1 = subprocess.run(["python3", metric_script, "-i", result_filename, "-a", annotations] + [str(x) for x in metric_params],
                       capture_output=True, text=True, cwd=DATA_DIR)
    if p1.returncode == 0:
       with open(main_result_path, 'a') as f:
           f.write(f'{ID_dataset}\t{ID_method}\t{ID_metric}\t{p1.stdout}')
    else:
       print(p1.stderr)
       print(f'{ID_dataset}_{ID_method}_{ID_metric} failed')
       # exit(1) # DO NOT EXIT because needs to run later combination
        

