import os
import sys
import pandas as pd
import argparse
import pickle as cp


if not os.path.exists("main_result.txt"):
    with open('main_result.txt', 'w') as f:
        f.write(f'ID_dataset\tID_method\tID_metric\tvalue\n')

if not os.path.exists("results/"):
    os.makedirs("results/")

# Get list of datasets, methods, and metrics; construct run list
def get_list_data_files(chosen_list=None):
    df = pd.read_json("1_datasets_log.ndjson",lines=True, dtype={'ID_dataset': str})
    if chosen_list is not None:
        df_chosen = df[df['ID_dataset'].isin(chosen_list)]
    else:
        df_chosen = df
    return df_chosen.to_dict(orient='records')
    

def get_list_methods(chosen_list=None):
    df = pd.read_json("2_methods_log.ndjson",lines=True, dtype={'ID_method': str})
    if chosen_list is not None:
        df_chosen = df[df['ID_method'].isin(chosen_list)]
    else:
        df_chosen = df
    return df_chosen.to_dict(orient='records')

def get_list_metrics(chosen_list=None):
    df = pd.read_json("3_metrics_log.ndjson",lines=True, dtype={'ID_metric': str})
    if chosen_list is not None:
        df_chosen = df[df['ID_metric'].isin(chosen_list)]
    else:
        df_chosen = df
    return df_chosen.to_dict(orient='records')


# Run the dataset, method, and metric if id matches
def run_dataset_method_metric(datasets, methods, metrics):
    need_run = []
    log = pd.read_csv("main_result.txt", sep='\t',dtype={'ID_dataset': str, 'ID_method': str, 'ID_metric': str})
    
    for dataset in datasets:
        ID_dataset = str(dataset["ID_dataset"])
        for method in methods:
            ID_method = str(method["ID_method"])
            # result_filename = "results/" + ID_dataset + "_" + ID_method + "_result.tsv"
            for metric in metrics:
                ID_metric = str(metric["ID_metric"])
                lookup = log[(log['ID_dataset'] == ID_dataset) & (log['ID_method'] == ID_method) & (log['ID_metric'] == ID_metric)]
                # run the method
                if lookup.empty:
                    need_run.append([dataset, method, metric])
                    if len(need_run) < 30:
                        print(f"{ID_dataset}_{ID_method}_{ID_metric}",end=", ")
                #else
                    #print(f"Dataset {ID_dataset} with method {ID_method} and metric {ID_metric} exists, skipping run")
    return need_run


def parse_args():
    parser = argparse.ArgumentParser(description='Process datasets and methods.')
    parser.add_argument('--datasets', type=str, help='Dataset indices (e.g., 1,2,3 or 1:5)')
    parser.add_argument('--methods', type=str, help='Method indices (e.g., 1,2,3 or 1:5)')
    parser.add_argument('--metrics', type=str, help='Metric indices (e.g., 1,2,3 or 1:5)')
    return parser.parse_args()


def get_chosen_arguments(arg_str):
    """
    Parse argument string into a list of integers.
    Handles both comma-separated values and range notation.
    
    Args:
        arg_str: String like "1,2,3" or "1:5" or "36-1:10,37-1:10", but not "36:38,39-1:10
    Returns:
        List of integers or None if arg_str is None
    """
    if arg_str is None:
        return None
    
    try:
        # Handle range notation
        if '-' in arg_str: # hyphenated for repeats of one method, treat as strings
            if ':' in arg_str: # 36-1:10,37-1:10
                list_parts = arg_str.split(',')
                result = []
                for part in list_parts:
                    if ':' not in part:
                        result.append(part)
                        continue
                    prefix, range_part = part.split('-')
                    start, end = map(int, range_part.split(':'))
                    result.extend([f"{prefix}-{i}" for i in range(start, end + 1)])
                return result
            else: # single hyphenated value
                return [str(x) for x in arg_str.split(',')]
        elif ':' in arg_str: # only range notation
            start, end = map(int, arg_str.split(':'))
            range_list = [str(x) for x in range(start, end + 1)]
            return range_list
        else: # simple comma-separated values
            return [str(x) for x in arg_str.split(',')]
    except (ValueError, SyntaxError) as e:
        print(f"Error parsing argument {arg_str}: {e}")
        sys.exit(1)


def main():
    args = parse_args()

    chosen_datasets = get_chosen_arguments(args.datasets)
    chosen_methods = get_chosen_arguments(args.methods)
    chosen_metrics = get_chosen_arguments(args.metrics)

    # Print results for verification
    print(f"Chosen datasets: {chosen_datasets}")
    print(f"Chosen methods: {chosen_methods}")
    print(f"Chosen metrics: {chosen_metrics}") 
    
    run_list = run_dataset_method_metric(get_list_data_files(chosen_datasets), get_list_methods(chosen_methods), get_list_metrics(chosen_metrics))

    cp.dump(run_list, open("4_run_list", "wb"))
    print(f"\nConstructed run list: {len(run_list)} jobs to run")

if __name__ == "__main__":
    main()

