import sys
import ast
import argparse
import scipy.stats as stats
import numpy as np
import pandas as pd
import pickle as cp
from sklearn.metrics import auc
from matplotlib import pyplot as plt
import re


def extract_disease_code_from_dataset(datasetid_file):
    id_to_code = {}
    line_re = re.compile(r'^(\S+)\s+(.*?)\s*\[([^\]]+)\]\s*$')  # group1=id, group2=full name, group3=code
    with open(datasetid_file, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            m = line_re.match(line)
            if m:
                gid, name, code = m.group(1), m.group(2).strip(), m.group(3).strip()
                id_to_code[gid] = code
            else:
                # fallback: try to find last [CODE] and first token
                if '[' in line and ']' in line:
                    try:
                        code = line[line.rfind('[')+1:line.rfind(']')].strip()
                        gid = line.split(None, 1)[0]
                        # disease name is everything between first whitespace and the bracket
                        name = line[len(gid):line.rfind('[')].strip()
                        id_to_code[gid] = code
                    except Exception:
                        # skip malformed line
                        continue
    return id_to_code
                
def extract_phenotype_rankings(pickle_file_path, cancer_id, microarray_mapping):
    with open(pickle_file_path, 'rb') as f:
        data_dict = cp.load(f)
    if cancer_id in data_dict.keys():
        phenotype_df = data_dict[cancer_id]
        phenotype_df.columns = ['GO_ID','GO_name','REL_SCORE','MATCHED_GENES','TOTAL_GENES']
        return phenotype_df
    elif cancer_id in microarray_mapping:
        phenotype_df = data_dict[microarray_mapping[cancer_id]]
        phenotype_df.columns = ['GO_ID','GO_name','REL_SCORE','MATCHED_GENES','TOTAL_GENES']
        return phenotype_df
    else:
        print("Cancer ID not found in the phenotype relevance file")
        return None   
    
def calculate_phenotype_relevance(input_file, phenotype_df):
    """
    Compute PR-AUC with phenotype relevance gene sets as positives
    
    :param input_file: Path to enriched terms file
    :param phenotype_df: DataFrame with phenotype relevance scores
    :return: Dictionary with relevance ratio, p-value, n_common, n_total
    """
    try:
        table = pd.read_csv(input_file, header=0,sep="\t")
        if 'pvalue' in table.columns:
            table.columns = ["term", "value", "pvalue"]
            table = table.sort_values(by=['pvalue','term'], ascending=[True,True])
            higher_is_better = False
        else:
            table.columns = ["term","value"]
            table = table.sort_values(by=['value','term'], ascending=[False,True]) # sorting by on_frequency and estimate
            higher_is_better = True
    except pd.errors.EmptyDataError:
        return {"error_fewterms": 0}
    
    # Retrieve all terms in combined term, compare with phenotype_df and keep only those that are in phenotype_df['GO_ID']
    valid_terms = set(phenotype_df['GO_ID'])
    def find_matching_term(x):
        terms = x.split(" & ")
        matching_terms = [t for t in terms if t in valid_terms]
        return matching_terms[0] if matching_terms else terms[0]
    table["term"] = table["term"].apply(find_matching_term)
    
    terms = np.array(table['term'])
    if higher_is_better:
        scores = np.array(table['value'], dtype=float)
    else:
        scores = np.array(table['pvalue'], dtype=float)
    
    valid_phenotype_terms = set.intersection(set(table['term']), set(phenotype_df['GO_ID']))
    N_GS = table.shape[0]
    
    # Precision-recall
    y_true = [1 if x in valid_phenotype_terms else 0 for x in terms]
    precision, recall = compute_precision_recall(y_true, scores, higher_is_better=higher_is_better)
    pheno_pr_auc = auc(recall, precision)
    prevalence = len(valid_phenotype_terms) / N_GS

    common = pd.merge(table, phenotype_df, left_on='term', right_on='GO_ID', how='inner')
    # Filter common terms with score > 0 or p-value < 1 (if not filtered, every method return all N_GS terms)
    if higher_is_better:
        common = common[common['value'] > 0]
    else:
        common = common[common['pvalue'] < 1]
    # res = {"n_common": len(valid_phenotype_terms), "n_common_selected": int(common.shape[0]), "n_total": int(N_GS),
    #        "max_relevance": round(common['REL_SCORE'].max(), 6), "avg_relevance": round(common['REL_SCORE'].mean(), 6), 
    #        "relevance_ratio": round(relevance_ratio, 6),
    #        "pheno_pr_auc": round(pheno_pr_auc, 6), "pheno_prevalence": round(prevalence, 4)}
    res = {"n_common": len(valid_phenotype_terms), "n_total": int(N_GS),
           "pheno_pr_auc": round(pheno_pr_auc, 6), "pheno_prevalence": round(prevalence, 4)}
    return res


def compute_precision_recall(y_true, y_scores, higher_is_better=True):
    precision = []
    recall = []
    thresholds = sorted(set(y_scores), reverse=higher_is_better) # sorted in ascending order by default, reverse=True when scores are posterior probs
    # plotting from higher to lower scores (when reverse=True) and lower to higher p-values (when reverse=False)
    for idx, thresh in enumerate(thresholds): # iterate through unique score values from high to low
        y_pred = [1 if (score >= thresh and higher_is_better==True) or (score <= thresh and higher_is_better==False) else 0 for score in y_scores]
        tp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
        fp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
        fn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)
        
        prec = tp / (tp + fp) if (tp + fp) > 0 else 1.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        if idx == 0 and rec != 0.0: # idx=0 is always the best score threshold (or lowest p-value)
            recall.append(0.0)
            precision.append(prec) # project to y-axis for left-boundary
        precision.append(prec)
        recall.append(rec)
    return precision, recall

def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Compute area under the precision-recall curve using phenotype relevance gene sets as positives.')
    
    parser.add_argument('--annot', '-a', required=True, 
                        help='Path to annotation file')
    parser.add_argument('--gene_file', '-g', required=True,
                        help='Path to gene file')
    parser.add_argument('--phenorelevance', '-p', required=True,
                        help="Path to dataset-specific phenotype relevance file")
    parser.add_argument('--input', '-i', required=True,
                        help='Path to list of enriched terms')
    return parser.parse_args(argv)


if __name__== "__main__":
    args=parse_inputs(sys.argv[1:])
    cancer_id = args.gene_file.split("/")[-1].split(".")[0].split("_de_limma")[0]
    datasetid_file = "GseId2Disease.txt" # for mapping GEO dataset IDs to disease codes
    microarray_mapping = extract_disease_code_from_dataset(datasetid_file)
    phenotype_df = extract_phenotype_rankings(args.phenorelevance, cancer_id, microarray_mapping=microarray_mapping)
    res = calculate_phenotype_relevance(args.input, phenotype_df)
    print(res)
