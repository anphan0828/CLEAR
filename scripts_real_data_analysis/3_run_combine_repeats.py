"""
This script is used to combine results from multiple repeats of CLEAR.
It reads in the result files from multiple repeats (with the same dataset and method), combines them,
and computes performance metrics such as phenotype precision-recall AUC, average term size, and average overlap coefficient.
The script should be as a standalone script (after all dataset-method repeat runs are done).
After running this script, the results will be appended to 'main_result.txt' with the same ID_metric as in 3_metrics_log.ndjson,
and method-repeat suffix removed.
"""

import os
import sys
import ast
import argparse
import scipy.stats as stats
import numpy as np
import pandas as pd
import pickle as cp
import random
from sklearn.metrics import auc
from matplotlib import pyplot as plt
import re
from sklearn.decomposition import PCA

def extract_disease_code_from_dataset(datasetid_file):
    id_to_code = {}
    # id_to_name = {}
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
                # id_to_name[gid] = name
            else:
                # fallback: try to find last [CODE] and first token
                if '[' in line and ']' in line:
                    try:
                        code = line[line.rfind('[')+1:line.rfind(']')].strip()
                        gid = line.split(None, 1)[0]
                        # disease name is everything between first whitespace and the bracket
                        name = line[len(gid):line.rfind('[')].strip()
                        id_to_code[gid] = code
                        # id_to_name[gid] = name
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

def combine_repeats(input_files):
    combined_table = pd.DataFrame()
    for repeat, input_file in enumerate(input_files):
        table = pd.read_csv(input_file, header=0,sep="\t")
        if 'pvalue' in table.columns:
            table.columns = ["term", "value", "pvalue"]
            table = table.sort_values(by=['pvalue','term'], ascending=[True,True])
            higher_is_better = False
        else:
            table.columns = ["term","value"]
            table = table.sort_values(by=['value','term'], ascending=[False,True]) # sorting by on_frequency and estimate
            higher_is_better = True
        if combined_table.empty:
            combined_table = table
        else:
            combined_table = pd.merge(combined_table, table, on="term", how="outer", suffixes=('', f'_{repeat}'))
    if higher_is_better:
        combined_table['sum_value'] = combined_table[[col for col in combined_table.columns if col.startswith('value')]].sum(axis=1)
    else:
        combined_table['sum_value'] = combined_table[[col for col in combined_table.columns if col.startswith('pvalue')]].sum(axis=1)
    return combined_table, higher_is_better

def pheno_pr(combined_table, higher_is_better, phenotype_df):    
    # Retrieve all terms in combined term, compare with phenotype_df and keep only those that are in phenotype_df['GO_ID']
    valid_terms = set(phenotype_df['GO_ID'])
    def find_matching_term(x):
        terms = x.split(" & ")
        matching_terms = [t for t in terms if t in valid_terms]
        return matching_terms[0] if matching_terms else terms[0]
    combined_table["term"] = combined_table["term"].apply(find_matching_term)
    N_GS = combined_table.shape[0]
    valid_phenotype_terms = set.intersection(set(combined_table['term']), set(phenotype_df['GO_ID']))
    prevalence = len(valid_phenotype_terms) / len(combined_table)
    
    # Calculate sum of values across repeats
    # combined_table['sum_value'] = combined_table[[col for col in combined_table.columns if col.startswith('value')]].sum(axis=1)
    combined_table['true_value'] = combined_table['term'].apply(lambda x: 1 if x in valid_phenotype_terms else 0)    
    
    # Return result across repeats
    precision, recall = compute_precision_recall(combined_table['true_value'], combined_table['sum_value'], higher_is_better=higher_is_better)
    pr_auc = auc(recall, precision)
    res = {"n_common": len(valid_phenotype_terms), "n_total": int(N_GS),
           "pheno_pr_auc": round(pr_auc, 6), "pheno_prevalence": round(prevalence, 4)}
    # Plotting all repeats (different from ranking by sum_value)
    # pr_auc_repeats = []
    # for col in combined_table.columns:
    #     if col.startswith('value'):
    #         col_scores = np.array(combined_table[col], dtype=float)
    #         if col.startswith('value'):
    #             col_higher_is_better = True
    #         else:
    #             col_higher_is_better = False
    #         precision, recall = compute_precision_recall(combined_table['true_value'], col_scores, higher_is_better=col_higher_is_better)
    #         pheno_pr_auc = auc(recall, precision)
    #         pr_auc_repeats.append(pheno_pr_auc)
    #         plt.plot(recall, precision, label=f'{col} (AUC={pheno_pr_auc:.4f})')
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.title('Precision-Recall Curve for Combined Repeats')
    # plt.legend()
    # plt.savefig('combined_repeats_pr_curve.png')
    # plt.close()
    return res

def calculate_termsize(combined_table, annotation_dict):
    names = np.array(combined_table['term'])
    T = annotation_dict
    term_sizes = []
    for term in names:
        term_size = len(T[term])
        # print(f"Term: {term}, Size: {term_size}")
        term_sizes.append(term_size)
    # print(term_sizes)
    result = np.mean(term_sizes)
    return {"score_termsize": round(result, 4)}

def calculate_overlap(combined_table, annotation_dict):
    names = np.array(combined_table['term'])
    T = annotation_dict
    overlap_coefficients = []
    for i, term in enumerate(names):
        for j, next_term in enumerate(names):
            if i != j:
                overlap = len(set(T[term]) & set(T[next_term])) / min(len(set(T[term])),len(set(T[next_term])))
                overlap_coefficients.append(overlap)
    res = np.mean(overlap_coefficients) if overlap_coefficients else 0.0
    return {"score_overlap": round(res, 4)}

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
       
# For running interactively (reading from 4_run_list)
# python3 4_prerun.py --datasets 1:39 --methods 63-1,64-1,65-1,66-1 --metrics 4
datasetid_file = "ext_data/GseId2Disease.txt"
microarray_mapping = extract_disease_code_from_dataset(datasetid_file)
run_list = cp.load(open("4_run_list", "rb"))
ID_pheno_metric = 1 # phenotype relevance metric ID from 3_metrics_log.ndjson
with open('main_result.txt', 'a') as f:
    for each in run_list:
        (dataset, method, metric) = each
        ID_dataset = str(dataset["ID_dataset"])
        gene_file = str(dataset["gene_list"])
        ID_method = str(method["ID_method"]).split("-")[0]  # remove repeat suffix
        cancer_id = gene_file.split("/")[-1].split(".")[0].split("_de_limma")[0]
        input_files = [os.path.join('results', x) for x in os.listdir('results') if x.startswith(f'{ID_dataset}_{ID_method}-') and x.endswith('_result.tsv')]
        phenotype_df = extract_phenotype_rankings("ext_data/GO_BP.cp", cancer_id, microarray_mapping=microarray_mapping)
        if phenotype_df is not None:
            combined_table, higher_is_better = combine_repeats(input_files)
            res = pheno_pr(combined_table, higher_is_better, phenotype_df)
            print(f'{ID_dataset}_{ID_method}: {res}')
            f.write(f'{ID_dataset}\t{ID_method}\t{ID_pheno_metric}\t{res}\n')
            
num=20
with open('main_result.txt', 'a') as f:
    for each in run_list:
        (dataset, method, metric) = each
        ID_dataset = str(dataset["ID_dataset"])
        gene_file = str(dataset["gene_list"])
        annotation_file = 'saved_data/' + str(dataset["annotations"])
        ID_method = str(method["ID_method"]).split("-")[0]  # remove repeat suffix
        cancer_id = gene_file.split("/")[-1].split(".")[0].split("_de_limma")[0]
        input_files = [os.path.join('results', x) for x in os.listdir('results') if x.startswith(f'{ID_dataset}_{ID_method}-') and x.endswith('_result.tsv')]
        
        # Resolve ties and get chosen terms
        combined_table, higher_is_better = combine_repeats(input_files)
        if higher_is_better:
            ties = combined_table.groupby('sum_value').apply(lambda x: len(x)).sort_index(ascending=False).to_dict()
            table_score = combined_table[['term','sum_value']].rename(columns={"sum_value":"score"})
        else:
            ties = combined_table.groupby('sum_value').apply(lambda x: len(x)).sort_index(ascending=True).to_dict()
            table_score = combined_table[['term','sum_value']].rename(columns={"sum_value":"score"})
        chosen_score = pd.DataFrame(columns=table_score.columns)
        for score, count in ties.items():
            if count >  num - len(chosen_score):
                chosen_score = pd.concat([chosen_score, table_score[table_score['score'] == score].sample(num-len(chosen_score),replace=False, random_state=42)])
            else:
                chosen_score = pd.concat([chosen_score, table_score[table_score['score'] == score]])
        
        # Get annotation_dict
        annotation_df = pd.read_csv(annotation_file,sep="\t",header=0)
        annotation_df.columns=["term","T"]
        annotation_df['T'] = annotation_df['T'].apply(lambda x: ast.literal_eval(x))
        dict_T = annotation_df.set_index('term')['T'].to_dict()
        
        ID_metric_termsize = 3 # term size metric ID from 3_metrics_log.ndjson
        res_termsize = calculate_termsize(chosen_score, dict_T)
        print(f'{ID_dataset}_{ID_method}: {res_termsize}')
        f.write(f'{ID_dataset}\t{ID_method}\t{ID_metric_termsize}\t{res_termsize}\n')
        
        ID_metric_overlap = 2 # overlap metric ID from 3_metrics_log.ndjson
        res_overlap = calculate_overlap(chosen_score, dict_T)
        print(f'{ID_dataset}_{ID_method}: {res_overlap}')
        f.write(f'{ID_dataset}\t{ID_method}\t{ID_metric_overlap}\t{res_overlap}\n')

# Plotting PCA for repeats
fig, ax = plt.subplots(ncols=4, nrows=6, figsize=(20, 30))        
idx=0
for each in run_list:
    (dataset, method, metric) = each
    ID_dataset = str(dataset["ID_dataset"])
    gene_file = str(dataset["gene_list"])
    annotation_file = '/work/idoerg/ahphan/benchmark_GSEA/saved_data/' + str(dataset["annotations"])
    ID_method = str(method["ID_method"]).split("-")[0]  # remove repeat suffix
    if ID_method != '65':
        continue
    idx+=1
    cancer_id = gene_file.split("/")[-1].split(".")[0].split("_de_limma")[0]
    input_files = [os.path.join('/work/idoerg/ahphan/benchmark_GSEA/results', x) for x in os.listdir('/work/idoerg/ahphan/benchmark_GSEA/results') if x.startswith(f'{ID_dataset}_{ID_method}-') and x.endswith('_result.tsv')]
    
    # Resolve ties and get chosen terms
    combined_table, higher_is_better = combine_repeats(input_files)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(combined_table[[col for col in combined_table.columns if col.startswith('value')]].T)
    df_pca = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
    
    ax_idx = (idx-1)//4, (idx-1)%4
    ax[ax_idx].scatter(df_pca['PC1'], df_pca['PC2'], color='blue', label='All gene sets', alpha=0.7)
    
    # Get rows where at least one value* column is >= 0.5
    filtered_df = combined_table[combined_table[[col for col in combined_table.columns if col.startswith('value')]].ge(0.5).any(axis=1)] 
    if len(filtered_df) >= 2:
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(filtered_df[[col for col in filtered_df.columns if col.startswith('value')]].T)
        filtered_df_pca = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
        
        ax[ax_idx].scatter(filtered_df_pca['PC1'], filtered_df_pca['PC2'], color='red', label='Gene sets >=0.5', alpha=0.7)
    ax[ax_idx].set_title(f'{cancer_id}_{ID_method}')
    ax[ax_idx].set_xlabel('PC1')
    ax[ax_idx].set_ylabel('PC2')

plt.tight_layout()
plt.legend()
plt.savefig('combined_repeats_pca_gamma.png')