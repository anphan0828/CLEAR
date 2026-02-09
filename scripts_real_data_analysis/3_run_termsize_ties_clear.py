import sys
import ast
import argparse
import numpy as np
import pandas as pd


def average_termsize(table, annotation_dict):
    names = np.array(table['term'])
    T = annotation_dict
    term_sizes = []
    for term in names:
        term_size = len(T[term])
        term_sizes.append(term_size)
    result = np.mean(term_sizes)
    return result


def calculate_from_input(input,annotations,num):
    try:
        table = pd.read_csv(input, header=0,sep="\t")
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
    
    # If there are ties in table_nes or table_pval that is more than num, we sample from the ties
    if higher_is_better:
        ties = table.groupby('value').apply(lambda x: len(x)).sort_index(ascending=False).to_dict()
        table_score = table[['term','value']].rename(columns={"value":"score"})
    else:
        ties = table.groupby('pvalue').apply(lambda x: len(x)).sort_index(ascending=True).to_dict()
        table_score = table[['term','pvalue']].rename(columns={"pvalue":"score"})
    chosen_score = pd.DataFrame(columns=table_score.columns)
    for score, count in ties.items():
        if count >  num - len(chosen_score):
            chosen_score = pd.concat([chosen_score, table_score[table_score['score'] == score].sample(num-len(chosen_score),replace=False, random_state=42)])
        else:
            chosen_score = pd.concat([chosen_score, table_score[table_score['score'] == score]])
            
    annotation_df = pd.read_csv(annotations,sep="\t",header=0)
    annotation_df.columns=["term","T"]
    annotation_df['T'] = annotation_df['T'].apply(lambda x: ast.literal_eval(x))
    
    dict_T = annotation_df.set_index('term')['T'].to_dict()

    score_termsize = average_termsize(chosen_score, dict_T)
    res = {"score_termsize": round(score_termsize,4)}
    return res
        
        
def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description='Compute average gene set size of top n terms using p-value ranking or posterior probability')    
    
    parser.add_argument('--annotations', '-a', required=True,
                        help="Path to dataset-specific annotation file")
    parser.add_argument('--input', '-i', required=True,
                        help='Path to list of enriched terms')
    parser.add_argument('--num', '-n', default=20,
                        help="Top n terms included in the metric")
    return parser.parse_args(argv)


if __name__== "__main__":
    args = parse_inputs(sys.argv[1:])
    value = calculate_from_input(args.input,args.annotations, int(args.num))
    print(value)

