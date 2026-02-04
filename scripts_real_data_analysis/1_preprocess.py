import ast
import sys
import pickle as cp
import argparse
import os
import time
import obonet
import Bio
from Bio.UniProt import GOA
from datetime import datetime
import numpy as np
import pandas as pd
import ia as ia

class annotation_data:
    ##############################
    ###  Initializing Methods  ###
    ##############################
    
    def __init__(self, lower_cutoff, upper_cutoff,
                gene_file, annotation_file, obo, identifier='test', annt_type='gaf'):
        self.lower_cutoff = int(lower_cutoff)
        self.upper_cutoff = int(upper_cutoff)
        self.gene_file = gene_file
        self.annotation_file = annotation_file
        self.identifier = str(identifier)
        self.annt_type = annt_type
        self.obo = obo
        
        
    def get_expression_data(self):
        """
        Reads gene expression data from a file and stores it in the object.
        Returns:
            None
        """
        data = pd.read_csv(self.gene_file, sep='\t', header=0)
        if data.shape[1] != 2:
            if 'external_gene_name' in data.columns:
                data = data.loc[:, ['external_gene_name', "log2FoldChange"]]
            elif 'GeneSymbol' in data.columns:
                data = data.loc[:, ['GeneSymbol', "log2FoldChange"]]
            # data = data.loc[:, ["external_gene_name", "log2FoldChange"]]
        self.genes_list = np.array(data.iloc[:, 0])
        
        
    def combine_equal_terms(self,verbose=0): 
        #print "merge equal GO terms into one..."
        nr_terms = len(self.term_names)
        
        ind=self.n_genes_per_term.argsort()
        cc_ord=self.n_genes_per_term[ind]
        keep_these = np.zeros(0,dtype=int)
        delete_these = np.zeros(0,dtype=int)
        for i in range(nr_terms):
            for j in range(i+1,nr_terms):
                if cc_ord[i]<cc_ord[j]:
                    break
                elif set(self.T[ind[i]])==set(self.T[ind[j]]):
                    keep_these = np.append(keep_these,ind[i])
                    delete_these = np.append(delete_these,ind[j])
        
        term_should_be_kept = np.ones(nr_terms,dtype=bool)
        term_should_be_kept[delete_these] = False
        
        for i in range(keep_these.shape[0]):
            self.term_names[keep_these[i]]+=' & '+self.term_names[delete_these[i]]
        
        self.T = self.T[term_should_be_kept]
        self.term_names = self.term_names[term_should_be_kept]
        self.n_genes_per_term = self.n_genes_per_term[term_should_be_kept]
        
        # should order the terms in the combined term to ensure reproducibility
        for i, name in enumerate(self.term_names):
            if name.find('&')!=-1:
                names = name.split(' & ')
                sorted_names = sorted(names)
                new_name = ' & '.join(sorted_names)
                self.term_names[i] = new_name
                    
        if verbose>0:
            print("New combined nodes:\n")
            for name in self.term_names:
                if name.find('&')!=-1:
                    print(name)
    
    def delete_too_large_and_too_small_terms(self):
        if self.upper_cutoff>=self.lower_cutoff:
            indices_to_keep = np.bitwise_and(self.n_genes_per_term>=self.lower_cutoff,self.n_genes_per_term<=self.upper_cutoff)
        else:
            indices_to_keep = self.n_genes_per_term>=self.lower_cutoff
        self.T = self.T[indices_to_keep]
        self.term_names = self.term_names[indices_to_keep]
        self.n_genes_per_term = self.n_genes_per_term[indices_to_keep]
    
    def remove_genes_not_annotated_to_any_term(self,verbose=0):
        toberemoved=list(set(range(self.n_genes))-set.union(*[set([])] + [set(t) for t in self.T]))
        if toberemoved==[]:
            return
        toberemoved = np.array(toberemoved,dtype=int)
        delete_these_genes = np.zeros(self.n_genes,dtype=bool)
        delete_these_genes[toberemoved] = True
        if verbose>0:
            print("Deleted genes: "+', '.join(self.unique_genes[delete_these_genes]))
        cumsum = np.cumsum(delete_these_genes)
        for j in range(len(self.T)):
            for i in range(len(self.T[j])):
                self.T[j][i]-=cumsum[self.T[j][i]]
        
        self.unique_genes = self.unique_genes[~delete_these_genes]
        self.n_genes -= len(toberemoved)
    
    def implement_cutoffs(self, verbose=0):
        '''
        Optional function for size filtering step
        Output is cut data-specfic annotation file
        '''
        # remove terms that have too many or too little annotations
        self.delete_too_large_and_too_small_terms()
        
        # remove genes that are no longer annotated by any term
        self.remove_genes_not_annotated_to_any_term(verbose)
        
        
    def create_annotation_data(self,verbose=0):
        """creates all the important variables from the annotation_file."""
        if self.annt_type=="tsv":
            data = pd.read_csv(self.annotation_file,sep='\t',header=0)
            if data.shape[1]!=2:
                raise ValueError("TSV annotation file must have exactly two columns: EntryID and term")
            else:
                data.columns=["EntryID","term"]
        else:
            dt = self.annotation_file.split('/')[-1].split('goa')[0]
            if not os.path.exists('saved_data/save_'+dt+'_ancestors'):
                print("Propagating from GAF file")
                gaf = GOA.gafiterator(open(self.annotation_file,'r'))
                # with open(self.annotation_file, 'r') as f:
                #     for line in f: 
                #         if line.startswith('!date-generated:'):
                #             datestr = line.strip().split(': ')[1]
                #             print(datestr)
                #             dt = datetime.strptime(datestr, '%Y-%m-%dT%H:%M')
                #             dt = dt.strftime('%b%y')
                #             break
                print("GAF date generated:", dt)
                
                annotation_df = pd.DataFrame.from_dict(gaf, orient="columns").loc[:,["DB_Object_Symbol","GO_ID","Aspect"]]
                annotation_df.columns=["EntryID","term","aspect"]
                annotation_df=annotation_df.loc[annotation_df["aspect"] == "P",:]
                aspect_map={"P":"BPO"}
                annotation_df["aspect"] = annotation_df["aspect"].apply(lambda x: aspect_map[x])
                
                # Propagate gene ontology
                ontology_graph = ia.clean_ontology_edges(obonet.read_obo(self.obo))
                roots = {'BPO': 'GO:0008150'}
                subontologies = {aspect: ia.fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots} 
                data, ancestor_lookup = ia.propagate_terms(annotation_df, subontologies)
                if not os.path.exists('saved_data'):
                    os.makedirs('saved_data')
                cp.dump((data, ancestor_lookup), open('saved_data/save_'+dt+'_ancestors', 'wb')) # not dataset-specific
            else:
                print("Loading propagated data")
                data, ancestor_lookup = cp.load(open('saved_data/save_'+dt+'_ancestors', 'rb'))
                # Note: saved ancestor_lookup are single terms, not combined

        
        # Group by term, get list of genes for each term stored in list of lists, only for genes in dict_genes_list
        data_filtered = data[data['EntryID'].isin(self.genes_list)]
        data_by_term = data_filtered[['term','EntryID']].groupby('term')['EntryID'].apply(lambda x: ' '.join(x))
        genes_by_term = data_by_term.tolist()
        self.term_names = np.array(data_by_term.index.drop_duplicates())
        self.unique_genes = np.array(list(set(data_filtered['EntryID'])))
        
        T = []
        for genes in genes_by_term:
            T.append([gene for gene in genes.strip().split()])
        self.T=np.array(T,dtype=object)
        
        self.n_genes_per_term = np.array(list(map(len,self.T)))
        self.n_genes = len(self.unique_genes)
        
        # combine terms that share gene annotations
        self.combine_equal_terms(verbose)
        
    def filter_existing_annotation_data(self,verbose=0):
        """If annotation data already in ['term', 'T'] format, load it and filter out all genes not in self.genes_list"""
        
        annotation_df = pd.read_csv(self.annotation_file, sep='\t',header=0)
        annotation_df.columns=["term","T"]
        annotation_df['T'] = annotation_df['T'].apply(lambda x: ast.literal_eval(x))
        
        # Keep only genes in self.genes_list
        data = annotation_df.explode('T')
        data = data[data['T'].isin(self.genes_list)]
        data_by_term = data.groupby('term')['T'].apply(lambda x: ' '.join(x))
        genes_by_term = data_by_term.tolist()
        self.term_names = np.array(data_by_term.index.drop_duplicates())
        self.unique_genes = np.array(list(set(data['T'])))
       
        T = []
        for genes in genes_by_term:
            T.append([gene for gene in genes.strip().split()])
        self.T=np.array(T,dtype=object)
        
        self.n_genes_per_term = np.array(list(map(len,self.T)))
        self.n_genes = len(self.unique_genes)
        
        # combine terms that share gene annotations
        self.combine_equal_terms(verbose)
        
        
    def store_files(self):
        '''
        Store annotation files after modification and cutoffs implemented
        '''
        df = pd.DataFrame({'term': self.term_names, 'T': self.T})
        # write a tsv file for cut data-specific annotation data
        if not os.path.exists('saved_data'):
            os.makedirs('saved_data')
        if df.shape[0]>0:
            df.to_csv('saved_data/save_'+self.identifier+'_'+str(self.lower_cutoff)+'to'
                +str(self.upper_cutoff)+'_annotations.tsv', sep='\t', index=False)
        else:
            print("No terms to save")
    
    #####################################
    ###  Main Function of this Class  ###
    #####################################
    
    def runMe(self, verbose=False):
        start=time.time()
        
        # load the expression data
        self.get_expression_data()
        
        if verbose:
            print('Checkpoint A, passed time',time.time()-start)
        
        if self.annt_type != "tsv_dict":
            # load/create the annotation data, initialize self.term_names, self.unique_genes, self.T, self.n_genes_per_term
            self.create_annotation_data(verbose)
        else:
            # load existing annotation data in ['term', 'T'] format, filter out genes not in self.genes_list
            self.filter_existing_annotation_data(verbose)
        
        if verbose:
            print('Checkpoint B, passed time',time.time()-start)
        
        # implement cutoffs if lower_cutoff or upper_cutoff are not 0 or 100000
        if self.lower_cutoff!=0 or self.upper_cutoff!=100000:
            self.implement_cutoffs(verbose)
            
        self.store_files()
    

def parse_args():
    parser = argparse.ArgumentParser(prog="preprocess.py", description="Pre-process dataset-specific annotation file",add_help=True)
    required_arguments = parser.add_argument_group("Required arguments")
    required_arguments.add_argument('--gene_file', '-g', required=True,
                                    help="Path to a tab-separated file of a list of genes ranked by any metric")
    required_arguments.add_argument('--annt', '-a',
                                    help="Path to a species-specific GAF file or a tab-separated file with GO terms and genes annotated to them.")
    required_arguments.add_argument('--identifier', '-i',
                                    help="Identifier for dataset") 
    required_arguments.add_argument('--annt_type', '-t', choices=['gaf','tsv','tsv_dict'], default='gaf',
                                help="Type of annotation file, either gaf or tsv") 
    # Add argument for implementing cutoff
    optional_arguments = parser.add_argument_group("Optional arguments")
    optional_arguments.add_argument('--lower_cutoff','-l',default=0,
                                    help="Minimum number of genes in term")
    optional_arguments.add_argument('--upper_cutoff','-u',default=100000,
                                    help="Maximum number of genes in term")

    # Add optional argument if type of annotation file is gaf
    optional_arguments.add_argument('--obo', '-o', 
                                    help='Path to OBO ontology graph file if GAF was provided.')                                    
    args = parser.parse_args()
    return args


if __name__=="__main__":
    args=parse_args()
    if args.annt_type=="gaf" and args.obo is None:
        raise ValueError("If annotation file type is gaf, OBO file must be provided")
    m=annotation_data(args.lower_cutoff,args.upper_cutoff,args.gene_file,
                      args.annt,args.obo,args.identifier,args.annt_type)
    filename = 'saved_data/save_'+args.identifier+'_'+str(args.lower_cutoff)+'to'+str(args.upper_cutoff)+'_annotations.tsv'
    if not os.path.exists(filename): 
        m.runMe(verbose=0)
    else:
        print("File already exists, not overwriting")
