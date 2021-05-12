"""
Modules object - essentially stores info and houses useful functions for accessing the data in different formats
"""

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from sklearn.metrics.pairwise import cosine_similarity
import pickle as pkl
from os import path
from collections import defaultdict
import itertools

class Modules(object):
    def __init__(self, calc_modules=False, num_samples=1000):
        self.nash = 'Nonalcoholic Steatohepatitis'

        # all modules = all non disease modules + nash module
        # columns: gene, module
        self.all_mods_genes = pd.read_csv(r'../data/gene_maps/all_modules.csv')
        self.gene_embeddings = pd.read_csv(r'../data/gene_maps/embedding.csv', index_col=0)
        # dicts
        self.module_to_genes, self.gene_to_modules = self.map_modules_genes()
        
        # module embeddings
        if calc_modules:
            self.module_vectors = self.get_module_vectors()
            self.module_similarities = self.get_module_similarities()
        else:
            self.module_vectors = pd.read_csv('../data/module_vectors_summed.csv', index_col=0)
            self.module_similarities = pd.read_csv('../data/module_similarities.csv', index_col=0)

        self.num_samples = num_samples
        self.similarity_dists = self.get_similarity_distributions() # key is a sorted (ascending) tuple
        self.pvalues = pd.read_csv(r'../data/similarity_pvalues_sum.csv', index_col=0)
        # TODO: integrate pvalue functions into this and re-run if you have time
        # TODO: calculate individual gene effects
        self.num_modules = len(self.module_vectors.index)

        self.curated_genes = sorted(self.module_to_genes[self.nash])
        self.befree_genes = sorted(self.load_befree_genes(0))
        self.sven_genes = sorted(self.load_svensson_genes(0))


    # Function to create dicts to access genes in a module, module labels
    # Gets called on initialization to store dicts for quicker use
    ##########
    def map_modules_genes(self):
        """
        Return dicts:
         1. map module name to list of genes in module {module: [gene, ...]}
         2. gene name to list of modules gene is in {gene: [module, ...]}
        :return:
        """
        only_embedded = self.all_mods_genes[self.all_mods_genes['gene'].isin(self.gene_embeddings.index)]
        gene_modules = dict(only_embedded.groupby('gene')['module'].apply(list))
        module_genes = dict(only_embedded.groupby('module')['gene'].apply(list))
        return module_genes, gene_modules

    # Functions to get similarities between modules
    ######
    def sum_embeddings(self, mod_genes):
        """
    	Takes a list of genes and returns a summed vector to represent the genes in the module
    	"""
        mod_embeds = self.gene_embeddings.loc[mod_genes,:].to_numpy()
        if len(mod_genes) == 1:
            return mod_embeds.flatten()
        return np.sum(mod_embeds, axis=0)

    def get_module_vectors(self, summed=True):
        """
        Returns a dataframe with a vector for each module in the "all modules" list
        :param summed: True if using summed vectors, false if using averaged vectors
        :return: module_vectors
        """
        module_vecs= pd.DataFrame()
        module_names = list(self.module_to_genes.keys())
        for mod_name in module_names:
            mod_genes = list(set(self.module_to_genes[mod_name]) & set(self.gene_embeddings.index))
            mod_vector = self.sum_embeddings(mod_genes)
            module_vecs[mod_name] = mod_vector
        module_vecs = module_vecs.T
        module_vecs.to_csv('../data/module_vectors_summed.csv')
        return module_vecs

    def get_module_similarities(self):
        """
        Returns a dataframe with pairwise cosine similarities between all module vectors ands saves it as a csv
        :return: dataframe with index and columns = module names
        """
        similarities = pd.DataFrame(cosine_similarity(self.module_vectors), columns=self.module_vectors.index,
                                    index=self.module_vectors.index)
        similarities.to_csv('../data/module_similarities.csv')
        return similarities

    def get_similarity_distributions(self):
        """
        Imports a pickle file with a dict to each distribution if such a file already exists
        :return:
        """
        if path.exists('similarity_distributions.pkl'):
            return pkl.load(open('similarity_distributions.pkl', 'rb'))
        else:
            return {}

    def add_sim_distribution(self, size1, size2):
        """
        Adds a distribution to the dict of similarity distributions based on a combination of sizes
        :param size1: number of genes in first module
        :param size2: number of genes in second module
        :return:
        """
        genes = list(self.gene_embeddings.index)
        dist = []
        for i in range(self.num_samples):
            genes1 = np.random.choice(genes, size1)
            genes2 = np.random.choice(genes, size2)
            embed1 = self.sum_embeddings(genes1)
            embed2 = self.sum_embeddings(genes2)
            dist.append(float(cosine_similarity([embed1], [embed2])))
        key = tuple(sorted([size1, size2]))
        self.similarity_dists[key] = dist
        pkl.dump(self.similarity_dists, open('similarity_distributions.pkl', 'wb'))

    def get_pvalue(self, size1, size2, sim):
        """
        Calulates the pvalue for a similarity based on the number of genes in each module used to calculate
        :param size1: number of genes in 1st module
        :param size2: number of genes in 2nd module
        :param sim: similarity between modules
        :return: pvalue
        """
        key = tuple(sorted([size1, size2]))
        if key not in self.similarity_dists.keys():
            self.add_sim_distribution(size1, size2)
        pval = sum(np.array(self.similarity_dists[key]) - sim >= 0)/float(self.num_samples)
        return pval

    def calculate_module_pvalues(self):
        pvals = pd.DataFrame(index=self.module_vectors.index, columns=self.module_vectors.index)
        for mod1 in self.module_to_genes.keys():
            for mod2 in self.module_to_genes.keys():
                sim = self.module_similarities.loc[mod1, mod2]
                size1 = len(self.module_to_genes[mod1])
                size2 = len(self.module_to_genes[mod2])
                pvals.loc[mod1, mod2] = self.get_pvalue(size1, size2, sim)
        return pvals

    # Functions to get info about module similarities
    #######
    def get_modules_with_pval(self, mod_of_interest):
        """
        Returns the names of all modules and their similarities with a cosine similarity
        to mod_of_interest pvalue < bonferroni pval
        :param mod_of_interest: module name
        :return: dataframe of similarities (index = name_
        """
        pval = .05 / self.num_modules
        mod_pvalues = pd.DataFrame(self.pvalues[mod_of_interest])
        sig_mods = mod_pvalues[mod_pvalues[mod_of_interest] < pval].index
        sims = self.module_similarities.loc[sig_mods, mod_of_interest]
        return sims

    def prioritize_genes(self, genes, modules):
        """
        Returns a df with gene and score for cosine similarity of gene embedding to module vector, sorted by similarity
        :param genes: list of gene names
        :param module: module name to use for module vector
        :return: pandas df: cols = 'gene', 'cosine similarity'
        """
        # gene_scores = pd.DataFrame(columns=['gene', 'cosine similarity'])
        scores = pd.DataFrame(index=genes)
        for module in modules:
            module_vector = self.module_vectors.loc[module,].to_numpy()
            module_vector = module_vector.reshape(1, len(module_vector))
            gene_embeds = self.gene_embeddings.loc[genes,:].to_numpy()
            similarities = pd.DataFrame(cosine_similarity(gene_embeds, module_vector),
                                        columns=[module], index=genes).sort_values(by=module, ascending=False)
            scores = scores.join(similarities, how = 'left')
        return scores

    def filter_similarity_matrix(self):
        """
        Return matrix of cosine similarities filtered by p_value (bonferroni corrected)
        :return:
        """
        sig = self.pvalues > .05 / (self.num_modules ** 2)
        sims = self.module_similarities
        sims[sig] = 0
        return sims
    
    def enrichment_in_modules(self, gene_set):
        """
        Computes hypergeometric pvalue for set of genes enriched in each of the modules
        Returns pandas dataframe
        :param gene_set: set of genes to see if 
        """
        M = len(self.gene_embeddings.index)
        pvals = {}
        for module, genes in self.module_to_genes.items():
            n = len(genes)
            N = len(gene_set)
            x = len(set(genes) & set(gene_set))
            pvals[module] = hypergeom.sf(x-1, M, n, N)
        return pd.DataFrame.from_dict(pvals, orient='index')
            
    def get_hypergeom_pval(self, mod1, mod2):
        """
        Gets the hypergeometric p-value for enrichment of mod2 genes in mod_1
        :param mod1 name of module1 1
        :param mod2 name of module 2 (calculating enrichment of this in mod1)
        :return: p-value


        https://blog.alexlenail.me/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
        M is the population size (previously N) - total genes
        n is the number of successes in the population (previously K) - cluster count
        N is the sample size (previously n) - module count
        X is still the number of drawn “successes”. - num in cluster and module
        pval = hypergeom.sf(x-1, M, n, N)
        """
        M = len(self.gene_embeddings.index)
        mod1_genes = set(self.module_to_genes[mod1])
        mod2_genes = set(self.module_to_genes[mod2])
        n = len(mod1_genes)
        N = len(mod2_genes)
        x = len(mod1_genes & mod2_genes)
        return hypergeom.sf(x - 1, M, n, N)

    def save_hypergeom_pvals(self, mod1_list, mod2_list, file):
        """
        Gets a dataframe of p values for enrichment of modules in list 2 in modules in list 1, saves to file
        """
        vhyper = np.vectorize(self.get_hypergeom_pval)
        pvals = vhyper(mod1_list, mod2_list)
        pv_df = pd.DataFrame(pvals, index=mod1_list, columns=mod2_list)
        pv_df.to_csv(file)
        return pv_df

    # Functions for SVM use
    def get_module_features(self):
        """
        Removes the drug and NASH modules and returns a set of module scores to use as features for svm
        :return:
        """
        drug_modules = pd.read_csv('../data/gene_maps/drug.csv')
        drugs = set(drug_modules['module'])
        rem = list(drugs) + [self.nash]
        mods = list(set(self.module_to_genes.keys()) - set(rem))
        return self.prioritize_genes(self.gene_embeddings.index, mods).sort_index(axis=1).sort_index()

    def load_befree_genes(self, threshold):
        """
        Loads befree genes that meet specified threshold
        :param threshold:
        :return:
        """
        befree = pd.read_csv('../data/befree_genes.csv')
        befree = befree[befree['score'] > threshold]
        return list(set(befree['gene']) & set(self.gene_embeddings.index) - set(self.curated_genes))

    def load_svensson_genes(self, threshold):
        """"
        Loads genes from top 200 svensson genes that meet threshold
        :param threshold
        :return
        """
        sven_data = pd.read_csv('../data/srebp1c high genes_top200_high.csv').dropna()
        logfc = 'high Log2 Fold Change'
        sven_data.index = [x.upper() for x in list(sven_data['FeatureName'])]
        sven_data = sven_data[sven_data[logfc] > threshold]
        return list(set(sven_data.index) & set(self.gene_embeddings.index) - set(self.curated_genes) - set(self.befree_genes))

    def gene_module_normalized_scores(self):
        mod_scores = self.get_module_features()

        dct = defaultdict(list)
        for mod in mod_scores.columns:
            mod_genes = set(self.module_to_genes[mod]) & set(self.gene_embeddings.index)
            dct['num_genes'].append(len(mod_genes))

            embeds = self.gene_embeddings.loc[mod_genes,]
            gene_sims = cosine_similarity(embeds)
            dct['mean: gene-gene within module'].append(np.mean(gene_sims))
            dct['std: gene-gene score within module'].append(np.std(gene_sims))

            gene_mod_sims = mod_scores.loc[mod_genes, mod]
            dct['mean: gene-module vector'].append(np.mean(gene_mod_sims))
            dct['std: gene-module vector'].append(np.std(gene_mod_sims))

        output = pd.DataFrame(dct, index=mod_scores.columns)
        output.to_csv('../results/module_normalization.csv')