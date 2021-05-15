# nash_network

## General Info:
### Gene sets:
* **Embedded genes:** 14,707 genes embedded from String PPI ('data/gene_maps/embeddings.csv')

* **Curated genes:** 70 genes identified as Nash related from the curated list from Disgenet. Derived from UNIPROT, CGI, ClinGen, Genomics England, CTD (human subset), PsyGeNET, and Orphanet. 

* **Befree genes:** set of 308 genes identified as Nash related through the Befree system by Disgenet (genes in curated list removed). Each gene is accompanied by a confidence score.

* **Svensson genes:** 118 experimentally derived Nash related genes from Dr. Svensson's lab at Stanford, mapped from the list of the top 200 differentially expressed genes from scRNA seqencing of liver cells by diet induced mice. Each gene is accompanied by a score for the log Fold Change

* **Modules:** sets of genes with biological relevance that are used as features by the model. Overall there are 237 modules from different sources:
    * 47 immune response modules: Proteomics data from ImmProt – modules derived from 6982 proteins enriched in immune cells (PMID: 28263321)
    * 50 signaling pathway modules: Hallmark dataset from Molecular Signatures Database hallmark gene set collection
    * 137 metabolic modules: Genome scale metabolic pathways derived from Human Metabolic Reaction Database (Human GEM)
    * 3 cytokine related modules

### Embeddings and Module Vectors
**Gene embeddings:** 64 dimensional representations of genes embedded from the String Human PPI using node2vec. There were a total of 9,250,034 edges found in STRING between 14,707 genes. After filtering for a confidence score of above 800, 728,090 edges were included in training the embedding space.

node2vec hyperparameters: length of walks=30, number of walks=10, min count=1, batch word=6, window=10.

**Module vectors:** Module vectors are used to represent a set of genes in the embedding space. Vectors were created by summing the gene embeddings of all genes represented in a module.

**Distance metric:** Cosine similarity was used to calculate similarity scores between different gene embeddings and module vectors.

## Scripts
**modules.py:** 

Modules class to encapsulate gene sets and embeddings and do various calculations with them.
Creating a Modules instance loads the embeddings, module vectors, dicts to map genes to modules, befree and svensson genes, and summed module vectors. Setting calc_modules=True when initializing will calculate the module vectors and similarities from scratch, otherwise load from filed in 'data' folder. The list of curated genes, befree genes, and svensson genes can be accessed through the class's variables directly or by thresholding using class methods.

Example:
```
m = Modules() # initialize to load module info
m.gene_embeddings # pandas dataframe with genes as index, embedding vectors as columns
m.module_to_genes['cytokines'] # gets list of genes in cytokine module
m.gene_to_modules['ACSL1'] # gets list of modules ACSL1 is in
```

The most important methods for getting features to feed into the predictive model are
``` 
m.prioritize_genes(genes, modules) # returns matrix of cosine similarity between specified gene embeddings and module vectors
m.get_module_features # returns matrix of cosine similarity between all 14,707 embedded genes and all 237 modules
```

Note: this file contains a number of methods that were used for previous analysis but are not part of the predictive modeling, so most of the methods can be ignored


**refined_model.py:** 

Contains a Nash_Model class and methods to be used for training and testing predictive models using the modules and genesets described above. 

To train model:
1. Instantiate a Nash_Model object with desired parameters
2. Call train_all_curated
3. Call score_all_genes if desired

When the model is trained, the entire object is stored in the 'save_path' specified upon instantiation for future use. Any analysis based on a Nash_Model that has already been trained should use the 'dataset'and 'neg_test_genes' variables within the object to maintain continuity. 

Examples of how to instantiate the model with different configurations can be found in 'benchmarking.ipynb' and 'train_final_model.ipynb'


**benchmarking.ipynb** 

Does benchmarking of model parameters, including
1. Features: gene embeddings vs module scores vs module scores with feature selection
2. Sampling: SMOTE vs 1:2 undersampling
3. Model type: linear svc vs logistic regression (svc only included in benchmarking tables)
    
Also does benchmarking of how different score cuttofs for Svensson genes and Befree genes impact ROC and AP scores for the final trained model (detailed in 'train_final_model.ipynb')


**train_final_model.ipynb**

Trains and saves the final "refined" model that performed best based on benchmarking and Svensson/Befree genes as external test sets.

The final model used to score all genes included a Linear SVC, training using all 70 curated genes as positives and a random set of 200 held out as negative set, random undersampling to achieve 1:2 positive to negative ratio, gene-module cosine similarity scores as features, and feature selection to pick top 64 modules based on ANOVA F-value.

## Figures:

## To alter figures or models:

1. Install required packages: pip install -r requirements.txt
2. To alter figures or create results using trained model, load trained model using:
```
sys.path.append(r'../../scripts') # for unpickling
nash_svc = pkl.load(open('../../results/final_model_svc/nash_model_trained.pkl', 'rb'))
feature_selector = nash_svc.skb
model = nash_svc.clf
```
Note: figures rely on model output, so the path to the final model must have the most up to date trained model before running code to generate figures. 

To re-train model or use different configurations, change Nash_Model input in train_final_model and run all cells.


