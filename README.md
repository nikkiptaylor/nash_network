# nash_network

## General Info:
### Gene sets:
* **Embedded genes:** 14,707 genes embedded from String PPI ('data/gene_maps/embeddings.csv')

* **Curated genes:** 70 genes identified as Nash related from the curated list from DisGeNET. Derived from UNIPROT, CGI, ClinGen, Genomics England, CTD (human subset), PsyGeNET, and Orphanet. 

* **Befree genes:** set of 308 genes identified as Nash related through the Befree system by DisGeNET (genes in curated list removed). Each gene is accompanied by a confidence score.

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

**Distance metric:** Cosine similarity was used to calculate similarity scores between different gene embeddings and module vectors. Gene-module scores are indicated by the cosine similarity of a gene's embedding vector and a module's summed module vector.

## Data files
**gene_maps:**
Folder that contains csv files that map gene symbols to module names or embeddings. 
* all_modules.csv - combined list that maps genes to immune response modules ('me.csv'), hallmark signaling modules ('hallmark.csv'), metabolic modules ('ms.csv'), and drug modules (drug.csv)
* disease.csv - maps gene symbol to disease name, from the curated list at DisGeNET
* embedding.csv - maps gene symbol to embedding vector
* snps-tf.csv - maps gene symbols to label as a Nash related SNP, PGC1 transcription factor, or gene downstream of PGC1s


**befree_genes.csv:** list of befree genes from DisGeNET and their scores 

**module_vectors_summed.csv:** maps module name to module vector. This file gets re-calculated when initializing a Modules object with 'calc_modules=True'

**module_similarities.csv:** matrix with pairwise cosine similarity of all module vectors. Rows and columns are module names, values are cosine similarities. 

**svensson_200.csv:** file with top 200 genes identified by Dr. Svensson's lab through scRNA sequencing of mice. These are the genes (118 out of 200 after mapping and removing known genes) that are used as an external validation set in the predictive modeling. genes are listed under column 'FeatureName', log fold change is listed under 'high Log2 Fold Change'

**svensson_909.csv:** file with all 909 genes sent by Dr. Svensson. These are NOT used as ground truth Nash genes, they are used for general analysis of model scoring vs log fold change. 607 of these 909 genes can be mapped to the set of 14,707 embedded genes.

**similarity_pvalues_sum:** old file, to be ignored


## Scripts
**modules.py:** 

Modules class to load gene sets and embeddings and do various calculations with them.
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
# return matrix of cosine similarity between specified gene embeddings and module vectors
m.prioritize_genes(genes, modules) 
# returns matrix of cosine similarity between all 14,707 embedded genes and all 237 modules
m.get_module_features 
```

Note: this file contains a number of methods that were used for previous analysis but are not part of the predictive modeling, so most of the methods can be ignored


**refined_model.py:** 

Contains a Nash_Model class and methods to be used for training and testing predictive models using the modules and genesets described above. Th model is used to train and test a binary classifier for genes to predict whether genes are Nash related ('positive') or not ('negative').

To train model:
1. Instantiate a Nash_Model object with desired parameters
2. Call train_all_curated
3. Call score_all_genes if desired

When the model is trained, the entire object is stored in the 'save_path' specified upon instantiation for future use. Any analysis based on a Nash_Model that has already been trained should use the 'dataset'and 'neg_test_genes' variables within the object to maintain continuity. 

Examples of how to instantiate the model with different configurations can be found in 'train_final_model.ipynb'


**train_final_model.ipynb**


Part 1: Benchmarking - 

Does benchmarking of model parameters, including
1. Features: gene embeddings vs module scores vs module scores with feature selection
2. Sampling: SMOTE vs 1:2 undersampling
3. Model type: linear svc vs logistic regression (svc only included in benchmarking tables)
    
Also does benchmarking of how different score cuttofs for Svensson genes and Befree genes impact ROC and AP scores for the final trained model (detailed in 'train_final_model.ipynb')

Part 2: Train refined model and score genes - 

Trains and saves the final "refined" model that performed best based on benchmarking and Svensson/Befree genes as external test sets.

The final model used to score all genes included a Linear SVC, training using all 70 curated genes as positives and a random set of 200 held out as negative set, random undersampling to achieve 1:2 positive to negative ratio, gene-module cosine similarity scores as features, and feature selection to pick top 64 modules based on ANOVA F-value. Scores are computed using the predict_proba function, which is calibrated using Platt's scaling by setting 'probability=True' when initializing the SVC.

**module_dispersion.ipynb**

Notebook to calculate differnt metrics to determine how dispersed genes in a module are within the embedding space


## Results:
 **benchmarking:** Folder that contains the results from benchmarking the model and different cutoffs. Results are generated from part 1 of 'train_final_model.ipynb'
 1. befree_cutoffs.csv - results from benchmarking the cutoff for befree genes (i.e. which befree genes are used as positives in model testing). Rows are AUROC and AP scores, columns are different cutoffs. Cutoff of -1 means that all genes are used
 2. svensson_cutoffs.csv - results from benchmarking the cutoff for svensson genes (i.e. which svensson genes are used as positives in model testing). Rows are AUROC and AP scores, columns are different cutoffs. Cutoff of 0 means that all genes are used
 3. full_summary.csv - Summary of results from benchmarking different model configurations (as described in benchmarking.ipynb). Each combination of different parameters has two rows: one for AUROC, one for AP. The befree column tests the model trained with that combination of parameters on the befree genes as positives and the held out negative set. The svensson column tests the model on the svensson genes as positives and the held out negative set. Columns labeled 1-10 are the results of 10-fold cross validation using curated genes as positives and all other training genes as negatives. 
 
Each of these models is saved along with its benchmarking output:
* embeddings = model was trained using embeddings as features
* mod_scores = model was trained using gene-module scores as features, no feature selection
* fs_mod_scores = model was trained using gene-module scores as features, feature selection
     * over = SMOTE oversampling was used
     * under = 1:2 undersampling was used
         * SVC = linear SVC
             * nash_model_trained.pkl = the full nash_model object that contains the trained svc and feature selector
             * model.pkl = the trained svc (same as contained in nash_model_trained.pkl)
             
**final_model_svc:** Directory that contains the results from the refined model as described in part 2 of train_final_model.ipynb. Contains the saved Nash_Model object, a csv of all 14,707 genes and their scores from the sklearn's predict_proba, a csv of binary predictions for each of those genes using sklearn's predict, and the feature importance scores for each of the 237 modules and their p-values. 


**drug_targets_scored.csv:** list of drugs that have been in phase 2-3 for Nash mapped to their targets and the model's probability score for those targets


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

To re-train model or use different configurations, change Nash_Model input in part 2 train_final_model and run all cells below.


