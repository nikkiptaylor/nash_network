from modules import *
from sklearn.feature_selection import SelectKBest, f_classif
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, plot_confusion_matrix, plot_roc_curve, plot_precision_recall_curve, roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE
from sklearn.svm import SVC
import os
import pickle as pkl

plt.style.use('ggplot')


class Nash_Model(object):
    def __init__(self, save_path, use_modules=True, feat_sel=True, sample=None, model_type='SVC'):
        np.random.seed(0)
        self.M = Modules()
        
        # module scores vs embeddings as features
        if use_modules:
            # use gene-module vector cosine similarities
            self.dataset = self.M.get_module_features()
        else:
            # use embeddings directly
            self.dataset = self.M.gene_embeddings

        self.feat_sel = feat_sel
        self.sample = sample
        self.skb = SelectKBest(f_classif, k=64)
        
        if model_type == 'LR':
            self.clf = LogisticRegression(class_weight='balanced')
        else:
            self.clf = SVC(kernel='linear', class_weight='balanced', probability=True)
        
        self.save_path = self.set_save_path(save_path)

        # leave out 200 genes to use as negatives for testing
        # make sure negative genes are not in known Nash gene sets
        negative_genes = sorted(list(set(self.dataset.index) - set(self.M.curated_genes)
                                     - set(self.M.befree_genes) - set(self.M.sven_genes)))
        self.neg_test_genes = sorted(list(np.random.choice(negative_genes, 200)))
        self.neg_train_genes = sorted(list(set(negative_genes) - set(self.neg_test_genes)))


#     def load_trained_model(self, model_path, skb_path):
#         """
#         Loads model and feature selector from specified paths. Alter dataset for selected features
#         """

#         self.clf = pkl.load(open(model_path, 'rb'))
#         self.skb = pkl.load(open(skb_path, 'rb'))
#         self.dataset = pd.DataFrame(self.skb.transform(self.dataset), index=self.dataset.index)

    def set_save_path(self, path):
        """
        Updates save path and creates files
        """
        self.save_path = path
        if not os.path.exists(path):
            os.makedirs(path)
        return path
    
    def format_input(self, pos_genes, neg_genes):
        """
        Returns a formatted test or train set using the model's features and specified genes.
        Format to be used for all methods below for any variation of X, y

        :param pos_genes: list of gene symbols to be used as positives for the model
        :param neg_genes: list of gene symbols to be used as negatives for the model

        :return X: dataframe of dataset only including the positive and negative genes
        list of labels for each genes
        :return y: labels for genes in X
        """
        X = self.dataset.loc[pos_genes + neg_genes, ]
        y = np.array([1] * len(pos_genes) + [0] * len(neg_genes))
        return X, y

    def do_feat_sel(self, train_X, train_y):
        """
        Does SelectKBest (sklearn) feature selection on the given training dataset and returns a transformed with top 64
        features selected. Saves csv with feature scores to self.save_path

        :return transformed train X with top 64 features selected
        """
        # feature selection based on training set only
        train_X_fs = pd.DataFrame(self.skb.fit_transform(train_X, train_y), index=train_X.index)
        scores = pd.DataFrame({'score': self.skb.scores_, 'pval': self.skb.pvalues_}).T
        scores.columns = train_X.columns
        scores.to_csv(self.save_path + '/feature_scores.csv')
        return train_X_fs

    def train(self, train_X, train_y):
        """
        Does feature selection, resampling, and trains model on train_X and train_y
        """
        if self.feat_sel:
            train_X = self.do_feat_sel(train_X, train_y)

        train_X, train_y = self.sample.fit_resample(train_X, train_y)
        self.clf.fit(train_X, train_y)

    def train_all_curated(self, bench=False):
        """
        Trains svm on all 70 curated genes and saves model
        :param bench: if true, the model will do benchmarking
        """
        train_X, train_y = self.format_input(self.M.curated_genes, self.neg_train_genes)
        self.train(train_X, train_y)
        pkl.dump(self, open(self.save_path + '/nash_model_trained.pkl', 'wb'))
        if bench:
            self.benchmark(train_X, train_y)

        # do feature selection on dataset as a whole so it is easier to be scored
        if self.feat_sel:
            self.dataset = pd.DataFrame(self.skb.transform(self.dataset), index=self.dataset.index)

    def benchmark(self, train_X, train_y):
        """
        Does cross validation and tests the trained model on both the befree and svensson gene sets.
        Saves output to dataframe in self.save_path
        """
        roc_ap = []
        test_X, test_y = self.format_input(self.M.befree_genes, self.neg_test_genes)
        roc_ap.append(self.test(test_X, test_y))
        test_X, test_y = self.format_input(self.M.sven_genes, self.neg_test_genes)
        roc_ap.append(self.test(test_X, test_y))
        roc_ap.extend(self.cross_validate(train_X, train_y))

        index = list(range(1, 11))
        index = ['befree', 'sven'] + index
        output = pd.DataFrame([index] + list(zip(*roc_ap)), index=['label', 'roc', 'ap'])
        output.to_csv(self.save_path + '/benchmarking output')
#         print(output)

    def cross_validate(self, X, y):
        """
        Does Stratified 10 Fold cross validation on the given X and y (training set)
        :return list of tuples with (roc, ap) for each CV fold
        """
        roc_ap = []
        kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=0)
        for train_ix, test_ix in kfold.split(X, y):
            train_X, test_X = X.iloc[train_ix, :], X.iloc[test_ix, :]
            train_y, test_y = y[train_ix], y[test_ix]
            self.train(train_X, train_y)
            roc_ap.append(self.test(test_X, test_y))
        return roc_ap

    def test(self, test_X, test_y):
        """
        Tests trained model on given test set and returns Area Under the ROC Curve and Average Precision Score
        """
        if self.feat_sel:
            test_X = self.skb.transform(test_X)
        predicted = self.clf.predict_proba(test_X)[:, 1]
        return roc_auc_score(test_y, predicted), average_precision_score(test_y, predicted)

    def score_all_genes(self):
        """
        Scores all genes in dataset using predic_proba and saves scores and predictions to save_path
        :return:
        """
        scores = pd.DataFrame(self.clf.predict_proba(self.dataset), index=self.dataset.index)[1]
        scores = pd.DataFrame(scores).sort_values(1, ascending=False)
        scores['known'] = [int(g in list(self.M.befree_genes + self.M.curated_genes + self.M.sven_genes)) for g in scores.index]
        scores.columns = ['score', 'known']
        scores.to_csv(self.save_path + '/all_gene_scores.csv')

        predictions = pd.DataFrame(self.clf.predict(self.dataset), index=self.dataset.index)
        predictions['known'] = [int(g in list(self.M.befree_genes + self.M.curated_genes + self.M.sven_genes)) for g in predictions.index]
        predictions.to_csv(self.save_path + '/predictions.csv')





