from numpy import *
from anndata import AnnData
import numpy as np
from sklearn.feature_selection import SelectKBest, chi2
import pandas as pd
import scanpy as sc


class DataProcess:
    def __init__(self, name, feature_num, name1, name2, normalize,fileform):
        super().__init__()

        self.name = name
        self.feature_num = feature_num
        self.name1 = name1
        self.name2 = name2
        self.normalize = normalize
        self.fileform=fileform


    def prefilter_specialgenes(self, adata, Gene1Pattern="ERCC", Gene2Pattern="MT-"):
        id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
        id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
        id_tmp = np.logical_and(id_tmp1, id_tmp2)
        adata._inplace_subset_var(id_tmp)

    def selectdata(self, data, name, n):
        label = np.loadtxt('./data/' + name + '/trainlabel' + str(n) + '.txt')
        model1 = SelectKBest(chi2, k=self.feature_num)
        model = model1.fit(data, label)
        dfscores = pd.DataFrame(model.scores_)
        dfcolumns = pd.DataFrame(range(0, data.shape[1]))
        featureScores = pd.concat([dfcolumns, dfscores], axis=1)
        featureScores.columns = ['feature', 'score']  # naming the dataframe columns
        featureScores.sort_values(by=["score"], inplace=True, ascending=False)
        a_dict = featureScores.nlargest(self.feature_num, 'score').index
        return a_dict

    def selectAll(self, data, name, n):
        lists = []
        for i in range(0, n):
            lists.append(self.selectdata(data, name, i))
        np.savetxt('./data/' + name + '/featureindex.txt', lists, fmt='%d')

    def split_train_test(self):
        #data process
        if self.fileform=='csv':
            train_data = sc.read_csv('./data/' + self.name1 + '/' + self.name1 + '.csv')
            test_data = sc.read_csv('./data/' + self.name2 + '/' + self.name2 + '.csv')
        else:
            train_data = sc.read_h5ad('./data/' + self.name1 + '/' + self.name1 + '.h5ad')
            test_data = sc.read_h5ad('./data/' + self.name2 + '/' + self.name2 + '.h5ad')
        train_label = pd.read_csv('./data/' + self.name1 + '/' + self.name1 + 'label.txt')
        test_label = pd.read_csv('./data/' + self.name2 + '/' + self.name2 + 'label.txt')
        train_label = np.array(train_label['x'])
        train_data.obs['x'] = train_label
        train_data.var_names_make_unique(join='-')
        train_data.obs_names_make_unique(join='-')
        sc.pp.filter_cells(train_data, min_genes=100)
        sc.pp.filter_genes(train_data, min_cells=10)
        self.prefilter_specialgenes(train_data)
        if self.normalize:
            sc.pp.normalize_per_cell(train_data,counts_per_cell_after=10e4)
            sc.pp.log1p(train_data)
        train_data.var_names = [i.upper() for i in list(train_data.var_names)]


        test_label = np.array(test_label['x'])
        test_data.obs['x'] = test_label
        test_data.var_names_make_unique(join='-')
        test_data.obs_names_make_unique(join='-')
        sc.pp.filter_cells(test_data, min_genes=100)
        sc.pp.filter_genes(test_data, min_cells=10)
        self.prefilter_specialgenes(test_data)
        if self.normalize:
            sc.pp.normalize_per_cell(test_data, counts_per_cell_after=10e4)
            sc.pp.log1p(test_data)
        test_data.var_names = [i.upper() for i in list(test_data.var_names)]

        adata_tmp = []
        adata_tmp.append(train_data)
        adata_tmp.append(test_data)
        full_adata = AnnData.concatenate(*adata_tmp, join='inner', batch_key="dataset_batch",
                                         batch_categories=["train", "test"])  # inner
        del adata_tmp
        del test_data
        del train_data
        ref_id = full_adata.obs["dataset_batch"] == "train"
        adata_test = full_adata[~ref_id, :].copy()
        adata_train = full_adata[ref_id, :].copy()
        X_test = adata_test.X
        X_train = adata_train.X
        trainlabel=adata_train.obs['x']
        testlabel=adata_test.obs['x']

        #fearute select process
        #1.convert trainlabel to num
        label = pd.Series(adata_train.obs['x'], dtype='category')
        label_name = label.cat.categories
        pd.Series(label_name).to_csv('./data/' + self.name + '/labelunique.txt')
        label = label.cat.rename_categories(range(len(label.cat.categories)))
        np.savetxt('./data/' + self.name + '/trainlabel.txt', label)
        classnum = len(np.unique(label))
        #2. Multi classes labels are converted into multiple two classes labels
        self.two_label(self.name, classnum)
        #3.Feature select
        self.selectAll(X_train, self.name, classnum)
        sc.pp.scale(X_train, zero_center=True, max_value=6)
        sc.pp.scale(X_test, zero_center=True, max_value=6)
        #4.read feature index
        feature_index = np.loadtxt('./data/' + self.name + '/featureindex.txt', dtype=int)
        feature_index = feature_index.reshape(feature_index.shape[0] * feature_index.shape[1])
        #5.producing train and test after feature selecting
        X_train = X_train[:, feature_index]
        X_test = X_test[:, feature_index]
        classnum = len(label.unique())

        return classnum, X_train, X_test,trainlabel,testlabel

    def two_label(self, dataname, n):
        for j in range(0, n):
            label = np.loadtxt('./data/' + dataname + '/trainlabel.txt')
            if j in label:
                for i in range(0, len(label)):
                    if label[i] == j:
                        label[i] = 1
                    else:
                        label[i] = 0
                np.savetxt('./data/' + dataname + '/trainlabel' + str(j) + '.txt', label)




