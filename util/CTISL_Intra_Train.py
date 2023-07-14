from numpy import *
from sklearn.model_selection import StratifiedKFold
from collections import Counter

from util import CTISL_model
import numpy as np
from sklearn.feature_selection import SelectKBest, chi2
import pandas as pd
import scanpy as sc



class DataProcess:
    def __init__(self, name, feture_num, normalize, fileform):
        super().__init__()

        self.name = name
        self.feture_num = feture_num
        self.normalize = normalize
        self.fileform=fileform

    def prefilter_specialgenes(self, adata, Gene1Pattern="ERCC", Gene2Pattern="MT-"):

        id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
        id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
        id_tmp = np.logical_and(id_tmp1, id_tmp2)
        adata._inplace_subset_var(id_tmp)

    def label_to_number(self, label):
        label = pd.Series(label, dtype='category')
        label = label.cat.rename_categories(range(len(label.cat.categories)))
        np.savetxt('data/' + self.name + '/label.txt', label)
        return len(unique(label))

    def selectdata(self, data, name, n):
        label = np.loadtxt('./data/' + name + '/trainlabel' + str(n) + '.txt')
        model1 = SelectKBest(chi2, k=self.feture_num)
        model = model1.fit(data, label)
        dfscores = pd.DataFrame(model.scores_)
        dfcolumns = pd.DataFrame(range(0, data.shape[1]))
        featureScores = pd.concat([dfcolumns, dfscores], axis=1)
        featureScores.columns = ['feature', 'score']  # naming the dataframe columns
        featureScores.sort_values(by=["score"], inplace=True, ascending=False)

        a_dict = featureScores.nlargest(self.feture_num, 'score').index

        return a_dict

    def selectAll(self, data, name, n, m):
        lists = []
        for i in range(0, n):
            lists.append(self.selectdata(data, name, i))
        np.savetxt('./data/' + name + '/fetureindex' + str(m) + '.txt', lists, fmt='%d')

    def split_5k(self):
        if self.fileform=='csv':
            data = sc.read_csv('./data/' + self.name + '/' + self.name + '.csv')
        else:
            data = sc.read_h5ad('./data/' + self.name + '/' + self.name + '.h5ad')

        label = pd.read_csv('./data/' + self.name + '/' + self.name + 'label.txt')

        label_num=label['x']
        label = np.array(label['x'])
        data.obs['x'] = label

        #convert label to num
        label_num = pd.Series(label_num, dtype='category')
        label_num = label_num.cat.rename_categories(range(len(label_num.cat.categories)))

        label_num=np.array(label_num)
        data.obs['x_num']=label_num
        data.obs['label']=label
        sc.pp.filter_cells(data, min_genes=100)
        sc.pp.filter_genes(data, min_cells=10)

        self.prefilter_specialgenes(data)


        kfold = StratifiedKFold(n_splits=5)
        index=np.arange(data.X.shape[0])
        np.random.shuffle(index)
        data=data[index,:]
        class_num=len(data.obs['x_num'].unique())
        m = 1

        for train, test in kfold.split(data.X, data.obs['label']):
            X_train, X_test = data[train, :], data[test, :]
            y_train=X_train.obs['x_num']
            y_test = X_test.obs['x_num']

            if self.normalize:
                sc.pp.normalize_per_cell(X_train, counts_per_cell_after=10e4)
                sc.pp.log1p(X_train)
                sc.pp.normalize_per_cell(X_test, counts_per_cell_after=10e4)
                sc.pp.log1p(X_test)
            X_train = X_train.X
            X_test = X_test.X
            np.savetxt('./data/' + self.name + '/trainlabel.txt', y_train)
            np.savetxt('./data/' + self.name + '/testlabel.txt', y_test)
            # np.savetxt('./results/' + self.name + '/truelabel' + str(m) + '.txt', y_test)
            self.two_label(self.name, class_num)
            self.selectAll(X_train, self.name, class_num, m)
            sc.pp.scale(X_train, zero_center=True, max_value=6)
            sc.pp.scale(X_test, zero_center=True, max_value=6)
            np.savetxt('./data/' + self.name + '/testdata.txt', X_test)
            np.savetxt('./data/' + self.name + '/traindata.txt', X_train)
            feture_index = np.loadtxt('./data/' + self.name + '/fetureindex' + str(m) + '.txt',
                                      dtype=int)

            feature_index = feture_index.reshape(feture_index.shape[0] * feture_index.shape[1])


            X_train = X_train[:, feature_index]
            X_test = X_test[:, feature_index]

            np.savetxt('./data/' + self.name + '/train_1600+pro6.txt', X_train)
            np.savetxt('./data/' + self.name + '/test_1600+pro6.txt', X_test)

            f = open('./data/' + self.name + '/' + self.name + 'labelcount.txt', 'wt')
            print('trainlabel', file=f)
            print(Counter(y_train), file=f)
            print('testlabel', file=f)
            print(Counter(y_test), file=f)
            CTISL_model.Model(name=self.name,class_num=class_num,m_k=m,feature_num=self.feture_num,
                              traindata=X_train,testdata=X_test,trainlabel=y_train,testlabel=y_test)


            m = m + 1
        return

    def two_label(self, dataname, n):
        for j in range(0, n):
            label = np.loadtxt('./data/' + dataname + '/trainlabel.txt')
            for i in range(0, len(label)):
                if label[i] == j:
                    label[i] = 1
                else:
                    label[i] = 0
            np.savetxt('./data/' + dataname + '/trainlabel' + str(j) + '.txt', label)

