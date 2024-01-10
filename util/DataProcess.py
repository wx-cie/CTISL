from numpy import *
import numpy as np
from anndata import AnnData
from sklearn.neural_network import MLPClassifier
from sklearn.feature_selection import SelectKBest, chi2
import pandas as pd
import scanpy as sc
from mlxtend.feature_selection import ColumnSelector
from mlxtend.classifier import StackingCVClassifier
from sklearn.pipeline import make_pipeline
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.multiclass import OneVsRestClassifier
import joblib
import os


class DataProcess:
    def __init__(self, traindata,trainlabel,testdata, predictlabel,modelname,feature_num,normalize,fileform):
        super().__init__()
        self.traindata=traindata
        self.trainlabel=trainlabel
        self.testdata=testdata
        self.predictlabel=predictlabel
        self.modelname=modelname

        self.name = 'temp'
        if os.path.exists('./temp/') == False:
            os.mkdir('./temp/')
        self.feature_num = feature_num

        self.normalize = normalize
        self.fileform=fileform


    def prefilter_specialgenes(self, adata, Gene1Pattern="ERCC", Gene2Pattern="MT-"):
        id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
        id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
        id_tmp = np.logical_and(id_tmp1, id_tmp2)
        adata._inplace_subset_var(id_tmp)

    def selectdata(self, data, name, n):
        label = np.loadtxt('./' + name +'/'+ '/trainlabel' + str(n) + '.txt')
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
        np.savetxt('./' + name +'/'+ '/featureindex.txt', lists, fmt='%d')

    def split_train_test(self):
        #data process
        if self.fileform=='csv':
            train_data = sc.read_csv(self.traindata)
            test_data = sc.read_csv(self.testdata)
        else:
            train_data = sc.read_h5ad(self.traindata)
            test_data = sc.read_h5ad(self.testdata)
        train_label = pd.read_csv(self.trainlabel)
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


        #fearute select process
        #1.convert trainlabel to num
        label = pd.Series(adata_train.obs['x'], dtype='category')
        label_name = label.cat.categories
        pd.Series(label_name).to_csv('./' + self.name +'/'+ '/labelunique.txt')
        label = label.cat.rename_categories(range(len(label.cat.categories)))
        np.savetxt('./' + self.name +'/'+ '/trainlabel.txt', label)
        classnum = len(np.unique(label))
        #2. Multi classes labels are converted into multiple two classes labels
        self.two_label(self.name, classnum)
        #3.Feature select
        self.selectAll(X_train, self.name, classnum)
        sc.pp.scale(X_train, zero_center=True, max_value=6)
        sc.pp.scale(X_test, zero_center=True, max_value=6)
        #4.read feature index
        feature_index = np.loadtxt('./' + self.name +'/'+ '/featureindex.txt', dtype=int)
        feature_index = feature_index.reshape(feature_index.shape[0] * feature_index.shape[1])
        #5.producing train and test after feature selecting
        X_train = X_train[:, feature_index]
        X_test = X_test[:, feature_index]
        classnum = len(label.unique())
        self.CTISLModel(self.name,classnum,1,self.feature_num,X_train,X_test,trainlabel)

        return

    def two_label(self, dataname, n):
        for j in range(0, n):
            label = np.loadtxt('./' + dataname +'/'+ '/trainlabel.txt')
            if j in label:
                for i in range(0, len(label)):
                    if label[i] == j:
                        label[i] = 1
                    else:
                        label[i] = 0
                np.savetxt('./' + dataname +'/'+ '/trainlabel' + str(j) + '.txt', label)
    def CTISLModel(self,name, class_num, m_k,feature_num,traindata,testdata,trainlabel):
        if self.modelname=='MLP':
            sclf = MLPClassifier(hidden_layer_sizes=(100, 50), max_iter=300)
        else:
            lists = []
            for i in range(class_num):
                a_dict = range(i * feature_num, (i + 1) * feature_num)
                pipe1 = make_pipeline(ColumnSelector(cols=a_dict),
                                      LogisticRegression(C=0.1, random_state=11))
                pipe2 = make_pipeline(ColumnSelector(cols=a_dict),
                                      OneVsRestClassifier(
                                          svm.SVC(kernel='rbf', probability=True, C=0.5, random_state=11)))
                lists.append(pipe1)
                lists.append(pipe2)
                sclf = StackingCVClassifier(classifiers=lists,
                                            meta_classifier=LogisticRegression(C=1.5),
                                            use_probas=True,
                                            cv=3,
                   
                                            random_state=11)
            sclf = StackingCVClassifier(classifiers=lists,
                                        meta_classifier=LogisticRegression(C=1.5),
                                        use_probas=True,
                                        cv=3,
                                        random_state=11)

      
        sclf.fit(traindata, trainlabel)
        joblib.dump(sclf, './' + name +'/'+ str(m_k) + '.model')
        model = joblib.load('./' + name +'/'+ str(m_k) + '.model')
        predict_classes = model.predict(testdata)
        np.savetxt(self.predictlabel,predict_classes,fmt='%s')




        return

