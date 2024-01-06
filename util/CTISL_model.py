from mlxtend.feature_selection import ColumnSelector
from mlxtend.classifier import StackingCVClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.multiclass import OneVsRestClassifier
import joblib
import numpy as np



def Model(name, modelname,class_num, m_k,feature_num,traindata,testdata,trainlabel,testlabel):
    if modelname=='MLP':
        sclf = MLPClassifier(hidden_layer_sizes=(100, 50), max_iter=300)
    else:
        lists = []
        for i in range(class_num):
            a_dict = range(i * feature_num, (i + 1) * feature_num)
            pipe1 = make_pipeline(ColumnSelector(cols=a_dict),
                                  LogisticRegression(C=0.1, random_state=11))
            pipe2 = make_pipeline(ColumnSelector(cols=a_dict),
                                  OneVsRestClassifier(svm.SVC(kernel='rbf', probability=True, C=0.5, random_state=11)))
            lists.append(pipe1)
            lists.append(pipe2)

        sclf = StackingCVClassifier(classifiers=lists,
                                    meta_classifier=LogisticRegression(C=1.5),
                                    use_probas=True,
                                    cv=3,
                                    random_state=11)


    sclf.fit(traindata, trainlabel)

    joblib.dump(sclf, './model/' + name + str(m_k) + '.model')
    model = joblib.load('./model/' + name + str(m_k) + '.model')
    predict_classes = model.predict(testdata)

    acc = accuracy_score(testlabel, predict_classes)
    p_class, r_class, f_class, support_micro = precision_recall_fscore_support(
        testlabel, predict_classes)

    f = open('./results/' + name +'_'+str(m_k)+ 'result.txt', 'wt')
    print("acc : ", acc,file=f)
    print("macroF1score : ", f_class.mean(),file=f)
    print("medianF1score : ", np.median(f_class),file=f)
    f.close()

    return
