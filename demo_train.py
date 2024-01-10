from util import DataProcess
import argparse
#intra
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parameters to be used for Stacking')
    parser.add_argument('-Train', '--train', required=True, help='Path to the training dataset.')
    parser.add_argument('-Trainlabel', '--trainlabel', required=True, help='Path to the training dataset labels ')
    parser.add_argument('-Test', '--test', required=True, help='Path to the testing dataset.')
    parser.add_argument('-Predictlabel', '--predictlabel', required=True, help='Path to the predict label.')
    parser.add_argument('-Modelname', '--modelname', required=True, help='The input name can be either CTISL or MLP ')
    parser.add_argument('-Fileform', '--form', required=True, help='The input file format can be either csv or h5ad ')
    parser.add_argument('-Norm', '--norm', required=True, help='If it is raw data, it needs to be normalized.')

    args = parser.parse_args()

    traindata=args.train
    trainlabel=args.trainlabel
    testdata=args.test
    predictlabel=args.predictlabel
    modelname=args.modelname
    fileform = args.form
    normalize = args.norm
    feature_num=300
    DataProcess.DataProcess(traindata,trainlabel,testdata,predictlabel,modelname, feature_num,normalize,fileform).split_train_test()

