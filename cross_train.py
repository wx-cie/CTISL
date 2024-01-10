from util import CTISL_model,CTISL_Inter_Train
import os
import argparse
#inter/cross_species/cross_batches
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parameters to be used for Stacking')
    parser.add_argument('-Sourcename', '--source', required=True, help='input sourcename ')
    parser.add_argument('-Targetname', '--target', required=True, help='input targetname')
    parser.add_argument('-ResultPath', '--result', required=True, help='The folder path for the predicted results')
    parser.add_argument('-Modelname', '--modelname', required=True, help='The input name can be either CTISL or MLP')
    parser.add_argument('-Fileform', '--form', required=True, help='The input file format can be either csv or h5ad ')
    parser.add_argument('-Norm', '--norm', required=True, help='If it is raw data, it needs to be normalized.')

    args =parser.parse_args()
    source_name = args.source
    target_name = args.target
    result = args.result
    modelname=args.modelname
    fileform =args.form
    normalize=args.norm


    feature_num=300
    name = source_name+'_'+target_name
    data_path='./data/'+name
    if os.path.exists(data_path)==False:
        os.mkdir(data_path)
    dataprocess=CTISL_Inter_Train.DataProcess(name,feature_num,source_name,target_name,normalize,fileform)
    class_num,train,test,trainlabel,testlabel=dataprocess.split_train_test()
    CTISL_model.Model(name=name, modelname=modelname, result=result,class_num=class_num,m_k=1,feature_num=feature_num,traindata=train,
                      testdata=test,trainlabel=trainlabel,testlabel=testlabel)
