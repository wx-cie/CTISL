from util import CTISL_model,CTISL_Intra_Train
import argparse
#intra
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parameters to be used for Stacking')
    parser.add_argument('-Name', '--name', required=True, help='input name ')
    parser.add_argument('-Fileform', '--form', required=True, help='The input file format can be either csv or h5ad ')
    parser.add_argument('-Norm', '--norm', required=True, help='If it is raw data, it needs to be normalized.')

    args = parser.parse_args()
    name = args.name
    fileform = args.form
    normalize = args.norm
    feature_num=300
    dataprocess=CTISL_Intra_Train.DataProcess(name,feature_num,normalize,fileform)
    dataprocess.split_5k()