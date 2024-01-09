# CTISL: A Dynamic Stacking Multi-class Classification Approach for Identifying Cell Types from Single-cell RNA-sequencing Data

## Requirements
- Python 3.7.11
- scikit-learn 1.0.2
- scanpy 1.8.2
- numpy 1.19.1
- pandas 1.3.5

## Train CTISL on Intra-dataset

To train CTISL on the Intra-dataset, use the following command:

```shell
python Intra_train.py -Name '10Xv2' -Modelname 'CTISL' -Fileform 'h5ad' -Norm True
```

**Parameters:**
- Name: The Intra-dataset name (e.g., 10Xv2, 10Xv3, etc.)
- Modelname:The model name, which can be either 'CTISL' or 'MLP'.
- Fileform: The gene expression matrix file format, which can be either 'csv' or 'h5ad'.
- Norm: Specify whether the raw data needs to be normalized. Use 'True' for normalization.

## Train CTISL on Inter-data, Cross-batch, Cross-species

To train CTISL on Inter-data, Cross-batch, or Cross-species, use the following command:

```shell
python cross_train.py -Sourcename 'dentritic_batch_1' -Targetname 'dentritic_batch_2' -Modelname 'CTISL' -Fileform 'csv' -Norm False
```

**Parameters:**
- Sourcename: The name of the training set. For example, 'dentritic_batch_1' is used as the training set in the cross-batch experiment.
- Targetname: The name of the testing set. For example, 'dentritic_batch_2' is used as the testing set in the cross-batch experiment.
- Modelname:The model name, which can be either 'CTISL' or 'MLP'.
- Fileform: The gene expression matrix file format, which can be either 'csv' or 'h5ad'.
- Norm: Specify whether the raw data needs to be normalized. Use 'True' for normalization.

## Data Source
The following datasets were used:
- PBMC datasets: Downloaded from the Broad Institute Single Cell portal [link](https://portals.broadinstitute.org/single_cell/study/SCP424/single-cellcomparisonpbmc-data).
- Pancreas datasets: Downloaded from [link](https://hemberg-lab.github.io/scRNA.seq.datasets/).
- Airway datasets: Downloaded from [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102580).
- Dendritic datasets: Downloaded from [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820).
- Retina datasets: Downloaded from [link](https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/retina/  ).
  
## Run the code on your data
If you want to train and predict with your data, please use the following command.
```shell
python demo_train.py  -Train [*The path to your training dataset*] -Trainlabel [*The path to your training dataset labels*] -Test [*The path to your testing dataset*] -Predictlabel [*Location for storing predicted labels*]  -Modelname [*The model name] -Fileform [*Your data format*] -Norm [*Is the data standardized?*]
```
