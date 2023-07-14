# CTISL:a dynamic stacking multi-class classification approach for identifying cell types from single-cell RNA-sequencing data
## Requirements
python 3.7.11
scikit-learn 1.0.2
scanpy 1.8.2
numpy 1.19.1
pandas 1.3.5

## Train CTISL on Intra-dataset as follows:
 python Intra_train.py -Name '10Xv2' -Fileform 'h5ad' -Norm True
 Parameters and description are as follows:
 -Name: Intra-dataset(for example:10Xv2,10Xv3...)
 -Fileform: The gene expression matrix file can be in either 'csv' or 'h5ad' format.
 -Norm: If it is raw data, it needs to be normalized with 'norm=True'

## Train CTISL on Inter-data,Cross-batch,Cross-species as follows;
 python cross_train.py -Sourcename 'dentritic_batch_1' -Targetname 'dentritic_batch_2' -Fileform 'csv' -Norm False
 Parameters and description are as follows:
 -Sourcename: Training set name.In the example above, "dentritic_batch_1" is used as the training set in the cross-batch experiment.
 -Targetname: Testing set name.In the example above, "dentritic_batch_2" is used as the testing set in the cross-batch experiment.
 -Fileform: The gene expression matrix file can be in either 'csv' or 'h5ad' format. 
 -Norm: If it is raw data, it needs to be normalized with 'norm=True'
## Data Source

The PBMC datasets were downloaded from the Broad Institute Single Cell portal https://portals.broadinstitute.org/single_cell/study/SCP424/single-cellcomparisonpbmc-data
The Pancreas datasets were downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/
The Airway datasets were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102580
The Dendritic datasets were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820
The Retina datasets were downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/mouse/retina/
