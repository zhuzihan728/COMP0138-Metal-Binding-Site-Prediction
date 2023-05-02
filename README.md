# COMP0138-Metal-Binding-Site-Prediction

## Challenge
[Metal Binding Site Challenge](https://drive.google.com/drive/folders/1wQWuywtJPw70nzqjhN5r11F1Fk5DYawK) held by [Uniprot](https://www.uniprot.org/).

## Data Source
- [Source](https://ftp.ebi.ac.uk/pub/contrib/UniProt/prediction_challenges/1_metal_binding/)
- Format 
  
    The data is in the format of [FASTA](https://en.wikipedia.org/wiki/FASTA_format). 

    The annotation is in the format: Accession\<TAB\>Evidence\<TAB\>ChEBI-ID\<TAB\>Position. 
 
 ## Directory tree
 

```
.
├── data/                   zipped training and test data in FASTA.
├── hyper-files/            the hyperparameter gird search results in csv.
├── hyper_tune/             the hyperparameter tuning pipeline.
├── label_encode/           containing class encodings in json, and the ChEBI-ID of the metal classes and metal-binding annotation file provided by UNiProt.
├── models/                 trained models
├── result_analysis/        notebooks for model performance analysis and visualization
├── thres_tune/             threshold tuning results in csv
├── embed.ipynb             the sequence truncation and embedding pipelines
├── helper_fn_short_val.py  helper functions for metric calculation and others
├── inference.ipynb         an inference demo using TFE-11
├── train.ipynb             model training pipeline
├── *.yaml                  configuration files
└── *.md                    README
```


## Data Preparation
 - unzip the train and test data in fasta under the `data/` folder.
 - download the labels via [this link](https://drive.google.com/drive/folders/17YSJRBH2Dx0wZo21pxW-LsG6NJQkiWcx?usp=sharing)


## Sequence Truncation and Embedding 


