# COMP0138-Metal-Binding-Site-Prediction

## Challenge

[Metal Binding Site Challenge](https://drive.google.com/drive/folders/1wQWuywtJPw70nzqjhN5r11F1Fk5DYawK) held by [Uniprot](https://www.uniprot.org/).

## Data Source

- [Source](https://ftp.ebi.ac.uk/pub/contrib/UniProt/prediction_challenges/1_metal_binding/)
- Format

  The data is in the format of [FASTA](https://en.wikipedia.org/wiki/FASTA_format).

  The annotation is in the format: Accession\<TAB\>Evidence\<TAB\>ChEBI-ID\<TAB\>Position.

## Directory tree

```txt
.
├── data/                   zipped training and test data in FASTA.
├── hyper-files/            the hyperparameter gird search results in csv.
├── hyper_tune/             the hyperparameter tuning pipeline.
├── label_encode/           containing class encodings in json, and the ChEBI-ID of the metal classes and metal-binding annotation file provided by UNiProt.
├── labels/                 the label files in npz in correspondence with the data files in data/
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
- download the labels via [this link](https://drive.google.com/drive/folders/17YSJRBH2Dx0wZo21pxW-LsG6NJQkiWcx?usp=sharing).

Each data file and its corresponding label file should share the same file name.

## Sequence Truncation and Embedding

- fill in the `embed.yaml` your configuration.

  ```yaml
  pLM: Ankh # or ProtT5
  data: # a list of file names of the data to be embedded and the label (use same name)
    - file_name1
    - file_name2
    ...

  data_dir: /path/to/data/ # the directory where the data is stored

  label_dir: /path/to/labels/ # the directory where the labels are stored

  embed_save_dir: /path/to/embeds/ # the directory where the model is saved

  label_save_dir: /path/to/truncated_labels/ # the directory where the labels are saved

  truncate: 512 # the maximum length of the input sequence

  ```

- run the notebook `embed.ipynb`.

## Model Training

- make a json file for your metal class encoding, the format is as follows:

  ```json
  {
      "0": [CHEBI-ID1, METAL-NAME1],
      "1": [CHEBI-ID2, METAL-NAME2],
      ...
  }

  ```

- fill in the `train.yaml` your configuration.

  ```yaml
      model: TFE # or CNN2L

      class_encode_path: /path/to/json # put your class encoding in json here

      truncated_label_path: /path/to/truncated_labels/ # put your truncated labels here, please split train positive and negative labels into separate files, name the pos one containing the keyword "pos" and "train" and the neg one containing the keyword "neg" and "train". The test label file should contain "test".

      truncated_embed_path: /path/to/truncated_embeds/ # put your truncated embeddings here, please split train positive and negative embeddings into separate files, name the pos one containing the keyword "pos" and "train" and the neg one containing the keyword "neg" and "train". The test embedding file should contain "test".

      CNN2L:
      hidden_channel: 128
      hidden_layer_num: 2
      kernel_size: 17
      lr: 0.001
      label_weight: [0.228, 5.802]
      batch_size: 16

      TFE:
      hidden_dim: 128
      num_encoder_layers: 2
      num_heads: 4
      dropout: 0.2
      lr: 0.0007585775750291836
      label_weight: [0.78324, 8.46187]
      batch_size: 16

  ```

- run the notebook `train.ipynb`.
- the last cell performs a threshold tuning, please suggest the model ckpt and the file path to save the tuning results.

## Hyperparameter Tuning

- fill in the `hyper_tune\hypertune.yaml` your configuration.

  ```yaml

      model: TFE # or CNN2L

      class_encode_path: /path/to/json # put your class encoding in json here

      truncated_label_path: /path/to/truncated_labels/

      truncated_embed_path: /path/to/truncated_embeds/

      batch_size: 16

      CNN2L:
      hidden_channel: [64, 128]
      hidden_layer_num: [2]
      kernel_size: [13, 15]
      lr: [0.001, 0.0005]
      label_weight: [[0.228, 5.802], [0.78324, 8.46187]]

      TFE:
      hidden_dim: [64, 128]
      num_encoder_layers: [2, 3]
      num_heads: [2, 4]
      dropout: [0.1, 0.2]
      lr: [0.0007585775750291836, 0.001]
      label_weight: [[0.228, 5.802], [0.78324, 8.46187]]


  ```

- run the notebook `hyper_tune\hyper_tune.ipynb`.

> for now the tuning results are stored in txt and requires a further parsing. We will update the notebook to save the results in csv.

## Model Inference

- run the notebook `inference.ipynb`, specifying the sequence in str and probability thresholds in list of float in the second last cell.

## Dependencies

python in dev: 3.10.10

Libs:

- scikit-learn
- numpy
- pandas
- biopython
- torch
- torch lightning
- h5py
- tqdm
- ankh
- transformers
- matplotlib
- seaborn (additionally for visualization)
