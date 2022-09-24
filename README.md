# Is-PTS1
A tool for detecting the PTS1 signal in a fasta sequence and validate the putative peroxisomal location

To use the PTS_validate.py script, you first need to download and install:

#### 1) seqvec 
Instructions are available here: https://github.com/Rostlab/SeqVec.

or 

`pip install seqvec`

https://pypi.org/project/seqvec/

#### 2) UniRep and the related weight files, in this case we used the 1900_weights.

https://github.com/churchlab/UniRep

Make sure you download the 1900_weights directory and place it together with the other files in this repository.

For example you can first install awscli: 

`pip install awscli`

Then download the weights with

`aws s3 sync --no-sign-request --quiet s3://unirep-public/1900_weights/ 1900_weights`


#### 3) A pre-computed model 'LR_model2.sav'/'SVM_model.sav' is also required (uncomment the selected model in the script). 

#### 4) Additional requirements

Suggested packages:

- `numpy 1.17.2`
- `biopython 1.77`
- `tensorflow 1.14`
- `pandas 0.25.1`
- `scikit-learn 0.22`
- `seqvec 0.4.1`
- `scipy 1.4.1`
- `overrides 3.1.0`

