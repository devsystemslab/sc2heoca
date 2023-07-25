# sc2heoca
A method to map new scRNA-seq data to HEOCA

## Install
```
git clone 
pip install 
```

## Download the model from Zenodo
```
mkdir -p heoca_scpoli_model.v1.0
wget https://zenodo.org/record/8181496/files/heoca_scpoli_model.v1.0.zip
tar xvzf heoca_scpoli_model.v1.0.zip -C heoca_scpoli_model.v1.0
rm heoca_scpoli_model.v1.0.zip
```

## Run model

```
from sc2heoca import loadmodel, mapping
```
