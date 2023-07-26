# sc2heoca
A method to map new scRNA-seq data to HEOCA

## Install
```
git clone git@github.com:devsystemslab/sc2heoca.git
pip install sc2heoca
```

## Download the model from Zenodo
```
mkdir -p heoca_scpoli_model.v1.0
wget https://zenodo.org/record/8185826/files/heoca_scpoli_model.v1.0.zip
tar xvzf heoca_scpoli_model.v1.0.zip -C heoca_scpoli_model.v1.0
rm heoca_scpoli_model.v1.0.zip
```
### All available models from HEOCA project
* [HEOCA model (all organoids)](https://zenodo.org/record/8185826/files/heoca_scpoli_model.v1.0.zip)
* [HIOCA model (intestine organoid)](https://zenodo.org/record/8185826/files/hioca_scpoli_model.v1.0.zip)
* [HLOCA model (lung organoid)](https://zenodo.org/record/8185826/files/hioca_scpoli_model.v1.0.zip)
* [HICA model (intestine tissues)](https://zenodo.org/record/8185826/files/hioca_scpoli_model.v1.0.zip)

## Run model

```
import scanpy as sc
from sc2heoca.sc2heoca import Query

# read sample
adata = sc.read_10x_mtx('Chan_NatCommun_2022/', prefix = 'GSM5628936_SCNPO2-')

model_dir = "heoca_scpoli_model.v1.0"
heoca_query = Query(model_dir, adata, 'GSM5628936_SCNPO2')
adata2 = heoca_query.run_scpoli()

sc.pl.umap(adata2, color='predict_level_2', frameon=False, size=5)

```
![](figures/GSM5628936_SCNPO2.png)
