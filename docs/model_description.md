
## Prepare reference atlas and model

* Download the reference atlas from Zenodo
```
mkdir -p heoca_atlas
cd heoca_atlas
wget https://zenodo.org/records/10977447/files/gut_scpoli_integration.h5ad
```

* Download the reference model from Zenodo
```
wget https://zenodo.org/record/8186773/files/heoca_scpoli_model.v1.0.zip
tar xvzf heoca_scpoli_model.v1.0.zip -C heoca_atlas
rm heoca_scpoli_model.v1.0.zip
```

* Reference atlas and available reference models from HEOCA project

    - [HEOCA atlas (all organoids)](https://zenodo.org/records/10977447/files/gut_scpoli_integration.h5ad)
    - [HEOCA model (all organoids)](https://zenodo.org/record/8186773/files/heoca_scpoli_model.v1.0.zip)
    - [HIOCA model (intestine organoid)](https://zenodo.org/record/8186773/files/hioca_scpoli_model.v1.0.zip)
    - [HLOCA model (lung organoid)](https://zenodo.org/record/8186773/files/hioca_scpoli_model.v1.0.zip)
    - [HICA model (intestine tissue)](https://zenodo.org/record/8186773/files/hica_scpoli_model.v1.0.zip)
