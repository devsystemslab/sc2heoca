# Examples

## Prerequisites

1. Please install sc2heoca if you haven't done so already, using the [installation documentation](installation.md).

2. If you installed sc2heoca via conda, which is recommended, make sure to activate the environment before you run it.

```
conda activate sc2heoca
```

3. Download the sc2heoca models.

```
mkdir -p sc2heoca
```

5. Download and unpack the example data:

```
tar xvzf sc2heoca_example_data.tgz
rm sc2heoca_example_data.tgz
```


##  Data

To run a full sc2heoca analysis you will need:

* reference model files.
* single-cell RNA-seq expression row counts.

These files are present in the example data for fibroblasts and for primary heart tissue:

```
$ tree sc2heoca_example_data/

sc2heoca_example_data/
├── heoce_mocel
│   ├── ...
├── scRNAseq
│   ├── mtx
│   ├── .
│   └── .
└── README.txt

3 directories, 17 files
```

Details for the different steps are described below.

### Load new sample
```
adata_query = sc.read_10x_mtx(path='Chan_NatCommun_2022', 
                              prefix='GSM5628936_SCNPO2-')
```
### Query new sample
```
model_dir = "heoca_atlas"
query = Query(model_dir=model_dir, load_ref=True)
adata_query = query.run_scpoli(adata_query=adata_query, 
                               sample_name='GSM5628936_SCNPO2')
```

### Plot query result UMAP
```
adata4plot = query.merge4plot(adata_query)
sc.pl.umap(adata4plot, color=['predict_level_2'], palette=query.colorpalette,
           frameon=False, size=5)
```

### Find DE genes to HECOA

```
de_res = query.find_de_genes(adata_query, tissue='lung')

```
