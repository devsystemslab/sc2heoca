# Examples

### Prerequisites

1. Please install sc2heoca if you haven't done so already, using the [installation documentation](installation.md).

2. If you installed sc2heoca via conda, which is recommended, make sure to activate the environment before you run it.

```
conda activate sc2heoca
```

3. Download the sc2heoca models (updated models can be found in the model description section).

```
mkdir -p sc2heoca
```

4. Download and unpack the example data (e.g. example counts of lung COPD data, GSE186017):

```
tar xvzf sc2heoca_example_data.tgz
rm sc2heoca_example_data.tgz
```

###  Data

To run a full sc2heoca analysis you will need:

* reference model files (updated models can be found in the model description section)
* single-cell RNA-seq expression row counts (e.g. example counts of lung COPD data; GEO database - GSE186017)

This file structure should be present in the example data:

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
sc.pp.filter_cells(adata_query, min_genes=1000)
sc.pp.filter_cells(adata_query, max_genes=8000)
```
### Query new sample
```
from sc2heoca.query import Query

model_dir = "heoca_atlas"
query = Query(model_dir=model_dir, load_ref=True)
adata_query = query.run_scpoli(adata_query=adata_query, 
                               sample_name='GSM5628936_SCNPO2')
```
### Perform the label transfer
```
adata_query, adata_knn_graph = query.scpoli_label_transfer(adata_query)
```
### Plot query result UMAP
```
adata4plot = query.merge4plot(adata_query)
sc.pl.umap(adata4plot, color=['predict_level_2'], 
            palette=query.colorpalette,
            frameon=False, size=5)
```
![COPD query to reference plot](../figures/GSM5628936_SCNPO2.png)
(UMAP might differ slightly when newest model is used)

### Find DE genes to HEOCA

```
de_res = query.find_de_genes(adata_query)

```

### Off-target analysis
```
from sc2heoca.target import Target

target = Target('fetal', model_dir)
target.adata_latent_source.obs["celltype"] = 'na'
target_res = target.get_target(adata_query=adata_query, 
                            sample_name='sample_name')
```

### Maturation analysis

To run the maturation analysis, ensure that `R` is installed along with the R package [miloR](https://github.com/MarioniLab/miloR). The path to the `R` executable should also be provided.

```
path2r = '~/miniconda3/envs/sc2heoca/'
mat_res = Maturation(path2r, 'fetal', f'{model_dir}/yu_atlas_milo.rds')
```
