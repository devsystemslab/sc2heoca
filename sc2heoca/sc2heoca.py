import os
from copy import copy
import sys
import math
import pickle

import numpy as np
import pandas as pd
from sklearn import preprocessing
from scipy.io import mmread
from scipy import sparse
from scarches.models.scpoli import scPoli
from sklearn.neighbors import KNeighborsClassifier

import scanpy as sc
import scanpy as sc
import anndata


def init_sample(adata0, name=None):
    
    malat1 = adata0.var_names.str.startswith('MALAT1')
    mito_genes = adata0.var_names.str.startswith('MT-')
    rb_genes = adata0.var_names.str.startswith(("RPS","RPL"))

    remove = np.add(mito_genes, malat1)
    remove = np.add(remove, rb_genes)
    keep = np.invert(remove)
    adata0 = adata0[:,keep]

    adata0 = anndata.AnnData.concatenate(*[adata0, empty_adata], join='outer', fill_value=0)

    if name is not None:
        adata0.obs['sample_id']=name
    else:
        adata0.obs['sample_id']='sample2query'

    adata0.layers['counts']=adata0.X
    sc.pp.normalize_total(adata0, target_sum=1e4)
    sc.pp.log1p(adata0)
    adata0.raw = adata0

    adata0 = adata0[:,[i for i in empty_adata.var.index if i in adata0.var.index]]

    adata0.obs['level_1']='na'
    adata0.obs['level_2']='na'
    adata0.obs['level_3']='na'

    return adata0



def run_scpoli(adata0, ref_path=scpoli_model, umap_model=umap_model):

    scpoli_query = scPoli.load_query_data(
        adata=adata0,
        reference_model=ref_path,
        labeled_indices=[],
    )

    early_stopping_kwargs = {
        "early_stopping_metric": "val_prototype_loss",
        "mode": "min",
        "threshold": 0,
        "patience": 20,
        "reduce_lr": True,
        "lr_patience": 13,
        "lr_factor": 0.1,
    }

    scpoli_query.train(
        n_epochs=5,
        pretraining_epochs=4,
        early_stopping_kwargs=early_stopping_kwargs,
        eta=10,
        alpha_epoch_anneal=100
    )

    results_dict = scpoli_query.classify(adata.X, adata0.obs["sample_id"].values)

    #get latent representation of query data
    data_latent= scpoli_query.get_latent(
        adata0.X, 
        adata0.obs["sample_id"].values,
        mean=True
    )

    adata_latent = sc.AnnData(data_latent)
    adata_latent.obs = adata0.obs.copy()

    #get label annotations
    # adata_latent.obs['level_2_pred'] = results_dict['level_2']['preds'].tolist()
    # adata_latent.obs['level_2_uncert'] = results_dict['level_2']['uncert'].tolist()
    # adata_latent.obs['classifier_outcome'] = (
    #     adata_latent.obs['level_2_pred'] == adata_latent.obs['level_2']
    # )

    #get prototypes
    labeled_prototypes = scpoli_query.get_prototypes_info()
    labeled_prototypes.obs['study'] = 'labeled prototype'
    unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
    unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

    que_embedding = umap_model.transform(adata_latent.X)
    adata_latent.obsm['X_umap'] = que_embedding
    # adata_all = adata_ref.concatenate(adata_que, batch_key='query')

    # predict tissue
    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source.to_df(), 
        adata_latent_source.obs.tissue)
    adata_latent.obs['predict_tissue'] = knn.predict(adata_latent.to_df())

    knn_res = pd.DataFrame(knn.predict_proba(adata_latent.to_df()))
    knn_res.columns=['prob_'+i for i in knn.classes_]
    knn_res.index=adata_latent.to_df().index

    adata_latent.obs['predict_tissue'] = knn.predict(adata_latent.to_df())
    adata_latent.obs = pd.merge(adata_latent.obs, knn_res, left_index=True, right_index=True)

    # predict level1
    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source.to_df(), 
        adata_latent_source.obs.level_1_late)
    adata_latent.obs['predict_level_1_late'] = knn.predict(adata_latent.to_df())

    # predict level2
    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source.to_df(), 
        adata_latent_source.obs.level_2_late)
    adata_latent.obs['predict_level_2_late'] = knn.predict(adata_latent.to_df())

    # predict level3
    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source.to_df(), 
        adata_latent_source.obs.level_3_late)
    adata_latent.obs['predict_level_3_late'] = knn.predict(adata_latent.to_df())

    return adata_latent

def load_model(model_dir):
    global scpoli_model, adata_latent_source, umap_model, empty_file

    scpoli_model = f"{model_dir}/scpoli_model/"
    adata_latent_source = sc.read_h5ad(f"{model_dir}/adata_latent_source.h5ad")
    umap_model = pickle.load(open(f"{model_dir}/umap_model.sav", 'rb'))
    empty_adata = sc.read_h5ad("{model_dir}/empty.h5ad")

def scquery(adata):
    adata0 = init_sample()
    adata_latent = run_scpoli(adata0, ref_path=scpoli_model, umap_model=umap_model)

    return adata_latent
    
