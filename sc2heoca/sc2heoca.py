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


model_dir = "new_20230705"


def read_sample(input_dir):
    adata = sc.read_10x_mtx(input_dir)

    return adata

def init_sample(adata0, name=None, empty_file=None):
    
    malat1 = adata0.var_names.str.startswith('MALAT1')
    mito_genes = adata0.var_names.str.startswith('MT-')
    rb_genes = adata0.var_names.str.startswith(("RPS","RPL"))

    remove = np.add(mito_genes, malat1)
    remove = np.add(remove, rb_genes)
    keep = np.invert(remove)

    adata0 = adata0[:,keep]

    empty_adata = sc.read_h5ad(empty_file)

    adata0 = anndata.AnnData.concatenate(*[adata0, empty_adata], join='outer', fill_value=0)

    if name is not None:
        adata0.obs['sample_id']=name
        
    adata0.layers['counts']=adata0.X
    sc.pp.normalize_total(adata0, target_sum=1e4)
    sc.pp.log1p(adata0)
    adata0.raw = adata0

    adata0 = adata0[:,[i for i in empty_adata.var.index if i in adata0.var.index]]

    adata0.obs['level_1']='na'
    adata0.obs['level_2']='na'
    adata0.obs['level_3']='na'
    
    if name is not None:
        adata0.obs['publication']='FMI_UP_2023'
        adata0.obs['time']= [i[1] for i in adata0.obs.sample_id.str.split('_')]
        adata0.obs['protocol'] = 'na'
        adata0.obs['replicate']= [i[2] for i in adata0.obs.sample_id.str.split('_')]
        adata0.obs['condition0']= ['_'.join(i[:2]) for i in adata0.obs.sample_id.str.split('_')]
        adata0.obs['condition']= [i[0] for i in adata0.obs.sample_id.str.split('_')]
    return adata0



def run_scpoli(adata, ref_path):
    ref_path = f'{project_dir}/results/{results_dir}/gut_scpoli_model/scpoli_model/'

    scpoli_query = scPoli.load_query_data(
        adata=adata,
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

    results_dict = scpoli_query.classify(adata.X, adata000.obs["sample_id"].values)

    return results_dict

def load_model(model_dir):

    adata_latent_source3 = sc.read_h5ad(f"{project_dir}/results/{results_dir}/gut_scpoli_model/adata_latent_source.h5ad")

    late_annot = pd.read_csv(f"{project_dir}/results/{results_dir}/gut_scpoli_model/gut_scpoli_late_annot.txt",
                index_col=0, sep='\t')

    adata_latent_source3.obs = pd.merge(adata_latent_source3.obs, late_annot, 
                                        left_index=True, right_index=True, how='left')


def query_sample():
    #get latent representation of query data
    data_latent= scpoli_query.get_latent(
        adata000.X, 
        adata000.obs["sample_id"].values,
        mean=True
    )

    adata_latent3 = sc.AnnData(data_latent)
    adata_latent3.obs = adata000.obs.copy()

    #get label annotations
    adata_latent3.obs['level_2_pred'] = results_dict['level_2']['preds'].tolist()
    adata_latent3.obs['level_2_uncert'] = results_dict['level_2']['uncert'].tolist()
    adata_latent3.obs['classifier_outcome'] = (
        adata_latent3.obs['level_2_pred'] == adata_latent3.obs['level_2']
    )

    #get prototypes
    labeled_prototypes = scpoli_query.get_prototypes_info()
    labeled_prototypes.obs['study'] = 'labeled prototype'
    unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
    unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

    umap_model = pickle.load(open(f"{project_dir}/results/{results_dir}/gut_scpoli_model/umap_model.sav", 'rb'))
    que_embedding = umap_model.transform(adata_latent3.X)
    adata_latent3.obsm['X_umap'] = que_embedding
    # adata_all = adata_ref.concatenate(adata_que, batch_key='query')

    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source3.to_df(), 
        adata_latent_source3.obs.level_2_late)

    adata_latent3.obs['predict_level_2_late'] = knn.predict(adata_latent3.to_df())

    knn = KNeighborsClassifier(n_neighbors=100)
    knn.fit(adata_latent_source3.to_df(), 
        adata_latent_source3.obs.tissue)

    # adata_latent30 = adata_latent3[adata_latent3.obs.Batch==sample]
    knn_res = pd.DataFrame(knn.predict_proba(adata_latent3.to_df()))
    knn_res.columns=['prob_'+i for i in knn.classes_]
    knn_res.index=adata_latent3.to_df().index

    adata_latent3.obs['predict_tissue'] = knn.predict(adata_latent3.to_df())
    adata_latent3.obs = pd.merge(adata_latent3.obs, knn_res, left_index=True, right_index=True)




