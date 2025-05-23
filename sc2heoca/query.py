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
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestNeighbors

import anndata
import scanpy as sc
from scarches.models.scpoli import scPoli

from . import PACKAGE_DIR
from sc2heoca.utils import seed_everything, load_colorpalette
from sc2heoca.de import get_matched_transcriptome, test_de_paired

def clear_genes(adata):
    gene_file = os.path.join(PACKAGE_DIR, "db", "hg_genes_clear.txt")

    clear_genes = pd.read_csv(gene_file, header=None)[0].tolist()
    sub_clear_genes = [i for i in clear_genes if i in adata.var.index.tolist()]
    adata = adata[:, sub_clear_genes]
    
    return adata

def init_sample(adata, empty_adata, sample_name=None):
    
    malat1 = adata.var_names.str.startswith('MALAT1')
    mito_genes = adata.var_names.str.startswith('MT-')
    rb_genes = adata.var_names.str.startswith(("RPS","RPL"))
    hb_genes = adata.var_names.str.contains('^HB[^(P)]')

    remove = np.add(mito_genes, malat1)
    remove = np.add(remove, rb_genes)
    remove = np.add(remove, hb_genes)
    keep = np.invert(remove)
    adata = adata[:,keep]

    adata = clear_genes(adata)

    adata = anndata.AnnData.concatenate(*[adata, empty_adata], join='outer', fill_value=0)

    adata.layers['counts']=adata.X
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    adata = adata[:,[i for i in empty_adata.var.index if i in adata.var.index]]

    adata.obs['tissue']='na'
    adata.obs['level_1']='na'
    adata.obs['level_2']='na'
    adata.obs['level_3']='na'

    if sample_name is not None:
        adata.obs['sample_id'] = sample_name
    else:
        adata.obs['sample_id']='sample2query'

    return adata

def find_uni_genes(auc_de_res_exp, auc_de_res_ctrl, cutoff):
    exp_genes = auc_de_res_exp[(auc_de_res_exp.auc>cutoff)&
                     (auc_de_res_exp.pvals_adj<0.01)&
                     (auc_de_res_exp.group=='query')].names.tolist()

    ctrl_genes = auc_de_res_ctrl[(auc_de_res_ctrl.auc>cutoff)&
                     (auc_de_res_ctrl.pvals_adj<0.01)&
                     (auc_de_res_ctrl.group=='query')].names.tolist()

    uni_genes = [i for i in exp_genes if i not in ctrl_genes]
    
    return uni_genes

class Query:
    def __init__(self, model_dir, load_ref=False):

        self.scpoli_model = f"{model_dir}/scpoli_model/"
        self.adata_latent_source = sc.read_h5ad(f"{model_dir}/adata_latent_source.h5ad")
        self.umap_model = pickle.load(open(f"{model_dir}/umap_model.sav", 'rb'))
        self.empty_adata = sc.read_h5ad(f"{model_dir}/empty.h5ad")
        self.colorpalette = load_colorpalette()

        if load_ref:
            self.adata = sc.read_h5ad(f"{model_dir}/gut_scpoli_integration.h5ad")

    def run_scpoli(self, adata_query, sample_name=None):
        seed_everything(0)
        
        adata_query = init_sample(adata_query, self.empty_adata, sample_name)

        scpoli_query = scPoli.load_query_data(
            adata=adata_query,
            reference_model=self.scpoli_model,
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

        #get latent representation of query data
        data_latent= scpoli_query.get_latent(
            adata_query, 
            mean=True
        )

        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = adata_query.obs.copy()

        #get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'
        
        adata_query.obsm['scpoli_latent']=adata_latent.to_df().to_numpy()
        
        return adata_query

    def scpoli_label_transfer(self, adata_query):

        que_embedding = self.umap_model.transform(adata_query.obsm['scpoli_latent'])
        adata_query.obsm['X_umap'] = que_embedding
        # adata_all = adata_ref.concatenate(adata_que, batch_key='query')

        # # predict tissue
        # knn = KNeighborsClassifier(n_neighbors=10)
        # knn.fit(self.adata_latent_source.to_df(), 
        #         self.adata_latent_source.obs.tissue)
        # adata_query.obs['predict_tissue'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # knn_res = pd.DataFrame(knn.predict_proba(adata_query.obsm['scpoli_latent']))
        # knn_res.columns=['prob_'+i for i in knn.classes_]
        # knn_res.index=adata_query.obs.index

        # adata_query.obs['predict_tissue'] = knn.predict(adata_query.obsm['scpoli_latent'])
        # adata_query.obs = pd.merge(adata_query.obs, knn_res, left_index=True, right_index=True)

        # predict level1
        knn = KNeighborsClassifier(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.level_1_late)
        adata_query.obs['predict_level_1'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # predict level2
        knn = KNeighborsClassifier(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.level_2_late)
        adata_query.obs['predict_level_2'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # # predict leiden
        # knn = KNeighborsClassifier(n_neighbors=10)
        # knn.fit(self.adata_latent_source.to_df(), 
        #         self.adata_latent_source.obs['leiden_10.0'])
        # adata_query.obs['predict_leiden'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # # predict leiden
        # knn = KNeighborsClassifier(n_neighbors=10)
        # knn.fit(self.adata_latent_source.to_df(), 
        #         self.adata_latent_source.obs['leiden_20.0'])
        # adata_query.obs['predict_leiden2'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # predict level3
        knn = KNeighborsClassifier(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.level_3_late)
        adata_query.obs['predict_level_3'] = knn.predict(adata_query.obsm['scpoli_latent'])

        # predict dist
        knn = NearestNeighbors(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df())
        knn_res = knn.kneighbors(adata_query.obsm['scpoli_latent'])
        mydist = pd.DataFrame(knn_res[0]).mean(1)
        adata_query.obs['mean_dist'] = mydist.tolist()

        # get neighbors
        knn = NearestNeighbors(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df())
        adata_knn_graph = knn.kneighbors_graph(adata_query.obsm['scpoli_latent'])

        # predict detail_tissue
        if "detail_tissue" in self.adata_latent_source.obs.columns:
            self.adata_latent_source.obs.detail_tissue = self.adata_latent_source.obs.detail_tissue.astype('str')
            knn = KNeighborsClassifier(n_neighbors=10)
            knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.detail_tissue)
            adata_query.obs['predict_detail_tissue'] = knn.predict(adata_query.obsm['scpoli_latent'])

        ## all str columns become object type TODO
        for i in adata_query.obs.columns:
            if adata_query.obs[i].dtype == 'object':
                adata_query.obs[i] = adata_query.obs[i].astype('category')

        return adata_query, adata_knn_graph

    def merge4plot(self, adata_query):
        merged_adata = anndata.AnnData.concatenate(*[adata_query, 
                                                     self.adata_latent_source], join='outer', fill_value=0)
        
        ##TODO magic way to remove var
        merged_adata = merged_adata[:,~merged_adata.var.index.isin([str(i) for i in range(10)])]

        return merged_adata

    def find_de_genes(self, adata_query):
        
        adata_query, adata_knn_graph = self.scpoli_label_transfer(adata_query)
        adata_bg = get_matched_transcriptome(adata_query, self.adata, adata_knn_graph)
        atlas_de_res = test_de_paired(adata_query.raw.to_adata(), adata_bg, num_threads=20)

        return atlas_de_res
