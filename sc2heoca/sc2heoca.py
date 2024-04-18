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
from wilcoxauc import wilcoxauc,top_markers

import anndata
import scanpy as sc
from scarches.models.scpoli import scPoli

from . import PACKAGE_DIR

def load_colorplate():
    color_file = os.path.join(PACKAGE_DIR, "db", "gut_scpoli_color.txt")

    plate_level_all = pd.read_csv(color_file, sep='\t', header=None, index_col=0)
    plate_level_all = plate_level_all.to_dict()[1]

    return plate_level_all

def clear_genes(adata):
    gene_file = os.path.join(PACKAGE_DIR, "db", "hg_genes_clear.txt")

    clear_genes = pd.read_csv(gene_file, header=None)[0].tolist()
    sub_clear_genes = [i for i in clear_genes if i in adata.var.index.tolist()]
    adata = adata[:, sub_clear_genes]
    
    return adata

def init_sample(adata, empty_adata):
    
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
        

        self.colorplate = load_colorplate()

        if load_ref:
            self.adata = sc.read_h5ad(f"{model_dir}/gut_scpoli_integration.h5ad")
    
    def run_scpoli(self, adata_query, sample_name):

        if sample_name is not None:
            adata_query.obs['sample_id'] = sample_name
        else:
            adata_query.obs['sample_id']='sample2query'

        adata = init_sample(adata_query, self.empty_adata)

        scpoli_query = scPoli.load_query_data(
            adata=adata,
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

        results_dict = scpoli_query.classify(adata.X, adata.obs["sample_id"].values)

        #get latent representation of query data
        data_latent= scpoli_query.get_latent(
            adata.X, 
            adata.obs["sample_id"].values,
            mean=True
        )

        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = adata.obs.copy()

        #get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

        que_embedding = self.umap_model.transform(adata_latent.X)
        adata_latent.obsm['X_umap'] = que_embedding
        # adata_all = adata_ref.concatenate(adata_que, batch_key='query')

        # predict tissue
        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
            self.adata_latent_source.obs.tissue)
        adata_latent.obs['predict_tissue'] = knn.predict(adata_latent.to_df())

        knn_res = pd.DataFrame(knn.predict_proba(adata_latent.to_df()))
        knn_res.columns=['prob_'+i for i in knn.classes_]
        knn_res.index=adata_latent.to_df().index

        adata_latent.obs['predict_tissue'] = knn.predict(adata_latent.to_df())
        adata_latent.obs = pd.merge(adata_latent.obs, knn_res, left_index=True, right_index=True)

        # predict level1
        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
            self.adata_latent_source.obs.level_1_late)
        adata_latent.obs['predict_level_1'] = knn.predict(adata_latent.to_df())

        # predict level2
        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
            self.adata_latent_source.obs.level_2_late)
        adata_latent.obs['predict_level_2'] = knn.predict(adata_latent.to_df())

        # predict level3
        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
            self.adata_latent_source.obs.level_3_late)
        adata_latent.obs['predict_level_3'] = knn.predict(adata_latent.to_df())

        # predict dist
        knn = NearestNeighbors(n_neighbors=10)
        knn.fit(self.adata_latent_source.to_df())
        knn_res = knn.kneighbors(adata_latent.to_df())
        mydist = pd.DataFrame(knn_res[0]).mean(1)

        adata.obs['predict_tissue']=adata_latent.obs.predict_tissue
        adata.obs['predict_level_1']=adata_latent.obs.predict_level_1
        adata.obs['predict_level_2']=adata_latent.obs.predict_level_2
        adata.obs['predict_level_3']=adata_latent.obs.predict_level_3
        adata.obsm['X_umap'] = adata_latent.obsm['X_umap']
        adata.obs['mean_dist'] = mydist.tolist()

        return adata
    
    def merge4plot(self, adata_query):
        merged_adata = anndata.AnnData.concatenate(*[adata_query, 
                                                     self.adata_latent_source], join='outer', fill_value=0)
        
        ##TODO magic way to remove var
        merged_adata = merged_adata[:,~merged_adata.var.index.isin([str(i) for i in range(10)])]

        return merged_adata

    def __find_de_genes(self, tissue, adata_query, detail_tissue=None, method='wilcoxauc'):
        if detail_tissue is None:
            adata_subset = self.adata[(self.adata.obs.tissue==tissue)].copy()
        else:
            adata_subset = self.adata[(self.adata.obs.tissue==tissue)&(self.adata.obs.detail_tissue==detail_tissue)].copy()
            
        adata_subset.obs['sample_state'] = 'atlas'
        adata_query.obs['sample_state'] = 'query'
        
        adata_merged = anndata.AnnData.concatenate(*[adata_subset, adata_query], join='outer', fill_value=0)
        
        de_res = {}

        # adata_merged.obs.predict_level_2 = adata_merged.obs.predict_level_2.astype('category')

        for celltype in adata_subset.obs.level_2.astype('str').unique():
            adata_query_counts = adata_query.obs.predict_level_2.value_counts()

            if celltype in adata_query_counts and adata_query_counts[celltype]>=5:

                adata_merged_sub = adata_merged[(adata_merged.obs.level_2_late==celltype)|
                                            (adata_merged.obs.predict_level_2==celltype)].copy()

                if method == 'wilcoxauc':
                    auc_res = wilcoxauc(adata_merged_sub, group_name='sample_state')
                    de_res[f'query_{celltype}'] = auc_res[(auc_res.group=='query')&(auc_res.auc>0.6)&(auc_res.pvals_adj<0.01)]['names'].tolist()

                else:
                    sc.tl.rank_genes_groups(adata_merged_sub, 'sample_state', method='wilcoxon',key_added = "wilcoxon")
                    # de_res[f'atlas_{celltype}'] = sc.get.rank_genes_groups_df(adata_merged_sub, group='atlas', key='wilcoxon')['names']
                    # de_res[f'query_{celltype}'] = sc.get.rank_genes_groups_df(adata_merged_sub, group='query', key='wilcoxon')['names']
                    de_res[f'query_{celltype}'] = sc.get.rank_genes_groups_df(adata_merged_sub, group='query', 
                                                                            key='wilcoxon', log2fc_min=1, 
                                                                            pval_cutoff=0.01)['names'].tolist()

        return de_res

    def find_de_genes(self, tissue, adata_query, adata_ctrl=None, method='wilcoxauc', detail_tissue=None):
        if adata_ctrl is None:
            de_res = self.__find_de_genes(tissue, adata_query, method=method, detail_tissue=detail_tissue)
            de_res_final = pd.DataFrame.from_dict(de_res, orient='index').T
        else:
            de_res1 = self.__find_de_genes(tissue, adata_query, method=method, detail_tissue=detail_tissue)
            de_res2 = self.__find_de_genes(tissue, adata_ctrl, method=method, detail_tissue=detail_tissue)

            uni_genes={}
            for i in de_res1:
                if i in de_res2:
                    uni_genes[i] = [j for j in de_res1[i] if j not in de_res2[i]]
                    
            de_res_final = pd.DataFrame.from_dict(uni_genes, orient='index').T
    
        return de_res_final
    