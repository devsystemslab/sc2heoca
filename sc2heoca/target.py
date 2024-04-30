import os
import pickle

import numpy as np
import pandas as pd

from sklearn.neighbors import KNeighborsClassifier

import anndata
import scanpy as sc
from scarches.models.scpoli import scPoli

from . import PACKAGE_DIR

def load_colorpalette():
    color_file = os.path.join(PACKAGE_DIR, "db", "gut_scpoli_color.txt")

    plate_level_all = pd.read_csv(color_file, sep='\t', header=None, index_col=0)
    plate_level_all = plate_level_all.to_dict()[1]

    return plate_level_all

class Target:
    def __init__(self, model, model_dir):
        
        self.model = model
        self.scpoli_model = f"{model_dir}/scpoli_model/"
        self.adata_latent_source = sc.read_h5ad(f"{model_dir}/adata_latent_source.h5ad")
        self.umap_model = pickle.load(open(f"{model_dir}/umap_model.sav", 'rb'))
        self.empty_adata = sc.read_h5ad(f"{model_dir}/empty.h5ad")
        self.colorpalette = load_colorpalette()

        if self.model == 'adult':
            self.adata_latent_source.obs.organ_tissue = self.adata_latent_source.obs.organ_tissue.astype('str')
            self.adata_latent_source.obs.loc[self.adata_latent_source.obs.organ_tissue.isin(['Small_Intestine','Large_Intestine']), 'organ_tissue'] = 'Intestine'

    def __correct_tissue(self, map_res, on_tissue, cutoff=0.3):
        ref_cellnum = self.adata_latent_source.obs.groupby(['celltype','tissue']).count().reset_index().set_index('tissue').pivot(columns='celltype', values='n_genes').T
        ref_cellpect = ref_cellnum.apply(lambda x: x/ref_cellnum.sum(1))

        correct_celltypes = {}
        for tissue in ref_cellpect:
            correct_celltypes[tissue] = ref_cellpect[tissue][ref_cellpect[tissue]>cutoff].index.tolist()

        if on_tissue in correct_celltypes:
            for celltype in correct_celltypes[on_tissue]:
                # map_res.loc[map_res.predict_celltype==celltype, 
                #             f'prob_{on_tissue}'] = map_res.loc[map_res.predict_celltype==celltype, 
                #                                             f'predict_celltype_prob']
                map_res.loc[map_res.predict_celltype==celltype, 
                            f'prob_{on_tissue}'] = map_res.loc[map_res.predict_celltype==celltype, 
                                         [f'predict_celltype_prob', f'prob_{on_tissue}']].max(1)
                
        map_res['predict_tissue'] = [i[5:] for i in map_res[[i for i in map_res.columns if i.startswith('prob_')]].idxmax(1)]

        return map_res

    def __run_scpoli(self, adata_query, sample_name, on_tissue=None):
        adata_query = adata_query.raw.to_adata()
        adata_query = anndata.AnnData.concatenate(*[adata_query, self.empty_adata], join='outer', fill_value=0)
        adata_query = adata_query[:,[i for i in self.empty_adata.var.index if i in adata_query.var.index]]
        adata_query.obs['celltype'] = 'na'

        if self.model == 'adult':
            adata_query.obs['Batch'] = sample_name
        elif self.model == 'fetal':
            adata_query.obs['Sample'] = sample_name

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
            n_epochs=10,
            pretraining_epochs=8,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=10,
            alpha_epoch_anneal=100
        )

        if self.model == 'adult':
            results_dict = scpoli_query.classify(adata_query.X, adata_query.obs["Batch"].values)
        elif self.model == 'fetal':
            results_dict = scpoli_query.classify(adata_query.X, adata_query.obs["Sample"].values)

        #get latent representation of query data
        if self.model == 'adult':
            data_latent= scpoli_query.get_latent(
                adata_query.X, 
                adata_query.obs["Batch"].values,
                mean=True
            )
        elif self.model == 'fetal':
            data_latent= scpoli_query.get_latent(
                adata_query.X, 
                adata_query.obs["Sample"].values,
                mean=True
            )
        adata_latent = sc.AnnData(data_latent)
        adata_latent.obs = adata_query.obs.copy()

        #get label annotations
        adata_latent.obs['celltype_pred'] = results_dict['celltype']['preds'].tolist()
        adata_latent.obs['celltype_uncert'] = results_dict['celltype']['uncert'].tolist()
        adata_latent.obs['classifier_outcome'] = (
            adata_latent.obs['celltype_pred'] == adata_latent.obs['celltype']
        )

        #get prototypes
        labeled_prototypes = scpoli_query.get_prototypes_info()
        labeled_prototypes.obs['study'] = 'labeled prototype'
        unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
        unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

        que_embedding = self.umap_model.transform(adata_latent.X)
        adata_latent.obsm['X_umap'] = que_embedding

        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.celltype)

        adata_latent.obs['predict_celltype'] = knn.predict(adata_latent.to_df())
        knn_res = pd.DataFrame(knn.predict_proba(adata_latent.to_df()))
        adata_latent.obs['predict_celltype_prob'] = knn_res.max(1).tolist()

        knn = KNeighborsClassifier(n_neighbors=100)
        knn.fit(self.adata_latent_source.to_df(), 
                self.adata_latent_source.obs.tissue)

        knn_res = pd.DataFrame(knn.predict_proba(adata_latent.to_df()))
        knn_res.columns=['prob_'+i for i in knn.classes_]
        knn_res.index=adata_latent.to_df().index

        adata_latent.obs['predict_tissue'] = knn.predict(adata_latent.to_df())
        adata_latent.obs = pd.merge(adata_latent.obs, knn_res, left_index=True, right_index=True)

        if on_tissue != None:
            map_res = self.__correct_tissue(adata_latent.obs, on_tissue)
        else:
            map_res = adata_latent.obs
        return map_res
    
    def __get_confidence(self, map_res, on_tissue=None):
        tissue_order = ["intestine", "lung", "liver",  "pancreas", "prostate", 
                        "salivarygland", "stomach", "biliarysystem",  "thyroid"]

        conf_res = map_res.groupby(['predict_tissue']).count()['prob_liver'].reset_index()

        raw_res = pd.DataFrame({'predict_tissue':tissue_order})
        conf_res = pd.merge(raw_res, conf_res, on='predict_tissue', how='left').fillna(0)
        conf_res = conf_res.rename(columns={'prob_liver':'cell_num'})
        conf_res['pect'] = conf_res.cell_num/conf_res.cell_num.sum()
        if on_tissue != None:
            conf_res['target'] = ['on' if i==on_tissue else 'off' for i in conf_res.predict_tissue]
        
        sample_prob=[]
        for index, row in map_res.iterrows():
            sample_prob.append(row[f"prob_{row['predict_tissue']}"])

        map_res['mean_confidence'] = sample_prob

        map_res = map_res.groupby(['predict_tissue']).mean()[['mean_confidence']].reset_index()

        conf_res = pd.merge(conf_res, map_res, on='predict_tissue', how='left')
        
        conf_res = conf_res.rename(columns={'predict_tissue':'target_tissue'}).set_index('target_tissue')
        
        return conf_res

    def get_target(self, adata_query, sample_name, on_tissue=None):
        map_res = self.__run_scpoli(adata_query, sample_name, on_tissue)
        conf_res = self.__get_confidence(map_res, on_tissue)

        return map_res, conf_res
    