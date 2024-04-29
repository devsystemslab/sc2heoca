import os
from copy import copy

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from . import PACKAGE_DIR

import torch
import random



def seed_everything(TORCH_SEED):
    random.seed(TORCH_SEED)
    os.environ['PYTHONHASHSEED'] = str(TORCH_SEED)
    np.random.seed(TORCH_SEED)
    torch.manual_seed(TORCH_SEED)
    torch.cuda.manual_seed_all(TORCH_SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

class Similarity:
    def __init__(self, r_path, model_file):
        self.r_path = r_path

        os.environ["R_HOME"] = self.r_path
        from rpy2.robjects.packages import importr
        base = importr('base')
        miloR = importr('miloR')
        
        self.milo_model = base.readRDS(f"{model_file}")

    def run_milo(self, adata, sample_name):

        os.environ["R_HOME"] = self.r_path
        import rpy2.robjects as robjects
        import anndata2ri
        anndata2ri.activate()

        adata_sce = anndata2ri.py2rpy(adata)

        r_source = robjects.r['source']

        milo_file = os.path.join(PACKAGE_DIR, "r", "run_milo.R")
        r_source(milo_file)

        r_getname = robjects.globalenv['run_milo']

        res = r_getname(self.milo_model, adata_sce, 'cell_ontology_class')
        res['sample']=sample_name
        
        return res
    
    