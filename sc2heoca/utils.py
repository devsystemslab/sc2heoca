import os
import numpy as np
import pandas as pd

import torch
import random
from . import PACKAGE_DIR


def seed_everything(TORCH_SEED):
    random.seed(TORCH_SEED)
    os.environ['PYTHONHASHSEED'] = str(TORCH_SEED)
    np.random.seed(TORCH_SEED)
    torch.manual_seed(TORCH_SEED)
    torch.cuda.manual_seed_all(TORCH_SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def load_colorpalette():
    color_file = os.path.join(PACKAGE_DIR, "db", "gut_scpoli_color.txt")

    plate_level_all = pd.read_csv(color_file, sep='\t', header=None, index_col=0)
    plate_level_all = plate_level_all.to_dict()[1]

    return plate_level_all
