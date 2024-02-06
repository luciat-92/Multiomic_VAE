import os
from pathlib import Path
from utils.utils import safe_create_dir, str2path

if os.path.exists('/group/iorio/'):
    GROUP_DIR = str2path('/group/iorio')
elif os.path.exists('/Volumes/iorio/'):
    GROUP_DIR = str2path('/Volumes/iorio')
    
ROOT_DIR = GROUP_DIR / 'lucia' / 'Multiomic_VAE'
DATA_DIR = ROOT_DIR / 'data'
RESULTS_DIR = ROOT_DIR / 'experiments'

RAW_DATA_DIR = GROUP_DIR / 'Datasets'
PERSONAL_DATA_DIR = GROUP_DIR / 'lucia' / 'datasets'
safe_create_dir(DATA_DIR)
safe_create_dir(RESULTS_DIR)