import os
import errno
import pandas as pd
from IPython.display import display, HTML
import json

def load_config(config_filename):
    with open('config.json') as f:
        config = json.load(f)
    return config

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def sub_dir_path(d):
    return filter(os.path.isdir, [os.path.join(d, f) for f in os.listdir(d)])