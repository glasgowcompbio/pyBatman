import os
import errno
import pandas as pd
from IPython.display import display, HTML
import json

from constants import WD, DEFAULT_CONFIG, USER_CONFIG, DEFAULT_DATABASE, USER_DATABASE

def load_config():
    # the default config file is located in the same directory as this script
    current_dirname = os.path.dirname(os.path.realpath(__file__))
    default_config = os.path.join(current_dirname, DEFAULT_CONFIG)
    selected = load_or_default(USER_CONFIG, default_config)
    print 'Using configuration %s' % selected
    with open(selected) as f:
        config = json.load(f)
        return config

def get_db_path():
    selected = load_or_default(USER_DATABASE, DEFAULT_DATABASE)
    print 'Using database %s' % selected
    return selected

def load_or_default(user_file, default_path):
    user_path = os.path.join(WD, user_file)
    if os.path.isfile(user_path):
        return user_path
    else:
        return default_path

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