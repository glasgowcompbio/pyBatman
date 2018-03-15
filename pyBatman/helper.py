import os
import errno
import pandas as pd
from IPython.display import display, HTML
import json
import pprint

from constants import WD, CONFIG, DATABASE, PATTERN

def load_config():
    user_path = os.path.join(WD, CONFIG)
    print 'Using configuration %s' % user_path
    with open(user_path) as f:
        config = json.load(f)
        config['pattern'] = config['pattern'] if 'pattern' in config else PATTERN
        config['verbose'] = True if config['verbose'] == 'True' else False
        config['correct_spectra'] = True if config['correct_spectra'] == 'True' else False

        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(config)

        return config

def get_db_path():
    db_file = os.path.join(WD, DATABASE)
    print 'Using database %s' % db_file
    return db_file


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