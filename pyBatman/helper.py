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


# https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_intervals(intervals):
    """
    A simple algorithm can be used:
    1. Sort the intervals in increasing order
    2. Push the first interval on the stack
    3. Iterate through intervals and for each one compare current interval
       with the top of the stack and:
       A. If current interval does not overlap, push on to stack
       B. If current interval does overlap, merge both intervals in to one
          and push on to stack
    4. At the end return stack
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged