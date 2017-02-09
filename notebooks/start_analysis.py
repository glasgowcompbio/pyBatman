# -*- coding: utf-8 -*-

import os
import sys
import glob
import re
import nbconvert

import pandas as pd

sys.path.append('/home/rstudio/codes')
from pyBatman.nbrun import run_notebook
from pyBatman import load_config, get_db_path, sub_dir_path, mkdir_p
from pyBatman.constants import WD, SPECTRA_DIR, BACKGROUND_DIR, OUTPUT_DIR

def process_spectra(nb_name_input):

    if os.path.isdir(WD): # if wd exists (has been mapped)

        spectra_found = os.path.isdir(SPECTRA_DIR)
        background_found = os.path.isdir(BACKGROUND_DIR)
        if spectra_found and background_found:
            input_spectra = sub_dir_path(SPECTRA_DIR)
            input_spectra = natural_sort(input_spectra)
            mkdir_p(OUTPUT_DIR)

            for i in range(len(input_spectra)):

                sd = input_spectra[i]
                base_name = os.path.basename(os.path.normpath(sd))
                print 'Now processing %d/%d: %s' % (i+1, len(input_spectra), base_name)

                nb_name_output = os.path.join(OUTPUT_DIR, '%s.html' % base_name)
                timeout = -1 # never times out
                nb_kwargs = {'spectra_dir': sd}
                try:
                    run_notebook(nb_name_input, nb_name_output, timeout=timeout, nb_kwargs=nb_kwargs)
                except nbconvert.preprocessors.execute.CellExecutionError, e:
                    print 'Failed to process this spectra: %s' % str(e)

def get_file_id(csv_filename):
    basename = os.path.basename(os.path.normpath(csv_filename))
    filename_no_extension = os.path.splitext(basename)[0]
    tokens = filename_no_extension.split('.')
    try: # Only for curve studies, extract e.g. 2001 from 'Curve_study.ID_2001.fP-EDTA.130510'
        file_id = tokens[1]
        file_id = file_id.replace('ID_', '')
    except IndexError:
        file_id = tokens[0]
    return file_id

def get_df(csv_filename):
    df = pd.read_csv(csv_filename, header=0, index_col=0)
    file_id = get_file_id(csv_filename)
    df.rename(columns={'Concentration (mM)': file_id}, inplace=True)
    return df

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def process_output(csv_files):
    dfs = []
    csv_files = natural_sort(csv_files)
    for csv_filename in csv_files:
        df = get_df(csv_filename)
        dfs.append(df)
    combined_df = pd.concat(dfs, axis=1)
    combined_df = combined_df.transpose()
    combined_df = combined_df.sort_index()
    combined_df.to_csv(os.path.join(OUTPUT_DIR, 'results.csv'))

def main():

    # process all the spectra
    nb_name_input = '/home/rstudio/codes/notebooks/predict_spectra_concentrations.ipynb'
    process_spectra(nb_name_input)

    # combined all the individual results into a single csv file
    csv_dir = os.path.join(OUTPUT_DIR, 'csv')
    csv_files = glob.glob(os.path.join(csv_dir, '*.csv'))
    process_output(csv_files)

if __name__ == "__main__": main()