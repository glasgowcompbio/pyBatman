#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import datetime
import errno
import tempfile
import shutil
import copy
import cPickle
import gzip
from collections import OrderedDict
from distutils.dir_util import copy_tree
from random import randint

import multiprocessing
import json
import numpy as np
import nmrglue as ng
import pandas as pd
import scipy.interpolate as interpolate

import matplotlib.pyplot as plt
import pylab as plt
import seaborn
from IPython.display import display

from .parallel_calls import par_run_bm
from .models import Spectra, Database
from .helper import mkdir_p, load_config, sub_dir_path

def start_analysis(config_filename):

    config = load_config(config_filename)

    background_dir = config['background_dir']
    pattern = config['pattern']
    working_dir = config['working_dir']
    db = config['database']
    spectra_dir = config['spectra_dir']
    n_burnin = config['n_burnin']
    n_sample = config['n_sample']
    n_iter = config['n_iter']
    tsp_concentration = config['tsp_concentration']
    output_dir = config['output_dir']
    verbose = True if config['verbose'] == 'True' else False

    input_spectra = sub_dir_path(spectra_dir)
    input_backgrounds = sub_dir_path(background_dir)

    pipeline = PyBatmanPipeline(input_backgrounds, pattern, working_dir, db)
    for sd in input_spectra:
        df, fit_results = pipeline.predict_conc(sd, n_burnin, n_sample, n_iter,
                                                tsp_concentration, verbose=verbose)
        pipeline.save_results(sd, df, fit_results, output_dir)

class PyBatmanPipeline(object):

    def __init__(self, background_dir, pattern, working_dir, database, make_plot=True):

        self.background_dir = background_dir
        self.pattern = pattern
        self.working_dir = working_dir
        self.db = Database(from_file=database)
        self.make_plot = make_plot
        self.spiked_bm = {}

        for m in self.db.metabolites:
            if m == 'TSP':
                tsp_multiplet = self.db.metabolites[m].multiplets[0]
                self.tsp_rel_intensity = tsp_multiplet.rel_intensity
                self.tsp_range = tsp_multiplet.ppm_range

    def load_spiked(self, metabolite_name, concentration, spectra_dir, verbose=False):

        print 'Loading spectra for %s (%d)' % (metabolite_name, concentration)
        key = (metabolite_name, concentration)
        bm = PyBatman([spectra_dir], self.background_dir, self.pattern,
            self.working_dir, self.db, verbose=verbose)
        self.spiked_bm[key] = bm

        # get default parameters
        metabolites = [metabolite_name]
        default = bm.get_default_params(metabolites)

        # perform background correction using default parameters
        exclude_regions = [self.tsp_range]
        bm.background_correct(default, make_plot=self.make_plot, exclude=exclude_regions)

        # perform baseline correction using default parameters
        # bm.baseline_correct(default)

    def update_rel_intensities(self, db, metabolite_name, concentrations, tsp_concentration):

        print 'Updating relative intensities for %s' % metabolite_name

        copy_db = copy.deepcopy(db)
        m = copy_db.find(metabolite_name)
        for u in m.multiplets:

            rel_intensities = []
            for concentration in concentrations: # for all std concentrations

                # get the spectral data
                key = (metabolite_name, concentration)
                bm = self.spiked_bm[key]
                spectra_data = bm.spectra_data

                # compute the corrected intensities from peak area ratio to TSP
                corrected = self._get_corrected_intensity(spectra_data, u.ppm_range,
                    tsp_concentration, concentration, self.tsp_rel_intensity)
                rel_intensities.append(corrected)

            corrected = np.median(rel_intensities)
            print 'Rel. intensities: initial = %f, corrected = %f\n' % (u.rel_intensity, corrected)

            # update the relative intensity of the multiplet in copy_db
            u.rel_intensity = corrected

            # update in the dataframe representation as well
            copy_db.df.rel_intensity[copy_db.df.name == metabolite_name] = corrected

        return copy_db

    def fit_single_metabolite(self, spectra_dir, metabolite_name, n_burnin, n_sample, n_iter,
                              verbose=False, n_plots=3, correct_background=True):

        bm = PyBatman([spectra_dir], self.background_dir, self.pattern,
            self.working_dir, self.db, verbose=verbose)

        # get default parameters
        metabolites = [metabolite_name]
        default = bm.get_default_params(metabolites)
        options = default.set('nItBurnin', n_burnin).set('nItPostBurnin', n_sample)

        if correct_background:
            bm.background_correct(default, make_plot=self.make_plot)

        print 'Fitting %s' % metabolite_name
        meta_fits = self._iterate(bm, options, n_iter, n_plots=n_plots)
        # mean_rmse = np.array([out.rmse() for out in meta_fits]).mean()
        return meta_fits

    def predict_conc(self, spectra_dir, n_burnin, n_sample, n_iter, tsp_concentration, verbose=False):

        multiplets = self.db.get_names()
        if verbose:
            print multiplets

        fit_results = {}
        for name in multiplets:
            print
            print '================================================================='
            print 'Now fitting %s for %s' % (name, spectra_dir)
            print '================================================================='
            print
            if name == 'TSP':
                correct_background = False
            else:
                correct_background = True
            fit_results[name] = self.fit_single_metabolite(spectra_dir, name, n_burnin, n_sample, n_iter,
                                                           correct_background=correct_background)

        # predict the concentrations of metabolites
        print
        print '================================================================='
        print 'Results'
        print '================================================================='
        print
        df = self.get_results(multiplets, fit_results, tsp_concentration)
        return df, fit_results

    def save_results(self, sd, df, fit_results, output_dir):
        base_name = os.path.basename(os.path.normpath(sd))
        # timestr = time.strftime("%Y%m%d_%H%M%S")
        output_csv = os.path.join(output_dir, '%s.csv' % base_name)
        output_model = os.path.join(output_dir, '%s.p' % base_name)

        # save the dataframe to csv and also the model as pickled objects
        display(df)
        print 'Saving CSV output to %s' % output_csv
        df.to_csv(output_csv, index=False)

        with gzip.GzipFile(output_model, 'wb') as f:
            print 'Saving model to %s' % output_model
            cPickle.dump(fit_results, f, protocol=cPickle.HIGHEST_PROTOCOL)

    def get_betas(self, bm_out_list):
        all_betas = []
        for bm_out in bm_out_list:
            df = bm_out.beta_df
            names = df.index.values
            first_col = df.ix[:,0]
            betas = first_col.values
            all_betas.append(betas)
        return names, np.array(all_betas)

    def get_results(self, multiplets, fit_results, tsp_concentration):
        betas = []
        for name in multiplets:
            if name == 'TSP':
                continue
            fit_names, fit_betas = self.get_betas(fit_results[name])
            b = fit_betas.flatten()
            betas.append(b)

        betas = np.array(betas)
        betas = betas.transpose()

        plt.boxplot(betas)
        ticks = np.arange(len(multiplets)) + 1
        plt.xticks(ticks, multiplets, rotation='vertical')
        plt.title('Betas')
        plt.show()

        metabolite_names = multiplets
        metabolite_betas = betas
        fit_names, fit_betas = self.get_betas(fit_results['TSP'])
        tsp_betas = fit_betas.flatten()

        beta_m = np.median(metabolite_betas, axis=0)
        beta_tsp = np.median(tsp_betas, axis=0)
        predicted = beta_m / beta_tsp * tsp_concentration
        rows = zip(metabolite_names, predicted)
        df = pd.DataFrame(rows, columns=['Metabolite', 'Concentration (μM)'])
        return df

    def _get_corrected_intensity(self, spectra_data, integrate_range, tsp_conc, metabo_conc, tsp_protons):

        ppm = spectra_data[:, 0]
        intensity = spectra_data[:, 1]

        # integrate over the peak area for the metabolite
        idx = (ppm > integrate_range[0]) & (ppm < integrate_range[1])
        selected_intensity = intensity[idx]
        metabo_area = np.trapz(selected_intensity)

        # integrate the TSP area
        lower = -0.05
        upper = 0.05
        idx = (ppm > lower) & (ppm < upper)
        selected_intensity = intensity[idx]
        tsp_area = np.trapz(selected_intensity)

        metabo_conc = float(metabo_conc)

        # compute the corrected relative intensity
        correction = tsp_protons * (metabo_area/tsp_area) * (tsp_conc/metabo_conc)
        print 'tsp_area=%f, metabo_area=%f, correction=%f' % (tsp_area, metabo_area, correction)

        return correction

    def _iterate(self, bm, options, n_iter, n_plots=3):

        if n_iter > 1: # map each mcmc chain to a parallel worker in the pool

            # prepare the parameters for each chain
            params = []
            for i in range(n_iter):
                param = (i, bm, options)
                params.append(param)

            try:
                num_cores = multiprocessing.cpu_count()
                print 'n_iter=%d num_cores=%d' % (n_iter, num_cores)
                pool = multiprocessing.Pool(num_cores)
                iter_results = pool.map(par_run_bm, params)
            finally:
                pool.close()
                pool.join()

        else: # just a single chain
            print 'START MCMC iteration %d' % 0
            bm_r = bm.run(options)
            print 'FINISH MCMC iteration %d' % 0
            iter_results = [bm_r]

        # plot the results if necessary
        bm_outputs = []
        assert len(iter_results) == n_iter
        for i in range(n_iter):
            bmo = PyBatmanOutput(iter_results[i], options)
            bm_outputs.append(bmo)
            if self.make_plot and i < n_plots:
                bmo.plot_fit()

        return bm_outputs

class PyBatman(object):

    def __init__(self, input_dirs, background_dirs, pattern, working_dir, db, verbose=False):

        self.db = db
        self.input_dirs = input_dirs
        self.background_dirs = background_dirs
        self.pattern = pattern
        self.working_dir = working_dir
        self.verbose = verbose

        # see http://stackoverflow.com/questions/3925096/how-to-get-only-the-last-part-of-a-path-in-python
        self.labels = []
        for input_dir in self.input_dirs:
            last = os.path.basename(os.path.normpath(input_dir))
            self.labels.append(last)

        # extract spectra and background
        self.spectra = self._load_data(self.input_dirs, self.pattern)
        self.background = self._load_data(self.background_dirs, self.pattern)

        # combine spectra to have the same ppm scales
        self.spectra_data, self.mean_bg = self._combine_spectra(self.spectra, self.background)

    def plot_spectra(self, name):

        meta_ranges = self.db.metabolites[name].ppm_range()

        for lower, upper in meta_ranges:

            ppm = self.spectra_data[:, 0]
            idx = (ppm > lower) & (ppm < upper)
            n_row, n_col = self.spectra_data.shape

            plt.figure()
            x = ppm[idx]
            for i in range(1, n_col):
                intensity = self.spectra_data[:, i]
                y = intensity[idx]
                plt.plot(x, y)
            plt.gca().invert_xaxis()
            plt.xlabel('ppm')
            plt.ylabel('intensity')
            plt.title('%s at (%.4f-%.4f)' % (name, lower, upper))
            plt.show()

    def get_default_params(self, names):
        # check if it's TSP
        is_tsp = False
        if 'TSP' in names:
            if len(names) > 1:
                raise ValueError('Not recommended to fit TSP with other metabolites')
            is_tsp = True

        # select only the metabolites we need
        selected = self._select_metabolites(names)
        options = PyBatmanOptions(selected)

        # and their ranges
        ranges = []
        for m in selected:
            if m.ppm_range() > 0:
                ranges.extend(m.ppm_range()) # a list of tuples
        ranges = ['(%s, %s)' % r for r in ranges] # convert to a list of string
        ranges_str = ' '.join(ranges)

        # set no. of spectra and no. of processors to use
        spec_no = len(self.spectra)
        para_proc = multiprocessing.cpu_count()

        # set common parameters
        options = options.set('ppmRange', ranges_str)             \
                             .set('specNo', '1-%d' % spec_no)     \
                             .set('paraProc', para_proc)          \
                             .set('nItBurnin', 19000)             \
                             .set('nItPostBurnin', 1000)          \
                             .set('thinning', 5)                  \
                             .set('tauMean', -0.01)               \
                             .set('tauPrec', 2)                   \
                             .set('rdelta', 0.005)                \
                             .set('csFlag', 0)

        # special parameters for TSP since it's so different from the rest
        if is_tsp:
            default = options.set('muMean', 1.5)                  \
                                 .set('muVar', 0.1)               \
                                 .set('muVar_prop', 0.01)         \
                                 .set('nuMVar', 0)                \
                                 .set('nuMVarProp', 0.1)
        else:
            default = options.set('muMean', 0)                     \
                                 .set('muVar', 0.01)               \
                                 .set('muVar_prop', 0.0002)        \
                                 .set('nuMVar', 0.0025)            \
                                 .set('nuMVarProp', 0.01)


        return default

    def run(self, options, plot=True):

        # create input dir for batman
        temp_dir, batman_input, batman_output = self._create_batman_dirs(self.working_dir)
        spectra_file = os.path.join(batman_input, 'NMRdata_temp.txt')
        self._write_spectra_data(spectra_file, self.spectra_data, self.labels)
        self._write_parameters(options, batman_input)

        # run batman using rpy2
        import rpy2.robjects as robjects
        sink_r = robjects.r['sink']
        if self.verbose:
            sink_r()
        else:
            sink_r('/dev/null')

        from rpy2.robjects.packages import importr
        importr('batman')
        batman_r = robjects.r['batman']
        bm = batman_r(runBATMANDir=temp_dir, txtFile=spectra_file, figBatmanFit=False)

        # TODO: make traceplots etc
        plot_batman_fit_r = robjects.r['plotBatmanFit']
        plot_batman_fit_stack_r = robjects.r['plotBatmanFitStack']
        plot_rel_con_r = robjects.r['plotRelCon']
        plot_meta_fit_r = robjects.r['plotMetaFit']
        plot_batman_fit_r(bm, showPlot=False)
        plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
        plot_rel_con_r(bm, showPlot=False)
        plot_meta_fit_r(bm, showPlot=False)

        if not self.verbose:
            sink_r()

        # cleaning up ..
        shutil.rmtree(temp_dir)
        if self.verbose:
            print 'Deleted', temp_dir

        return bm

    def baseline_correct(self, options):

        ppm = self.spectra_data[:, 0]
        n_row, n_col = self.spectra_data.shape

        for i in range(1, n_col):

            intensity = self.spectra_data[:, i]
            selected = options.selected
            label = self.labels[i-1]

            for metabolite in selected:
                ranges = metabolite.ppm_range()
                if len(ranges) > 0:
                    for lower, upper in ranges:
                        idx = (ppm > lower) & (ppm < upper)
                        selected_intensity = intensity[idx]
                        selected_ppm = ppm[idx]

                        # divide the array into halves and find the mininum indices from each portion
                        left, right = np.array_split(selected_intensity, 2)
                        to_find = [left.min(), right.min()]
                        to_find_idx = np.in1d(selected_intensity, to_find)
                        min_pos = np.where(to_find_idx)[0]

                        # node_list will be the indices of ...
                        first = 0
                        last = len(selected_intensity)-1
                        # node_list = [first, last]
                        node_list = [first]         # the first element
                        node_list.extend(min_pos)   # the min elements from left and right
                        node_list.append(last)      # the last element
                        corrected_intensity = ng.proc_bl.base(selected_intensity, node_list) # piece-wise baseline correction

                        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
                        title = '(%.4f-%.4f)' % (lower, upper)
                        ax1.plot(selected_ppm, selected_intensity, 'b', label='Spectra')
                        ax1.plot(selected_ppm, corrected_intensity, 'g', label='Baseline Corrected')
                        ax1.set_xlabel('ppm')
                        ax1.set_ylabel('intensity')
                        ax1.legend(loc='best')
                        ax1.invert_xaxis()
                        ax1.set_title(title)

                        intensity[idx] = corrected_intensity
                        idx = (ppm > lower-0.05) & (ppm < upper+0.05)
                        title = 'Corrected (%.4f-%.4f)' % (lower-0.05, upper+0.05)
                        ax2.plot(ppm[idx], intensity[idx], 'b')
                        ax2.set_xlabel('ppm')
                        ax2.invert_xaxis()
                        ax2.set_title(title)

                        plt.suptitle('%s %s' % (label, metabolite.name), fontsize=16, y=1.08)
                        plt.tight_layout()
                        plt.show()
                        assert np.array_equal(self.spectra_data[:, i], intensity)

    def background_correct(self, options, make_plot=True, exclude=[]):

        ppm = self.spectra_data[:, 0]
        background = self.mean_bg
        n_row, n_col = self.spectra_data.shape

        for i in range(1, n_col):

            intensity = np.copy(self.spectra_data[:, i])
            corrected_intensity = intensity - background
            self.spectra_data[:, i] = corrected_intensity

            selected = options.selected
            label = self.labels[i-1]

            for metabolite in selected:
                ranges = metabolite.ppm_range()
                if len(ranges) > 0:
                    for lower, upper in ranges:
                        idx = (ppm > lower) & (ppm < upper)
                        if make_plot:
                            f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
                            title = '(%.4f-%.4f)' % (lower, upper)
                            ax1.plot(ppm[idx], intensity[idx], 'b', label='Spectra')
                            ax1.plot(ppm[idx], background[idx], 'k--', label='Background')
                            ax1.plot(ppm[idx], corrected_intensity[idx], 'g', label='Spectra - Background')
                            ax1.set_xlabel('ppm')
                            ax1.set_ylabel('intensity')
                            ax1.legend(loc='best')
                            ax1.invert_xaxis()
                            ax1.set_title(title)

                            idx = (ppm > lower-0.05) & (ppm < upper+0.05)
                            title = 'Corrected (%.4f-%.4f)' % (lower-0.05, upper+0.05)
                            ax2.plot(ppm[idx], corrected_intensity[idx], 'b')
                            ax2.set_xlabel('ppm')
                            ax2.invert_xaxis()
                            ax2.set_title(title)

                            plt.suptitle('%s %s' % (label, metabolite.name), fontsize=16, y=1.08)
                            plt.tight_layout()
                            plt.show()

            if len(exclude)> 0:
                for lower, upper in exclude:
                    idx = (ppm > lower) & (ppm < upper)
                    self.spectra_data[idx, i] = intensity[idx]

    def _select_metabolites(self, names):
        selected = []
        for name in names:
            selected.append(self.db.find(name))
        return selected

    def _create_batman_dirs(self, working_dir):

        # create the main working directory
        mkdir_p(working_dir)
        prefix = datetime.datetime.now().strftime("%G%m%d_%H%M%S_")
        temp_dir = tempfile.mkdtemp(prefix=prefix, dir=working_dir)
        batman_dir = os.path.join(temp_dir, 'runBATMAN')

        # create batman input and output dirs
        batman_input = os.path.join(batman_dir, 'BatmanInput')
        batman_output = os.path.join(batman_dir, 'BatmanOutput')
        mkdir_p(batman_input)
        mkdir_p(batman_output)
        if self.verbose:
            print 'Working directory =', working_dir
            print '- batman_input =', batman_input
            print '- batman_output =', batman_output

        return temp_dir, batman_input, batman_output

    def _write_spectra_data(self, spectra_file, spectra_data, labels):
        columns = ['ppm']
        columns.extend(labels)
        df = pd.DataFrame(spectra_data, columns=columns)
        df.to_csv(spectra_file, index=False, sep='\t')
        if self.verbose:
            print 'Spectra written to', spectra_file

    def _write_parameters(self, options, batman_input, parallel=False, seed=None):

        selected = options.selected
        metabolites_list_df = self._write_metabolites_list(selected, batman_input, 'metabolitesList.csv')
        multi_data_user_df = self._write_multiplet_data(selected, batman_input, 'multi_data_user.csv')
        chem_shift_per_spec_df = self._write_chem_shift(selected, batman_input, 'chemShiftPerSpec.csv')
        if self.verbose:
            display(metabolites_list_df)
            display(multi_data_user_df)
            display(chem_shift_per_spec_df)

        if not parallel:
            options = options.set('paraProc', 1)

        if seed is not None:
            random_seed = seed
        else:
            random_seed = randint(0, 1e6)
        options = options.set('randSeed', random_seed)

        options_lines = self._write_batman_options(options, batman_input, 'batmanOptions.txt')
        if self.verbose:
            print
            print '----------------------------------------------------------------'
            print 'Parameters'
            print '----------------------------------------------------------------'
            for line in options_lines:
                print line
            print

    def _write_metabolites_list(self, metabolites, batman_input, file_name):
        columns = ['Metabolite']
        data = [m.name for m in metabolites]
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(batman_input, file_name)
        if self.verbose:
            print 'metabolites list = %s' % out_file
        df.to_csv(out_file, header=False, index=False, sep=',')

        return df

    # create multi_data_user.csv
    def _write_multiplet_data(self, metabolites, batman_input, file_name):
        columns = ['Metabolite', 'pos_in_ppm', 'couple_code', 'J_constant',
                   'relative_intensity', 'overwrite_pos', 'overwrite_truncation',
                   'Include_multiplet']
        data = []
        OVERWRITE_POS = 'n'
        OVERWRITE_TRUNCATION = 'n'
        INCLUDE_MULTIPLET = 1
        for m in metabolites:
            for u in m.multiplets:
                row = (u.parent.name, u.ppm, u.couple_code, u.j_constant, u.rel_intensity,
                        OVERWRITE_POS, OVERWRITE_TRUNCATION, INCLUDE_MULTIPLET)
                data.append(row)
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(batman_input, file_name)
        if self.verbose:
            print 'multiplet data = %s' % out_file
        df.to_csv(out_file, index=False, sep=',')

        return df

    # create chemShiftPerSpec.csv
    def _write_chem_shift(self, metabolites, batman_input, file_name):
        columns = ['multiplets', 'pos_in_ppm']
        columns.extend(self.labels)
        data = []
        ns = ['n'] * len(self.labels)
        for m in metabolites:
            for u in m.multiplets:
                row = [u.parent.name, u.ppm]
                row.extend(ns)
                data.append(row)
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(batman_input, file_name)
        if self.verbose:
            print 'chem shift = %s' % out_file
        df.to_csv(out_file, index=False, sep=',')

        return df

    # create batmanOptions.txt
    def _write_batman_options(self, options, batman_input, file_name):
        out_file = os.path.join(batman_input, file_name)
        if self.verbose:
            print 'batman options = %s' % out_file

        option_lines = []
        with open(out_file, 'w') as f:
            params = options.params
            params_desc = options.params_desc
            for key in params:
                value = str(params[key])
                description = params_desc[key]
                line = key + ' - ' + description + ': ' + value
                option_lines.append(line)
                f.write('%s\n' % line)

        return option_lines

    def _get_matching_paths(self, input_dirs, pattern):

        matching = []
        for input_dir in input_dirs:

            # find all the child directories
            sub_dirs = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]

            for path in sub_dirs:

                # check the pulse program if it exists
                if os.path.isdir(path):
                    pp = os.path.join(path, 'pulseprogram')
                    if os.path.isfile(pp): # if exists

                        # if it contains the pattern then store this path
                        with open(pp, 'r') as f:
                            head = [next(f) for x in xrange(2)]
                            if pattern in head[1]:
                                matching.append(path)

        # special check for a Bruker folder contains multiple PROC that match the pattern
        if len(input_dirs) == 1 and len(matching) > 1:
            print 'WARNING: multiple matches'
            for match in matching:
                print '- %s' % match
            print 'found for the pattern "%s" in %s' % (pattern, input_dirs[0])
            print 'The first match %s is used' % (matching[0])
            matching = [matching[0]]

        return matching

    def _mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _prepare_data(self, matching, data_dir):
        out_dirs = []
        for f in range(len(matching)):
            path = matching[f]
            out_dir = os.path.join(data_dir, 'pdata', str(f))
            mkdir_p(out_dir)
            if self.verbose:
                print 'Copied spectra to', out_dir
            copy_tree(path, out_dir)
            out_dirs.append(out_dir)
        return out_dirs

    def _load_single_spectra(self, spectra_dir):
        p_data = os.path.join(spectra_dir, 'pdata/1')
        if self.verbose:
            print 'Processing', p_data

        dic, data = ng.bruker.read_pdata(p_data)
        udic = ng.bruker.guess_udic(dic, data)
        uc = ng.fileiobase.uc_from_udic(udic, 0)

        x = []
        y = []
        for ppm in uc.ppm_scale():
            x.append(ppm)
            y.append(data[uc(ppm, 'ppm')])
        x = np.array(x)
        y = np.array(y)

        return Spectra(x, y)

    def _load_data(self, input_dirs, pattern):
        matching = self._get_matching_paths(input_dirs, pattern)
        spectra_list = [self._load_single_spectra(input_dir) for input_dir in matching]
        return spectra_list

    def _resample(self, xs, ys, ppm):
        new_ys = []
        for x, y in zip(xs, ys):
            new_y = interpolate.interp1d(x, y)(ppm)
            new_ys.append(new_y)
        return new_ys

    def _get_ppms(self, spectra_list):
        xs = [spectra.ppm for spectra in spectra_list]
        minx = min([x[0] for x in xs])
        maxx = max([x[-1] for x in xs])
        ppm = np.linspace(minx, maxx, num=len(xs[0]))
        return ppm

    def _combine_spectra(self, spectra_list, background_list):

        ppm = self._get_ppms(spectra_list + background_list)

        # process the sample spectra first
        xs = [spectra.ppm for spectra in spectra_list]
        ys = [spectra.intensity for spectra in spectra_list]
        new_ys = self._resample(xs, ys, ppm)

        combined = [ppm]
        combined.extend(new_ys)
        combined = np.array(combined).transpose()
        if self.verbose:
            print 'Loaded', combined.shape

        # process the background too if available
        if len(background_list) > 0:

            xs = [spectra.ppm for spectra in background_list]
            ys = [spectra.intensity for spectra in background_list]
            new_ys = self._resample(xs, ys, ppm)

            # get the average of all the background spectra
            all_bg = np.array(new_ys)
            mean_bg = all_bg.mean(axis=0)

        else:
            mean_bg = None

        return combined, mean_bg

class PyBatmanOptions(object):

    def __init__(self, selected, params=None):

        self.selected = selected
        if params is None: # default parameters

            # These parameters will be written to the batmanOptions.txt file, so
            # DO NOT CHANGE THE ORDER OF PARAMETERS HERE!!
            # since the parser in batman seems to break easily
            self.params = OrderedDict()
            self.params['ppmRange']      = None         # set inside PyBatman.get_default_params(), e.g. '(-0.05, 0.05)'
            self.params['specNo']        = None         # set inside PyBatman.get_default_params(), e.g. '1-4'
            self.params['paraProc']      = None         # set inside PyBatman.get_default_params(), e.g. 4
            self.params['negThresh']     = -0.5         # probably don't change
            self.params['scaleFac']      = 900000       # probably don't change
            self.params['downSamp']      = 1            # probably don't change
            self.params['hiresFlag']     = 1            # probably don't change
            self.params['randSeed']      = None         # set inside PyBatman.get_default_params(), e.g. 25
            self.params['nItBurnin']     = None         # set inside PyBatman.get_default_params(), e.g. 9000
            self.params['nItPostBurnin'] = None         # set inside PyBatman.get_default_params(), e.g.
            self.params['multFile']      = 2            # probably don't change
            self.params['thinning']      = None         # set inside PyBatman.get_default_params(), e.g. 5
            self.params['cfeFlag']       = 0            # probably don't change
            self.params['nItRerun']      = 5000         # unused
            self.params['startTemp']     = 1000         # probably don't change
            self.params['specFreq']      = 600          # probably don't change
            self.params['a']             = 0.00001      # probably don't change
            self.params['b']             = 0.000000001  # probably don't change
            self.params['muMean']        = None         # set inside PyBatman.get_default_params(), e.g. 0
            self.params['muVar']         = None         # set inside PyBatman.get_default_params(), e.g. 0.1
            self.params['muVar_prop']    = None         # set inside PyBatman.get_default_params(), e.g. 0.002
            self.params['nuMVar']        = None         # set inside PyBatman.get_default_params(), e.g. 0.0025
            self.params['nuMVarProp']    = None         # set inside PyBatman.get_default_params(), e.g. 0.1
            self.params['tauMean']       = None         # set inside PyBatman.get_default_params(), e.g. -0.05
            self.params['tauPrec']       = None         # set inside PyBatman.get_default_params(), e.g. 2
            self.params['rdelta']        = None         # set inside PyBatman.get_default_params(), e.g. 0.002
            self.params['csFlag']        = None         # set inside PyBatman.get_default_params(), e.g. 0

        else:
            self.params = params

        # load parameter descriptions
        current_dirname = os.path.dirname(os.path.realpath(__file__))
        desc_filename = os.path.join(current_dirname, 'params_desc.json')
        with open(desc_filename) as f:
            data = json.load(f)
        self.params_desc = {}
        for key, value in data:
            self.params_desc[key] = value

    # returns a new copy each time
    def set(self, key, val):
        copy = self.params.copy()
        copy[key] = val
        return PyBatmanOptions(self.selected, params=copy)

class PyBatmanOutput(object):

    def __init__(self, bm, options):

        self.output_dir = bm[bm.names.index('outputDir')][0]
        self.bm = bm
        self.options = options

        beta = bm[bm.names.index('beta')]
        # beta_sam = bm[bm.names.index('betaSam')]
        specfit = bm[bm.names.index('sFit')]

        from rpy2.robjects import pandas2ri
        pandas2ri.activate() # to convert between R <-> pandas dataframes

        self.beta_df = pandas2ri.ri2py(beta)
        # self.beta_sam_df = pandas2ri.ri2py(beta_sam).transpose()

        self.specfit_df = pandas2ri.ri2py(specfit)
        self.ppm = self.specfit_df['ppm'].values
        self.original_spectrum = self.specfit_df['Original.spectrum'].values
        self.metabolites_fit = self.specfit_df['Metabolites.fit'].values
        self.wavelet_fit = self.specfit_df['Wavelet.fit'].values
        self.overall_fit = self.specfit_df['Overall.fit'].values

    def plot_fit(self):

        for metabolite in self.options.selected:

            for lower, upper in metabolite.ppm_range():

                idx = (self.ppm > lower) & (self.ppm < upper)
                ppm = self.ppm[idx]
                original_spectrum = self.original_spectrum[idx]
                metabolites_fit = self.metabolites_fit[idx]
                wavelet_fit = self.wavelet_fit[idx]

                plt.figure()
                plt.plot(ppm, original_spectrum, color='blue', label='Original Spectra')
                plt.plot(ppm, metabolites_fit, color='green', label='Metabolites Fit')
                plt.plot(ppm, wavelet_fit, color='red', label='Wavelet Fit')
                # plt.plot(self.ppm, self.overall_fit, '--', color='black', label='Combined Fit')
                plt.legend(loc='best')
                plt.gca().invert_xaxis()
                title = 'Fit Results -- %s (%.4f-%.4f)' % (metabolite.name, lower, upper)
                plt.title('%s' % title)
                plt.show()

    # def plot_beta_sam(self):
    #     self.beta_sam_df.boxplot()

    def rmse(self):
        error = np.sqrt((self.error_vect() ** 2).mean())
        return error

    def error_vect(self):
        return self.original_spectrum - self.metabolites_fit