#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyBatman import PyBatmanPipeline

db = PyBatmanPipeline.load_db('my_db.p')
background_dir = [
    '/home/rstudio/NMR/Hannah_NMR/HR_standards.ID_Hannah_Dialysed-serum--PBS-1.fS.160928',
    '/home/rstudio/NMR/Hannah_NMR/HR_standards.ID_Hannah_Dialysed-serum--PBS-2.fS.160928',
    '/home/rstudio/NMR/Hannah_NMR/HR_standards.ID_Hannah_Dialysed-serum--PBS-3.fS.160928'
]
pattern = 'cpmg'
working_dir = 'test/temp'

pipeline = PyBatmanPipeline(background_dir, pattern, working_dir, db, make_plot=False)
tsp_concentration = 2320
burn_in = 100
samples = 100
iterations = 4

metabolite = 'Acetate'
spectra_dir = '/home/rstudio/NMR/Hannah_NMR/HR_standards.ID_Hannah_100-_M-acetate.fS.160928'
results = pipeline.predict_conc(spectra_dir, metabolite, burn_in, samples, iterations, tsp_concentration,
                                verbose=False)