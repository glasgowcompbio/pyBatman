# pyBatman
A pipeline to infer metabolite concentrations from 1D NMR data using Bayesian AuTomated Metabolite Analyser for NMR data (BATMAN)

This pipeline can be used to infer metabolite concentrations from 1D NMR spectra. It loads data in Bruker format using [nmrglue](https://www.nmrglue.com/), performs baseline and background correction, and uses BATMAN [1] for inference. A multiplet database optimised for human serum is also provided (my_db.p) in this repo.

[A docker container is provided](https://github.com/joewandy/docker-batman), which includes all the dependencies.

[1] Hao, Jie, et al. "BATMAN—an R package for the automated quantification of metabolites from nuclear magnetic resonance spectra using a Bayesian model." Bioinformatics 28.15 (2012): 2088-2090


