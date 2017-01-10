from random import randint
from time import sleep

# needed to fix 'undefined symbol: PC' when importing robjects on mac
import readline

# def process_bm(temp_dir, spectra_file, verbose):
#
#     import rpy2.robjects as robjects
#     sink_r = robjects.r['sink']
#     if verbose:
#         sink_r()
#     else:
#         sink_r('/dev/null')
#
#     # actually runs batman here using rpy2
#     from rpy2.robjects.packages import importr
#     importr('batman')
#     batman_r = robjects.r['batman']
#     print temp_dir, spectra_file
#     print type(temp_dir), type(spectra_file)
#     bm = batman_r(runBATMANDir=temp_dir, txtFile=spectra_file, figBatmanFit=False)
#
#     # TODO: make traceplots etc
#     plot_batman_fit_r = robjects.r['plotBatmanFit']
#     plot_batman_fit_stack_r = robjects.r['plotBatmanFitStack']
#     plot_rel_con_r = robjects.r['plotRelCon']
#     plot_meta_fit_r = robjects.r['plotMetaFit']
#     plot_batman_fit_r(bm, showPlot=False)
#     plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
#     plot_rel_con_r(bm, showPlot=False)
#     plot_meta_fit_r(bm, showPlot=False)
#
#     if not verbose:
#         sink_r()
#
#     return bm

def par_run_bm(params):

    i, bm, options = params

    print 'START MCMC iteration %d' % i
    bm_r = bm.run(options)
    print 'FINISH MCMC iteration %d' % i

    return bm_r