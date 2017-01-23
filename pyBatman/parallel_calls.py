from random import randint
from time import sleep

# needed to fix 'undefined symbol: PC' when importing robjects on mac
import readline

def par_run_bm(params):
    i, bm, options = params
    print 'START MCMC iteration %d' % i
    bm_r = bm.run(options)
    print 'FINISH MCMC iteration %d' % i
    return bm_r