#!/usr/bin/env python
"""
---------------------------------------------------------------------------------------------------------------------
    DESCRIPTION
        Optimization of overburden pressure formulation

        In this program we seek to optimize the overburden pressure formulation of Eric Anderson 1976:
        'A Point Energy and Mass Balance Model of a Snow Cover', formula 3.29.
        Starting with parameter values given by the paper as initial guess, can we find a different set of 
        parameters that fits better to Antarctic firn density profiles?

        We pose this question as a non-linear optimization problem with a single-valued cost function (the misfit) 
        with multiple degrees of freedom (the parameters from the equation). A single computation of the cost function 
        requires us to compute density profiles of Antarctic point locations for which we have measurements. For these
        computations we use a simple dry firn densification model (DFDM) with meteorological forcing from RACMO2.3. 
        The cost function can be defined by the user, but generally it would involve the shape of the density profile 
        or certain characteristic points such as the z550 and z830 depth. 

        An optimization algorithm from the scipy.optimize toolbox is used: BFGS. 

    DATE
        4-FEB-2016
   
    AUTHOR
        L.vankampenhout@uu.nl
---------------------------------------------------------------------------------------------------------------------
"""

from __future__ import print_function
import sys, os, logging, datetime
from subprocess import check_call
import scipy.optimize
import numpy as np

logger = logging.getLogger('optimization')
fdmexe = './dfdm/build/dfdm.exe'

class FdmSettings:
    heat = "true"
    compaction = "overburden"
    
    def __init__(self, x):
        self.eta0 = x[0]
        self.c5 = x[1]
        self.c6 = x[2]

    def writeIniFileWithParameters(self, inifile):
        f = open(inifile, 'w')
        print('[general]',file=f)
        print('heat = '+self.heat,file=f)
        print('compaction = '+self.compaction,file=f)
        print('',file=f)
        print('[overburden]',file=f)
        print('eta0 = '+str(self.eta0),file=f)
        print('c5 = '+str(self.c5),file=f)
        print('c6 = '+str(self.c6),file=f)
        f.close()

kEval = 0 # number of cost function evalutions

def cost_function(x):
    """ Single valued cost function as called by scipy """
    global kEval
    rootdir = os.getcwd()
    state = FdmSettings(x)
    
    rundir = 'run'+str(kEval).zfill(3)
    if os.path.isdir(rundir):
        logger.error('Directory '+rundir+' already exists') 
        sys.exit(1)
    else:
        args = ['mkdir', rundir]
        check_call(args)

    os.chdir(rundir)
    inifile = 'settings.ini'
    state.writeIniFileWithParameters(inifile)

    # run DFDM program in subdirectory
    args = [rootdir+'/'+fdmexe]
    check_call(args)

    # compute cost function

    os.chdir(rootdir)
    kEval += 1
    return abs(x[1])


def optimize(): 
    initialize_logging()

    x0 = np.asarray([9e5, 0.08,0.023]) # initial guess

    cost_function(x0)
    cost_function(x0)
    #res = scipy.optimize.minimize(cost_function, x0, method='BFGS', tol=1e-3)
    #print(res)
    print('kEval = '+str(kEval))


def initialize_logging():
    """
    Initialize logging. Prints to both the console and a log file, at configurable levels
    """
    now = datetime.datetime.now().strftime("%Y%m%d-%H%M")

    #set root logger to debug level        
    logging.getLogger().setLevel(logging.DEBUG)

    #formatter = logging.Formatter('%(asctime)s %(name)s:%(process)d %(levelname)s %(message)s')
    formatter = logging.Formatter('%(levelname)s %(message)s')

    # create logfile specific for coupler
    #log_filename = 'optimize.log.'+now
    #fh = logging.FileHandler(log_filename)
    #fh.setLevel(logging.DEBUG)
    #fh.setFormatter(formatter)
    #logger.addHandler(fh)
    #logger.info('Logging output to %s', log_filename)
    
    # create console handler and attach to root logger
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)
    
 

if __name__ == "__main__":
    optimize()

