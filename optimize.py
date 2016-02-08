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
import math

logger = logging.getLogger('optimization')
fdmexe = './dfdm/build/dfdm.exe'
rootdir = os.getcwd()

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
        print('',file=f)

        print('[forcing]',file=f)
        print('dt = 21600 ; 6 hourly data ',file=f)
        print('f_acc = '+rootdir+'/forcing/acc.core03.6H.nc',file=f)
        print('f_wind10m = '+rootdir+'/forcing/w10m.core03.6H.nc',file=f)
        print('f_tskin = '+rootdir+'/forcing/tskin.core03.6H.nc',file=f)
        f.close()

    def toString(self):
        print('eta0 = '+str(self.eta0))
        print('c5 = '+str(self.c5))
        print('c6 = '+str(self.c6))


kEval = 0 # number of cost function evalutions

class CacheLast(object):
    """ Decorator function to workout bug that causes
        twice the same function call 
            https://github.com/scipy/scipy/issues/4076
    """
    def __init__(self, f):
        self.f = f
        self.last_x = None
        self.last_f = None
        self.ncall = 0

    def __call__(self, x):
        if np.all(x == self.last_x):
            return self.last_f
        else:
            self.last_x = x
            self.last_f = self.f(x)
            self.ncall += 1
            return self.last_f


def cost_function(x):
    """ Single valued cost function as called by scipy """
    global kEval
    print(repr(x))
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
    z550_target = 7.00
    z830_target = 55.
    f = open('z550.txt', 'r')
    z550 = float(f.readline())
    f.close()
    f = open('z830.txt', 'r')
    z830 = float(f.readline())
    f.close()

    if (z550 < 0 or z830 < 0):
        # simulation reached maxYear without attaining target density
        cost = 1e9
    else:
        cost = abs(z550-z550_target) + abs(z830-z830_target)

    os.chdir(rootdir)

    kEval += 1
    print('kEval = '+str(kEval)+', cost = '+str(cost))
    return cost

def optimize(): 
    initialize_logging()

    x0 = np.asarray([9e5, 0.08,0.023]) # initial guess
    func = CacheLast(cost_function) 
    res = scipy.optimize.minimize(func, x0, bounds=[(8e5,10e5), (0.02, 0.16), (0.02,0.025)], method='L-BFGS-B')
    #res = scipy.optimize.minimize(func, x0, method='BFGS')

    print(res)


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

