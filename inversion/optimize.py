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
now = datetime.datetime.now().strftime("%Y%m%d-%H%M")
rootdir = os.getcwd()

class FdmSettings:
    heat = "false"
    compaction = "overburden"
    
    def __init__(self, x):
        #x0 = np.asarray([9e5, 0.08,0.023]) # initial guess
        self.eta0   = (x[0] + 1) * 9e5
        self.c5     = (x[1] + 1) * 0.08
        self.c6     = (x[2] + 1) * 0.023
        self.rho_s  = (x[3] + 1) * 150

        #for i, xi in enumerate(x):
        #    logger.info('x['+str(i)+'] = '+str(xi))

        logger.info('x[0] = ' + str(x[0])+', x[1] = ' + str(x[1])+', x[2] = ' + str(x[2])+ ', x[3] = '+ str(x[3]))
        logger.info('eta0 = ' + str(self.eta0) + ', c5 = '+str(self.c5)+', c6 = '+str(self.c6)+', rho_s = '+str(self.rho_s))

    def writeIniFileWithParameters(self, inifile):
        global max_depth
        global max_year
        global problem
        f = open(inifile, 'w')
        print('[general]',file=f)
        print('heat = '+self.heat,file=f)
        print('compaction = '+self.compaction,file=f)
        print('max_depth = '+str(max_depth),file=f)
        print('max_year = '+str(max_year),file=f)
        print('',file=f)

        print('[density]',file=f)
        print('rho_s = '+str(self.rho_s),file=f)
        print('',file=f)

        print('[overburden]',file=f)
        print('eta0 = '+str(self.eta0),file=f)
        print('c5 = '+str(self.c5),file=f)
        print('c6 = '+str(self.c6),file=f)
        print('',file=f)

        print('[forcing]',file=f)
        print('dt = 21600 ; 6 hourly data ',file=f)
        print('f_acc = '+rootdir+'/forcing/acc.'+problem.tag+'.6H.nc',file=f)
        print('f_wind10m = '+rootdir+'/forcing/w10m.'+problem.tag+'.6H.nc',file=f)
        print('f_tskin = '+rootdir+'/forcing/tskin.'+problem.tag+'.6H.nc',file=f)
        f.close()

    def toString(self):
        logger.info('eta0 = '+str(self.eta0))
        logger.info('c5 = '+str(self.c5))
        logger.info('c6 = '+str(self.c6))


class Problem:
    def __init__(self, tag, z550, z830):
        self.tag = tag
        self.z550 = z550
        self.z830 = z830


kEval = 0 # number of cost function evalutions

def cost_function(x):
    """ Single valued cost function as called by scipy """
    global kEval
    global max_depth, max_year
    global problem

    max_depth = 125. # maximum depth of 1D firn model at which to consider the simulation as failed (prevents excessive runtimes)
    max_year = 3500 # maximum number of years at which to consider the simulation as failed

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

    f = open('z550.txt', 'r')
    z550 = float(f.readline())
    f.close()
    f = open('z830.txt', 'r')
    z830 = float(f.readline())
    f.close()

    if (z830 > 0 and z550 > 0):
        # normal termination
        cost = abs(z550-problem.z550) + abs(z830-problem.z830)
    elif (z830 < 0 and z550 > 0):
        # simulation exited by maxYear or maxDepth without attaining z830
        cost = abs(z550-problem.z550) + 400
    else:
        # simulation reached maxYear without attaining any target density
        cost = 1e4 * (1+sum(abs(x)))

    os.chdir(outdir)

    kEval += 1
    logger.info('kEval = '+str(kEval)+', cost = '+str(cost))
    return cost


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


def optimize(): 
    global problem

    core01 = Problem("core01", z550=25., z830=100.)
    core03 = Problem("core03", z550=7., z830=55.)

    # choose core site here
    problem = core01

    initialize_outdir()
    initialize_logging()

    func = CacheLast(cost_function) 
    """
        List of constrained optimization methods
            from: http://docs.scipy.org/doc/scipy/reference/optimize.html
        L-BFGS-B    tends to jump to extreme values alot, so needs tight bounds.
                    Also, in optimization of eta0 alone, with bounds (-1e-3, 1e-3)
                    and with self.eta0 = 1e9 * x[0] + 9e5 formulation, appears to get 
                    stuck in alternating sequence and never finishes.
        TNC         very slow to come out of initial guess
        SLSQP       stops after 2 iterations,  message: 'Positive directional derivative for linesearch'
    """
    #x0 = np.asarray([0, 0.08,0.023]) # initial guess
    x0 = np.asarray([0,0,0,0]) # initial guess

    #INFO x[0] = 0.135104056169, x[1] = 0.135183621789, x[2] = 0.135920060328, x[3] = -0.0489701430606
    # warm start
    x0 = np.asarray([0.135104056169,0.135183621789,0.135920060328, -0.0489701430606])

    #res = scipy.optimize.minimize(func, x0, bounds=[(-1e-3,1e-3)], method='L-BFGS-B')
    #res = scipy.optimize.minimize(func, x0, bounds=[(8e5,10e5), (0.02, 0.16), (0.02,0.025)], method='L-BFGS-B')
    #res = scipy.optimize.minimize(func, x0, bounds=[(-1e-3,1e-3), (0.02, 0.16), (0.02, 0.026)], method='L-BFGS-B')
    minbound = -1.
    maxbound = 1.
    eps = 1e-3
    res = scipy.optimize.minimize(func, x0, bounds=[(minbound,maxbound), 
                                                    (minbound,maxbound), 
                                                    (minbound,maxbound), 
                                                    (-0.5,0.5)], method='L-BFGS-B', options={'eps':eps})
    logger.info(res)


def initialize_outdir():
    """ 
    Create output directory and CD
    """
    global now, outdir
    outdir = 'out.'+problem.tag+"."+now
    args = ['mkdir', outdir]
    check_call(args)
    os.chdir(outdir)
    outdir = os.getcwd() # full path


def initialize_logging():
    """
    Initialize logging. Prints to both the console and a log file, at configurable levels
    """
    global now

    #set root logger to debug level        
    logging.getLogger().setLevel(logging.DEBUG)

    #formatter = logging.Formatter('%(asctime)s %(name)s:%(process)d %(levelname)s %(message)s')
    formatter = logging.Formatter('%(levelname)s %(message)s')

    # create console handler and attach to root logger
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)

    # create logfile and attach to local logger
    log_filename = 'optimize.log.'+now
    fh = logging.FileHandler(log_filename)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info('Logging output to %s', log_filename)
 

if __name__ == "__main__":
    optimize()

