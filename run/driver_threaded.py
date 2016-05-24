#!/usr/bin/env python

# Threading based on:
# http://chriskiehl.com/article/parallelism-in-one-line/

import ConfigParser
import os, sys
from subprocess import check_call
from shutil import move
from time import strftime
from multiprocessing.dummy import Pool as ThreadPool 

###### BEGIN USER SETTINGS ######

num_threads = 4
model_exe='dfdm.exe'
config_name = 'settings.ini'

###### END USER SETTINGS ######


def main():
   global model_exe
   model_exe = os.path.abspath(model_exe)

   out_dir = "out."+strftime("%Y%m%d-%H%M%S")
   os.mkdir(out_dir)
   os.chdir(out_dir)

   pool = ThreadPool(num_threads) 
   pool.map(runModel, createRunDirs())
   
def runModel(rundir):
   print(rundir+" starting")
   args=[model_exe, config_name]
   check_call(args,cwd=rundir)

def createRunDirs():
   rundirs = [ ]

   for physics in ['CROCUS','Anderson','CLM45','HL','Pimienta']:
      for case in ['core01', 'core03', 'SipleDome']:
      #for case in ['core03']:
         config = ConfigParser.RawConfigParser()
      
         config.add_section('general')
         config.add_section('physics')
         config.add_section('stopping_criteria')
         config.add_section('fresh_snow_density')
         config.add_section('forcing')
      
         config.set('general', 'dt', '21600')
         
         config.set('stopping_criteria', 'which_stop', '3')
         config.set('stopping_criteria', 'years', '4000')
         config.set('stopping_criteria', 'dens', '850')
         
         #config.set('fresh_snow_density', 'density', '200.')
         
         config.set('forcing', 'which_forcing', '1')
         config.set('forcing', 'netcdf_dt', '21600')
         config.set('forcing', 'netcdf_acc', '/Users/leo/workspace/github/compaction/forcing/'+case+'/acc.6H.nc')
         config.set('forcing', 'netcdf_w10m', '/Users/leo/workspace/github/compaction/forcing/'+case+'/w10m.6H.nc')
         config.set('forcing', 'netcdf_tskin', '/Users/leo/workspace/github/compaction/forcing/'+case+'/tskin.6H.nc')

         config.set('physics', 'which_diffusion', '2')
      
         if (physics == 'CROCUS'):
            # CROCUS
            config.set('physics', 'which_compaction', '5')
            config.set('physics', 'which_metamorphism', '2')
            config.set('fresh_snow_density', 'which_fsd', '3')
         elif (physics == 'Anderson'):
            # Anderson 1979 (no wind density)
            config.set('physics', 'which_compaction', '2')
            config.set('physics', 'which_metamorphism', '1')
            config.set('fresh_snow_density', 'which_fsd', '4')
         elif (physics == 'CLM45'):
            # CLM 4.5 == Anderson 1979 and wind dependence Liston 2007
            config.set('physics', 'which_compaction', '2')
            config.set('physics', 'which_metamorphism', '1')
            config.set('fresh_snow_density', 'which_fsd', '5')
         elif (physics == 'HL'):
            # Herron&Langway with Helsen2008 density
            config.set('physics', 'which_compaction', '1')
            config.set('physics', 'which_metamorphism', '0')
            config.set('fresh_snow_density', 'which_fsd', '1')
         elif (physics == 'Pimienta'):
            # Pimienta with Helsen2008 density
            config.set('physics', 'which_compaction', '3')
            config.set('physics', 'which_metamorphism', '0')
            config.set('fresh_snow_density', 'which_fsd', '1')
         else:
            print('undefined value for physics : '+physics)
            sys.exit(-1)
      
         caseDir = os.path.join(physics,case)
         os.makedirs(caseDir)
         rundirs.append(caseDir)

         # Writing our configuration file
         with open(os.path.join(caseDir,config_name), 'wb') as configfile:
             config.write(configfile)

   return rundirs

if __name__ == "__main__":
   main()
