#!/usr/bin/env python

# Threading based on:
# http://www.boschmans.net/2010/02/03/simple-python-threading-thread-example/

import ConfigParser
import os, sys
from subprocess import check_call
from shutil import move
from time import strftime
from threading import Thread, Lock
from Queue import Queue

###### BEGIN USER SETTINGS ######

num_threads = 4
model_exe='dfdm.exe'
config_name = 'settings.ini'

###### END USER SETTINGS ######

jobs = Queue()
singlelock = Lock()

def main():
   global model_exe
   model_exe = os.path.abspath(model_exe)

   out_dir = "out."+strftime("%Y%m%d-%H%M%S")
   os.mkdir(out_dir)
   os.chdir(out_dir)

   # spawn threads
   for i in range(num_threads):
      t = Thread(target=worker)
      t.deamon = True
      print("Thread "+str(i)+" starting")
      t.start()

   # built queue
   for rundir in createRunDirs():
      print(rundir)
      jobs.put(rundir)
   
def worker():
   HERE = os.getcwd()
   while True:
      try: 
         job = jobs.get(True, 1) # timeout of 1 second
         singlelock.acquire()
         print(job)
         singlelock.release()
         os.chdir(job)
         args=[model_exe, config_name]
         check_call(args)
         os.chdir(HERE)
         jobs.task_done()
         singlelock.acquire()
         print(str(jobs.qsize())+" jobs remaining")
         singlelock.release()
      except:
         break

def createRunDirs():
   rundirs = [ ]
   HERE = os.getcwd()
   for physics in ['CROCUS','Anderson','CLM45','HL','Pimienta']:
      out_dir = physics
      os.mkdir(out_dir)
      os.chdir(out_dir)
      
      #for case in ['core01', 'core03', 'SipleDome']:
      for case in ['core03']:
         rundirs.append(os.path.join(HERE,out_dir,case))
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
            check_call(['touch','CROCUS'])
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
      
         os.mkdir(case)
         # Writing our configuration file
         with open(os.path.join(case,config_name), 'wb') as configfile:
             config.write(configfile)

      os.chdir(HERE)
   return rundirs

if __name__ == "__main__":
   main()
   jobs.join()
