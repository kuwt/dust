import subprocess as sbprc
import numpy as np
import sys
import fileinput
import argparse
import os
import math as mth
import h5py

pars = argparse.ArgumentParser()
pars.add_argument("exe_path",help='path to the executables')
pars.add_argument('-t','--tolerance', nargs='?', help='tolerance for the test',
                  required=False, default=1e-10, type=float)
args = pars.parse_args()
exe_path = args.exe_path
exe_path = os.path.abspath(exe_path)
pre    = exe_path+'/dust_pre'
solver = exe_path+'/dust'
post   = exe_path+'/dust_post'
tol = args.tolerance

pwd = os.getcwd()
sets_descr = ['Rotor with a cgns nacelle']

#run the first set of simulations, set a
os.chdir("./reg_test_b")
pre_ret = sbprc.call([pre])
runs = ['dust_base.in']
descr = ['basic configuration']
runs_names = ['base']
sol_dsets = ['/ParticleWake/WakePoints', '/ParticleWake/WakeVort', 
      '/Components/Comp001/Solution/Vort', '/Components/Comp001/Solution/Pres',
      '/Components/Comp005/Solution/Vort', '/Components/Comp005/Solution/Pres']
sol_descr = ['Particles position ', 'Particles Intensity', 
             'Blade Intensity    ', 'Blade Pressure     ',
             'Nacelle Intensity  ', 'Nacelle Pressure   ']
errors = np.zeros([len(runs),len(sol_dsets)])
for i, run in enumerate(runs):
  solv_ret = sbprc.call([solver,run])

  ref_name = 'ref_'+runs_names[i]
  check_name = 'test_'+runs_names[i]
  suffix = '.h5'
  folder = './Output/'
  res = '_res_0011'

  ref_filename = folder+ref_name+res+suffix
  check_filename = folder+check_name+res+suffix
  ref_file = h5py.File(ref_filename, 'r')
  check_file = h5py.File(check_filename, 'r')
  for j, dset_name in enumerate(sol_dsets):
    ref_data   = np.array(ref_file[dset_name])
    check_data = np.array(check_file[dset_name])
    err = np.linalg.norm(check_data-ref_data)/np.linalg.norm(ref_data)
    errors[i,j] = err


#erease all the produced files
files = os.listdir('Output/')
for file in files:
  if file.startswith('test'):
    os.remove(os.path.join('Output/',file))
os.remove('geo_input.h5')

#print the errors:
print('Difference w.r.t. reference:')
for i, run in enumerate(descr):
  print('In run ',run,':')
  for j, dset in enumerate(sol_descr):
    print('Difference on ',dset,': ',errors[i,j])


#whatever happened, get back to the previous folder
os.chdir(pwd)
