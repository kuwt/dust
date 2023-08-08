import os
import subprocess as sp

path2dust = os.path.join('dust')
if os.path.isfile(os.path.join(path2dust,'precice-dust-convergence.log')): 
    convergence = os.path.join(path2dust,'precice-dust-convergence.log')
    events = os.path.join(path2dust,'precice-dust-events.json')
    iteration =  os.path.join(path2dust,'precice-dust-iterations.log')

    sp.run('rm -fv ' + convergence, shell=True)
    sp.run('rm -fv ' + events,      shell=True)
    sp.run('rm -fv ' + iteration,   shell=True)
###################################################
# run dust
####################################################
oldpwd = os.getcwd()
dust_input = 'dust.in'
os.chdir(path2dust)
sp.run('dust ' + os.path.join(dust_input) + ' &', shell=True) 
os.chdir(oldpwd)