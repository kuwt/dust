import os

import tempfile 
import os, sys
import precice
from precice import *

import subprocess as sp
import numpy as np
#import pdb
import time
from dataclasses import dataclass
import string
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy') 
from mbc_py_interface import mbcNodal


##############################################
# Input parameters
###############################################
dt_input = 0.001
path2mbdyn = os.path.join('structure') 
path2precice = os.path.join('')

class CouplingTemporaryDataStorage :
    def __init__(self):
      self.namelist = [] 
      self.typelist = []  # currently only vector supported
      self.iolist = []
      self.dataDictionary = {}

      self.addField('Position','vector','write')
      self.addField('Velocity','vector','write')
      self.addField('Rotation','vector','write')
      self.addField('AngularVelocity','vector','write')
      self.addField('Force','vector','read')
      self.addField('Moment','vector','read')

    def addField(self, name, data_type, io_type):
      self.namelist.append(name)
      self.typelist.append(data_type)
      self.iolist.append(io_type)
      
    def setData(self, name, data):
      if name in self.namelist:
        self.dataDictionary[name] = data
      else:
        print("Error setting {}. No such item".format(name))
        
    def getData(self, name):
      if name in self.namelist:
        return self.dataDictionary.get(name)
      else:
        print("Error getting {}. No such item".format(name))
        return []
        
###################################################
# cleanup previous .log and cached files
####################################################
if os.path.isfile(os.path.join(path2mbdyn,'precice-MBDyn-iterations.json')):
    events = os.path.join(path2mbdyn,'precice-MBDyn-events.json')
    iteration = os.path.join(path2mbdyn,'precice-MBDyn-iterations.log') 
    sp.run('rm -fv ' + events,      shell=True)
    sp.run('rm -fv ' + iteration,   shell=True)

if os.path.isdir(os.path.join(path2mbdyn,'__pycache__')):      
    pycache = os.path.join(path2mbdyn,'__pycache__')
    sp.run('rm -rfv ' + pycache, shell=True)

#remove temporary folder 
sp.run('rm -rfv ' + os.path.join(path2mbdyn,'mbdyn','.mbdyn_*'),        shell=True)

# remove dummy precice-run folder
if os.path.isdir(os.path.join(path2precice,'precice-run')):
    sp.run('rm -rfv '+ os.path.join(path2precice,'precice-run'), shell=True)
    sp.run('rm -fv ' + os.path.join(path2precice,'nnodes.dat'),  shell=True)


###################################################
# set up socket
####################################################
tmpdir = tempfile.mkdtemp('', '.mbdyn_')
socket_path = tmpdir + '/mbdyn.sock'
print(' socket_path: ', socket_path)
os.environ['MBSOCK'] = socket_path

###################################################
# run mbdyn
####################################################
mbdyn_model = 'wing.mbd'
os.chdir(path2mbdyn)  
sp.run('mbdyn ' + mbdyn_model + ' -o ' + 'Output' + ' 2>&1 &' , shell=True)

###################################################
# precice setup
####################################################

########### load coupling nodes #############3
rr = np.loadtxt('./refConfigNodes.in')
if ( len(rr.shape) == 1 ):
  rr = rr.reshape( 1,3 )
nnodes = np.size(rr,0)

############## setup mbcNodal, the socket communication to mbdyn ##############3
@dataclass
class MbdynSocketData:
    path: string = "",
    host: string = ""
    port: int = 0
    timeout: int = -1 
    verbose: int = 0      
    data_and_next: int = 1
    refnode: int = 0
    nnodes: int = 1
    labels: int = 0
    rot: int = 0x100 #!!!!!!!!!!!!!!!!
    accels:int = 1

socketData = MbdynSocketData(path = socket_path, nnodes = nnodes)
nodal = mbcNodal(socketData.path, \
                socketData.host, \
                socketData.port, \
                socketData.timeout, \
                socketData.verbose, \
                socketData.data_and_next, \
                socketData.refnode, \
                socketData.nnodes, \
                socketData.labels, \
                socketData.rot, \
                socketData.accels )

nodal.negotiate()
nodal.recv()
nodal.n_f[:] = 0.   # initialize something? 
nodal.n_m[:] = 0.

fid = open('./../nnodes.dat', "w")
fid.write('%d' % ( socketData.nnodes ) )
fid.close()

############# init Temporary data storage ###################
couplingData = CouplingTemporaryDataStorage()

n = socketData.nnodes; nd = 3
couplingData.setData("Position",np.reshape( nodal.n_x , (n, nd) ))
couplingData.setData("Velocity",np.reshape( nodal.n_xp , (n, nd) ))
couplingData.setData("Rotation",np.reshape( nodal.n_theta , (n, nd) ))
couplingData.setData("AngularVelocity",np.reshape( nodal.n_omega , (n, nd) ))
couplingData.setData("Force",np.zeros((n, nd)))
couplingData.setData("Moment",np.zeros((n, nd)))

############## setup precice interface #################
participant_name = 'MBDyn'
precice_config_file_name = './../precice-config.xml'
comm_rank = 0
comm_size = 1   
precice_interface = precice.Interface( participant_name, precice_config_file_name, comm_rank, comm_size )

mesh_id = precice_interface.get_mesh_id( 'MBDynNodes')
vertex_ids = precice_interface.set_mesh_vertices( mesh_id, rr )

dt_precice = precice_interface.initialize()
precice_interface.initialize_data()

################################################
###         Main loop                     ########
###################################################
cowic = precice.action_write_iteration_checkpoint()
coric = precice.action_read_iteration_checkpoint()
force = np.zeros((n, nd))

t = 0.
niter = 0.
is_ongoing = precice_interface.is_coupling_ongoing()
while ( is_ongoing ):
  niter = niter + 1
  if (niter == 1):
    start = time.time()
  print('\033[0m   Iteration ▶ ', niter)
  
  #> save the data for reloading purpose
  if ( precice_interface.is_action_required( cowic ) ):
    pos_t = couplingData.getData("Position")
    vel_t = couplingData.getData("Velocity")
    rot_t = couplingData.getData("Rotation")
    ome_t = couplingData.getData("AngularVelocity")
    precice_interface.mark_action_fulfilled( cowic )

  #> Set MBDyn nodal values of forces and moments
  for i in np.arange(n):
    nodal.n_f[i*nd:(i+1)*nd] = couplingData.getData("Force")[i,:] 
    nodal.n_m[i*nd:(i+1)*nd] = couplingData.getData("Moment")[i,:] 
  
  dt = min( dt_input, dt_precice )

  #> Communication with MBDyn 
  if ( nodal.send(False) ):
    break

  # Receive data from MBDyn
  if ( nodal.recv() ):
    print('**** break, after nodal.recv() ****'); break

  #> Read position and velocity from MBDyn
  couplingData.setData("Position",np.reshape( nodal.n_x , (n, nd) ))
  couplingData.setData("Velocity",np.reshape( nodal.n_xp , (n, nd) ))
  couplingData.setData("Rotation",np.reshape( nodal.n_theta , (n, nd) ))
  couplingData.setData("AngularVelocity",np.reshape( nodal.n_omega , (n, nd) ))

  #> Write to precice (i.e. AeroSolver)
  data_id = precice_interface.get_data_id( 'Position', mesh_id )
  precice_interface.write_block_vector_data( data_id, vertex_ids,couplingData.getData("Position"))
  data_id = precice_interface.get_data_id( 'Velocity', mesh_id )
  precice_interface.write_block_vector_data( data_id, vertex_ids,couplingData.getData("Velocity"))
  data_id = precice_interface.get_data_id( 'Rotation', mesh_id )
  precice_interface.write_block_vector_data( data_id, vertex_ids,couplingData.getData("Rotation"))
  data_id = precice_interface.get_data_id( 'AngularVelocity', mesh_id )
  precice_interface.write_block_vector_data( data_id, vertex_ids,couplingData.getData("AngularVelocity"))

  #> advance step
  is_ongoing = precice_interface.is_coupling_ongoing()
  dt_precice = precice_interface.advance(dt)

  #> Receive data from precice (i.e. AeroSolver) and then set nodal.n_f field
  data_id = precice_interface.get_data_id( 'Force', mesh_id )
  data = precice_interface.read_block_vector_data( data_id, vertex_ids)
  couplingData.setData("Force",data)
  data_id = precice_interface.get_data_id( 'Moment', mesh_id )
  data = precice_interface.read_block_vector_data( data_id, vertex_ids)
  couplingData.setData("Moment",data)

  #> Check convergence: iterate or finalize the timestep
  if ( precice_interface.is_action_required( coric ) ): # dt not converged, reload
    couplingData.setData("Position",pos_t)
    couplingData.setData("Velocity",vel_t)
    couplingData.setData("Rotation",rot_t)
    couplingData.setData("AngularVelocity",ome_t)

    precice_interface.mark_action_fulfilled( coric )
  else: # dt converged
    if (nodal.send(True) ):
      break
    # Receive data from MBDyn
    if (nodal.recv() ):
      print('**** break, after nodal.recv() ****'); break
    t = t + dt
    niter = 0
    end = time.time()
    print('\033[0m   Elapsed Time:', round((end - start),8), 'sec')
    print('\033[0;32m ------------------------------- ')
    print('\033[0;32m ◀ Simulation Time ▶ ', round(t, 6))
    print('\033[0;32m ------------------------------- \033[0m ')
        

precice_interface.finalize()
nodal.destroy()
os.environ.pop("MBSOCK")    