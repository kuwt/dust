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


#> MBDyn model parameters 
dt_input = 0.001
path2mbdyn = os.path.join('structure') 
path2precice = os.path.join('')

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
# run precice
####################################################

########### load coupling nodes #############3
rr = np.loadtxt('./refConfigNodes.in')
if ( len(rr.shape) == 1 ):
  rr = rr.reshape( 1,3 )
nnodes = np.size(rr,0)

############## setup mbcNodal, the socket communication to mbdyn ##############3
@dataclass
class SocketData:
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

socketData = SocketData(path = socket_path, nnodes = nnodes)
nodal = mbcNodal(socketData.path, \
                socketData.host, \
                socketData.port, \
                socketData.timeout, \
                socketData.verbose, \
                socketData.data_and_next, \
                socketData.refnode, \
                socketData.nnodes, \
                socketData.labels, \
                0x100, \
                socketData.accels )

nodal.negotiate()
nodal.recv()
nodal.n_f[:] = 0.   # initialize something? 
nodal.n_m[:] = 0.

fid = open('./../nnodes.dat', "w")
fid.write('%d' % ( socketData.nnodes ) )
fid.close()

############ setup data dictionary, not sure its exact purpose #################
nodal_data = { "Position":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Velocity":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Rotation":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "AngularVelocity":{ \
                    "type":"vector", "io":"write", "data":[] }, \
                  "Force":{ \
                    "type":"vector", "io":"read" , "data":[] }, \
                  "Moment":{ \
                    "type":"vector", "io":"read" , "data":[] } \
                }

n = socketData.nnodes; nd = 3
nodal_data["Position"       ]["data"] = np.reshape( nodal.n_x    , (n, nd) )
nodal_data["Velocity"       ]["data"] = np.reshape( nodal.n_xp   , (n, nd) )
nodal_data["Rotation"       ]["data"] = np.reshape( nodal.n_theta, (n, nd) )
nodal_data["AngularVelocity"]["data"] = np.reshape( nodal.n_omega, (n, nd) )
nodal_data["Force" ]["data"] = np.zeros((n, nd))
nodal_data["Moment"]["data"] = np.zeros((n, nd))


############  setup precice data ???? #########################33
partic = { \
          'name':'MBDyn', \
          'mesh':{ \
              'name':'MBDynNodes', 'id':[], 'node_ids':[], \
              'nodes':[], 'nnodes':[], 'dim':[] }, \
          'fields':{}  \
          }
field_dict = { 'name':'field_name', 'id':[], 'data':[],
            'type':'scalar/vector', 'io':'read/write' }

fieldlist = [ 'Position', 'Velocity', 'Rotation', \
              'AngularVelocity', 'Force', 'Moment' ]
typelist  = [ 'vector', 'vector', 'vector', \
              'vector', 'vector', 'vector'  ]
iolist    = [ 'write', 'write', 'write', 'write', 'read', \
              'read' ]

############## setup precice interface #################
precice_config_file_name = './../precice-config.xml'
comm_rank = 0
comm_size = 1   

precice_interface = precice.Interface( partic['name'], \
                                    precice_config_file_name, \
                                    comm_rank, comm_size )


############## fill in information to partic #################
partic['mesh']['id'] = \
         precice_interface.get_mesh_id( partic['mesh']['name'] )

for i in np.arange( len(fieldlist) ):
  fieldn = fieldlist[i]
  if ( precice_interface.has_data( fieldn,partic['mesh']['id'] ) ):
    field = field_dict.copy()
    field['name'] = fieldn
    field['id']   = precice_interface.get_data_id( field['name'], \
                                        partic['mesh']['id'] )
    field['type'] = nodal_data[fieldn]['type']
    field['io']   = nodal_data[fieldn]['io']
    partic['fields'][field['name']] = field

partic['mesh']['nodes'] = rr 
partic['mesh']['nnodes'] = np.size( nodal_data['Position']['data'], 0 )
partic['mesh']['dim']    = np.size( nodal_data['Position']['data'], 1 )
partic['mesh']['node_id'] = precice_interface.set_mesh_vertices( \
                partic['mesh']['id'], partic['mesh']['nodes'] )

############# initialize precice #######################
dt_precice = precice_interface.initialize()
precice_interface.initialize_data()

########### run precice #########################
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
  
  if ( precice_interface.is_action_required( cowic ) ):
    pos_t = nodal_data['Position']['data'][:,:] 
    vel_t = nodal_data['Velocity']['data'][:,:] 
    rot_t = nodal_data['Rotation']['data'][:,:] 
    ome_t = nodal_data['AngularVelocity']['data'][:,:]
    precice_interface.mark_action_fulfilled( cowic )

  #> Set MBDyn nodal values of forces and moments
  for i in np.arange(n):
    nodal.n_f[i*nd:(i+1)*nd] = nodal_data['Force' ]['data'][i,:] 
    nodal.n_m[i*nd:(i+1)*nd] = nodal_data['Moment']['data'][i,:] 
  
  dt = min( dt_input, dt_precice )

  #> === Communication with MBDyn ===========================
  if ( nodal.send(False) ):
    break

  # Receive data from MBDyn
  if ( nodal.recv() ):
    print('**** break, after nodal.recv() ****'); break

  #> Read position and velocity from MBDyn
  nodal_data['Position'       ]['data'] = np.reshape( nodal.n_x    , (n, nd)) 
  nodal_data['Velocity'       ]['data'] = np.reshape( nodal.n_xp   , (n, nd))
  nodal_data['Rotation'       ]['data'] = np.reshape( nodal.n_theta, (n, nd))
  nodal_data['AngularVelocity']['data'] = np.reshape( nodal.n_omega, (n, nd))

  #> Write to AeroSolver
  for fie in partic['fields']:
    if ( partic['fields'][fie]['io'] == 'write' ):
      if ( partic['fields'][fie]['type'] == 'scalar' ):
        precice_interface.write_block_scalar_data( \
                                partic['fields'][fie]['id'], \
                                partic['mesh']['node_id'],   \
                                nodal_data[fie]['data'] )
      if ( partic['fields'][fie]['type'] == 'vector' ):
        precice_interface.write_block_vector_data( \
                                partic['fields'][fie]['id'], \
                                partic['mesh']['node_id'],   \
                                nodal_data[fie]['data'] )

  #> advance steo
  is_ongoing = precice_interface.is_coupling_ongoing()
  dt_precice = precice_interface.advance(dt)

  #> Receive data from AeroSolver and set nodal.n_f field
  for fie in partic['fields']:
    if ( partic['fields'][fie]['io'] == 'read' ):
      if ( partic['fields'][fie]['type'] == 'scalar' ):
        nodal_data[fie]['data'] = \
          precice_interface.read_block_scalar_data( \
                                partic['fields'][fie]['id'], \
                                partic['mesh']['node_id']    )
      if ( partic['fields'][fie]['type'] == 'vector' ):
        nodal_data[fie]['data'] = \
          precice_interface.read_block_vector_data( \
                                partic['fields'][fie]['id'], \
                                partic['mesh']['node_id']    )

  #> Check convergence: iterate or finalize the timestep
  if ( precice_interface.is_action_required( coric ) ): # dt not converged
    nodal_data['Position'       ]['data'] = pos_t
    nodal_data['Velocity'       ]['data'] = vel_t
    nodal_data['Rotation'       ]['data'] = rot_t
    nodal_data['AngularVelocity']['data'] = ome_t
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