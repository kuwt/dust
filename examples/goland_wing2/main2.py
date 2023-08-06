import tempfile 
import os, sys
# set to path of MBDyn support for python communication
sys.path.append('/usr/local/mbdyn/libexec/mbpy') 

from mbc_py_interface import mbcNodal

import precice
from precice import *
 
import subprocess as sp
import numpy as np
import pdb
import time

class MBDynInterface:
  """ MBDyn interface collecting kinematic and forces """
  #> Constructor ------------------------------------------------
  def __init__(self):
    self.initialized = False
    self.data = { "Position":{ \
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

  #> Inner class ------------------------------------------------
  class Socket:
    """ Class containing the socket parameters """
    """ for comm between MBDyn and mbc_py      """
    def __init__( self, \
                  path="", host="", port=0, \
                  timeout=-1, verbose=0, data_and_next=1, \
                  refnode=0, nnodes=1, labels = 0, rot=0x100, \
                  accels=0 ):
      self.path          = path
      self.host          = host          
      self.port          = port          
      self.timeout       = timeout       
      self.verbose       = verbose       
      self.data_and_next = data_and_next 
      self.refnode       = refnode       
      self.nnodes        = nnodes        
      self.labels        = labels        
      self.rot           = rot           
      self.accels        = accels

      print(' in Socket __init__. self.nnodes: ', self.nnodes)

  #> Methods ----------------------------------------------------
  #> initialize -------------------------------------------------
  def initialize(self,path='.mbdyn/mbdyn.sock', \
                      verbose=0, \
                      nnodes=1, \
                      accels=0, \
                      forceType='Nodal', \
                      dumpAuxFile=False, \
                      dumpAuxFilen='./../nnodes.dat'):

    #> Initialize communication
    self.socket = self.Socket(path=path, verbose=0, \
                              nnodes=nnodes, accels=accels )

    if ( forceType == 'Nodal' ):
      self.nodal = mbcNodal(self.socket.path, \
                            self.socket.host, \
                            self.socket.port, \
                            self.socket.timeout, \
                            self.socket.verbose, \
                            self.socket.data_and_next, \
                            self.socket.refnode, \
                            self.socket.nnodes, \
                            self.socket.labels, \
                            self.socket.rot, \
                            self.socket.accels )
      print(' self.socket.nnodes: ', self.socket.nnodes)
      #> Negotiate communication and receive kinematics from MBDyn
      self.nodal.negotiate()
      self.nodal.recv()

      if ( dumpAuxFile ):
        fid = open(dumpAuxFilen, "w")
        fid.write('%d' % ( self.socket.nnodes ) )
        fid.close()

      #> Save coordinates of the MBDyn nodes
      n = self.socket.nnodes; nd = 3
      self.data["Position"       ]["data"] = np.reshape( self.nodal.n_x    , (n, nd) )
      self.data["Velocity"       ]["data"] = np.reshape( self.nodal.n_xp   , (n, nd) )
      self.data["Rotation"       ]["data"] = np.reshape( self.nodal.n_theta, (n, nd) )
      self.data["AngularVelocity"]["data"] = np.reshape( self.nodal.n_omega, (n, nd) )

      #> Initialize read fields
      self.data["Force" ]["data"] = np.zeros((n, nd))
      self.data["Moment"]["data"] = np.zeros((n, nd))

      #> Initialize nodal.n_f, n_m
      self.nodal.n_f[:] = 0.;  self.nodal.n_m[:] = 0.
      # self.nodal.send(True)

    else:
      sys.exit(" So far, only forceType='Nodal' is implemented")

  #> finalize ---------------------------------------------------
  def finalize(self):
    self.nodal.destroy()
    
  #> refConfigNodes ---------------------------------------------
  def refConfigNodes(self, filen='./refConfigNodes.in'):
    rr = np.loadtxt( fname = filen )

    if ( len(rr.shape) == 1 ):
      rr = rr.reshape( 1,3 )
    self.rr = rr 
    return rr


class MBDynAdapter:
  class Participant:
    def __init__(self, name='MBDyn'):
      """ Do nothing """
      self.name  = name
      self.field = [] 

    class Mesh:
      def __init__(self, name='MBDynNodes', id=0):
        self.name = 'MBDynNodes'
        self.id   = id
    class Field:
      def __init__(self, name='', id=0, data=[]):
        self.name = name 
        self.id   = id   
        self.data = data 


  def __init__(self, mbdInterface, \
              config_file_name = './../precice-config.xml'):

    self.debug = True

    #> Initialize PreCICE participant dictionary, p for participant
    self.p = { \
              'name':'MBDyn', \
              'mesh':{ \
                  'name':'MBDynNodes', 'id':[], 'node_ids':[], \
                  'nodes':[], 'nnodes':[], 'dim':[] }, \
              'fields':{}  \
              }
    field_dict = { 'name':'field_name', 'id':[], 'data':[],
                'type':'scalar/vector', 'io':'read/write' }
  
    #> Use MBDyn interface, containing info about MBDyn-mbc_py
    self.mbd = mbdInterface

    # All the data that can be exchanged with MBDyn
    fieldlist = [ 'Position', 'Velocity', 'Rotation', \
                  'AngularVelocity', 'Force', 'Moment' ]
    typelist  = [ 'vector', 'vector', 'vector', \
                  'vector', 'vector', 'vector'  ]
    iolist    = [ 'write', 'write', 'write', 'write', 'read', \
                  'read' ]
    
    #> === PreCICE interface ====================================
    comm_rank = 0; comm_size = 1   
    
    self.interface = precice.Interface( self.p['name'], \
                                        config_file_name, \
                                        comm_rank, comm_size )

    self.dim = self.interface.get_dimensions()

    self.p['mesh']['id'] = \
              self.interface.get_mesh_id( self.p['mesh']['name'] )

    #> === Fields ===============================================
    # get_data_id
    for i in np.arange( len(fieldlist) ):
      fieldn = fieldlist[i]
      if ( self.interface.has_data( fieldn, self.p['mesh']['id'] ) ):
        field = field_dict.copy()
        field['name'] = fieldn
        field['id']   = self.interface.get_data_id( field['name'], \
                                            self.p['mesh']['id'] )
        field['type'] = self.mbd.data[fieldn]['type']
        field['io']   = self.mbd.data[fieldn]['io']
        self.p['fields'][field['name']] = field
    
    #> === Mesh =================================================
    #> Nodes: set_mesh_vertices
    self.p['mesh']['nodes'] = self.mbd.rr 
    
    # old: initial configuration as the reference configuration
    self.p['mesh']['nnodes'] = np.size( self.mbd.data['Position']['data'], 0 )
    self.p['mesh']['dim']    = np.size( self.mbd.data['Position']['data'], 1 )
    self.p['mesh']['node_id'] = self.interface.set_mesh_vertices( \
                    self.p['mesh']['id'], self.p['mesh']['nodes'] )
    
    print('self.mbd.socket.nnodes',self.mbd.socket.nnodes)
    if ( self.debug ): 
      for i in np.arange(self.mbd.socket.nnodes):
        print('  ', i,': ', self.p['mesh']['node_id'][i],' , ', \
                            self.p['mesh']['nodes'][i,:] )
    
    #> === Initialize ===========================================
    
    self.dt_precice = self.interface.initialize()
    
    self.is_ongoing = self.interface.is_coupling_ongoing()
    if ( self.debug ): 
      print(' interface.is_coupling_ongoing: \033[0;32m ▶ %s' % self.is_ongoing )

    #> === Initialize data ======================================
    cowid = precice.action_write_initial_data()
    if ( self.interface.is_action_required( cowid ) ):
      if ( self.p['fields'][fie]['io'] == 'write' ):
        if ( self.p['fields'][fie]['type'] == 'scalar' ):
          self.interface.write_block_scalar_data( \
                                  self.p['fields'][fie]['id'], \
                                  self.p['mesh']['node_id'],   \
                                  self.mbd.data[fie]['data'] )
        if ( self.p['fields'][fie]['type'] == 'vector' ):
          self.interface.write_block_vector_data( \
                                  self.p['fields'][fie]['id'], \
                                  self.p['mesh']['node_id'],   \
                                  self.mbd.data[fie]['data'] )
      
      self.interface.mark_action_fulfilled( cowic )

    self.interface.initialize_data()


  def runPreCICE(self, dt_set):

    n = self.mbd.socket.nnodes; nd = 3
    
    cowic = precice.action_write_iteration_checkpoint()
    coric = precice.action_read_iteration_checkpoint()

    dt_precice = self.dt_precice

    force = np.zeros((n, nd))

    t = 0.
    niter = 0.
    is_ongoing = self.interface.is_coupling_ongoing()
    while ( is_ongoing ):
      niter = niter + 1
      if (niter == 1):
        start = time.time()
      print('\033[0m   Iteration ▶ ', niter)
      
      if ( self.interface.is_action_required( cowic ) ):
        pos_t = self.mbd.data['Position']['data'][:,:] 
        vel_t = self.mbd.data['Velocity']['data'][:,:] 
        rot_t = self.mbd.data['Rotation']['data'][:,:] 
        ome_t = self.mbd.data['AngularVelocity']['data'][:,:]
        self.interface.mark_action_fulfilled( cowic )

      #> Set MBDyn nodal values of forces and moments
      for i in np.arange(n):
        self.mbd.nodal.n_f[i*nd:(i+1)*nd] = self.mbd.data['Force' ]['data'][i,:] 
        self.mbd.nodal.n_m[i*nd:(i+1)*nd] = self.mbd.data['Moment']['data'][i,:] 
      
      dt = min( dt_set, dt_precice )
    
      #> === Communication with MBDyn ===========================
      if ( self.mbd.nodal.send(False) ):
        break

      # Receive data from MBDyn
      if ( self.mbd.nodal.recv() ):
        print('**** break, after nodal.recv() ****'); break
    
      #> Read position and velocity from MBDyn
      self.mbd.data['Position'       ]['data'] = np.reshape( self.mbd.nodal.n_x    , (n, nd)) 
      self.mbd.data['Velocity'       ]['data'] = np.reshape( self.mbd.nodal.n_xp   , (n, nd))
      self.mbd.data['Rotation'       ]['data'] = np.reshape( self.mbd.nodal.n_theta, (n, nd))
      self.mbd.data['AngularVelocity']['data'] = np.reshape( self.mbd.nodal.n_omega, (n, nd))

      #> Write to AeroSolver
      for fie in self.p['fields']:
        if ( self.p['fields'][fie]['io'] == 'write' ):
          if ( self.p['fields'][fie]['type'] == 'scalar' ):
            self.interface.write_block_scalar_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id'],   \
                                    self.mbd.data[fie]['data'] )
          if ( self.p['fields'][fie]['type'] == 'vector' ):
            self.interface.write_block_vector_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id'],   \
                                    self.mbd.data[fie]['data'] )
    
      is_ongoing = self.interface.is_coupling_ongoing()
      dt_precice = self.interface.advance(dt)
    
      #> Receive data from AeroSolver and set nodal.n_f field
      for fie in self.p['fields']:
        if ( self.p['fields'][fie]['io'] == 'read' ):
          if ( self.p['fields'][fie]['type'] == 'scalar' ):
            self.mbd.data[fie]['data'] = \
              self.interface.read_block_scalar_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id']    )
          if ( self.p['fields'][fie]['type'] == 'vector' ):
            self.mbd.data[fie]['data'] = \
              self.interface.read_block_vector_data( \
                                    self.p['fields'][fie]['id'], \
                                    self.p['mesh']['node_id']    )

      #> Check convergence: iterate or finalize the timestep
      if ( self.interface.is_action_required( coric ) ): # dt not converged
        self.mbd.data['Position'       ]['data'] = pos_t
        self.mbd.data['Velocity'       ]['data'] = vel_t
        self.mbd.data['Rotation'       ]['data'] = rot_t
        self.mbd.data['AngularVelocity']['data'] = ome_t
        self.interface.mark_action_fulfilled( coric )
      else: # dt converged
        if ( self.mbd.nodal.send(True) ):
          break
        # Receive data from MBDyn
        if ( self.mbd.nodal.recv() ):
          print('**** break, after nodal.recv() ****'); break
        t = t + dt
        niter = 0
        end = time.time()
        print('\033[0m   Elapsed Time:', round((end - start),8), 'sec')
        print('\033[0;32m ------------------------------- ')
        print('\033[0;32m ◀ Simulation Time ▶ ', round(t, 6))
        print('\033[0;32m ------------------------------- \033[0m ')
            
    print(' before finalize() ')
    self.interface.finalize()
    print(' Finalize auxiliary coupling ')
    
    #> ~~~ PreCICE coupling between mbdyn and aero solver ~~~~~~~~~~~   
    self.mbd.nodal.destroy()




















#> MBDyn model parameters 
dt = 0.001
path2mbdyn = os.path.join('structure') 
path2dust = os.path.join('dust')
path2precice = os.path.join('')
mbdyn_model = 'wing.mbd'
################################3
# cleanup previous .log and cached files
###############################
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

if os.path.isfile(os.path.join(path2dust,'precice-dust-convergence.log')): 
    convergence = os.path.join(path2dust,'precice-dust-convergence.log')
    events = os.path.join(path2dust,'precice-dust-events.json')
    iteration =  os.path.join(path2dust,'precice-dust-iterations.log')

    sp.run('rm -fv ' + convergence, shell=True)
    sp.run('rm -fv ' + events,      shell=True)
    sp.run('rm -fv ' + iteration,   shell=True)
    
# remove dummy precice-run folder
if os.path.isdir(os.path.join(path2precice,'precice-run')):
    sp.run('rm -rfv '+ os.path.join(path2precice,'precice-run'), shell=True)
    sp.run('rm -fv ' + os.path.join(path2precice,'nnodes.dat'),  shell=True)


################################3
# set up socket
###############################
tmpdir = tempfile.mkdtemp('', '.mbdyn_')
path = tmpdir + '/mbdyn.sock'
print(' path: ', path)
os.environ['MBSOCK'] = path

###########################3
# run dust
######################
oldpwd = os.getcwd()
os.chdir(path2dust)
sp.run('dust ' + os.path.join('dust.in') + ' &', shell=True) 
os.chdir(oldpwd)

###########################3
# run mbdy
######################
oldpwd = os.getcwd()
try:
    os.chdir(path2mbdyn)  
    sp.run('mbdyn ' + mbdyn_model + ' -o ' + 'Output' + ' 2>&1 &' , shell=True)
    os.chdir(oldpwd)
    cwd = os.getcwd()
except:
    pass

###########################3
# run precice
######################
oldpwd = os.getcwd()
os.chdir(path2mbdyn)
#> Does the other solver needs MBDyn nodes to build its own mesh?
writeNodes = False
#> Construct MBDyn/mbc_py interface
#> Initialize MBDyn/mbc_py interface: negotiate and recv()
mbd = MBDynInterface()
# get node number 
rr = mbd.refConfigNodes() 
nnodes = np.size(rr,0)
# initialize socket for data exchange
mbd.initialize( path=path, verbose=1, nnodes=nnodes, accels=1, \
                dumpAuxFile=True )

#> Send MBDyn exposed nodes to the other solver, if needed
n = mbd.socket.nnodes
print(' n: ', n)
print(" Initialize MBDyn adapter ")
try: 
    adapter = MBDynAdapter(mbd)
except:
    pass
if ( adapter.debug ):
    print(' participant: ', adapter.p["name"] )
    print(' solver     : ', adapter.p["mesh"]["name"] )
    print(' fields     : ', adapter.p["fields"] )
#> Start coupled simulation with PreCICE
try: 
    adapter.runPreCICE(dt) 
except:
    pass
os.environ.pop("MBSOCK")    
os.chdir(oldpwd)
cwd = os.getcwd() 

