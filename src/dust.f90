!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!!
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!!
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!!
!! Authors:
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!          Andrea Colli
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================

!> This is the main file of the DUST solver

program dust

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_basename, &
  check_file_exists

use mod_geometry, only: &
  t_geo, &
  create_geometry, update_geometry, &
  t_tedge,  destroy_geometry, destroy_elements

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_doublet, only: &
  initialize_doublet

use mod_surfpan, only: &
  t_surfpan , initialize_surfpan

use mod_vortlatt, only: &
  t_vortlatt 

use mod_stripe, only: &
  t_stripe

use mod_liftlin, only: &
  update_liftlin, t_liftlin_p, &
  build_ll_kernel, &
  solve_liftlin, solve_liftlin_piszkin

use mod_actuatordisk, only: &
  update_actdisk

use mod_vortline, only: &
  initialize_vortline

use mod_vortpart, only: &
  initialize_vortpart

use mod_c81, only: &
  t_aero_tab

use mod_linsys_vars, only: &
  t_linsys

use mod_linsys, only: &
  initialize_linsys, assemble_linsys, solve_linsys, destroy_linsys, &
  dump_linsys

use mod_pressure_equation, only: &
  dump_linsys_pres, press_normvel_der, initialize_pressure_sys, &
  assemble_pressure_sys, solve_pressure_sys

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

use mod_wake, only: &
  t_wake, initialize_wake, update_wake, &
  prepare_wake, load_wake, complete_wake, destroy_wake

use mod_vtk_out, only: &
  vtk_out_bin

use mod_tecplot_out, only: &
  tec_out_sol_bin

use mod_hdf5_io, only: &
  h5loc, initialize_hdf5, destroy_hdf5, new_hdf5_file, open_hdf5_file, &
  close_hdf5_file, new_hdf5_group, open_hdf5_group, close_hdf5_group, &
  write_hdf5, write_hdf5_attr, read_hdf5, read_hdf5_al, append_hdf5

use mod_dust_io, only: &
  save_status, load_solution, load_time

use mod_viscosity, only: &
  viscosity_effects

use mod_octree, only: &
  initialize_octree, destroy_octree, sort_particles, t_octree, &
  apply_multipole_panels

use mod_math, only: & 
  cross, dot

#if USE_PRECICE
  use mod_precice, only: &
    t_precice
#endif

use mod_wind, only: &
  variable_wind

implicit none

!> Run-id
integer                           :: run_id(10)

!> Input
!> Main parameters parser
type(t_parse)                     :: prms
character(len=*), parameter       :: input_file_name_def = 'dust.in'
character(len=max_char_len)       :: input_file_name
character(len=max_char_len)       :: target_file
character(len=extended_char_len)  :: message

!> Time parameters
real(wp)                          :: time
integer                           :: it, nstep, nout
real(wp)                          :: t_last_out, t_last_debug_out
real(wp)                          :: time_no_out, time_no_out_debug
logical                           :: time_2_out, time_2_debug_out
logical                           :: output_start
real(wp)                          :: dt_debug_out

!> Main variables
!> All the implicit elements, sorted first static then moving
type(t_impl_elem_p), allocatable  :: elems(:)
!> All the explicit elements
type(t_expl_elem_p), allocatable  :: elems_expl(:)
!> Only the lifting line elements
type(t_liftlin_p), allocatable    :: elems_ll(:)
!> Only the actuator disk elements
type(t_expl_elem_p), allocatable  :: elems_ad(:)
!> All the elements (panels+ll)
type(t_pot_elem_p), allocatable   :: elems_tot(:)
!> All the non corrected elements (panels+vl+ll+a-vl_nl)
type(t_pot_elem_p), allocatable   :: elems_non_corr(:)
!> All the corrected vortex lattice elements 
type(t_pot_elem_p), allocatable   :: elems_corr(:) 

!> Geometry
type(t_geo)                       :: geo
!> Trailing edge
type(t_tedge)                     :: te
!> Airfoil table data
type(t_aero_tab), allocatable     :: airfoil_data(:)
!> Linear system
type(t_linsys)                    :: linsys
!> Wake
type(t_wake)                      :: wake
!> Timing vars
real(t_realtime)                  :: t1 , t0, t00, t11, t22
!> I/O prefixes
character(len=max_char_len)       :: frmt, frmt_vl
character(len=max_char_len)       :: basename_debug

real(wp), allocatable             :: res_old(:)
real(wp), allocatable             :: surf_vel_SurfPan_old(:,:)
real(wp), allocatable             :: nor_SurfPan_old(:,:)
real(wp), allocatable             :: al_kernel(:,:), al_v(:)

!> VL viscous correction
integer                           :: i_el, i_c, i_s, i, sel, i_p, i_c2, i_s2
integer                           :: it_vl
real(wp)                          :: tol, diff, max_diff 
real(wp)                          :: d_cd(3), vel(3), v(3), a_v, area_stripe, dforce_stripe(3), e_d(3), e_l(3)
real(wp)                          :: nor(3), tang_cen(3), u_v, q_inf

!> relaxation 
real(wp), allocatable             :: residual_vl(:), residual_vl_old(:), residual_vl_delta(:)
real(wp)                          :: rel_aitken

!> octree parameters
type(t_octree)                    :: octree

#if USE_PRECICE
  !> PreCICE data 
  type(t_precice)                 :: precice
  integer                         :: bool
  logical                         :: precice_convergence
  integer                         :: j
  real(wp)                        :: sum_force(3)
#endif


call printout(nl//'>>>>>> DUST beginning >>>>>>'//nl)

t00 = dust_time()

call get_run_id(run_id)

!> Modules initialization 
call initialize_hdf5()

!> Input reading 
if(command_argument_count().gt.0) then
  call get_command_argument(1, value = input_file_name)
else
  input_file_name = input_file_name_def
endif

call printout(nl//'Reading input parameters from file "'//&
                trim(input_file_name)//'"'//nl)

!> Define the parameters to be read
!> Time
call prms%CreateRealOption('tstart', "Starting time")
call prms%CreateRealOption('tend',   "Ending time")
call prms%CreateRealOption('dt',     "time step")
call prms%CreateIntOption ('timesteps', "number of timesteps")
call prms%CreateRealOption('dt_out', "output time interval")
call prms%CreateRealOption('dt_debug_out', "debug output time interval")
call prms%CreateIntOption ('ndt_update_wake', 'n. dt between two wake updates', '1')

!> Input
call prms%CreateStringOption('geometry_file','Main geometry definition file')
call prms%CreateStringOption('reference_file','Reference frames file','no_set')

!> Output
call prms%CreateStringOption('basename','oputput basename','./')
call prms%CreateStringOption('basename_debug','oputput basename for debug','./')
call prms%CreateLogicalOption('output_start', "output values at starting &
                                                          & iteration", 'F')
call prms%CreateLogicalOption('output_detailed_geo', "output at each &
                    &timestep the detailed geometry in the results file", 'F')
call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', '1')

!> Restart
call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
call prms%CreateStringOption('restart_file','restart file name')
call prms%CreateLogicalOption('reset_time','reset the time from previous execution?','F')

!> Parameters: reference conditions 
call prms%CreateRealArrayOption('u_inf', "free stream velocity", '(/1.0, 0.0, 0.0/)')
call prms%CreateRealOption('u_ref', "reference velocity")             
call prms%CreateRealOption('P_inf', "free stream pressure", '101325')    
call prms%CreateRealOption('rho_inf', "free stream density", '1.225')   
call prms%CreateRealOption('a_inf', "Speed of sound", '340.0')        ! m/s   for dimensional sim
call prms%CreateRealOption('mu_inf', "Dynamic viscosity", '0.000018') ! kg/ms

!> Wake 
call prms%CreateIntOption('n_wake_panels', 'number of wake panels','4')
call prms%CreateIntOption('n_wake_particles', 'number of wake particles', '10000')
call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
                                &particles bounding box', '(/-10.0, -10.0, -10.0/)')
call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
                                &particles bounding box', '(/10.0, 10.0, 10.0/)')
call prms%CreateRealOption('implicit_panel_scale', "Scaling of the first implicit wake panel", '0.3')
call prms%CreateRealOption('implicit_panel_min_vel', "Minimum velocity at the trailing edge", '1.0e-8')
call prms%CreateLogicalOption('rigid_wake','rigid wake?','F')
call prms%CreateRealArrayOption('rigid_wake_vel', "rigid wake velocity" )
call prms%CreateLogicalOption('join_te','join trailing edge','F')
call prms%CreateRealOption('join_te_factor', "join the trailing edges when closer than factor*te element size",'1.0' )

!> Regularisation 
call prms%CreateRealOption('far_field_ratio_doublet', &
      "Multiplier for far field threshold computation on doublet", '10.0')
call prms%CreateRealOption('far_field_ratio_source', &
      "Multiplier for far field threshold computation on sources", '10.0')
call prms%CreateRealOption('doublet_threshold', &
      "Thresold for considering the point in plane in doublets", '1.0e-6')
call prms%CreateRealOption('rankine_rad', &
      "Radius of Rankine correction for vortex induction near core", '0.1')
call prms%CreateRealOption('vortex_rad', &
      "Radius of vortex core, for particles", '0.1')
call prms%CreateRealOption('k_vortex_rad', &
      "Radius coefficient of vortex core, for particles", '1.0') ! default is ON
call prms%CreateRealOption('cutoff_rad', &
      "Radius of complete cutoff  for vortex induction near core", '0.001')

!> Lifting line elements
call prms%CreateStringOption('ll_solver','Solver for the LL elements', &
                          &'GammaMethod')
call prms%CreateLogicalOption('ll_reynolds_corrections', &
                          &'Use Reynolds corrections for the .c81 tables?', 'F')
call prms%CreateRealOption('ll_reynolds_corrections_nfact', &
                          &'Exponent in (Re/Re_T)^n correction', '0.2')
call prms%CreateIntOption('ll_max_iter', &
                          &'Maximum number of iteration in LL algorithm', '100')
call prms%CreateRealOption('ll_tol', 'Tolerance for the relative error in &
                          &fixed point iteration for LL','1.0e-6' )
call prms%CreateRealOption('ll_damp', 'Damping param in fixed point iteration &
                          &for LL used to avoid oscillations','25.0')
call prms%CreateLogicalOption('ll_stall_regularisation', &
                          &'Avoid "unphysical" separations in inner sections of LL?','T')
call prms%CreateIntOption('ll_stall_regularisation_nelems', &
                          &'Number of "unphysical" separations to be removed', '1' )
call prms%CreateIntOption('ll_stall_regularisation_niters', &
                          &'Number of timesteps between two regularisations', '1' )
call prms%CreateRealOption('ll_stall_regularisation_alpha_stall', &
                          &'Stall angle used as threshold for regularisation [deg]', '15.0' )
call prms%CreateRealOption('ll_artificial_viscosity', &
                          &'Constant artificial viscosity for regularisation', '0.0' )
call prms%CreateLogicalOption('ll_artificial_viscosity_adaptive', &
                          &'Adaptive artificial viscosity algorithm', 'F' )
call prms%CreateRealOption('ll_artificial_viscosity_adaptive_alpha', &
                          &'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
call prms%CreateRealOption('ll_artificial_viscosity_adaptive_dalpha', &
                          &'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
call prms%CreateLogicalOption('ll_loads_avl', &
                          &'Use AVL expression for inviscid load computation','T')

!> VL correction parameter 
call prms%CreateRealOption('vl_relax', 'Relaxation factor for rhs update','0.3')
call prms%CreateIntOption('vl_maxiter', &
                          &'Maximum number of iteration in VL algorithm', '100')
call prms%CreateRealOption('vl_tol', 'Tolerance for the absolute error on lift coefficient in &
                          &fixed point iteration for VL','1.0e-4' )
call prms%CreateIntOption('vl_start_step', &
                          &'Step in which the VL correction start', '0')
call prms%CreateLogicalOption('vl_dynstall', 'Dynamic stall on corrected VL', 'F')
call prms%CreateLogicalOption('aitken_relaxation', 'Employ aitken acceleration method during &   
                              &the fixed point iteration', 'T')  

!> Octree and multipole data 
call prms%CreateLogicalOption('fmm','Employ fast multipole method?','T')
call prms%CreateLogicalOption('fmm_panels','Employ fast multipole method &
                              &also for panels?','F')
call prms%CreateRealOption('box_length','length of the octree box')
call prms%CreateIntArrayOption('n_box','number of boxes in each direction')
call prms%CreateRealArrayOption( 'octree_origin', "rigid wake velocity" )
call prms%CreateIntOption('n_octree_levels','number of octree levels')
call prms%CreateIntOption('min_octree_part','minimum number of octree particles')
call prms%CreateIntOption('multipole_degree','multipole expansion degree')
call prms%CreateLogicalOption('dyn_layers','Use dynamic layers','F')
call prms%CreateIntOption('nmax_octree_levels','maximum number of octree levels')
call prms%CreateRealOption('leaves_time_ratio','Ratio that triggers the &
                                          &increase of the number of levels')

!> Models options
call prms%CreateLogicalOption('vortstretch','Employ vortex stretching','T')
call prms%CreateLogicalOption('vortstretch_from_elems','Employ vortex stretching&
                              & from geometry elements','F')
call prms%CreateLogicalOption('divergence_filtering','Employ divergence filtering','T')
call prms%CreateRealOption('filter_time_scale','Filter timescale','40.0')
call prms%CreateLogicalOption('diffusion','Employ vorticity diffusion','T')
call prms%CreateLogicalOption('turbulent_viscosity','Employ turbulent &
                              &viscosity','F')
call prms%CreateLogicalOption('penetration_avoidance','Employ penetration avoidance','F')
call prms%CreateRealOption('penetration_avoidance_check_radius', &
      'Check radius for penetration avoidance','5.0')
call prms%CreateRealOption('penetration_avoidance_element_radius', &
      'Element impact radius for penetration avoidance','1.5')
call prms%CreateLogicalOption('viscosity_effects','Simulate viscosity &
                                                              & effects','F')
call prms%CreateLogicalOption('particles_redistribution','Employ particles &
                                                        &redistribution','F')
call prms%CreateIntOption('octree_level_solid','Level at which the panels &
                          & are considered for particles redistribution')
call prms%CreateRealOption('particles_redistribution_ratio','How many times &
          &a particle need to be smaller than the average of the cell to be&
          & eliminated','3.0')

!> HCAS
call prms%CreateLogicalOption('HCAS','Hover Convergence Augmentation System', 'F')
call prms%CreateRealOption('HCAS_time','HCAS deployment time')
call prms%CreateRealArrayOption('HCAS_velocity','HCAS velocity')

!> Variable wind
call prms%CreateLogicalOption('gust','Gust perturbation','F')
call prms%CreateStringOption('gust_type','Gust model','AMC')
call prms%CreateRealArrayOption('gust_origin','Gust origin point')
call prms%CreateRealArrayOption('gust_front_direction','Gust front direction vector')
call prms%CreateRealArrayOption('gust_front_speed','Gust front speed')
call prms%CreateRealOption('gust_u_des','Design gust velocity')
call prms%CreateRealArrayOption('gust_perturbation_direction','Gust perturbation &
                              &direction vector','(/0.0, 0.0, 1.0/)')
call prms%CreateRealOption('gust_gradient','Gust gradient')
call prms%CreateRealOption('gust_start_time','Gust starting time','0.0')

!> preCICE
#if USE_PRECICE
call prms%CreateStringOption('precice_config','PreCICE configuration file','./../precice-config.xml')
#endif

!> Get the parameters and print them out
call printout(nl//'====== Input parameters: ======')
call check_file_exists(input_file_name,'dust main')
call prms%read_options(input_file_name, printout_val=.true.)

!> Initialize all the parameters
nout = 0  !> Reset the numbering for output files
output_start = getlogical(prms, 'output_start')
call init_sim_param(sim_param, prms, nout, output_start)

!> Remaining parameters
if(countoption(prms, 'dt_debug_out') .lt. 1) then
  dt_debug_out = sim_param%dt_out
else
  dt_debug_out = getreal(prms, 'dt_debug_out')
endif

if(countoption(prms, 'basename_debug') .lt. 1) then
  basename_debug = sim_param%basename
else
  basename_debug = getstr(prms,'basename_debug')
endif

sim_param%basename_debug = basename_debug

#if USE_PRECICE
!> Initialize PreCICE 
! Do it here because it needs to read precice_config path from dust.in
call precice%initialize()
#endif


!> Parameters Initializations 
call initialize_doublet()
call initialize_vortline()
call initialize_surfpan()

!> Check that tend .gt. tinit
if ( sim_param%tend .le. sim_param%t0 ) then
  write(message,*) 'The end time of the simulation',sim_param%tend,'is&
                  & lower or equal to the start time',sim_param%t0,'.&
                  & Remembver that when restarting without resetting the&
                  & time the start time is taken from the restart result!'
  call error('dust','',message)
end if


!> Printout the parameters
if (sim_param%debug_level .ge. 3) then
  write(message,*) 'Initial time tstart: ', sim_param%t0;             call printout(message)
  write(message,*) 'Final time tend:     ', sim_param%tend;           call printout(message)
  write(message,*) 'Time step dt:        ', sim_param%dt;             call printout(message)
  write(message,*) 'Output interval:     ', sim_param%dt_out;         call printout(message)
  write(message,*) 'Output first step:   ', output_start;             call printout(message)
  write(message,*) 'Debug level:         ', sim_param%debug_level;    call printout(message)
  write(message,*) 'Free stream velocity:', sim_param%u_inf;          call printout(message)
  write(message,*) 'Maximum wake panels: ', sim_param%n_wake_panels;  call printout(message)
  write(message,*) 'Results basename:    ', trim(sim_param%basename); call printout(message)
  write(message,*) 'Debug basename:      ', trim(basename_debug);     call printout(message)
endif


!> Check that the basenames are valid 
call check_basename(trim(sim_param%basename),'dust main')
if(sim_param%debug_level.ge.10)  call check_basename(trim(basename_debug),'dust main')

!> Simulation parameters 
nstep = sim_param%n_timesteps
allocate(sim_param%time_vec(sim_param%n_timesteps))
sim_param%time_vec = (/ ( sim_param%t0 + &
          real(i-1,wp)*sim_param%dt, i=1,sim_param%n_timesteps ) /)

!> Geometry creation 
call printout(nl//'====== Geometry Creation ======')

t0 = dust_time()
target_file = trim(sim_param%basename)//'_geo.h5'

call create_geometry(sim_param%GeometryFile, sim_param%ReferenceFile, &
                    input_file_name, geo, te, elems, elems_expl, elems_ad, &
                    elems_ll, elems_corr, elems_non_corr, elems_tot, airfoil_data, target_file, run_id)

t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Created geometry in: ' , t1 - t0,' s.'
  call printout(message)
endif

if(sim_param%debug_level .ge. 15) &
  call debug_printout_geometry_minimal(elems, geo, basename_debug, 0)
if(sim_param%debug_level .ge. 15) &
  call debug_ll_printout_geometry(elems_ll, geo, basename_debug, 0)

!> TODO: check whether to move these calls before, and precisely what they do
if(sim_param%debug_level .ge. 7) call ignoredParameters(prms)
call finalizeParameters(prms)

!> Initialize PreCICE mesh and fields 
#if USE_PRECICE
  te%t_hinged = te%t
  call precice%initialize_mesh( geo )
  call precice%initialize_fields()
  call precicef_initialize(precice%dt_precice)
  call precice%update_elems(geo, elems_tot, te ) ! TEST
#endif

!> Initialization 
if(sim_param%use_fmm) then
  call printout(nl//'====== Initializing Octree ======')
  t0 = dust_time()
  call initialize_octree(sim_param%BoxLength, sim_param%NBox, &
                        sim_param%OctreeOrigin, sim_param%NOctreeLevels, &
                        sim_param%MinOctreePart, sim_param%MultipoleDegree, &
                        sim_param%RankineRad, octree)
  t1 = dust_time()
  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,F9.3,A)') 'Initialized octree in: ' , t1 - t0,' s.'
    call printout(message)
  endif
endif

call printout(nl//'====== Initializing Wake ======')

call initialize_wake(wake, geo, te, sim_param%n_wake_panels, &
      sim_param%n_wake_panels, sim_param%n_wake_particles)

call printout(nl//'====== Initializing Linear System ======')
t0 = dust_time()
call initialize_linsys(linsys, geo, elems, elems_expl, wake ) 

t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Initialized linear system in: ' , t1 - t0,' s.'
  call printout(message)
endif

!> Restart 
if (sim_param%restart_from_file) then
  call load_solution(sim_param%restart_file, geo%components, geo%refs)
  call load_wake(sim_param%restart_file, wake, elems_tot)

else ! Set to zero the intensity of all the singularities

  do i_el = 1 , size(elems)      ! implicit elements (vr, sp)
      elems(i_el)%p%mag = 0.0_wp
  end do
  
  do i_el = 1 , size(elems_expl) ! explicit elements (ll, ad)
    elems_expl(i_el)%p%mag = 0.0_wp
  end do
  if (size(elems_ll) .gt. 0) then
    do i_el = 1, size(elems_ll)
      elems_ll(i_el)%p%Gamma_old = 0.0_wp
    end do
  endif
endif

t22 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') nl//'------ Completed all preliminary operations &
                              &in: ' , t22 - t00,' s.'
  call printout(message)
endif

!> Build kernel for regularisation (averaging) process for LL 
if (( size(elems_ll) .gt. 0 )) then
  allocate(al_v(size(elems_ll))) 
  al_v = 0.0_wp  
  call build_ll_kernel(elems_ll, al_v, al_kernel)
  deallocate(al_v)
end if

!> Initialize coupling 
#if USE_PRECICE
  call precicef_ongoing(precice%is_ongoing)  
  write(*,*) ' is coupling ongoing: ', precice%is_ongoing
  call precicef_ongoing(precice%is_ongoing)
  precice_convergence = .true.
  write(*,*) ' is coupling ongoing: ', precice%is_ongoing
  write(*,*) ' dt_precice         : ', precice%dt_precice
#endif


!=========================== Time Cycle ==============================
call printout(nl//'////////// Performing Computations //////////')
time = sim_param%t0
sim_param%time_old = sim_param%t0 + 1
t_last_out = time
t_last_debug_out = time
time_no_out = 0.0_wp
time_no_out_debug = 0.0_wp

allocate(surf_vel_SurfPan_old(geo%nSurfpan,3)) ; surf_vel_SurfPan_old = 0.0_wp
allocate(     nor_SurfPan_old(geo%nSurfpan,3)) ;      nor_SurfPan_old = 0.0_wp

allocate(res_old(size(elems)))
if ( sim_param % restart_from_file ) then
  do i_el = 1 , size(elems)
    res_old(i_el) = elems(i_el)%p%mag
  end do
else
  res_old = 0.0_wp
end if

!> Start time cycle 
t11 = dust_time()
it = 0
#if USE_PRECICE
it = 1
  do while ( ( it .lt. nstep ) .and. ( precice%is_ongoing .eq. 1 ) ) 
#else
  do while ( ( it .lt. nstep ) )
    it = it + 1
#endif

    sim_param%time_old = sim_param%time

    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,I5,A,I5,A,F9.4)') nl//'--> Step ',it,' of ', &
                                        nstep, ' simulation time: ', time
      call printout(message)
      t22 = dust_time()
      write(message,'(A,F9.3,A)') 'Elapsed wall time: ', t22 - t00
      call printout(message)
    endif

#if USE_PRECICE
    if ( precice_convergence ) then
#endif

    call init_timestep(time)

#if USE_PRECICE
      precice_convergence = .false.
    end if
#endif

    if ( mod( it-1, sim_param%ndt_update_wake ) .eq. 0 ) then
      call prepare_wake(wake, elems_tot, octree)
    end if

    call update_liftlin(elems_ll,linsys)
    call update_actdisk(elems_ad,linsys)

    !> Debug geometry printing
    if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
              call debug_printout_geometry(elems, geo, basename_debug, it)
    if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
              call debug_ll_printout_geometry(elems_ll, geo, basename_debug, it)

#if USE_PRECICE

    call precicef_action_required( precice%write_it_checkp , bool )
    if ( bool .eq. 1 ) then ! Save old state
      !> Save old state: forces and moments
      do j = 1, size(precice%fields)
        if ( trim(precice%fields(j)%fio) .eq. 'write' ) then
          precice%fields(j)%cdata = precice%fields(j)%fdata
        end if
      end do
      !> PreCICE action fulfilled
      call precicef_mark_action_fulfilled( precice%write_it_checkp )
    end if

    !> Read data from structural solver
    do i = 1, size(precice%fields)
      if ( trim(precice%fields(i)%fio) .eq. 'read' ) then
        if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
          call precicef_read_bsdata( precice%fields(i)%fid, &
                                    precice%mesh%nnodes  , &
                                    precice%mesh%node_ids, &
                                    precice%fields(i)%fdata(1,:) )
        elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
            call precicef_read_bvdata( precice%fields(i)%fid, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%fields(i)%fdata )
        endif
      end if
    end do

    !> Update dust geometry ( elems and first wake panels )
    call precice%update_elems( geo, elems_tot, te )

    !> Update geo_data()
    do i_el = 1, size(elems_tot)
      call elems_tot(i_el)%p%calc_geo_data( &
                            geo%points(:,elems_tot(i_el)%p%i_ver) )
    end do

    !> Update near-field wake
    call precice%update_near_field_wake( geo, wake, te )

    !> Update dt--> mbdyn should take care of the dt and send it to precice (TODO)
#else
#endif

    !> Calculate the normal velocity derivative for the pressure equation
    call press_normvel_der(geo, elems, surf_vel_SurfPan_old)

    !>-------------- Assemble the system ------
    t0 = dust_time()
    sel = size(elems) ! total number of elements

    call assemble_linsys(linsys, geo, elems, elems_expl, wake)
    call assemble_pressure_sys(linsys, geo, elems, wake)
    t1 = dust_time()

    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Assembled linear system in: ' , t1 - t0,' s.'
      call printout(message)
    endif


    !debug output of the system
    if ((sim_param%debug_level .ge. 50).and.time_2_debug_out) then
      write(frmt,'(I4.4)') it
      call dump_linsys(linsys, trim(basename_debug)//'A_'//trim(frmt)//'.dat', &
                              trim(basename_debug)//'b_'//trim(frmt)//'.dat' )
      call dump_linsys_pres(linsys, &
                              trim(basename_debug)//'Apres_'//trim(frmt)//'.dat', &
                              trim(basename_debug)//'bpres_'//trim(frmt)//'.dat')
    endif

    !> Solve the pressure system 
    if ( it .gt. 1 .and. geo%nSurfPan .gt. 0 ) then
      call solve_pressure_sys(linsys)
    end if

    !------ Solve the system ------
    t0 = dust_time()
    if (linsys%rank .gt. 0) then
      call solve_linsys(linsys)
    endif
    t1 = dust_time()

    sel = size(elems)

    !> compute dGamma_dt for unsteady contribution
!$omp parallel do private(i_el)
    do i_el = 1 , sel
      elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
    end do
!$omp end parallel do

if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)')  'Solved linear system in: ' , t1 - t0,' s.'
  call printout(message)
endif

!debug print of the results
if (sim_param%debug_level .ge. 20.and.time_2_debug_out) &
                      call debug_printout_result(linsys, basename_debug, it)

  !------ Update the explicit part ------  % v-----implicit elems: p,v
  if ( size(elems_ll) .gt. 0 ) then

    if ( trim(sim_param%llSolver) .eq. 'GammaMethod' ) then ! Gamma-method
      if (sim_param%time .gt. sim_param%time_old)  then
        do i_el = 1, size(elems_ll)
          elems_ll(i_el)%p%Gamma_old_old = elems_ll(i_el)%p%Gamma_old
          elems_ll(i_el)%p%Gamma_old = elems_ll(i_el)%p%mag
        enddo
      endif
      call solve_liftlin(elems_ll, elems_tot, elems , elems_ad , &
              (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it)

    elseif ( trim(sim_param%llSolver) .eq. 'AlphaMethod' ) then
      call solve_liftlin_piszkin(elems_ll, elems_tot, elems , elems_ad , &
            (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it,&
            al_kernel )
  else
    call error('dust','dust',' Wrong string for LLsolver. &
          &This parameter should have been set equal to "GammaMethod" (default) &
          &in init_sim_param() routine. Something went wrong. Stop')
    end if
  end if

!------ Compute loads -------
! Implicit elements: vortex rings and 3d-panels
! 2019-07-23: D.Isola suggested to implement AVL formula for VL elements
! so far, select type() to keep the old formulation for t_surfpan and
! use AVL formula for t_vortlatt
#if USE_PRECICE
  !$omp parallel do private(i_el)
    do i_el = 1 , sel
      ! ifort bugs workaround:
      ! apparently it is not possible to call polymorphic methods inside
      ! select cases for intel, need to call these for all elements and for the
      ! vortex lattices it is going to be a dummy empty function call
      call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
                geo%components(elems(i_el)%p%comp_id)%coupling_node_rot)
      call elems(i_el)%p%compute_dforce()
    end do
  !$omp end parallel do
#else
  !$omp parallel do private(i_el)
    do i_el = 1 , sel
      call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
              geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g)
      call elems(i_el)%p%compute_dforce()
    end do
  !$omp end parallel do
#endif


  ! ifort bugs workaround:
  ! since even if the following calls looks thread safe, they mess up with
  ! ifort and parallel runs, so the cycle is executed another time just for the
  ! vortex lattices
  if ( geo%nVortLatt .gt. 0) then
    !$omp parallel do private(i_el)
      do i_el = 1 , sel      
        select type(el => elems(i_el)%p)        
          class is(t_vortlatt)          
            call el%get_vel_ctr_pt( elems_tot, (/ wake%pan_p, wake%rin_p/), wake%vort_p)
            !> compute dforce using AVL formula
        end select
      end do 
    !$omp end parallel do
    do i_el = 1 , sel      
      select type(el => elems(i_el)%p)        
        class is(t_vortlatt)          
          call el%compute_dforce_jukowski(.true.) 
        ! update the pressure field, p = df.n / area
          ! compute vel at 1/4 chord (some approx, see the comments in the fcn)
          ! update the pressure field, p = df.n / area
          el%pres = sum(el%dforce * el%nor)/el%area
      end select
    end do
  end if 
  
  !> Vl correction for viscous forces 
  if (sim_param%vl_correction) then
    tol = sim_param%vl_tol
    diff = 1.0_wp 
    it_vl = 0
    
    max_diff = tol + 1e-6_wp
    linsys%skip = .true.
    t0 = dust_time()
    

    !> Select time step to start the vl correction 
    !  (avoid strange behaviour at the begining of simulation)
    if (it .gt. sim_param%vl_startstep) then 

      !> allocate residual terms for Aitken acceleration  
      allocate(residual_vl(size(linsys%b)));        residual_vl = 0.0_wp
      allocate(residual_vl_old(size(linsys%b)));    residual_vl_old = 0.0_wp  
      allocate(residual_vl_delta(size(linsys%b)));  residual_vl_delta = 0.0_wp
      
      do while (max_diff .gt. tol .and. it_vl .lt. sim_param%vl_maxiter)
        
        max_diff = 0.0_wp
        diff = 0.0_wp
        
        do i_c = 1, size(geo%components)
          if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
            trim(geo%components(i_c)%aero_correction) .eq. 'true') then 
              
            !> calculate geo data and initial correction 
            if (it_vl .eq. 0) then 
              !> calc geo quantities  
              !$omp parallel do private(i_s)
              do i_s = 1, size(geo%components(i_c)%stripe)
                call geo%components(i_c)%stripe(i_s)%calc_geo_data(geo%components(i_c)%stripe(i_s)%ver) 
              end do 
              !$omp end parallel do 
              !> calc induced velocity from all components: vel_w
              !$omp parallel do private(i_s)
              do  i_s = 1, size(geo%components(i_c)%stripe)
                call geo%components(i_c)%stripe(i_s)%get_vel_ctr_pt(elems_non_corr, (/ wake%pan_p, wake%rin_p /), wake%vort_p)
              end do 
              !$omp end parallel do 
            endif 

            !$omp parallel do private(i_s, i_s2, i_c2, vel, v) schedule(dynamic, 4) 
            do  i_s = 1, size(geo%components(i_c)%stripe)
              !> calc velocity induced by stripe component: vel 
              vel = 0.0_wp
              do i_c2 = 1, size(geo%components)
                if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
                    trim(geo%components(i_c)%aero_correction) .eq. 'true') then 
                  do i_s2 = 1, size(geo%components(i_c2)%stripe)
                    call geo%components(i_c2)%stripe(i_s2)%compute_vel_stripe(geo%components(i_c)%stripe(i_s)%cen , v)
                    vel = vel + v         
                  end do 
                endif 
              end do                                  
              geo%components(i_c)%stripe(i_s)%vel = vel
              call geo%components(i_c)%stripe(i_s)%correction_c81_vortlatt(airfoil_data, linsys, diff, residual_vl, it_vl, i_s)
            !$omp atomic
                max_diff = max(diff, max_diff) 
            !$omp end atomic
            end do
            !$omp end parallel do               
          end if 
        end do 

        !> Debug output of the system
        if ((sim_param%debug_level .ge. 50) .and. time_2_debug_out) then
          write(frmt,'(I4.4)') it
          write(frmt_vl,'(I4.4)') it_vl
          call dump_linsys(linsys,  &
                          trim(basename_debug)//'A_'//trim(frmt)//'_it_'//trim(frmt_vl)//'.dat', &
                          trim(basename_debug)//'b_'//trim(frmt)//'_it_'//trim(frmt_vl)//'.dat' )
        endif

        !> relaxation factor 
        residual_vl_delta = residual_vl - residual_vl_old 
        if (sim_param%rel_aitken .and. it_vl .gt. 2) then 
          rel_aitken = -rel_aitken* dot(residual_vl_old, residual_vl_delta) / &
                                    dot(residual_vl_delta, residual_vl_delta)
          linsys%b = linsys%b + rel_aitken*residual_vl 
        else
          linsys%b = linsys%b + sim_param%vl_relax*residual_vl 
          rel_aitken = -sim_param%vl_relax* dot(residual_vl_old, residual_vl_delta) / & 
                                            dot(residual_vl_delta, residual_vl_delta)
        endif 
        residual_vl_old = residual_vl

        !> Solve the factorized system
        if (linsys%rank .gt. 0) then
          call solve_linsys(linsys)     
        endif

        it_vl = it_vl + 1 
        
        !$omp parallel do private(i_el)
        do i_el = 1 , sel      
          elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
        enddo 
        !$omp end parallel do 

        do i_el = 1, size(elems_non_corr)
          select type(el => elems_non_corr(i_el)%p)        
            class is(t_vortlatt)   
              !> compute dforce using AVL formula with prandtl glauert 
              call el%compute_dforce_jukowski(.true.) 
          end select            
        end do
        do i_el = 1, size(elems_corr)
          select type(el => elems_corr(i_el)%p)
            class is(t_vortlatt)   
              !> compute dforce using AVL formula without prandtl glauert correction since it is 
              !  already contained in the .c81 table 
              call el%compute_dforce_jukowski(.false.) 
          end select
        end do 


      end do !(while)

      deallocate(residual_vl, residual_vl_old, residual_vl_delta)
      if(it_vl .eq. sim_param%vl_maxiter) then
        call warning('dust','dust','max iteration reached for non linear vl:&
                    &increase VLmaxiter!') 
        write(message,*) 'Last iteration error: ', max_diff
        call printout(message)
      endif
        
      !> Viscous and pressure drag correction 
      do i_c = 1, size(geo%components)
        if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
          trim(geo%components(i_c)%aero_correction) .eq. 'true') then 

          do i_s = 1, size(geo%components(i_c)%stripe) 
            
            nor = geo%components(i_c)%stripe(i_s)%nor
            tang_cen = geo%components(i_c)%stripe(i_s)%tang_cen
            a_v = geo%components(i_c)%stripe(i_s)%alpha*pi/180.0_wp
            area_stripe = geo%components(i_c)%stripe(i_s)%area 
            u_v = geo%components(i_c)%stripe(i_s)%vel_2d
            q_inf = 0.5_wp*sim_param%rho_inf * u_v ** 2.0_wp * area_stripe

            d_cd =  0.5_wp * sim_param%rho_inf * u_v**2.0_wp * &                 
                    geo%components(i_c)%stripe(i_s)%cd * &
                    (tang_cen * cos(a_v) + nor * sin(a_v))

            dforce_stripe = 0.0_wp  
            do i_p = 1, size(geo%components(i_c)%stripe(i_s)%panels)
              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce = &
                          geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce + &
                          d_cd * geo%components(i_c)%stripe(i_s)%panels(i_p)%p%area
              !> update pressure 
              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%pres = & 
                          sum(geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce * &
                              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%nor)/ & 
                              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%area 
              
              dforce_stripe = dforce_stripe + &
                              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce
            end do

            !> update cl and cd            
            e_l = nor * cos(a_v) - tang_cen * sin(a_v)
            e_d = nor * sin(a_v) + tang_cen * cos(a_v)
            e_l = e_l / norm2(e_l)
            e_d = e_d / norm2(e_d)

            geo%components(i_c)%stripe(i_s)%aero_coeff(1) = dot(dforce_stripe,e_l) / q_inf 
            geo%components(i_c)%stripe(i_s)%aero_coeff(2) = dot(dforce_stripe,e_d) / q_inf
            
          end do
        end if 
      end do         
    endif
    linsys%skip = .false.
    t1 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Solved nonlinear vortex lattice in: ' , t1 - t0,' s.'
      call printout(message)
    endif
  end if 

  ! Explicit elements:
  ! - liftlin: _pres and _dforce computed in solve_liftin()
  ! - actdisk: avg delta_pressure and force computed here,
  !            to include thier effects in postpro (e.g. integral loads)
!$omp parallel do private(i_el)
  do i_el = 1 , size(elems_ad)
    call elems_ad(i_el)%p%compute_pres( &     ! update surf_vel field too
          geo%refs( geo%components(elems_ad(i_el)%p%comp_id)%ref_id )%R_g)
    call elems_ad(i_el)%p%compute_dforce()
  end do
!$omp end parallel do

#if USE_PRECICE
    sum_force = 0.0_wp
    do i = 1, size(elems_tot)
      sum_force = sum_force + elems_tot(i)%p%dforce
    end do

    !> Update force and moments to be passed to the structural solver
    call precice%update_force(geo, elems_tot)

    call precicef_ongoing(precice%is_ongoing)
    if (precice%is_ongoing .eq. 1) then
      call precicef_advance( precice%dt_precice )
      !write(*,*) ' ++++ dt_precice: ', precice%dt_precice
    end if

    !> Write force and moments to structural solver
    do i = 1, size(precice%fields)
      if (trim(precice%fields(i)%fio) .eq. 'write') then
        if (trim(precice%fields(i)%ftype) .eq. 'scalar') then
          call precicef_write_bsdata( precice%fields(i)%fid, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%fields(i)%fdata(1,:) )
        elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
          call precicef_write_bvdata( precice%fields(i)%fid, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%fields(i)%fdata )
        endif
      end if
    end do

    call precicef_action_required( precice%read_it_checkp, bool )

    if ( bool .eq. 1 ) then ! timestep not converged
      !> Reload checkpoint state
      do j = 1, size(precice%fields)
        if ( trim(precice%fields(j)%fio) .eq. 'write' ) then
          precice%fields(j)%fdata = precice%fields(j)%cdata
        end if
      end do
      call precicef_mark_action_fulfilled( precice%read_it_checkp )
    else ! timestep converged
      precice_convergence = .true.
      !> Finalize timestep
      ! Do the same actions as a simulation w/o coupling
      ! *** to do *** check if something special is needed
#else
#endif

    !> Print the results 
    if(time_2_out)  then
      nout = nout+1
      call save_status(geo, wake, nout, time, run_id)
    endif

    !> Viscous Effects and Flow Separations 
    ! some computation of surface quantities and vorticity to be released.
    ! Free vortices will be introduced in prepare_wake(), some lines below
    if(sim_param%use_ve) call viscosity_effects( geo , elems , te )

    !> Treat the wake: this needs to be done after output, 
    !                  in practice the update is for the next iteration)
  
    t0 = dust_time()
    if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
      call update_wake(wake, elems_tot, octree)
    end if
  
    t1 = dust_time()

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Updated wake in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    !> Pressure integral equation 
    !> save old velocity on the surfpan (before updating the geom, few lines below)
    do i_el = 1 , geo%nSurfPan
      select type ( el => elems(i_el)%p ) ; class is ( t_surfpan )
        surf_vel_SurfPan_old( geo%idSurfPanG2L(i_el) , : ) = el%ub   ! el%surf_vel
            nor_SurfPan_old( geo%idSurfPanG2L(i_el) , : ) = el%nor   ! el%surf_vel
      end select
    end do
    t0 = dust_time()
    !> update geometry 
    if(it .lt. nstep) then
      time = min(sim_param%tend, sim_param%time_vec(it+1))
      !> Update geometry
      call update_geometry(geo, te, time, .false., .true.)
      if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
            call complete_wake(wake, geo, elems_tot, te)
      end if
    endif
    t1 = dust_time() 

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Updated geometry in: ' , t1 - t0,' s.'
      call printout(message)
    endif
    !> Save old solution (at the previous dt) of the linear system
    res_old = linsys%res


    !> Update nor_old (moved from geo/mod_geo.f90/update_geometry(), l.2220 approx
!$omp parallel do private(i_el)
    do i_el = 1 , size(elems_tot)
      elems_tot(i_el)%p%nor_old = elems_tot(i_el)%p%nor
    end do
!$omp end parallel do

    !> Save alpha_old for dynamic stall 
    do i_c = 1, size(geo%components)
      if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
          trim(geo%components(i_c)%aero_correction) .eq. 'true' .and. &
          sim_param%vl_dynstall) then 
!$omp parallel do private(i_s)
          do i_s = 1, size(geo%components(i_c)%stripe)
            geo%components(i_c)%stripe(i_s)%alpha_old = geo%components(i_c)%stripe(i_s)%alpha 
          end do 
!$omp end parallel do  
      end if
    end do 

#if USE_PRECICE
      ! *** to do *** dirty implementation. it-update moved into #if USE_PRECICE
      ! statement, to avoid double time counter update, when the code is not coupled
      ! with external softwares through PRECICE.
      ! #if .not. USE_PRECICE, it-update occurs at the beginning of the time loop,
      ! approximately at l.615.
      ! Is there a reasone why it-update should occur in two different places of the
      ! time loop (at the begin w/o precice, at the end w/ precice)?
      !> Update n. time step
      it = it + 1
    endif ! End of the if statement that check whether the timestep
          ! has converged or not
#endif

  enddo !> while do 
call printout(nl//'\\\\\\\\\\  Computations Finished \\\\\\\\\\')
deallocate(res_old)
!> End Time Cycle 

!> Finalize PreCICE 
#if USE_PRECICE
  call precicef_finalize()
#endif

!> Cleanup 
call destroy_wake(wake)
call destroy_octree(octree)
call destroy_linsys(linsys)
call destroy_elements(geo)
!call destroy_geometry(geo, elems_tot)
call destroy_hdf5()

t22 = dust_time()
!> Debug
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Completed all computations in ',t22-t00,' s'
  call printout(message)
endif
call printout(nl//'<<<<<< DUST end <<<<<<'//nl)

contains

!> Functions (maybe put everything in a dedicated module) 
subroutine get_run_id (run_id)
  integer, intent(out) :: run_id(10)
  real(wp)             :: randr
  integer              :: maxi, randi

  !> First 8 values are the date and time
  call date_and_time(VALUES=run_id(1:8))

  !> Last 3 values are 2 random integers
  maxi = huge(maxi)
  call random_number(randr)
  randi = int(randr*real(maxi,wp))
  run_id(9) = randi

  call random_number(randr)
  randi = int(randr*real(maxi,wp))
  run_id(10) = randi
end subroutine

!> Initialize all the parameters reading them from the the input file
subroutine init_sim_param(sim_param, prms, nout, output_start)
  class(t_sim_param)      :: sim_param
  type(t_parse)           :: prms
  integer, intent(inout)  :: nout
  logical, intent(inout)  :: output_start
  
  !> Timing
  sim_param%t0                  = getreal(prms, 'tstart')
  sim_param%tend                = getreal(prms, 'tend')
  sim_param%dt_out              = getreal(prms,'dt_out')
  sim_param%debug_level         = getint(prms, 'debug_level')
  sim_param%output_detailed_geo = getlogical(prms, 'output_detailed_geo')
  sim_param%ndt_update_wake     = getint(prms, 'ndt_update_wake')
  
  !> Reference environment values
  sim_param%P_inf               = getreal(prms,'P_inf')
  sim_param%rho_inf             = getreal(prms,'rho_inf')
  sim_param%a_inf               = getreal(prms,'a_inf')
  sim_param%mu_inf              = getreal(prms,'mu_inf')
  sim_param%nu_inf              = sim_param%mu_inf/sim_param%rho_inf
  sim_param%u_inf               = getrealarray(prms, 'u_inf', 3)
  !> Check on reference velocity
  if ( countoption(prms,'u_ref') .gt. 0 ) then
    sim_param%u_ref = getreal(prms, 'u_ref')
  else
    sim_param%u_ref = norm2(sim_param%u_inf)
    if (sim_param%u_ref .le. 0.0_wp) then
      call error('dust','dust','No reference velocity u_ref provided but &
      &zero free stream velocity. Provide a non-zero reference velocity. &
      &Stopping now before producing invalid results')
    endif
  end if
  
  !> Wake parameters
  sim_param%n_wake_panels         = getint(prms, 'n_wake_panels')
  sim_param%n_wake_particles      = getint(prms, 'n_wake_particles')
  sim_param%particles_box_min     = getrealarray(prms, 'particles_box_min',3)
  sim_param%particles_box_max     = getrealarray(prms, 'particles_box_max',3)
  sim_param%rigid_wake            = getlogical(prms, 'rigid_wake')
  sim_param%rigid_wake_vel        = sim_param%u_inf   !> initialisation
  !> Check on wake panels
  if(sim_param%n_wake_panels .lt. 1) then
    sim_param%n_wake_panels = 1
    call warning('dust','dust','imposed a number of wake panels rows &
                &LOWER THAN 1. At least one row of panels is mandatory, &
                &the simulation will proceed with "n_wake_panels = 1"')
  endif
  if ( sim_param%rigid_wake ) then
    if ( countoption(prms,'rigid_wake_vel') .eq. 1 ) then
      sim_param%rigid_wake_vel    = getrealarray(prms, 'rigid_wake_vel',3)
    else if ( countoption(prms,'rigid_wake_vel') .le. 0 ) then
      call warning('dust','dust','no rigid_wake_vel parameter set, &
            &with rigid_wake = T; rigid_wake_vel = u_inf')
      sim_param%rigid_wake_vel    = sim_param%u_inf
    end if
  end if

  !> Trailing edge 
  sim_param%join_te               = getlogical(prms, 'join_te')
  if (sim_param%join_te) sim_param%join_te_factor=getreal(prms,'join_te_factor')

  !> Names
  sim_param%basename              = getstr(prms, 'basename')
  sim_param%GeometryFile          = getstr(prms, 'geometry_file')
  sim_param%ReferenceFile         = getstr(prms, 'reference_file')

  !> Method parameters
  sim_param%FarFieldRatioDoublet  = getreal(prms, 'far_field_ratio_doublet')
  sim_param%FarFieldRatioSource   = getreal(prms, 'far_field_ratio_source')
  sim_param%DoubletThreshold      = getreal(prms, 'doublet_threshold')
  sim_param%RankineRad            = getreal(prms, 'rankine_rad')
  sim_param%VortexRad             = getreal(prms, 'vortex_rad')
  sim_param%KVortexRad             = getreal(prms, 'k_vortex_rad')
  sim_param%CutoffRad             = getreal(prms, 'cutoff_rad')
  sim_param%first_panel_scaling   = getreal(prms, 'implicit_panel_scale')
  sim_param%min_vel_at_te         = getreal(prms, 'implicit_panel_min_vel')
  sim_param%use_vs                = getlogical(prms, 'vortstretch')
  sim_param%vs_elems              = getlogical(prms, 'vortstretch_from_elems')
  sim_param%use_vd                = getlogical(prms, 'diffusion')
  sim_param%use_tv                = getlogical(prms, 'turbulent_viscosity')
  sim_param%use_ve                = getlogical(prms, 'viscosity_effects')
  sim_param%use_pa                = getlogical(prms, 'penetration_avoidance')
  !> Check on penetration avoidance 
  if(sim_param%use_pa) then
    sim_param%pa_rad_mult = getreal(prms, 'penetration_avoidance_check_radius')
    sim_param%pa_elrad_mult = getreal(prms,'penetration_avoidance_element_radius')
  endif

  !> Lifting line elements
  sim_param%llSolver                        = getstr(    prms, 'll_solver')
  sim_param%llReynoldsCorrections           = getlogical(prms, 'll_reynolds_corrections')
  sim_param%llReynoldsCorrectionsNfact      = getreal(   prms, 'll_reynolds_corrections_nfact')
  sim_param%llMaxIter                       = getint(    prms, 'll_max_iter'            )
  sim_param%llTol                           = getreal(   prms, 'll_tol'                )
  sim_param%llDamp                          = getreal(   prms, 'll_damp'               )
  sim_param%llStallRegularisation           = getlogical(prms, 'll_stall_regularisation')
  sim_param%llStallRegularisationNelems     = getint(    prms, 'll_stall_regularisation_nelems')
  sim_param%llStallRegularisationNiters     = getint(    prms, 'll_stall_regularisation_niters')
  sim_param%llStallRegularisationAlphaStall = getreal(   prms, 'll_stall_regularisation_alpha_stall')
  sim_param%llArtificialViscosity           = getreal(   prms, 'll_artificial_viscosity')
  sim_param%llArtificialViscosityAdaptive   = getlogical(prms, 'll_artificial_viscosity_adaptive')
  sim_param%llLoadsAVL                      = getlogical(prms, 'll_loads_avl')
  !> check LL inputs
  if ((trim(sim_param%llSolver) .ne. 'GammaMethod') .and. &
      (trim(sim_param%llSolver) .ne. 'AlphaMethod')) then
    write(*,*) ' sim_param%llSolver : ' , trim(sim_param%llSolver)
    call warning('dust','init_sim_param',' Wrong string for LLsolver. &
                &This parameter is set equal to "GammaMethod" (default) &
                &in init_sim_param() routine.')
    sim_param%llSolver = 'GammaMethod'
  end if

  if ( trim(sim_param%llSolver) .eq. 'GammaMethod' ) then
    if ( sim_param%llArtificialViscosity .gt. 0.0_wp ) then
      call warning('dust','init_sim_param','LLartificialViscoisty set as an input, &
          & different from zero, but ll regularisation available only if&
          & LLsolver = AlphaMethod. LLartificialViscosity = 0.0')
      sim_param%llArtificialViscosity = 0.0_wp
      sim_param%llArtificialViscosityAdaptive = .false.
      sim_param%llArtificialViscosityAdaptive_Alpha  = 0.0_wp
      sim_param%llArtificialViscosityAdaptive_dAlpha = 0.0_wp
    end if
    if ( sim_param%llArtificialViscosityAdaptive ) then
      call warning('dust','init_sim_param','LLartificialViscosityAdaptive set&
          & as an input, but ll adaptive regularisation available only if&
          & LLsolver = AlphaMethod')
    end if
  end if

  !> The user is required to set _Alpha and _dAlpha for ll adaptive regularisation
  if (trim(sim_param%llSolver) .eq. 'AlphaMethod') then
    if (sim_param%llArtificialViscosityAdaptive) then
      if ((countoption(prms,'LLartificialViscosityAdaptive_Alpha')  .eq. 0) .or. &
          (countoption(prms,'LLartificialViscosityAdaptive_dAlpha') .eq. 0)) then
        call error('dust','init_sim_param','LLartificialViscosityAdaptive_Alpha or&
          & LLartificialViscosity_dAlpha not set as an input, while LLartificialViscosityAdaptive&
          & is set equal to T. Set these parameters [deg].')
      else
        sim_param%llArtificialViscosityAdaptive_Alpha = getreal(prms,'ll_artificial_viscosity_adaptive_alpha')
        sim_param%llArtificialViscosityAdaptive_dAlpha= getreal(prms,'ll_artificial_viscosity_adaptive_dalpha')
      end if
    end if
  end if
  !write(*,*) ' sim_param%llSolver : ' , trim(sim_param%llSolver)
  !> VL correction 
  sim_param%vl_tol                        = getreal(prms, 'vl_tol')
  sim_param%vl_relax                      = getreal(prms, 'vl_relax')
  sim_param%vl_maxiter                    = getint(prms, 'vl_maxiter')
  sim_param%vl_startstep                  = getint(prms, 'vl_start_step')
  sim_param%rel_aitken                    = getlogical(prms, 'aitken_relaxation')
  !>  VL Dynamic stall
  sim_param%vl_dynstall                   = getlogical(prms, 'vl_dynstall')  
  
  !> Octree and FMM parameters
  sim_param%use_fmm                       = getlogical(prms, 'fmm')

  if(sim_param%use_fmm) then
    sim_param%use_fmm_pan                 = getlogical(prms, 'fmm_panels')
    sim_param%BoxLength                   = getreal(prms, 'box_length')
    sim_param%NBox                        = getintarray(prms, 'n_box',3)
    sim_param%OctreeOrigin                = getrealarray(prms, 'octree_origin',3)
    sim_param%NOctreeLevels               = getint(prms, 'n_octree_levels')
    sim_param%MinOctreePart               = getint(prms, 'min_octree_part')
    sim_param%MultipoleDegree             = getint(prms,'multipole_degree')
    sim_param%use_dyn_layers              = getlogical(prms,'dyn_layers')

    if(sim_param%use_dyn_layers) then
      sim_param%NMaxOctreeLevels          = getint(prms, 'nmax_octree_levels')
      sim_param%LeavesTimeRatio           = getreal(prms, 'leaves_time_ratio')
    else
      sim_param%NMaxOctreeLevels          = sim_param%NOctreeLevels
    endif

    sim_param%use_pr                      = getlogical(prms, 'particles_redistribution')

    if(sim_param%use_pr) then
      sim_param%part_redist_ratio         = getreal(prms,'particles_redistribution_ratio')
      if ( countoption(prms,'OctreeLevelSolid') .gt. 0 ) then
        sim_param%lvl_solid               = getint(prms, 'octree_level_solid')
      else
        sim_param%lvl_solid               = max(sim_param%NOctreeLevels-2,1)
      endif
    endif
  else
    sim_param%use_fmm_pan = .false.
  endif

  !> HCAS
  sim_param%hcas                          = getlogical(prms,'HCAS')
  if(sim_param%hcas) then
    sim_param%hcas_time                   = getreal(prms,'HCAS_time')
    sim_param%hcas_vel                    = getrealarray(prms,'HCAS_velocity',3)
  endif

  !> Variable_wind
  sim_param%use_gust                      = getlogical(prms, 'gust')

  if(sim_param%use_gust) then
    sim_param%GustType                    = getstr(prms,'gust_type')
    sim_param%gust_origin                 = getrealarray(prms, 'gust_origin',3)

    if(countoption(prms,'GustFrontDirection') .gt. 0) then
      sim_param%gust_front_direction      = getrealarray(prms, 'gust_front_direction',3)
    else
      sim_param%gust_front_direction      = sim_param%u_inf
    end if
    
    if(countoption(prms,'GustFrontSpeed') .gt. 0) then
      sim_param%gust_front_speed          = getreal(prms, 'gust_front_speed')
    else
      sim_param%gust_front_speed          = norm2(sim_param%u_inf)
    end if
    sim_param%gust_u_des                  = getreal(prms,'gust_u_des')
    sim_param%gust_perturb_direction      = getrealarray(prms,'gust_perturbation_direction',3)
    sim_param%gust_time                   = getreal(prms,'gust_start_time')
    sim_param%gust_gradient               = getreal(prms,'gust_gradient')
  end if
  
  !> PreCICE
#if USE_PRECICE
    sim_param%precice_config              = getstr(prms,'precice_config')
#endif
  
  !> Manage restart
  sim_param%restart_from_file             = getlogical(prms,'restart_from_file')
  if (sim_param%restart_from_file) then

    sim_param%reset_time                  = getlogical(prms,'reset_time')
    sim_param%restart_file                = getstr(prms,'restart_file')
    
    !> Removing leading "./" if present to avoid issues when restarting
    if(sim_param%basename(1:2) .eq. './') sim_param%basename = sim_param%basename(3:)
    
    if(sim_param%restart_file(1:2) .eq. './') sim_param%restart_file = sim_param%restart_file(3:)
    
    call printout('RESTART: restarting from file: '//trim(sim_param%restart_file))
    sim_param%GeometryFile = sim_param%restart_file(1:len(trim(sim_param%restart_file))-11) //'geo.h5'

    !restarting the same simulation, advance the numbers
    if(sim_param%restart_file(1:len(trim(sim_param%restart_file))-12).eq. &
                                                trim(sim_param%basename)) then
    read(sim_param%restart_file(len(trim(sim_param%restart_file))-6:len(trim(sim_param%restart_file))-3),*) nout
      call printout('Identified restart from the same simulation, keeping the&
                    & previous output numbering')
      !> avoid rewriting the same timestep
      output_start = .false.
    endif
    if(.not. sim_param%reset_time) call load_time(sim_param%restart_file, sim_param%t0)
  endif

  !> Check the number of timesteps
  if(CountOption(prms,'dt') .gt. 0) then
    if( CountOption(prms,'timesteps') .gt. 0) then
      call error('init_sim_param','dust','Both number of timesteps and dt are&
      & set, but only one of the two can be specified')
    else
      !> get dt and compute number of timesteps
      sim_param%dt     = getreal(prms, 'dt')
      sim_param%n_timesteps = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
                              !(+1 for the zero time step)
    endif
  else
    !> get number of steps, compute dt
    sim_param%n_timesteps = getint(prms, 'timesteps')
    sim_param%dt =  (sim_param%tend-sim_param%t0)/&
                      real(sim_param%n_timesteps,wp)
    sim_param%n_timesteps = sim_param%n_timesteps + 1
                            !add one for the first step
  endif

end subroutine init_sim_param

!> Perform preliminary procedures each timestep, mainly chech if it is time
!! to perform output or not
subroutine init_timestep(t)
  real(wp), intent(in) :: t

  sim_param%time = t

  !if (real(t-t_last_out) .ge. real(sim_param%dt_out)) then
  if (real(time_no_out) .ge. real(sim_param%dt_out)) then
    time_2_out = .true.
    t_last_out = t
    time_no_out = 0.0_wp
  else
    time_2_out = .false.
  endif

  !if (real(t-t_last_debug_out) .ge. real(dt_debug_out)) then
  if (real(time_no_out_debug) .ge. real(dt_debug_out)) then
    time_2_debug_out = .true.
    t_last_debug_out = t
    time_no_out_debug = 0.0_wp
  else
    time_2_debug_out = .false.
  endif

  !If it is the last timestep output the solution, unless dt_out is set
  !longer than the whole execution, declaring implicitly that no output is
  !required.
  if((it .eq. nstep) .and. (sim_param%dt_out .le. sim_param%tend)) then
    time_2_out = .true.
    time_2_debug_out = .true.
  endif

  !if requested, print the output also at the first step (t0)
  if ((t.eq.sim_param%t0) .and. output_start) then
    t_last_out = t
    t_last_debug_out = t
    time_2_out = .true.
    time_2_debug_out = .true.
  endif

  time_no_out = time_no_out + sim_param%dt
  time_no_out_debug = time_no_out_debug + sim_param%dt

end subroutine init_timestep

subroutine debug_printout_result(linsys, basename, it)
  type(t_linsys),   intent(in) :: linsys
  character(len=*), intent(in) :: basename
  integer,          intent(in) :: it
  real(wp), allocatable        :: res(:,:)
  character(len=max_char_len)  :: sit

  allocate(res(1,linsys%rank+linsys%n_expl))
  res(1,:) = (/linsys%res,linsys%res_expl(:,1)/)
  write(sit,'(I4.4)') it
  call write_basic(res,trim(basename)//'_result_'//trim(sit)//'.dat')
  deallocate(res)

end subroutine debug_printout_result

subroutine debug_printout_geometry(elems, geo, basename, it)
  type(t_impl_elem_p),   intent(in) :: elems(:)
  type(t_geo),      intent(in)      :: geo
  character(len=*), intent(in)      :: basename
  integer,          intent(in)      :: it
  real(wp), allocatable             :: norm(:,:), cent(:,:), velb(:,:)
  integer, allocatable              :: el(:,:), conn(:,:)
  character(len=max_char_len)       :: sit
  integer                           :: ie, iv
  !> surf_vel and vel_phi
  real(wp), allocatable             :: surf_vel(:,:), vel_phi(:,:)
  integer                           :: i_e
  !> chtls
  real(wp)                          :: f(5)
  integer                           :: n_neigh

  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  !> surf_vel and vel_phi
  allocate( surf_vel(3,size(elems)), vel_phi(3,size(elems)) )
  surf_vel = -666.6_wp ; vel_phi = -666.6_wp

  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    velb(:,ie) = elems(ie)%p%ub
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
    do iv=1,elems(ie)%p%n_ver
      if(associated(elems(ie)%p%neigh(iv)%p)) then
        conn(iv, ie) = elems(ie)%p%neigh(iv)%p%id
      else
        conn(iv, ie) = 0
      endif
    enddo

    !> surf_vel and vel_phi for surfpan only
    select type( el => elems(ie)%p )
      class is(t_surfpan)
        surf_vel(:,ie) = el%surf_vel   ! elems(ie)%p%surf_vel
        !>  2. CHTLS method ------
        vel_phi(:,ie) = 0.0_wp
        f = 0.0_wp ; n_neigh = 0
        do i_e = 1 , el%n_ver
          if ( associated(el%neigh(i_e)%p) ) then
            n_neigh = n_neigh + 1
            f(n_neigh) = - ( el%neigh(i_e)%p%mag - el%mag )
          end if
        end do
        f(n_neigh+1) = sum(el%nor * (-variable_wind(el%cen, sim_param%time) - el%uvort + el%ub) )
        vel_phi(:,ie) = matmul( el%chtls_stencil , f(1:n_neigh+1) )
    end select
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(velb,       trim(basename)//'_mesh_velb_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_mesh_elems_'  //trim(sit)//'.dat')
  call write_basic(conn,       trim(basename)//'_mesh_conn_'   //trim(sit)//'.dat')
  deallocate(norm, cent, el, conn, velb)

  ! surf_vel and vel_phi
  call write_basic(surf_vel,   trim(basename)//'_mesh_surfvel_'  //trim(sit)//'.dat')
  call write_basic(vel_phi ,   trim(basename)//'_mesh_velphi_'   //trim(sit)//'.dat')
  deallocate(surf_vel,vel_phi)

end subroutine debug_printout_geometry

!------------------------------------------------------------------------------

subroutine debug_printout_geometry_minimal(elems,geo,basename, it)
  type(t_impl_elem_p),   intent(in) :: elems(:)
  type(t_geo),      intent(in)      :: geo
  character(len=*), intent(in)      :: basename
  integer,          intent(in)      :: it
  real(wp), allocatable             :: norm(:,:), cent(:,:)
  integer, allocatable              :: el(:,:)
  character(len=max_char_len)       :: sit
  integer                           :: ie
  integer(h5loc)                    :: h5fid

  allocate(norm(3,size(elems)), cent(3,size(elems)))
  allocate(el(4,size(elems))); el = 0

  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_mesh_elems_'  //trim(sit)//'.dat')

  call new_hdf5_file(trim(basename)//'_geo_'  //trim(sit)//'.h5',h5fid)
  call write_hdf5(geo%points,'points',h5fid)
  call close_hdf5_file(h5fid)

  deallocate(norm, cent, el)

end subroutine debug_printout_geometry_minimal

! ---------------------------------------------------------------------

subroutine debug_ll_printout_geometry(elems, geo, basename, it)
  type(t_liftlin_p),   intent(in)  :: elems(:)
  type(t_geo),      intent(in)     :: geo
  character(len=*), intent(in)     :: basename
  integer,          intent(in)     :: it
  real(wp), allocatable            :: norm(:,:), cent(:,:), velb(:,:)
  integer, allocatable             :: el(:,:), conn(:,:)
  character(len=max_char_len)      :: sit
  integer                          :: ie, iv
  
  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  
  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    velb(:,ie) = elems(ie)%p%ub
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
    do iv=1,elems(ie)%p%n_ver
      if(associated(elems(ie)%p%neigh(iv)%p)) then
        conn(iv, ie) = elems(ie)%p%neigh(iv)%p%id
      else
        conn(iv, ie) = 0
      endif
    enddo
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_ll_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_ll_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(velb,       trim(basename)//'_ll_mesh_velb_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_ll_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_ll_mesh_elems_'  //trim(sit)//'.dat')
  call write_basic(conn,       trim(basename)//'_ll_mesh_conn_'   //trim(sit)//'.dat')
  deallocate(norm, cent, el, conn, velb)

end subroutine debug_ll_printout_geometry

end program dust




