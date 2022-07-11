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
!!=========================================================================


!> Module to define the structrues containing the simulation parameters
!!
module mod_sim_param

use mod_param, only: &
  wp, max_char_len

use mod_hdf5_io, only: &
  h5loc, &
  write_hdf5_attr

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

implicit none

public :: t_sim_param, sim_param, create_param_main, create_param_pre, &
          create_param_post

private

type t_sim_param

  !For debugging purpose
  character(max_char_len) :: basename_debug

  !Time:
  !> Start time
  real(wp) :: t0
  !> Time step
  real(wp) :: dt
  !> Final time
  real(wp) :: tend
  !> Number of timesteps
  integer  :: n_timesteps
  !> Vector of time instants
  real(wp) , allocatable :: time_vec(:)
  !> Actual time
  real(wp) :: time
  !> Previous time
  real(wp) :: time_old
  !> ndt between 2 wake updates
  integer :: ndt_update_wake

  !> Output detailed geometry each timestep
  logical :: output_detailed_geo

  !Physical parameters:
  !> Free stream pressure
  real(wp) :: P_inf
  !> Free stream density
  real(wp) :: rho_inf
  !> Free stream velocity
  real(wp) :: u_inf(3)
  !> Reference velocity (magnitude of u_inf unless specified)
  real(wp) :: u_ref
  !> Free stream speed of sound
  real(wp) :: a_inf
  !> Free stream dynamic viscosity
  real(wp) :: mu_inf
  !> Free stream kinematic viscosity
  real(wp) :: nu_inf

  !Wake
  !> Scaling of the first implicit panel
  real(wp) :: first_panel_scaling
  !> Minimum velocity at the trailing edge
  real(wp) :: min_vel_at_te
  !> Is the wake rigid?
  logical :: rigid_wake
    !> Velocity of the rigid wake
    real(wp)  :: rigid_wake_vel(3)
  !> Number of wake panels
  integer :: n_wake_panels
  !> Number of wake particles
  integer :: n_wake_particles
  !> Minimum and maximum of the particles box
  real(wp) :: particles_box_min(3)
  real(wp) :: particles_box_max(3)
  !> Join close trailing edges
  logical :: join_te
    !> All trailing edges closer than join_te_factor will be joined
    real(wp) :: join_te_factor

  !Method parameters
  !> Multiplier for far field threshold computation on doublet
  real(wp) :: FarFieldRatioDoublet
  !> Multiplier for far field threshold computation on sources
  real(wp) :: FarFieldRatioSource
  !> Thresold for considering the point in plane in doublets
  real(wp) :: DoubletThreshold
  !> Rankine Radius for vortices
  real(wp) :: RankineRad
  !> Vortex Radius for vortex particles
  real(wp) :: VortexRad
  !> Vortex Radius coefficient for vortex particles
  !> if too low or negative reverts to original behaviour (VortexRad)
  real(wp) :: KVortexRad  
  !> Complete cutoff radius
  real(wp) :: CutoffRad
  !> use the vortex stretching or not
  logical :: use_vs
  !> use the vortex stretching from elements
  logical :: vs_elems
  !> use the divergence filtering
  logical :: use_divfilt
  !> time scale of the divergence filter
  real(wp) :: filt_eta
  !> use the vorticity diffusion or not
  logical :: use_vd
  !> use turbulent viscosity or not
  logical :: use_tv
  !> use the penetration avoidance
  logical :: use_pa
  !> check radius for penetration avoidance
  real(wp) :: pa_rad_mult
  !> element impact radius for penetration avoidance
  real(wp) :: pa_elrad_mult
  !> simulate viscosity effects or not
  logical :: use_ve

  !Lifting Lines
  character(len=max_char_len) :: llSolver
  !> Reynolds corrections of .c81 tables
  logical  :: llReynoldsCorrections
  !> n factor for Reynolds corrections of .c81 tables: (Re/Re_T)^n
  real(wp) :: llReynoldsCorrectionsNfact
  !> Maximum number of iteration in LL algorithm
  integer  :: llMaxIter
  !> Tolerance for the relative error in fixed point iteration for LL
  real(wp) :: llTol
  !> Damping param in fixed point iteration for LL used to avoid oscillations
  real(wp) :: llDamp
  !> Avoid "unphysical" separations in inner sections of LL? :: llTol
  logical  :: llStallRegularisation
  !> Number of "unphysical" separations that can be removed
  integer  :: llStallRegularisationNelems
  !> Number of iterations between two regularisation processes
  integer  :: llStallRegularisationNiters
  !> Reference stall AOA for regularisation
  real(wp) :: llStallRegularisationAlphaStall
  !> Constant Artificial Viscosity for regularisation
  real(wp) :: llArtificialViscosity
  !> Adaptive Artificial Viscosity algorithm
  logical  :: llArtificialViscosityAdaptive
  !> Adaptive Artificial Viscosity algorithm, reference AOA
  real(wp) :: llArtificialViscosityAdaptive_Alpha
  !> Adaptive Artificial Viscosity algorithm, blending interval
  ! between full regularisation and no regularisation ( AOA )
  real(wp) :: llArtificialViscosityAdaptive_dAlpha
  !> Use AVL expression for inviscid load computation ( ~ VL )
  logical  :: llLoadsAVL

  
  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
    !> Employing the FMM method also for panels
    logical :: use_fmm_pan
    !> Size of the Octree box
    real(wp) :: BoxLength
    !> Number of boxes in each direction
    integer :: NBox(3)
    !> Origin of the Octree system of boxes
    real(wp) :: OctreeOrigin(3)
    !> Number of Octree levels
    integer :: NOctreeLevels
    !> Minimum number of particles for each box
    integer :: MinOctreePart
    !> Multipole expansion degree
    integer :: MultipoleDegree
    !> Use dynamic levels
    logical :: use_dyn_layers
      !> Maximum number of octree levels
      integer :: NMaxOctreeLevels
      !> Time ratio that triggers the increase of levels
      real(wp) :: LeavesTimeRatio
    !> use particles redistribution
    logical :: use_pr
      !> Level at which is checked the presence of panels
      integer :: lvl_solid
      real(wp) :: part_redist_ratio

  !HCAS parameters
  !> Use hcas
  logical :: hcas
    !> Time of deployment of the hcas
    real(wp) :: hcas_time
    !> Velocity of the hcas
    real(wp) :: hcas_vel(3)


  !Handling parameters:
  !> Debug level
  integer :: debug_level
  !> Output interval
  real(wp) :: dt_out
  !> Basename
  character(len=max_char_len) :: basename
  !> Geometry file
  character(len=max_char_len) :: GeometryFile
  !> References file
  character(len=max_char_len) :: ReferenceFile
  !> Restart from file
  logical :: restart_from_file
  !> Restart file
  character(len=max_char_len) :: restart_file
  !> Reset the time after restart
  logical :: reset_time

  !Variable wind
  !> Gust
  logical :: use_gust
  !> Gust type
  character(len=max_char_len) :: GustType
  !> Gust parameters
  real(wp) :: gust_origin(3)
  real(wp) :: gust_front_direction(3)
  real(wp) :: gust_front_speed
  real(wp) :: gust_u_des
  real(wp) :: gust_perturb_direction(3)
  real(wp) :: gust_gradient
  real(wp) :: gust_time

  !> Vl correction 
  logical   :: vl_correction = .false. 
  real(wp)  :: vl_tol
  real(wp)  :: vl_relax
  integer   :: vl_maxiter
  integer   :: vl_startstep
  integer   :: vl_iter_ave 
  logical   :: vl_dynstall = .false.
  logical   :: rel_aitken  = .false. 
  logical   :: vl_ave      = .false.

  !> PreCICE
#if USE_PRECICE
  character(len=max_char_len) :: precice_config
#endif
  
contains
  procedure, pass(this) :: save_param => save_sim_param
end type t_sim_param

type(t_sim_param) :: sim_param

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Subroutines for main dust, dust_pre and dust_post prms creation 
!  to avoid clutter in source files
!> Create parameter object for dust main parameters
subroutine create_param_main(prms)
  type(t_parse), intent(inout) :: prms
  
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
  call prms%CreateIntOption('n_wake_panels', 'number of wake panels','1')
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
                            &'Use AVL expression for inviscid load computation','F')
  
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
  call prms%CreateLogicalOption('vl_average', 'Average panel intensity between the last iterations', 'F')  
  call prms%CreateIntOption('vl_average_iter', 'Number of iterations to average', '10')  
  
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

end subroutine create_param_main

!> Create parameter object for dust preprocessor parameters
subroutine create_param_pre(prms)
  type(t_parse), intent(inout) :: prms
    
  call prms%CreateStringOption('comp_name','Component Name', multiple=.true.)
  call prms%CreateStringOption('geo_file','Geometry definition files', multiple=.true.)
  call prms%CreateStringOption('ref_tag','Reference Tag of the component', multiple=.true.)
  call prms%CreateRealOption('tol_se_wing','Global parameter for closing gaps','0.001')
  call prms%CreateRealOption('inner_product_te','Global parameter for edge identification','-0.5')
  call prms%CreateStringOption('file_name','Preprocessor output file')

end subroutine create_param_pre

!> Create parameter object for dust postprocessor parameters
subroutine create_param_post(prms, sbprms, bxprms)
  type(t_parse), intent(inout) :: prms
  type(t_parse), pointer, intent(inout) :: sbprms , bxprms

  call prms%CreateStringOption('basename','Base name of the processed data')
  call prms%CreateStringOption('data_basename','Base name of the data to be &
                                &processed')
  
  call prms%CreateRealOption( 'far_field_ratio_doublet', &
        "Multiplier for far field threshold computation on doublet", '10.0')
  call prms%CreateRealOption( 'far_field_ratio_source', &
        "Multiplier for far field threshold computation on sources", '10.0')
  call prms%CreateRealOption( 'doublet_threshold', &
        "Thresold for considering the point in plane in doublets", '1.0e-6')
  call prms%CreateRealOption( 'rankine_rad', &
        "Radius of Rankine correction for vortex induction near core", '0.1')
  call prms%CreateRealOption( 'vortex_rad', &
        "Radius of vortex core, for particles", '0.1')
  call prms%CreateRealOption( 'cutoff_rad', &
        "Radius of complete cutoff  for vortex induction near core", '0.001')
  
  call prms%CreateSubOption('analysis','Definition of the motion of a frame', &
                            sbprms, multiple=.true.)
  call sbprms%CreateStringOption('type','type of analysis')
  call sbprms%CreateStringOption('name','specification of the analysis')
  call sbprms%CreateIntOption('start_res', 'Starting result of the analysis')
  call sbprms%CreateIntOption('end_res', 'Final result of the analysis')
  call sbprms%CreateIntOption('step_res', 'Result stride of the analysis')
  call sbprms%CreateLogicalOption('average', 'Perform time averaging','F')
  call sbprms%CreateLogicalOption('wake', 'Output also the wake for &
                                  &visualization','T')
  call sbprms%CreateLogicalOption('separate_wake', 'Output the wake in a separate &
                                  &way','F')
  call sbprms%CreateStringOption('format','Output format')
  call sbprms%CreateStringOption('component','Component to analyse', &
                                multiple=.true.)
  call sbprms%CreateStringOption('hinge_tag','Hinge to analyse', &
                                multiple=.true.)                              
  call sbprms%CreateStringOption('variable','Variables to be saved: velocity, pressure or&
                                & vorticity', multiple=.true.)
  
  ! probe output -------------
  call sbprms%CreateStringOption('input_type','How to specify probe coordinates',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('point','Point coordinates in dust_post.in',&
                                multiple=.true.)
  call sbprms%CreateStringOption('file','File containing the coordinates of the probes',&
                                multiple=.true.)
  ! flow field output --------
  call sbprms%CreateIntArrayOption( 'n_xyz','number of points per coordinate',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('min_xyz','lower bounds of the box',&
                                multiple=.true.)
  call sbprms%CreateRealArrayOption('max_xyz','upper bounds of the box',&
                                multiple=.true.)
  ! loads --------------------
  call sbprms%CreateStringOption('comp_name','Components where loads are computed',&
                                multiple=.true.)
  call sbprms%CreateStringOption('reference_tag','Reference frame where loads&
                              & are computed',multiple=.true.)
  ! sectional loads ----------
  call sbprms%CreateRealArrayOption('axis_dir','Direction of the axis defined the reference&
                              & points for sectional loads analisys', multiple=.true.)
  call sbprms%CreateRealArrayOption('axis_nod','Node belonging to the axis used for sectional&
                              & loads analisys', multiple=.true.)
  call sbprms%CreateLogicalOption('lifting_line_data', 'Output lifting line data&
                              & alongside sectional loads','F')
  call sbprms%CreateLogicalOption('vortex_lattice_data', 'Output of corrected vortex lattice data&
                              & alongside sectional loads','F')
  ! chordwise loads
  call sbprms%CreateIntOption('n_station','Number of stations where the loads are extracted' &
                              , multiple=.true.)
  call sbprms%CreateRealArrayOption('span_station','Spanwise coordinates in the component reference frame&
                              & where they are extracted', multiple=.true.)
  
  
  ! sectional loads: box -----
  call sbprms%CreateSubOption('box_sect','Definition of the box for sectional loads', &
                            bxprms)
  call bxprms%CreateRealArrayOption('ref_node','reference node to build the box')
  call bxprms%CreateRealArrayOption('face_vec','vector identifying the direction of the base side &
                            &of the sections')
  call bxprms%CreateRealArrayOption('face_bas','dimension along faceVec of the first and last sections')
  call bxprms%CreateRealArrayOption('face_hei','dimension orthogonal to faceVec of the first and last sections')
  call bxprms%CreateRealArrayOption('span_vec','vector defining the out-of plane direction of the box')
  call bxprms%CreateRealOption('span_len','dimension along the spanVec direction')
  call bxprms%CreateIntOption('num_sect','number of sections')
  call bxprms%CreateLogicalOption('reshape_box','logical input to reshape the box if &
                            &it is "too large"')
  call sbprms%CreateRealArrayOption('axis_mom','axis for the computation of the moment. Perpendicular to sections')
  
  sbprms=>null()

end subroutine create_param_post


subroutine save_sim_param(this, loc)
  class(t_sim_param) :: this
  integer(h5loc), intent(in) :: loc

  call write_hdf5_attr(this%t0, 't0', loc)
  call write_hdf5_attr(this%dt, 'dt', loc)
  call write_hdf5_attr(this%tend, 'tend', loc)
  call write_hdf5_attr(this%output_detailed_geo, 'output_detailed_geo', loc)
  call write_hdf5_attr(this%P_inf, 'P_inf', loc)
  call write_hdf5_attr(this%rho_inf, 'rho_inf', loc)
  call write_hdf5_attr(this%u_inf, 'u_inf', loc)
  call write_hdf5_attr(this%u_ref, 'u_ref', loc)
  call write_hdf5_attr(this%a_inf, 'a_inf', loc)
  call write_hdf5_attr(this%mu_inf, 'mu_inf', loc)
  call write_hdf5_attr(this%first_panel_scaling, 'first_panel_scaling', loc)
  call write_hdf5_attr(this%min_vel_at_te, 'min_vel_at_te', loc)
  call write_hdf5_attr(this%rigid_wake, 'rigid_wake', loc)
  if(this%rigid_wake) &
    call write_hdf5_attr(this%rigid_wake_vel, 'rigid_wake_vel', loc)
  call write_hdf5_attr(this%n_wake_panels, 'n_wake_panels', loc)
  call write_hdf5_attr(this%n_wake_particles, 'n_wake_particles', loc)
  call write_hdf5_attr(this%particles_box_min, 'particles_box_min', loc)
  call write_hdf5_attr(this%particles_box_max, 'particles_box_max', loc)
  call write_hdf5_attr(this%join_te, 'join_te', loc)
  if(this%join_te) &
    call write_hdf5_attr(this%join_te_factor, 'join_te_factor', loc)


  call write_hdf5_attr(this%FarFieldRatioDoublet, 'FarFieldRatioDoublet', loc)
  call write_hdf5_attr(this%FarFieldRatioSource, 'FarFieldRatioSource', loc)
  call write_hdf5_attr(this%DoubletThreshold, 'DoubletThreshold', loc)
  call write_hdf5_attr(this%RankineRad, 'RankineRad', loc)
  call write_hdf5_attr(this%VortexRad, 'VortexRad', loc)
  call write_hdf5_attr(this%KVortexRad, 'KVortexRad', loc)
  call write_hdf5_attr(this%CutoffRad, 'CutoffRad', loc)
  call write_hdf5_attr(this%use_vs, 'Vortstretch', loc)
  if(this%use_vs) then
    call write_hdf5_attr(this%vs_elems, 'VortstretchFromElems', loc)
    call write_hdf5_attr(this%use_divfilt, 'DivergenceFiltering', loc)
    call write_hdf5_attr(1.0_wp/this%filt_eta*this%dt, 'FilterTimescale', loc)
  endif
  call write_hdf5_attr(this%use_vd, 'vortdiff', loc)
  call write_hdf5_attr(this%use_tv, 'turbvort', loc)
  call write_hdf5_attr(this%use_pa, 'PenetrationAvoidance', loc)
  if(this%use_pa) then
    call write_hdf5_attr(this%pa_rad_mult, 'PenetrationAvoidanceCheckRadius', loc)
    call write_hdf5_attr(this%pa_elrad_mult, 'PenetrationAvoidanceElementRadius', loc)
  endif
  call write_hdf5_attr(this%use_ve, 'ViscosityEffects', loc)
  call write_hdf5_attr(this%use_fmm, 'use_fmm', loc)
  if(this%use_fmm) then
    call write_hdf5_attr(this%use_fmm_pan, 'use_fmm_panels', loc)
    call write_hdf5_attr(this%BoxLength, 'BoxLength', loc)
    call write_hdf5_attr(this%Nbox, 'Nbox', loc)
    call write_hdf5_attr(this%OctreeOrigin, 'OctreeOrigin', loc)
    call write_hdf5_attr(this%NOctreeLevels, 'NOctreeLevels', loc)
    call write_hdf5_attr(this%MinOctreePart, 'MinOctreePart', loc)
    call write_hdf5_attr(this%MultipoleDegree, 'MultipoleDegree', loc)
    call write_hdf5_attr(this%use_dyn_layers, 'use_dyn_layers', loc)
    if(this%use_dyn_layers) then
      call write_hdf5_attr(this%NMaxOctreeLevels, 'NMaxOctreeLevels', loc)
      call write_hdf5_attr(this%LeavesTimeRatio, 'LeavesTimeRatio', loc)
    endif
    call write_hdf5_attr(this%use_pr, 'ParticlesRedistribution', loc)
    if(this%use_pr) then
      call write_hdf5_attr(this%lvl_solid, 'OctreeLevelSolid', loc)
      call write_hdf5_attr(this%part_redist_ratio, &
                                          'ParticlesRedistributionRatio', loc)
    endif
  endif
  call write_hdf5_attr(this%debug_level, 'debug_level', loc)
  call write_hdf5_attr(this%dt_out, 'dt_out', loc)
  call write_hdf5_attr(this%basename, 'basename', loc)
  call write_hdf5_attr(this%GeometryFile, 'GeometryFile', loc)
  call write_hdf5_attr(this%ReferenceFile, 'ReferenceFile', loc)
  call write_hdf5_attr(this%restart_from_file, 'restart_from_file', loc)
  if(this%restart_from_file) then
    call write_hdf5_attr(this%restart_file, 'restart_file', loc)
    call write_hdf5_attr(this%reset_time, 'reset_time', loc)
  endif
  call write_hdf5_attr(this%hcas, 'HCAS', loc)
  if(this%hcas) then
    call write_hdf5_attr(this%hcas_time, 'HCAS_time', loc)
    call write_hdf5_attr(this%hcas_vel, 'HCAS_velocity', loc)
  endif

end subroutine save_sim_param

end module mod_sim_param
