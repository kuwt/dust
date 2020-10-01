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
    !! Copyright (C) 2018-2020 Davide   Montagnani, 
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
    !!          Federico Fonte             <federico.fonte@outlook.com>
    !!          Davide Montagnani       <davide.montagnani@gmail.com>
    !!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
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

    !use mod_aero_elements, only: &
    !  c_elem, t_elem_p !, t_vp
    use mod_aeroel, only: &
      c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
      t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

    use mod_doublet, only: &
      initialize_doublet

    use mod_surfpan, only: &
      t_surfpan , initialize_surfpan

    use mod_vortlatt, only: &
      t_vortlatt

    use mod_liftlin, only: &
     update_liftlin, t_liftlin_p, &
     build_ll_kernel, &
     solve_liftlin, solve_liftlin_piszkin
    !solve_liftlin_optim , &
    !solve_liftlin_newton, &
    !solve_liftlin_optim_regul

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

#if USE_PRECICE
      use mod_precice, only: &
          t_precice
#endif

    implicit none

    !run-id
    integer :: run_id(10)

    !Input
    !> Main parameters parser
    type(t_parse) :: prms
    character(len=*), parameter :: input_file_name_def = 'dust.in'
    character(len=max_char_len) :: input_file_name
    character(len=max_char_len) :: target_file
    character(len=extended_char_len) :: message

    !Simulation parameters
    !type(t_sim_param) :: sim_param
    ! Asymptotic conditions

    !Time parameters
    real(wp) :: time
    integer  :: it, nstep, nout
    real(wp) :: t_last_out, t_last_debug_out
    real(wp) :: time_no_out, time_no_out_debug
    logical  :: time_2_out, time_2_debug_out
    logical  :: output_start
    real(wp) :: dt_debug_out

    !Main variables
    !> All the implicit elements, sorted first static then moving
    type(t_impl_elem_p), allocatable :: elems(:)
    !> All the explicit elements
    type(t_expl_elem_p), allocatable :: elems_expl(:)
    !> Only the lifting line elements
    !type(t_expl_elem_p), allocatable :: elems_ll(:)
    type(t_liftlin_p), allocatable :: elems_ll(:)
    !> Only the actuator disk elements
    type(t_expl_elem_p), allocatable :: elems_ad(:)
    !> All the elements (panels+ll)
    type(t_pot_elem_p), allocatable :: elems_tot(:)
    !> Geometry
    type(t_geo) :: geo
    !> Trailing edge
    type(t_tedge) :: te
    !> Airfoil table data
    type(t_aero_tab), allocatable :: airfoil_data(:)
    !> Linear system
    type(t_linsys) :: linsys
    !> Wake 
    type(t_wake) :: wake

    !> Timing vars
    real(t_realtime) :: t1 , t0, t00, t11, t22

    !I/O prefixes
    character(len=max_char_len) :: frmt
    character(len=max_char_len) :: basename_debug

    real(wp) , allocatable :: res_old(:)
    real(wp) , allocatable :: surf_vel_SurfPan_old(:,:)
    real(wp) , allocatable ::      nor_SurfPan_old(:,:)

    real(wp) , allocatable :: al_kernel(:,:) , al_v(:)

    integer :: i_el , i, j, sel

    !octree parameters
    type(t_octree) :: octree

#if USE_PRECICE
    !> --- PreCICE data ---------------------------------------------
    type(t_precice) :: precice
    integer :: bool
    logical :: precice_convergence
    !> --- PreCICE --------------------------------------------------
#endif

    !> debug --------
    real(wp) :: sum_force(3)
    !> debug --------


#if USE_PRECICE
    !> --- Initialize PreCICE ---------------------------------------
    call precice % initialize()
    !> --- Initialize PreCICE: done ---------------------------------
#endif

    call printout(nl//'>>>>>> DUST beginning >>>>>>'//nl)
    t00 = dust_time()

    call get_run_id(run_id)

    !------ Modules initialization ------
    call initialize_hdf5()

    !------ Input reading ------
    if(command_argument_count().gt.0) then
      call get_command_argument(1,value=input_file_name)
    else
      input_file_name = input_file_name_def
    endif

    ! define the parameters to be read
    ! time:
    call prms%CreateRealOption( 'tstart', "Starting time")
    call prms%CreateRealOption( 'tend',   "Ending time")
    call prms%CreateRealOption( 'dt',     "time step")
    call prms%CreateIntOption ( 'timesteps', "number of timesteps")
    call prms%CreateRealOption( 'dt_out', "output time interval")
    call prms%CreateRealOption( 'dt_debug_out', "debug output time interval")
    call prms%CreateIntOption ( 'ndt_update_wake', 'n. dt between two wake updates', '1')

    ! input:
    call prms%CreateStringOption('GeometryFile','Main geometry definition file')
    call prms%CreateStringOption('ReferenceFile','Reference frames file','no_set')

    ! output:
    call prms%CreateStringOption('basename','oputput basename','./')
    call prms%CreateStringOption('basename_debug','oputput basename for debug', &
                                                                             './')
    call prms%CreateLogicalOption( 'output_start', "output values at starting &
                                                              & iteration", 'F')
    call prms%CreateLogicalOption( 'output_detailed_geo', "output at each &
                        &timestep the detailed geometry in the results file", 'F')
    call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', &
                                                                             '1')

    ! restart
    call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
    call prms%CreateStringOption('restart_file','restart file name')
    call prms%CreateLogicalOption('reset_time','reset the time from previous &
                                   &execution?','F')

    ! parameters:
    ! --- reference conditions ---
    call prms%CreateRealArrayOption( 'u_inf', "free stream velocity", &
                                                           '(/1.0, 0.0, 0.0/)')
    call prms%CreateRealOption( 'u_ref', "reference velocity")             !
    call prms%CreateRealOption( 'P_inf', "free stream pressure", '0.0')    ! 
    call prms%CreateRealOption( 'rho_inf', "free stream density", '1.0')   !
    call prms%CreateRealOption( 'a_inf', "Speed of sound", '340.0')        ! m/s   for dimensional sim
    call prms%CreateRealOption( 'mu_inf', "Dynamic viscosity", '0.000018') ! kg/ms

    ! --- wake ---
    call prms%CreateIntOption('n_wake_panels', 'number of wake panels','4')
    call prms%CreateIntOption('n_wake_particles', 'number of wake particles', &
                                                                      '10000')
    call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
         &particles bounding box', '(/-10.0, -10.0, -10.0/)')
    call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
         &particles bounding box', '(/10.0, 10.0, 10.0/)')

    call prms%CreateRealOption( 'ImplicitPanelScale', &
                        "Scaling of the first implicit wake panel", '0.3')
    call prms%CreateRealOption( 'ImplicitPanelMinVel', &
                        "Minimum velocity at the trailing edge", '1.0e-8')
    call prms%CreateLogicalOption('rigid_wake','rigid wake?','F')
    call prms%CreateRealArrayOption( 'rigid_wake_vel', "rigid wake velocity" )
    call prms%CreateLogicalOption('join_te','join trailing edge','F')
    call prms%CreateRealOption( 'join_te_factor', "join the trailing edges &
                                 &when closer than factor*te element size",'1.0' )

    ! --- regularisation ---
    call prms%CreateRealOption( 'FarFieldRatioDoublet', &
          "Multiplier for far field threshold computation on doublet", '10.0')
    call prms%CreateRealOption( 'FarFieldRatioSource', &
          "Multiplier for far field threshold computation on sources", '10.0')
    call prms%CreateRealOption( 'DoubletThreshold', &
          "Thresold for considering the point in plane in doublets", '1.0e-6')
    call prms%CreateRealOption( 'RankineRad', &
          "Radius of Rankine correction for vortex induction near core", '0.1')
    call prms%CreateRealOption( 'VortexRad', &
          "Radius of vortex core, for particles", '0.1')
    call prms%CreateRealOption( 'CutoffRad', &
          "Radius of complete cutoff  for vortex induction near core", '0.001')

    ! --- lifting line elements ---
    call prms%CreateStringOption('LLsolver','Solver for the LL elements', &
            'GammaMethod')
    call prms%CreateLogicalOption('LLReynoldsCorrections', &
           'Use Reynolds corrections for the .c81 tables?', 'F')
    call prms%CreateRealOption('LLReynoldsCorrectionsNfact', &
           'Exponent in (Re/Re_T)^n correction', '0.2')
    call prms%CreateIntOption('LLmaxIter', & 
           'Maximum number of iteration in LL algorithm', '100')
    call prms%CreateRealOption('LLtol', 'Tolerance for the relative error in &
                               &fixed point iteration for LL','1.0e-6' )
    call prms%CreateRealOption('LLdamp', 'Damping param in fixed point iteration &
                               &for LL used to avoid oscillations','25.0')
    call prms%CreateLogicalOption('LLstallRegularisation', & 
           'Avoid "unphysical" separations in inner sections of LL?','T')
    call prms%CreateIntOption('LLstallRegularisationNelems', &
           'Number of "unphysical" separations to be removed', '1' )
    call prms%CreateIntOption('LLstallRegularisationNiters', &
           'Number of timesteps between two regularisations', '1' )
    call prms%CreateRealOption('LLstallRegularisationAlphaStall', &
           'Stall angle used as threshold for regularisation [deg]', '15.0' )
    call prms%CreateRealOption('LLartificialViscosity', &
           'Constant artificial viscosity for regularisation', '0.0' )
    call prms%CreateLogicalOption('LLartificialViscosityAdaptive', &
           'Adaptive artificial viscosity algorithm', 'F' )
    call prms%CreateRealOption('LLartificialViscosityAdaptive_Alpha', &
           'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
    call prms%CreateRealOption('LLartificialViscosityAdaptive_dAlpha', &
           'Adaptive Artificial Viscosity algorithm, reference AOA [deg]')
    call prms%CreateLogicalOption('LLloadsAVL', & 
           'Use AVL expression for inviscid load computation','T')

    !== Octree and multipole data == 
    call prms%CreateLogicalOption('FMM','Employ fast multipole method?','T')
    call prms%CreateLogicalOption('FMMPanels','Employ fast multipole method &
                                  &also for panels?','F')
    call prms%CreateRealOption('BoxLength','length of the octree box')
    call prms%CreateIntArrayOption('NBox','number of boxes in each direction')
    call prms%CreateRealArrayOption( 'OctreeOrigin', "rigid wake velocity" )
    call prms%CreateIntOption('NOctreeLevels','number of octree levels')
    call prms%CreateIntOption('MinOctreePart','minimum number of octree &
                                                                 &particles')
    call prms%CreateIntOption('MultipoleDegree','multipole expansion degree')
    call prms%CreateLogicalOption('DynLayers','Use dynamic layers','F')
    call prms%CreateIntOption('NMaxOctreeLevels','maximum number &
                                                              &of octree levels')
    call prms%CreateRealOption('LeavesTimeRatio','Ratio that triggers the &
                                              &increase of the number of levels')

    ! models options
    call prms%CreateLogicalOption('Vortstretch','Employ vortex stretching','T')
    call prms%CreateLogicalOption('VortstretchFromElems','Employ vortex stretching&
                                         & from geometry elements','F')
    call prms%CreateLogicalOption('Diffusion','Employ vorticity diffusion','T')
    call prms%CreateLogicalOption('TurbulentViscosity','Employ turbulent &
                                   &viscosity','F')
    call prms%CreateLogicalOption('PenetrationAvoidance','Employ penetration &
                                                                  & avoidance','F')
    call prms%CreateRealOption('PenetrationAvoidanceCheckRadius', &
         'Check radius for penetration avoidance','5.0')
    call prms%CreateRealOption('PenetrationAvoidanceElementRadius', &
         'Element impact radius for penetration avoidance','1.5')
    call prms%CreateLogicalOption('ViscosityEffects','Simulate viscosity &
                                                                  & effects','F')
    call prms%CreateLogicalOption('ParticlesRedistribution','Employ particles &
                                                            &redistribution','F')
    call prms%CreateIntOption('OctreeLevelSolid','Level at which the panels &
                              & are considered for particles redistribution')
    call prms%CreateRealOption('ParticlesRedistributionRatio','How many times &
             &a particle need to be smaller than the average of the cell to be&
             & eliminated','3.0')
    ! HCAS
    call prms%CreateLogicalOption('HCAS','Hover Convergence Augmentation System',&
                                   'F')
    call prms%CreateRealOption('HCAS_time','HCAS deployment time')
    call prms%CreateRealArrayOption('HCAS_velocity','HCAS velocity')

    ! get the parameters and print them out
    call printout(nl//'====== Input parameters: ======')
    call check_file_exists(input_file_name,'dust main')
    call prms%read_options(input_file_name, printout_val=.true.)

    ! initialize all the parameters
    nout = 0 !reset the numbering for output files
    output_start = getlogical(prms, 'output_start')
    call init_sim_param(sim_param, prms, nout, output_start)

    ! remaining parameters
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

    !-- Parameters Initializations --
    call initialize_doublet()
    call initialize_vortline()
    call initialize_vortpart()
    call initialize_surfpan()


    ! Check that tend .gt. tinit
    if ( sim_param%tend .le. sim_param%t0 ) then
        write(message,*) 'The end time of the simulation',sim_param%tend,'is&
                         & lower or equal to the start time',sim_param%t0,'.&
                         & Remembver that when restarting without resetting the&
                         & time the start time is taken from the restart result!'
        call error('dust','',message)
    end if


    !Printout the parameters
    if (sim_param%debug_level .ge. 3) then
      write(message,*) 'Initial time tstart: ', sim_param%t0; call printout(message)
      write(message,*) 'Final time tend:     ', sim_param%tend; call printout(message)
      write(message,*) 'Time step dt:        ', sim_param%dt; call printout(message)
      write(message,*) 'Output interval:     ', sim_param%dt_out; call printout(message)
      write(message,*) 'Output first step:   ', output_start; call printout(message)
      write(message,*) 'Debug level:', sim_param%debug_level; call printout(message)
      write(message,*) 'Free stream velocity:', sim_param%u_inf; call printout(message)
      write(message,*) 'Maximum wake panels:', sim_param%n_wake_panels; call printout(message)
      write(message,*) 'Results basename: ', trim(sim_param%basename); call printout(message)
      write(message,*) 'Debug basename: ', trim(basename_debug); call printout(message)
    endif


    !---- Check that the basenames are valid ----
      call check_basename(trim(sim_param%basename),'dust main')
      if(sim_param%debug_level.ge.10) &
        call check_basename(trim(basename_debug),'dust main')

    !---- Simulation parameters ----
    !!nstep = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1 
    !!       !(+1 for the zero time step)
    !!sim_param%n_timesteps = nstep
    nstep = sim_param%n_timesteps
    allocate(sim_param%time_vec(sim_param%n_timesteps))
    sim_param%time_vec = (/ ( sim_param%t0 + &
             real(i-1,wp)*sim_param%dt, i=1,sim_param%n_timesteps ) /)

    !------ Geometry creation ------
    call printout(nl//'====== Geometry Creation ======')
    t0 = dust_time()
    target_file = trim(sim_param%basename)//'_geo.h5'
    call create_geometry(sim_param%GeometryFile, sim_param%ReferenceFile, &
                         input_file_name, geo, te, elems, elems_expl, elems_ad, &
                         elems_ll, elems_tot, airfoil_data, target_file, run_id)

    t1 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Created geometry in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    if(sim_param%debug_level .ge. 15) &
          call debug_printout_geometry_minimal(elems, geo, basename_debug, 0)
    if(sim_param%debug_level .ge. 15) &
          call debug_ll_printout_geometry(elems_ll, geo, basename_debug, 0)

    !TODO: check whether to move these calls before, and precisely what they do
    if(sim_param%debug_level .ge. 7) call ignoredParameters(prms)
    call finalizeParameters(prms)

    !> --- Initialize PreCICE mesh and fields -----------------------
#if USE_PRECICE
      call precice % initialize_mesh( geo )
      call precice % initialize_fields()

!     !> check ---
!     do i = 1, size(precice%fields)
!       write(*,*) trim(precice%fields(i)%fname), precice%fields(i)%fid
!     end do
!     write(*,*) ' Using PreCICE: stop in dust.f90, after the initialization '
!     write(*,*) ' of the fields and the mesh, used for coupling, l.490.'; stop

#endif 
    !> --- Initialize PreCICE mesh and fields: done -----------------

    !------ Initialization ------

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
    call initialize_linsys(linsys, geo, elems, elems_expl, &
                           wake ) ! sim_param%u_inf)
    !call initialize_pressure_sys(linsys, geo, elems)

    t1 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Initialized linear system in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    ! Restart --------------
    !------ Reloading ------
    if (sim_param%restart_from_file) then
     call load_solution(sim_param%restart_file, geo%components, geo%refs)
     call load_wake(sim_param%restart_file, wake, elems_tot)

    ! Moved to mod_geo.f90/create_geometry()
    !! Update the initial relative position of the ref.sys.
    !call update_relative_initial_conditions(restart_file, sim_param%ReferenceFile, geo%refs) 

    else ! Set to zero the intensity of all the singularities

      do i_el = 1 , size(elems)      ! implicit elements (vr, sp)
         elems(i_el)%p%mag = 0.0_wp
      end do 
      do i_el = 1 , size(elems_expl) ! explicit elements (ll, ad)
         elems_expl(i_el)%p%mag = 0.0_wp
      end do 

    endif

    t22 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') nl//'------ Completed all preliminary operations &
                                 &in: ' , t22 - t00,' s.'
      call printout(message)
    endif

    ! === build kernel for regularisation (averaging) process for LL === 
    if ( ( size(elems_ll) .gt. 0 ) ) then
      allocate(al_v(size(elems_ll))) ; al_v = 0.0_wp
      call build_ll_kernel( elems_ll , al_v , al_kernel )
      deallocate(al_v)
    end if


    !> --- Initialize coupling --------------------------------------
#if USE_PRECICE
      call precicef_ongoing(   precice % is_ongoing)
      write(*,*) ' is coupling ongoing: ', precice % is_ongoing
      call precicef_initialize(precice % dt_precice) 
      call precicef_ongoing(   precice % is_ongoing)

      precice_convergence = .true.

      write(*,*) ' is coupling ongoing: ', precice % is_ongoing
      write(*,*) ' dt_precice         : ', precice % dt_precice
#endif
    !> --- Initialize PreCICE mesh and fields: done -----------------


    !====== Time Cycle ======
    call printout(nl//'////////// Performing Computations //////////')
    time = sim_param%t0
    t_last_out = time; t_last_debug_out = time
    time_no_out = 0.0_wp; time_no_out_debug = 0.0_wp

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


    t11 = dust_time()
    it = 0
#if USE_PRECICE
    it = 1
    do while ( ( it .lt. nstep ) .and. ( precice%is_ongoing .eq. 1 ) )
#else
    do while ( ( it .lt. nstep ) )
      it = it + 1
#endif

      if(sim_param%debug_level .ge. 1) then
        write(message,'(A,I5,A,I5,A,F9.4)') nl//'--> Step ',it,' of ', &
                                            nstep, ' simulation time: ', time
        call printout(message)
        t22 = dust_time()
        write(message,'(A,F9.3,A)') 'Elapsed wall time: ',t22-t00
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
    !   write(*,*) ' =================== call update_wake ==================== '
        call prepare_wake(wake, elems_tot, octree)
      end if
     
      call update_liftlin(elems_ll,linsys)
      call update_actdisk(elems_ad,linsys)

      !debug geometry printing
      if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
                call debug_printout_geometry(elems, geo, basename_debug, it)
      if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
                call debug_ll_printout_geometry(elems_ll, geo, basename_debug, it)

#if USE_PRECICE
      !> -------------------------------------------------------------
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

      ! check ---
      write(*,*) ' precice%mesh%nnodes : ', precice%mesh%nnodes
      ! check ---
      !> Read data from structural solver
      do i = 1, size(precice%fields)
        if ( trim(precice%fields(i)%fio) .eq. 'read' ) then
          if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
            ! check ---
            ! write(*,*) ' read scalar field: ', trim(precice%fields(i)%fname)
            ! check ---
            call precicef_read_bsdata( precice%fields(i)%fid, &
                                       precice%mesh%nnodes  , &
                                       precice%mesh%node_ids, &
                                       precice%fields(i)%fdata(1,:) )
          elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
            ! ! check ---
            ! write(*,*) ' read vector field: ', trim(precice%fields(i)%fname)
            ! write(*,*) '    mesh%nnodes  : ' , precice%mesh%nnodes  
            ! write(*,*) '    mesh%nodes_ids: ', precice%mesh%node_ids
            ! ! check ---
            call precicef_read_bvdata( precice%fields(i)%fid, &
                                       precice%mesh%nnodes  , &
                                       precice%mesh%node_ids, &
                                       precice%fields(i)%fdata )
          endif
          ! check ---
          write(*,*) ' precice%mesh%nnodes : ', precice%mesh%nnodes
          write(*,*) i, precice%fields(i)%fid, precice%fields(i)%fname
          do i_el = 1, size(precice%fields(i)%fdata,2)
            write(*,*) precice%fields(i)%fdata(:,i_el)
          end do
          write(*,*)
          ! check ---
        end if
      end do

      ! write(*,*) ' write something, before stop in dust.f90, l.700 '
      ! read(*,*)

      ! ! check ---
      ! stop
      ! ! check ---

      !> Update dust geometry ( elems and first wake panels )
      call precice % update_elems( geo, elems_tot )

      ! ! debug ---
      ! write(*,*) ' debug in dust.f90, l.715. geo%points:'
      ! do i_el = 1, size(geo%points,2)
      !   write(*,*) i_el, ' :   ', geo%points(:,i_el)
      ! end do
      ! write(*,*) ' stop in dust.f90, l.719. '; stop
      ! write(*,*) ' debug in dust.f90, l.715. elems_tot(1)%p%ver: '
      ! do i_el = 1, size(elems_tot(1)%p%ver,2)
      !   write(*,*) elems_tot(1)%p%i_ver(i_el), ': ', elems_tot(1)%p%ver(:,i_el)
      ! end do
      ! ! debug ---

      !> Update geo_data()
      do i_el = 1, size(elems_tot)
        call elems_tot(i_el)%p%calc_geo_data( &
                               geo%points(:,elems_tot(i_el)%p%i_ver) )
      end do

      !> Update near-field wake
      call precice % update_near_field_wake( geo, wake )

      !> Update dt
      ! *** to do *** : change the way of treating time integration.
      ! So far, only explicit fixed-dt integration is available.
      ! The law of motion of the rigid components is evaluated before
      ! the time loop starts
      ! dt = min ( dt_dust, dt_precice )

#else
#endif

      !Calculate the normal velocity derivative for the pressure equation
      call press_normvel_der(geo, elems, surf_vel_SurfPan_old)

      !------ Assemble the system ------
      !call prepare_wake(wake, geo, sim_param)
      t0 = dust_time()
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

      !------ Solve the pressure system ------
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
      ! compute time derivative of the result ( = i_vortex = -i_doublet ) -------
    !$omp parallel do private(i_el)
      do i_el = 1 , sel
        elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
      end do
    !$omp end parallel do
      ! res_old = linsys%res  ! moved to l.940 approx, when the dt is finalized

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

    !   call solve_liftlin_optim(elems_ll, elems_tot, elems , elems_ad , &
    !           (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it, &
    !           M_array )
    !   call solve_liftlin_optim_regul(elems_ll, elems_tot, elems , elems_ad , &
    !           (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it, &
    !           M_array )
    !   call solve_liftlin_newton(elems_ll, elems_tot, elems , elems_ad , &
    !           (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it, &
    !           M_array )
       end if
      end if

      !------ Compute loads -------
      ! Implicit elements: vortex rings and 3d-panels
      ! 2019-07-23: D.Isola suggested to implement AVL formula for VL elements
      ! so far, select type() to keep the old formulation for t_surfpan and 
      ! use AVL formula for t_vortlatt
    !$omp parallel do private(i_el)
      do i_el = 1 , sel

        ! ifort bugs workaround:
        ! apparently it is not possible to call polymorphic methods inside
        ! select cases for intel, need to call these for all elements and for the
        ! vortex lattices it is going to be a dummy empty function call

!           write(*,*) ' i_el : ', i_el
!           write(*,*) ' elems(i_el)%p%comp_id): ', elems(i_el)%p%comp_id
!           write(*,*) ' geo%components()%ref_id: ', &
!                        geo%components(elems(i_el)%p%comp_id)%ref_id
           
            call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
                    geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g)
            call elems(i_el)%p%compute_dforce()

      end do
    !$omp end parallel do

      ! ifort bugs workaround:
      ! since even if the following calls looks thread safe, they mess up with 
      ! ifort and parallel runs, so the cycle is executed another time just for the
      ! vortex lattices
      if ( geo%nVortLatt .gt. 0) then
        do i_el = 1 , sel
          select type( el => elems(i_el)%p ) 
            class is(t_vortlatt)
              ! compute vel at 1/4 chord (some approx, see the comments in the fcn)
              call el%get_vel_ctr_pt( elems_tot, (/ wake%pan_p, wake%rin_p/) )
              ! compute dforce using AVL formula
              call el%compute_dforce_jukowski()
              ! update the pressure field, p = df.n / area
              !el%pres = sum( el%dforce * el%nor ) / el%area
              !ifort bug workaround
              elems(i_el)%p%pres = sum( elems(i_el)%p%dforce * elems(i_el)%p%nor )&
                                   / elems(i_el)%p%area

              ! ! debug ---
              ! write(*,'(I4,A,3F10.5,A,3F10.5,A,3F10.5,A,3F10.5,A,F10.5)') &
              !                        i_el, '     ', el%vel_ctr_pt, &
              !                              '     ', el%dforce, &
              !                              '     ', el%nor,    &
              !                              '     ', el%dn_dt,  &
              !                              '     ', el%pres
              ! ! debug ---

          end select
        end do
      endif

      ! Explicit elements:
      ! - liftlin: _pres and _dforce computed in solve_liftin()
      ! - actdisk: avg delta_pressure and force computed here,
      !            to include thier effects in postpro (e.g. integral loads)
    !$omp parallel do private(i_el)
      do i_el = 1 , size(elems_ad)
    !   call elems_ad(i_el)%p%compute_pres(sim_param)
        call elems_ad(i_el)%p%compute_pres( &     ! update surf_vel field too
                 geo%refs( geo%components(elems_ad(i_el)%p%comp_id)%ref_id )%R_g)
        call elems_ad(i_el)%p%compute_dforce()
      end do
    !$omp end parallel do

#if USE_PRECICE

      write(*,*) ' ------------------------------------------------------------------- '
      write(*,*) ' debug in dust.f90, l.890: i, pres, dforce '
      sum_force = 0.0_wp
      do i = 1, size(elems_tot)
        sum_force = sum_force + elems_tot(i)%p%dforce
        write(*,*) i, ' : ', elems_tot(i)%p%pres, '      ', elems_tot(i)%p%dforce
      end do
      write(*,*) '                                                 ', sum_force
      write(*,*) ' ------------------------------------------------------------------- '
      write(*,*)

!     write(*,*) ' debug in dust.f90, l.890 '
!     do i = 1, 10
!     ! !> check mag, pres
!     ! write(*,*) i, ':  ', elems_ll(i)%p%mag , elems_ll(i+10)%p%mag , elems_ll(i+20)%p%mag, &
!     !               '   ', elems_ll(i)%p%pres, elems_ll(i+10)%p%pres, elems_ll(i+20)%p%pres
!     ! !> check nor, dforce ---> issues in the unsteady force contributions
!     ! write(*,*) i   , ':  ', elems_ll(i   )%p%nor , elems_ll(i   )%p%dforce, &
!     !                    sum( elems_ll(i   )%p%nor * elems_ll(i   )%p%dforce )
!     ! write(*,*) i+10, ':  ', elems_ll(i+10)%p%nor , elems_ll(i+10)%p%dforce, &
!     !                    sum( elems_ll(i+10)%p%nor * elems_ll(i+10)%p%dforce )
!     ! write(*,*) i+20, ':  ', elems_ll(i+20)%p%nor , elems_ll(i+20)%p%dforce, &
!     !                    sum( elems_ll(i+20)%p%nor * elems_ll(i+20)%p%dforce )
!     ! write(*,*)
!     ! write(*,*) i, ':  ', elems_ll(i   )%p%dGamma_dt , &
!     !                      elems_ll(i+10)%p%dGamma_dt , &
!     !                      elems_ll(i+20)%p%dGamma_dt
!     ! !> check nor, dforce ---> issues in the unsteady force contributions
!     ! write(*,*) i   , ':  ', elems_ll(i   )%p%dn_dt
!     ! write(*,*) i+10, ':  ', elems_ll(i+10)%p%dn_dt
!     ! write(*,*) i+20, ':  ', elems_ll(i+20)%p%dn_dt
!     ! write(*,*)
!     end do

      !> Update force and moments to be passed to the structural solver
      call precice % update_force( geo, elems_tot )

      call precicef_ongoing( precice%is_ongoing )
      if ( precice%is_ongoing .eq. 1 ) then
        call precicef_advance( precice%dt_precice )
        write(*,*) ' ++++ dt_precice: ', precice%dt_precice
      end if

      !> Write force and moments to structural solver
      do i = 1, size(precice%fields)
        if ( trim(precice%fields(i)%fio) .eq. 'write' ) then

          if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
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

          ! check ---
          write(*,*) ' precice%mesh%nnodes : ', precice%mesh%nnodes
          write(*,*) i, precice%fields(i)%fid, precice%fields(i)%fname
          do i_el = 1, size(precice%fields(i)%fdata,2)
            write(*,*) precice%fields(i)%fdata(:,i_el)
          end do
          write(*,*)
          ! check ---

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

      !----- Print the results -----
      if(time_2_out)  then
        nout = nout+1
        call save_status(geo, wake, nout, time, run_id)
      endif

      !------ Viscous Effects and Flow Separations ------
      ! some computation of surface quantities and vorticity to be released.
      ! Free vortices will be introduced in prepare_wake(), some lines below
      if(sim_param%use_ve) call viscosity_effects( geo , elems , te )

      !------ Treat the wake ------
      ! (this needs to be done after output, in practice the update is for the
      !  next iteration)
      t0 = dust_time()
      if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
    !   write(*,*) ' =================== call update_wake ==================== '
        call update_wake(wake, elems_tot, octree)
      end if
      t1 = dust_time()
      if(sim_param%debug_level .ge. 1) then
        write(message,'(A,F9.3,A)') 'Updated wake in: ' , t1 - t0,' s.'
        call printout(message)
      endif

      ! Pressure integral equation  +++++++
      ! save old velocity on the surfpan (before updating the geom, few lines below)
      do i_el = 1 , geo%nSurfPan
        select type ( el => elems(i_el)%p ) ; class is ( t_surfpan )
          surf_vel_SurfPan_old( geo%idSurfPanG2L(i_el) , : ) = el%ub   ! el%surf_vel
               nor_SurfPan_old( geo%idSurfPanG2L(i_el) , : ) = el%nor  ! el%surf_vel
        end select
      end do
      ! Pressure integral equation ++++++++

      if(it .lt. nstep) then 
        !time = min(sim_param%tend, time+sim_param%dt)
        time = min(sim_param%tend, sim_param%time_vec(it+1))
        !> Update geometry
        call update_geometry(geo, time, .false.)
        if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
    !     write(*,*) ' =================== call complete_wake ==================== '
          call complete_wake(wake, geo, elems_tot)
        end if
      endif

      !> Save old solution (at the previous dt) of the linear system
      res_old = linsys%res

      !> Update nor_old (moved from geo/mod_geo.f90/update_geometry(), l.2220 approx
    !$omp parallel do private(i_el)
      do i_el = 1 , size(elems_tot)
        elems_tot(i_el)%p%nor_old = elems_tot(i_el)%p%nor
      end do
    !$omp end parallel do


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

    enddo
    call printout(nl//'\\\\\\\\\\  Computations Finished \\\\\\\\\\')


    deallocate( res_old )
    !===== End Time Cycle ======

    !> --- Finalize PreCICE -----------------------------------------
#if USE_PRECICE
      call precicef_finalize()
#endif
    !> --- Finalize PreCICE: done -----------------------------------

    !------ Cleanup ------
    !call destroy_wake_panels(wake_panels)
    call destroy_wake(wake)
    call destroy_octree(octree)
    call destroy_linsys(linsys)
    call destroy_elements(geo)
    call destroy_geometry(geo, elems_tot)

    call destroy_hdf5()

    t22 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Completed all computations in ',t22-t00,' s'
      call printout(message)
    endif
    call printout(nl//'<<<<<< DUST end <<<<<<'//nl)



    !----------------------------------------------------------------------
    contains
    !----------------------------------------------------------------------

    subroutine get_run_id (run_id)
     integer, intent(out) :: run_id(10)

     real(wp) :: randr
     integer  :: maxi, randi

      !First 8 values are the date and time
      call date_and_time(VALUES=run_id(1:8))

      !Last 3 values are 2 random integers
      maxi = huge(maxi)
      call random_number(randr)
      randi = int(randr*real(maxi,wp))
      run_id(9) = randi

      call random_number(randr)
      randi = int(randr*real(maxi,wp))
      run_id(10) = randi

    end subroutine

    !----------------------------------------------------------------------

    !> Initialize all the parameters reading them from the the input file
subroutine init_sim_param(sim_param, prms, nout, output_start)
 class(t_sim_param) :: sim_param
 type(t_parse) :: prms
 integer, intent(inout) :: nout
 logical, intent(inout) :: output_start
 
  !timing
  sim_param%t0     = getreal(prms, 'tstart')
  sim_param%tend   = getreal(prms, 'tend')
  sim_param%dt_out = getreal(prms,'dt_out')
  sim_param%debug_level = getint(prms, 'debug_level')
  sim_param%output_detailed_geo = getlogical(prms, 'output_detailed_geo')
  sim_param%ndt_update_wake = getint(prms, 'ndt_update_wake')

  !Reference values
  sim_param%P_inf = getreal(prms,'P_inf')
  sim_param%rho_inf  = getreal(prms,'rho_inf')
  sim_param%u_inf = getrealarray(prms, 'u_inf', 3)
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
  sim_param%a_inf  = getreal(prms,'a_inf')
  sim_param%mu_inf = getreal(prms,'mu_inf')
  sim_param%nu_inf = sim_param%mu_inf/sim_param%rho_inf
  !Wake parameters
  sim_param%n_wake_panels = getint(prms, 'n_wake_panels')
  if(sim_param%n_wake_panels .lt. 1) then
    sim_param%n_wake_panels = 1
    call warning('dust','dust','imposed a number of wake panels rows &
                &LOWER THAN 1. At least one row of panels is mandatory, &
                &the simulation will proceed with "n_wake_panels = 1"')
  endif
  sim_param%n_wake_particles = getint(prms, 'n_wake_particles')
  sim_param%particles_box_min = getrealarray(prms, 'particles_box_min',3)
  sim_param%particles_box_max = getrealarray(prms, 'particles_box_max',3)
  sim_param%rigid_wake = getlogical(prms, 'rigid_wake')
  sim_param%rigid_wake_vel = sim_param%u_inf   ! initialisation
  if ( sim_param%rigid_wake ) then
    if ( countoption(prms,'rigid_wake_vel') .eq. 1 ) then
      sim_param%rigid_wake_vel = getrealarray(prms, 'rigid_wake_vel',3)
    else if ( countoption(prms,'rigid_wake_vel') .le. 0 ) then
      call warning('dust','dust','no rigid_wake_vel parameter set, &
           &with rigid_wake = T; rigid_wake_vel = u_inf')
      sim_param%rigid_wake_vel = sim_param%u_inf
    end if
  end if
  sim_param%join_te = getlogical(prms, 'join_te')
  if (sim_param%join_te) sim_param%join_te_factor=getreal(prms,'join_te_factor')
  !Names
  sim_param%basename = getstr(prms,'basename')
  sim_param%GeometryFile = getstr(prms,'GeometryFile')
  sim_param%ReferenceFile = getstr(prms,'ReferenceFile')
  !Method parameters
  sim_param%FarFieldRatioDoublet  = getreal(prms, 'FarFieldRatioDoublet')
  sim_param%FarFieldRatioSource  = getreal(prms, 'FarFieldRatioSource')
  sim_param%DoubletThreshold   = getreal(prms, 'DoubletThreshold')
  sim_param%RankineRad = getreal(prms, 'RankineRad')
  sim_param%VortexRad = getreal(prms, 'VortexRad')
  sim_param%CutoffRad  = getreal(prms, 'CutoffRad')
  sim_param%first_panel_scaling = getreal(prms,'ImplicitPanelScale')
  sim_param%min_vel_at_te  = getreal(prms,'ImplicitPanelMinVel')
  sim_param%use_vs = getlogical(prms, 'Vortstretch')
  sim_param%vs_elems = getlogical(prms, 'VortstretchFromElems')
  sim_param%use_vd = getlogical(prms, 'Diffusion')
  sim_param%use_tv = getlogical(prms, 'TurbulentViscosity')
  sim_param%use_pa = getlogical(prms, 'PenetrationAvoidance')
  if(sim_param%use_pa) then
    sim_param%pa_rad_mult = getreal(prms, 'PenetrationAvoidanceCheckRadius')
    sim_param%pa_elrad_mult = getreal(prms,'PenetrationAvoidanceElementRadius')
  endif
  sim_param%use_ve = getlogical(prms, 'ViscosityEffects')
  !Lifting line elements
  sim_param%llSolver                        = getstr(    prms, 'LLsolver')
  sim_param%llReynoldsCorrections           = getlogical(prms, 'LLreynoldsCorrections')
  sim_param%llReynoldsCorrectionsNfact      = getreal(   prms, 'LLreynoldsCorrectionsNfact')
  sim_param%llMaxIter                       = getint(    prms, 'LLmaxIter'            )
  sim_param%llTol                           = getreal(   prms, 'LLtol'                )
  sim_param%llDamp                          = getreal(   prms, 'LLdamp'               )
  sim_param%llStallRegularisation           = getlogical(prms, 'LLstallRegularisation')
  sim_param%llStallRegularisationNelems     = getint(    prms, 'LLstallRegularisationNelems')
  sim_param%llStallRegularisationNiters     = getint(    prms, 'LLstallRegularisationNiters')
  sim_param%llStallRegularisationAlphaStall = getreal(   prms, 'LLstallRegularisationAlphaStall')
  sim_param%llArtificialViscosity           = getreal(   prms, 'LLartificialViscosity')
  sim_param%llArtificialViscosityAdaptive   = getlogical(prms, 'LLartificialViscosityAdaptive')
  sim_param%llLoadsAVL                      = getlogical(prms, 'LLloadsAVL')
  !> check LL inputs
  if ( ( trim(sim_param%llSolver) .ne. 'GammaMethod' ) .and. &
       ( trim(sim_param%llSolver) .ne. 'AlphaMethod' ) ) then

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
  if ( trim(sim_param%llSolver) .eq. 'AlphaMethod' ) then
    if ( sim_param%llArtificialViscosityAdaptive) then
      if ( ( countoption(prms,'LLartificialViscosityAdaptive_Alpha')  .eq. 0 ) .or. &
           ( countoption(prms,'LLartificialViscosityAdaptive_dAlpha') .eq. 0 ) ) then

        call error('dust','init_sim_param','LLartificialViscosityAdaptive_Alpha or&
           & LLartificialViscosity_dAlpha not set as an input, while LLartificialViscosityAdaptive&
           & is set equal to T. Set these parameters [deg].')
  
      else
        sim_param%llArtificialViscosityAdaptive_Alpha = &
                  getreal(prms,'LLartificialViscosityAdaptive_Alpha')
        sim_param%llArtificialViscosityAdaptive_dAlpha= &
                  getreal(prms,'LLartificialViscosityAdaptive_dAlpha')
      end if 
    end if 
  end if
  write(*,*) ' sim_param%llSolver : ' , trim(sim_param%llSolver)
 
  !Octree and FMM parameters
  sim_param%use_fmm = getlogical(prms, 'FMM')
  if(sim_param%use_fmm) then
    sim_param%use_fmm_pan = getlogical(prms, 'FMMPanels')
    sim_param%BoxLength = getreal(prms, 'BoxLength')
    sim_param%NBox = getintarray(prms, 'NBox',3)
    sim_param%OctreeOrigin = getrealarray(prms, 'OctreeOrigin',3)
    sim_param%NOctreeLevels = getint(prms, 'NOctreeLevels')
    sim_param%MinOctreePart = getint(prms, 'MinOctreePart')
    sim_param%MultipoleDegree = getint(prms,'MultipoleDegree')
    
    sim_param%use_dyn_layers = getlogical(prms,'DynLayers')
    if(sim_param%use_dyn_layers) then
      sim_param%NMaxOctreeLevels = getint(prms, 'NMaxOctreeLevels')
      sim_param%LeavesTimeRatio = getreal(prms, 'LeavesTimeRatio')
    else
      sim_param%NMaxOctreeLevels = sim_param%NOctreeLevels
    endif
  
    sim_param%use_pr = getlogical(prms, 'ParticlesRedistribution')
    if(sim_param%use_pr) then
      sim_param%part_redist_ratio = getreal(prms,'ParticlesRedistributionRatio')
      if ( countoption(prms,'OctreeLevelSolid') .gt. 0 ) then
        sim_param%lvl_solid = getint(prms, 'OctreeLevelSolid')
      else
        sim_param%lvl_solid = max(sim_param%NOctreeLevels-2,1)
      endif
    endif
  else
    sim_param%use_fmm_pan = .false.
  endif

  !HCAS
  sim_param%hcas = getlogical(prms,'HCAS')
  if(sim_param%hcas) then
    sim_param%hcas_time = getreal(prms,'HCAS_time')
    sim_param%hcas_vel = getrealarray(prms,'HCAS_velocity',3)
  endif
  
  !Manage restart
  sim_param%restart_from_file = getlogical(prms,'restart_from_file')
  if (sim_param%restart_from_file) then
  
    sim_param%reset_time = getlogical(prms,'reset_time')
    sim_param%restart_file = getstr(prms,'restart_file')
    !Removing leading "./" if present to avoid issues when restarting
    if(sim_param%basename(1:2) .eq. './') sim_param%basename = sim_param%basename(3:)
    if(sim_param%restart_file(1:2) .eq. './') sim_param%restart_file = sim_param%restart_file(3:)
    call printout('RESTART: restarting from file: '//trim(sim_param%restart_file))
    sim_param%GeometryFile = sim_param%restart_file(1:len(trim(sim_param%restart_file))-11)&
                                                             &//'geo.h5'
  
    !restarting the same simulation, advance the numbers
    if(sim_param%restart_file(1:len(trim(sim_param%restart_file))-12).eq. &
                                                 trim(sim_param%basename)) then
    read(sim_param%restart_file(len(trim(sim_param%restart_file))-6:len(trim(sim_param%restart_file))-3),*) &
                                                                           nout 
      call printout('Identified restart from the same simulation, keeping the&
      & previous output numbering')
      !avoid rewriting the same timestep
      output_start = .false.
    endif
    if(.not. sim_param%reset_time) call load_time(sim_param%restart_file, sim_param%t0)
  
  endif


  !Check the number of timesteps
  if(CountOption(prms,'dt') .gt. 0) then
    if( CountOption(prms,'timesteps') .gt. 0) then
      !error
      call error('init_sim_param','dust','Both number of timesteps and dt are&
      & set, but only one of the two can be specified')
    else
      !get dt and compute number of timesteps
      sim_param%dt     = getreal(prms, 'dt')
      sim_param%n_timesteps = ceiling((sim_param%tend-sim_param%t0)/&
                                       sim_param%dt) + 1 
       !(+1 for the zero time step)
    endif
  else
    !get number of steps, compute dt
    sim_param%n_timesteps = getint(prms, 'timesteps')
    sim_param%dt =  (sim_param%tend-sim_param%t0)/&
                     real(sim_param%n_timesteps,wp)
    sim_param%n_timesteps = sim_param%n_timesteps + 1 
                            !add one for the first step
  endif
  
end subroutine init_sim_param

!----------------------------------------------------------------------

!DISCONTINUED: consider removing
!!subroutine copy_geo(sim_param, geo_file, run_id)
!! type(t_sim_param), intent(inout) :: sim_param
!! character(len=*), intent(inout)     :: geo_file
!! integer, intent(in)              :: run_id(10)
!!
!! character(len=max_char_len) :: target_file
!! integer :: estat, cstat
!! integer(h5loc) :: floc
!!
!!  !target file name: same as run basename with appendix
!!  target_file = trim(sim_param%basename)//'_geo.h5'
!!  
!!  if (trim(geo_file) .ne. trim(target_file)) then
!!  !Copy the geometry file
!!  call execute_command_line('cp '//trim(geo_file)//' '//trim(target_file), &
!!                                           exitstat=estat,cmdstat=cstat)
!!  if((cstat .ne. 0) .or. (estat .ne. 0)) &
!!    call error('dust','','System errors while trying to copy the geometry &
!!    &to the output path')
!!  endif
!!
!!
!!  !Attach the run_id to the file as an attribute
!!  call open_hdf5_file(trim(target_file), floc)
!!  call write_hdf5_attr(run_id, 'run_id', floc)
!!  call close_hdf5_file(floc)
!!
!!
!!  !Overwrite the geo file name, so that the copy is going to be
!!  !opened
!!  geo_file = trim(target_file)
!!
!!end subroutine copy_geo

!----------------------------------------------------------------------

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


!------------------------------------------------------------------------------

subroutine debug_printout_result(linsys, basename, it)
 type(t_linsys),   intent(in) :: linsys
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: it

 real(wp), allocatable :: res(:,:)
 character(len=max_char_len) :: sit

  allocate(res(1,linsys%rank+linsys%n_expl))
  !!res(1,:) = linsys%res
  res(1,:) = (/linsys%res,linsys%res_expl(:,1)/)
  write(sit,'(I4.4)') it
  call write_basic(res,trim(basename)//'_result_'//trim(sit)//'.dat')
  deallocate(res)

end subroutine debug_printout_result

!------------------------------------------------------------------------------

subroutine debug_printout_geometry(elems, geo, basename, it)
 type(t_impl_elem_p),   intent(in) :: elems(:)
 type(t_geo),      intent(in) :: geo
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: it

 real(wp), allocatable :: norm(:,:), cent(:,:), velb(:,:)
 integer, allocatable  :: el(:,:), conn(:,:)
 character(len=max_char_len) :: sit
 integer :: ie, iv

 ! surf_vel and vel_phi
 real(wp), allocatable :: surf_vel(:,:), vel_phi(:,:) 
 integer :: i_e

 ! chtls
 real(wp) :: f(5)
 integer :: n_neigh

  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  ! surf_vel and vel_phi
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

    ! surf_vel and vel_phi for surfpan only
    select type( el => elems(ie)%p )
     class is(t_surfpan) 
      surf_vel(:,ie) = el%surf_vel   ! elems(ie)%p%surf_vel

!     ! 1. FVM approx --------      
!     vel_phi(:,ie) = 0.0_wp
!     do i_e = 1 , el%n_ver    ! elems(ie)%p%n_ver
!       if ( associated(el%neigh(i_e)%p) ) then !  .and. &
!         vel_phi(:,ie) = vel_phi(:,ie) + &
!           el%pot_vel_stencil(:,i_e) * (el%neigh(i_e)%p%mag - el%mag)
!       end if
!     end do
!     vel_phi(:,ie)  = - vel_phi(:,ie)

      ! 2. CHTLS method ------      
      vel_phi(:,ie) = 0.0_wp
      f = 0.0_wp ; n_neigh = 0
      do i_e = 1 , el%n_ver
        if ( associated(el%neigh(i_e)%p) ) then
          n_neigh = n_neigh + 1
          f(n_neigh) = - ( el%neigh(i_e)%p%mag - el%mag )
        end if
      end do
      f(n_neigh+1) = sum(el%nor * (-sim_param%u_inf - el%uvort + el%ub) )
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
 type(t_geo),      intent(in) :: geo
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: it

 real(wp), allocatable :: norm(:,:), cent(:,:)
 integer, allocatable  :: el(:,:)
 character(len=max_char_len) :: sit
 integer :: ie
 integer(h5loc) :: h5fid

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
 !type(t_expl_elem_p),   intent(in) :: elems(:)
 type(t_liftlin_p),   intent(in) :: elems(:)
 type(t_geo),      intent(in) :: geo
 character(len=*), intent(in) :: basename
 integer,          intent(in) :: it

 real(wp), allocatable :: norm(:,:), cent(:,:), velb(:,:)
 integer, allocatable  :: el(:,:), conn(:,:)
 character(len=max_char_len) :: sit
 integer :: ie, iv


  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;

! only for surfpan !!!!!!
! ! surf_vel and vel_phi
! allocate( surf_vel(3,size(elems)), vel_phi(3,size(elems)) )
! surf_vel = -666.6 ; vel_phi = -666.6 
 
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

! only for surfpan !!!!!!
!   ! surf_vel and vel_phi for surfpan only
!   select type( el => elems(ie)%p )
!    class is(t_surfpan) 
!     surf_vel(:,ie) = el%surf_vel   ! elems(ie)%p%surf_vel
!     
!     vel_phi(:,ie) = 0.0_wp
!     do i_e = 1 , el%n_ver    ! elems(ie)%p%n_ver
!       if ( associated(el%neigh(i_e)%p) ) then !  .and. &
!         vel_phi(:,ie) = vel_phi(:,ie) + &
!           el%pot_vel_stencil(:,i_e) * (el%neigh(i_e)%p%mag - el%mag)
!       end if
!     end do
!     vel_phi(:,ie)  = - vel_phi(:,ie)
!
!   end select
  
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_ll_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_ll_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(velb,       trim(basename)//'_ll_mesh_velb_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_ll_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_ll_mesh_elems_'  //trim(sit)//'.dat')
  call write_basic(conn,       trim(basename)//'_ll_mesh_conn_'   //trim(sit)//'.dat')
  deallocate(norm, cent, el, conn, velb)

! only for surfpan !!!!!!
! ! surf_vel and vel_phi
! call write_basic(surf_vel,   trim(basename)//'_mesh_surfvel_'  //trim(sit)//'.dat')
! call write_basic(vel_phi ,   trim(basename)//'_mesh_velphi_'   //trim(sit)//'.dat')
! deallocate(surf_vel,vel_phi)

end subroutine debug_ll_printout_geometry

!------------------------------------------------------------------------------
!----------------------------------------------------------------------
!UNDER SCRUTINY: employs old stuff, to remove?
!subroutine debug_printout_loads(elems, basename_debug, it)
! type(t_elem_p),   intent(in) :: elems(:)
! character(len=*), intent(in) :: basename_debug
! integer,          intent(in) :: it
!
! real(wp), allocatable :: vel(:,:), cp(:,:), F_aero(:,:)
! integer :: i_el
!
!  allocate( vel(3,size(elems)) )
!  allocate(  cp(1,size(elems)) )
!  allocate(F_aero(3,1))
!  F_aero = 0.0_wp
!  do i_el = 1 , size(elems)
!    vel(:,i_el) = elems(i_el)%p%vel
!    cp (1,i_el) = elems(i_el)%p%cp
!    F_aero(:,1) = F_aero(:,1) - 0.5_wp * rho * norm2(uinf)**2.0_wp * cp(1,i_el) * &
!                     elems(i_el)%p%area * elems(i_el)%p%nor
!  end do
!  write(frmt,'(I4.4)') it
!  call write_basic(vel,trim(basename_debug)//'_velocity_'//trim(frmt)//'.dat')
!  call write_basic(cp ,trim(basename_debug)//'_cp_'//trim(frmt)//'.dat')
!  call write_basic(F_aero ,trim(basename_debug)//'_Faero_'//trim(frmt)//'.dat')
!  deallocate(vel,cp,F_aero)
!
!end subroutine debug_printout_loads

!----------------------------------------------------------------------

!Consider discontinuing
!DISCONTINUED
!subroutine debug_printout_wake(wake_panels, basename, it)
! type(t_wake), intent(in) :: wake_panels
! character(len=*), intent(in) :: basename
! integer,          intent(in) :: it
!
! real(wp), allocatable :: norm(:,:), cent(:,:), res(:,:)
! integer, allocatable :: el(:,:)
! character(len=max_char_len) :: sit
! integer :: ie, of, p1, p2
!
!  allocate(norm(3,size(wake_panels%pan_p)), cent(3,size(wake_panels%pan_p)))
!  allocate(el(4,size(wake_panels%pan_p))); el = 0
!  allocate(res(1,size(wake_panels%pan_p)))
!  do ie=1,size(wake_panels%pan_p)
!    norm(:,ie) = wake_panels%pan_p(ie)%p%nor
!    cent(:,ie) = wake_panels%pan_p(ie)%p%cen
!    p1 = wake_panels%i_start_points(1,mod(ie-1,wake_panels%n_wake_stripes)+1)
!    p2 = wake_panels%i_start_points(2,mod(ie-1,wake_panels%n_wake_stripes)+1)
!    of = wake_panels%n_wake_points*((ie-1)/wake_panels%n_wake_stripes)
!    el(1:4,ie) = (/of+p2, of+p1, of+p1+wake_panels%n_wake_points, &
!                     of+p2+wake_panels%n_wake_points/)
!    res(1,ie) = wake_panels%pan_p(ie)%p%mag
!  enddo
!  write(sit,'(I4.4)') it
!  call write_basic( &
!    reshape(wake_panels%w_points(:,:,1:wake_panels%wake_len+1),&
!    (/3,(wake_panels%n_wake_points)*(wake_panels%wake_len+1)/)), &
!    trim(basename)//'_wake_points_'//trim(sit)//'.dat')
!  call write_basic( &
!    reshape(wake_panels%w_vel(:,:,1:wake_panels%wake_len+1),&
!    (/3,(wake_panels%n_wake_points)*(wake_panels%wake_len+1)/)), &
!    trim(basename)//'_wake_vels_'//trim(sit)//'.dat')
!  call write_basic(norm, trim(basename)//'_wake_norm_'//trim(sit)//'.dat')
!  call write_basic(cent, trim(basename)//'_wake_cent_'//trim(sit)//'.dat')
!  call write_basic(el,   trim(basename)//'_wake_elems_'//trim(sit)//'.dat')
!  call write_basic(res,  trim(basename)//'_wake_result_'//trim(sit)//'.dat')
!  deallocate(norm, cent, el, res)
!end subroutine debug_printout_wake

!------------------------------------------------------------------------------

end program dust




