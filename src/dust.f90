!!=====================================================================
!!
!! Copyright (C) 2018 Politecnico di Milano
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
!!          Federico Fonte             <federico.fonte@polimi.it>
!!          Davide Montagnani       <davide.montagnani@polimi.it>
!!          Matteo Tugnoli             <matteo.tugnoli@polimi.it>
!!=====================================================================

!> This is a more structured version of the test code to build a sort of
!! architecture proof

program dust

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

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

use mod_liftlin, only: &
 update_liftlin, solve_liftlin

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

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

use mod_wake, only: &
  t_wake, initialize_wake, update_wake, &
  prepare_wake, load_wake, destroy_wake

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

use mod_octree, only: &
  initialize_octree, destroy_octree, sort_particles, t_octree
 

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
type(t_sim_param) :: sim_param
! Asymptotic conditions

!Time parameters
real(wp) :: time
integer  :: it, nstep, nout
real(wp) :: t_last_out, t_last_debug_out
logical  :: time_2_out, time_2_debug_out
logical  :: output_start
real(wp) :: dt_debug_out

!Main variables
!> All the implicit elements, sorted first static then moving
type(t_impl_elem_p), allocatable :: elems(:)
!> All the explicit elements
type(t_expl_elem_p), allocatable :: elems_expl(:)
!> Only the lifting line elements
type(t_expl_elem_p), allocatable :: elems_ll(:)
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

!Restart
logical :: restart
character(len=max_char_len) :: restart_file
logical :: reset_time


!I/O prefixes
character(len=max_char_len) :: frmt
character(len=max_char_len) :: basename_debug

real(wp) , allocatable :: res_old(:)

integer :: i_el , i

!octree parameters
type(t_octree) :: octree


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
call prms%CreateRealOption( 'dt_out', "output time interval")
call prms%CreateRealOption( 'dt_debug_out', "debug output time interval")

! input:
call prms%CreateStringOption('GeometryFile','Main geometry definition file')
call prms%CreateStringOption('ReferenceFile','Reference frames file','no_set')

! output:
call prms%CreateStringOption('basename','oputput basename','./')
call prms%CreateStringOption('basename_debug','oputput basename for debug', &
                                                                         './')
call prms%CreateLogicalOption( 'output_start', "output values at starting &
                                                          & iteration", 'F')
call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', &
                                                                         '0')

! restart
call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
call prms%CreatestringOption('restart_file','restart file name')
call prms%CreateLogicalOption('reset_time','reset the time from previous &
                               &execution?','F')

! parameters:
call prms%CreateRealArrayOption( 'u_inf', "free stream velocity", &
                                                       '(/1.0, 0.0, 0.0/)')
call prms%CreateRealOption( 'P_inf', "free stream pressure", '1.0')
call prms%CreateRealOption( 'rho_inf', "free stream density", '1.0')
call prms%CreateRealOption( 'a_inf', "Speed of sound", '340.0')  ! m/s
call prms%CreateRealOption( 'mu_inf', "Dynamic viscosity", '0.00001') ! kg/ms
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

call prms%CreateRealOption( 'FarFieldRatioDoublet', &
      "Multiplier for far field threshold computation on doublet", '10.0')
call prms%CreateRealOption( 'FarFieldRatioSource', &
      "Multiplier for far field threshold computation on sources", '10.0')
call prms%CreateRealOption( 'DoubletThreshold', &
      "Thresold for considering the point in plane in doublets", '1.0e-6')
call prms%CreateRealOption( 'RankineRad', &
      "Radius of Rankine correction for vortex induction near core", '0.1')
call prms%CreateRealOption( 'CutoffRad', &
      "Radius of complete cutoff  for vortex induction near core", '0.001')

call prms%CreateLogicalOption('rigid_wake','rigid wake?','F')
call prms%CreateRealArrayOption( 'rigid_wake_vel', "rigid wake velocity" )


!== Octree and multipole data == 
call prms%CreateLogicalOption('FMM','Employ fast multipole method?','T')
call prms%CreateRealOption('BoxLength','length of the octree box')
call prms%CreateIntArrayOption('NBox','number of boxes in each direction')
call prms%CreateRealArrayOption( 'OctreeOrigin', "rigid wake velocity" )
call prms%CreateIntOption('NOctreeLevels','number of octree levels')
call prms%CreateIntOption('MinOctreePart','minimum number of octree &
                                                             &particles')
call prms%CreateIntOption('MultipoleDegree','multipole expansion degree')
call prms%CreateLogicalOption('Vortstretch','Employ vortex stretching','T')


! get the parameters and print them out
call printout(nl//'====== Input parameters: ======')
call prms%read_options(input_file_name, printout_val=.true.)
!timing
sim_param%t0 = getreal(prms, 'tstart')
sim_param%tend   = getreal(prms, 'tend')
sim_param%dt     = getreal(prms, 'dt')
sim_param%dt_out = getreal(prms,'dt_out')
dt_debug_out = getreal(prms, 'dt_debug_out')
output_start = getlogical(prms, 'output_start')
sim_param%debug_level = getint(prms, 'debug_level')
!Reference values
sim_param%P_inf = getreal(prms,'P_inf')
sim_param%rho_inf  = getreal(prms,'rho_inf')
sim_param%u_inf = getrealarray(prms, 'u_inf', 3)
sim_param%a_inf  = getreal(prms,'a_inf')
sim_param%mu_inf = getreal(prms,'mu_inf')
!Wake parameters
sim_param%n_wake_panels = getint(prms, 'n_wake_panels')
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
!Names
sim_param%basename = getstr(prms,'basename')
basename_debug = getstr(prms,'basename_debug')
sim_param%GeometryFile = getstr(prms,'GeometryFile')
sim_param%ReferenceFile = getstr(prms,'ReferenceFile')
!Method parameters
sim_param%FarFieldRatioDoublet  = getreal(prms, 'FarFieldRatioDoublet')
sim_param%FarFieldRatioSource  = getreal(prms, 'FarFieldRatioSource')
sim_param%DoubletThreshold   = getreal(prms, 'DoubletThreshold')
sim_param%RankineRad = getreal(prms, 'RankineRad')
sim_param%CutoffRad  = getreal(prms, 'CutoffRad')
sim_param%first_panel_scaling = getreal(prms,'ImplicitPanelScale')
sim_param%min_vel_at_te  = getreal(prms,'ImplicitPanelMinVel')
sim_param%use_vs = getlogical(prms, 'Vortstretch')
!Octree and FMM parameters
sim_param%use_fmm = getlogical(prms, 'FMM')
sim_param%BoxLength = getreal(prms, 'BoxLength')
sim_param%NBox = getintarray(prms, 'NBox',3)
sim_param%OctreeOrigin = getrealarray(prms, 'OctreeOrigin',3)
sim_param%NOctreeLevels = getint(prms, 'NOctreeLevels')
sim_param%MinOctreePart = getint(prms, 'MinOctreePart')
sim_param%MultipoleDegree = getint(prms,'MultipoleDegree')




!-- Parameters Initializations --
call initialize_doublet(sim_param%FarFieldRatioDoublet, &
                        sim_param%DoubletThreshold, sim_param%RankineRad, &
                        sim_param%CutoffRad);
call initialize_vortline(sim_param%RankineRad, sim_param%CutoffRad);
call initialize_vortpart(sim_param%RankineRad, sim_param%CutoffRad);
call initialize_surfpan(sim_param%FarFieldRatioSource);
!reset the numbering for output files
nout = 0
!Manage restart
restart = getlogical(prms,'restart_from_file')
if (restart) then
  reset_time = getlogical(prms,'reset_time')
  restart_file = getstr(prms,'restart_file')
  call printout('RESTART: restarting from file: '//trim(restart_file))
  sim_param%GeometryFile = restart_file(1:len(trim(restart_file))-11)&
                                                           &//'geo.h5'

  !restarting the same simulation, advance the numbers
  if(restart_file(1:len(trim(restart_file))-12).eq. &
                                               trim(sim_param%basename)) then
  read(restart_file(len(trim(restart_file))-6:len(trim(restart_file))-3),*) &
                                                                         nout 
    call printout('Identified restart from the same simulation, keeping the&
    & previous output numbering')
    !avoid rewriting the same timestep
    output_start = .false.
  endif
  if(.not. reset_time) call load_time(restart_file, sim_param%t0)
endif

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


!---- Simulation parameters ----
nstep = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1 
       !(+1 for the zero time step)
sim_param%n_timesteps = nstep
allocate(sim_param%time_vec(sim_param%n_timesteps))
sim_param%time_vec = (/ ( sim_param%t0 + &
         dble(i-1)*sim_param%dt, i=1,sim_param%n_timesteps ) /)

!------ Geometry creation ------
call printout(nl//'====== Geometry Creation ======')
t0 = dust_time()
target_file = trim(sim_param%basename)//'_geo.h5'
call create_geometry(sim_param%GeometryFile, sim_param%ReferenceFile, &
                     input_file_name, geo, te, elems, elems_expl, elems_ad, &
                     elems_ll, elems_tot, airfoil_data, sim_param, &
                     target_file, run_id)

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
call ignoredParameters(prms)
call finalizeParameters(prms)


!------ Initialization ------
call printout(nl//'====== Initializing Wake ======')

call initialize_wake(wake, geo, te, sim_param%n_wake_panels, &
       sim_param%n_wake_panels, sim_param%n_wake_particles, &
       sim_param%particles_box_min, &
       sim_param%particles_box_max,  sim_param)

call printout(nl//'====== Initializing Linear System ======')
t0 = dust_time()
call initialize_linsys(linsys, geo, elems, elems_expl, &
                       wake,  sim_param%u_inf)
t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Initialized linear system in: ' , t1 - t0,' s.'
  call printout(message)
endif

!===========EXPERIMENTAL PART, OCTREE========
call printout(nl//'====== Initializing Octree ======')
t0 = dust_time()
call initialize_octree(sim_param%BoxLength, sim_param%NBox, &
                       sim_param%OctreeOrigin, sim_param%NOctreeLevels, &
                       sim_param%MinOctreePart, sim_param%MultipoleDegree, &
                       sim_param%RankineRad, sim_param, octree)
t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Initialized octree in: ' , t1 - t0,' s.'
  call printout(message)
endif
!============================================

!------ Reloading ------
if (restart) then
 call load_solution(restart_file, geo%components, geo%refs)
 call load_wake(restart_file, wake)
endif

t22 = dust_time()
write(message,'(A,F9.3,A)') nl//'------ Completed all preliminary operations &
                             &in: ' , t22 - t00,' s.'
call printout(message)


!====== Time Cycle ======
time = sim_param%t0
t_last_out = time; t_last_debug_out = time

allocate(res_old(size(elems))) ; res_old = 0.0_wp
t11 = dust_time()
do it = 1,nstep

  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,I5,A,I5,A,F7.2)') nl//'--> Step ',it,' of ', &
                                                 nstep, ' time: ', time
    call printout(message)
    t22 = dust_time()
    write(message,'(A,F9.3,A)') 'Elapsed time: ',t22-t00
    call printout(message)
  endif

  call init_timestep(time)

  !call update_geometry(geo, time, .false.)
  !call prepare_wake(wake, geo, sim_param, it)

  call update_liftlin(elems_ll,linsys)
  call update_actdisk(elems_ad,linsys,sim_param)


  if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
            call debug_printout_geometry(elems, geo, basename_debug, it)
  if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
            call debug_ll_printout_geometry(elems_ll, geo, basename_debug, it)


  !------ Assemble the system ------
  !call prepare_wake(wake, geo, sim_param)
  t0 = dust_time()

  call assemble_linsys(linsys, elems, elems_expl,  &
                       wake, sim_param%u_inf)
  t1 = dust_time()

  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,F9.3,A)') 'Assembled linear system in: ' , t1 - t0,' s.'
    call printout(message)
  endif

  if ((sim_param%debug_level .ge. 50).and.time_2_debug_out) then
    write(frmt,'(I4.4)') it
    call dump_linsys(linsys, trim(basename_debug)//'A_'//trim(frmt)//'.dat', &
                             trim(basename_debug)//'b_'//trim(frmt)//'.dat' )
  endif

  !------ Solve the system ------
  t0 = dust_time()
  if (linsys%rank .gt. 0) then
    call solve_linsys(linsys)
  endif
  t1 = dust_time()

  ! compute time derivative of the result ( = i_vortex = -i_doublet ) -------
  do i_el = 1 , size(elems)
    elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
  end do
  res_old = linsys%res

  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,F9.3,A)')  'Solved linear system in: ' , t1 - t0,' s.'
    call printout(message)
  endif

  if (sim_param%debug_level .ge. 20.and.time_2_debug_out) &
                      call debug_printout_result(linsys, basename_debug, it)

  !------ Update the explicit part ------  % v-----implicit elems: p,v
  if ( size(elems_ll) .gt. 0 ) then
    call solve_liftlin(elems_ll, elems_tot, elems , elems_ad , &
            (/ wake%pan_p, wake%rin_p/), wake%vort_p, sim_param, airfoil_data)
  end if

  !------ Compute loads -------
  ! Implicit elements: vortex rings and 3d-panels
  do i_el = 1 , size(elems)
    call elems(i_el)%p%compute_pres(sim_param)
    call elems(i_el)%p%compute_dforce(sim_param)
  end do
  ! Explicit elements:
  ! - liftlin: _pres and _dforce computed in solve_liftin()
  ! - actdisk: avg delta_pressure and force computed here,
  !            to include thier effects in postpro (e.g. integral loads)
  do i_el = 1 , size(elems_ad)
    call elems_ad(i_el)%p%compute_pres(sim_param)
    call elems_ad(i_el)%p%compute_dforce(sim_param)
  end do

  !if ((sim_param%debug_level .ge. 10).and.time_2_debug_out) then
  !  call debug_printout_loads(elems, basename_debug, it)
  !endif

  !------ Output the results  ------
  !Printout the wake
  !DISCONTINUED
  !if((sim_param%debug_level .ge. 17).and.time_2_debug_out)  &
  !                    call debug_printout_wake(wake, basename_debug, it)

  !Print the results
  if(time_2_out)  then
    nout = nout+1
    call save_status(geo, wake, sim_param, nout, time, run_id)
  endif

  !------ Treat the wake ------
  ! (this needs to be done after output, in practice the update is for the
  !  next iteration)
  t0 = dust_time()
  call update_wake(wake, elems_tot, octree, sim_param)
  t1 = dust_time()
  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,F9.3,A)') 'Updated wake in: ' , t1 - t0,' s.'
    call printout(message)
  endif

  time = min(sim_param%tend, time+sim_param%dt)
  call update_geometry(geo, time, .false.)
  call prepare_wake(wake, geo, sim_param)

enddo


deallocate( res_old )
!===== End Time Cycle ======


!------ Cleanup ------
!call destroy_wake_panels(wake_panels)
call destroy_wake(wake)
call destroy_octree(octree)
call destroy_linsys(linsys)
call destroy_elements(geo)
call destroy_geometry(geo, elems_tot)

call destroy_hdf5()

t22 = dust_time()
write(message,'(A,F9.3,A)') 'Completed all computations in ',t22-t00,' s'
call printout(message)
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

!---------------------------------------------------------------------_
!DISCONTINUED: consider removing
subroutine copy_geo(sim_param, geo_file, run_id)
 type(t_sim_param), intent(inout) :: sim_param
 character(len=*), intent(inout)     :: geo_file
 integer, intent(in)              :: run_id(10)

 character(len=max_char_len) :: target_file
 integer :: estat, cstat
 integer(h5loc) :: floc

  !target file name: same as run basename with appendix
  target_file = trim(sim_param%basename)//'_geo.h5'
  
  if (trim(geo_file) .ne. trim(target_file)) then
  !Copy the geometry file
  call execute_command_line('cp '//trim(geo_file)//' '//trim(target_file), &
                                           exitstat=estat,cmdstat=cstat)
  if((cstat .ne. 0) .or. (estat .ne. 0)) &
    call error('dust','','System errors while trying to copy the geometry &
    &to the output path')
  endif


  !Attach the run_id to the file as an attribute
  call open_hdf5_file(trim(target_file), floc)
  call write_hdf5_attr(run_id, 'run_id', floc)
  call close_hdf5_file(floc)


  !Overwrite the geo file name, so that the copy is going to be
  !opened
  geo_file = trim(target_file)

end subroutine copy_geo

!----------------------------------------------------------------------

!> Perform preliminary procedures each timestep, mainly chech if it is time
!! to perform output or not
subroutine init_timestep(t)
 real(wp), intent(in) :: t

  if (real(t-t_last_out) .ge. real(sim_param%dt_out)) then
    time_2_out = .true.
    t_last_out = t
  else
    time_2_out = .false.
  endif

  if (real(t-t_last_debug_out) .ge. real(dt_debug_out)) then
    time_2_debug_out = .true.
    t_last_debug_out = t
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

end subroutine init_timestep

!----------------------------------------------------------------------
!DISCONTINUED: substituted by the routines in dust_io and dust_post
!subroutine output_status(elems_tot, geo, wake_panels, basename, it, t)
! type(t_elem_p),   intent(in) :: elems_tot(:)
! type(t_geo),      intent(in) :: geo
! type(t_wake_panels), intent(in) :: wake_panels
! character(len=*), intent(in) :: basename
! integer,          intent(in) :: it
! real(wp), intent(in)         :: t
!
! integer, allocatable :: el(:,:), w_el(:,:)
! real(wp), allocatable :: w_points(:,:), w_res(:)
! integer :: ie, of, p1, p2
! character(len=max_char_len) :: sit
!
!  allocate(el(4,size(elems_tot))); el = 0
!  allocate(w_el(4,size(wake_panels%pan_p))); w_el = 0
!  allocate(w_points(3,(wake_panels%n_wake_points)*(wake_panels%wake_len+1)))
!  allocate(w_res(size(wake_panels%pan_p)))
!
!  !=== VTK output ===
!  do ie=1,size(elems_tot)
!    el(1:elems_tot(ie)%p%n_ver,ie) = elems_tot(ie)%p%i_ver
!  enddo
!  do ie=1,size(wake_panels%pan_p)
!    p1 = wake_panels%i_start_points(1,mod(ie-1,wake_panels%n_wake_stripes)+1)
!    p2 = wake_panels%i_start_points(2,mod(ie-1,wake_panels%n_wake_stripes)+1)
!    !of = ie-mod(ie,wake_panels%n_wake_stripes-1)
!    of = wake_panels%n_wake_points*((ie-1)/wake_panels%n_wake_stripes)
!    w_el(1:4,ie) = (/of+p2, of+p1, of+p1+wake_panels%n_wake_points, &
!                     of+p2+wake_panels%n_wake_points/)
!    w_res(ie) = wake_panels%pan_p(ie)%p%idou
!  enddo
!  w_points = reshape(wake_panels%w_points(:,:,1:wake_panels%wake_len+1),&
!    (/3,(wake_panels%n_wake_points)*(wake_panels%wake_len+1)/))
!  write(sit,'(I4.4)') it
!  call vtk_out_bin (geo%points, el, (/linsys%res,linsys%res_expl(:,1)/),  &
!                    w_points, w_el, w_res,  &
!                    trim(basename)//'_res_'//trim(sit)//'.vtu')
!  call tec_out_sol_bin(geo%points, el, (/linsys%res,linsys%res_expl(:,1)/),  &
!                    w_points, w_el, w_res, t,  &
!                    trim(basename)//'_res_'//trim(sit)//'.plt')
!
!
!  deallocate(el,w_el,w_points,w_res)
!
!end subroutine output_status

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
      
      vel_phi(:,ie) = 0.0_wp
      do i_e = 1 , el%n_ver    ! elems(ie)%p%n_ver
        if ( associated(el%neigh(i_e)%p) ) then !  .and. &
          vel_phi(:,ie) = vel_phi(:,ie) + &
            el%pot_vel_stencil(:,i_e) * (el%neigh(i_e)%p%mag - el%mag)
        end if
      end do
      vel_phi(:,ie)  = - vel_phi(:,ie)
 
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

! ------------------------------------------------------------------------------

subroutine debug_ll_printout_geometry(elems, geo, basename, it)
 type(t_expl_elem_p),   intent(in) :: elems(:)
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




