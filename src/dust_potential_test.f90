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
  wp, nl, max_char_len, pi, extended_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_math, only: &
  cross

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
  initialize_octree, sort_particles, t_octree, perform_multipole
  
use mod_vortpart, only: &
  t_vortpart, t_vortpart_p
implicit none

!run-id
integer :: run_id(10)

!Input
!> Main parameters parser
type(t_parse) :: prms
character(len=*), parameter :: input_file_name_def = 'dust.in'
character(len=max_char_len) :: input_file_name
character(len=max_char_len) :: geo_file_name
character(len=max_char_len) :: target_file
character(len=max_char_len) :: ref_file_name
character(len=extended_char_len) :: message

!Simulation parameters
type(t_sim_param) :: sim_param
! Asymptotic conditions
real(wp) :: uinf(3)
real(wp) :: rho , Pinf, a_inf , mu_inf  ! Re, Mach

!Time parameters
real(wp) :: tstart, tend, dt, time
integer  :: it, nstep, nout
real(wp) :: t_last_out, t_last_debug_out
logical  :: time_2_out, time_2_debug_out
logical  :: output_start
real(wp) :: dt_out, dt_debug_out

!Wake parameters
!> Number of wake panels(/rings)
integer :: n_wake_panels
!> Number of wake particles
integer :: n_parts
real(wp) :: part_box_min(3), part_box_max(3)
real(wp) :: wake_pan_scaling , wake_pan_minvel
logical  :: rigid_wake
real(wp) :: rigid_wake_vel(3)
character(len=max_char_len) :: rigid_wake_vel_str

!doublet parameters
real(wp) :: ff_ratio_dou, ff_ratio_sou, eps_dou, r_Rankine, r_cutoff

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
!> Level of debug output
integer :: debug_level

!Restart
logical :: restart
character(len=max_char_len) :: restart_file
logical :: reset_time


!I/O prefixes
character(len=max_char_len) :: frmt
character(len=max_char_len) :: basename
character(len=max_char_len) :: basename_debug

real(wp) , allocatable :: res_old(:)

integer :: i_el , i, j, k
integer :: iq, jq, kq
integer :: ip, npart_tot


!octree stuff
type(t_octree) :: octree
real(wp) :: BoxLength
integer :: NBox(3)
real(wp) :: OctreeOrigin(3)
integer :: NOctreeLevels, MinOctreePart
integer :: MultipoleDegree

!particles
type(t_vortpart), allocatable, target :: parts(:,:,:)
real(wp),target,  allocatable :: prt_mag(:,:,:)
type(t_vortpart_p), allocatable :: part_p(:)

real(wp), allocatable :: velocity_bf(:,:,:,:)
real(wp), allocatable, target :: velocity_fm(:,:,:,:)
real(wp), allocatable :: velocity_fm_long(:,:)

real(wp) :: partmag

real(wp) :: err, rel_err, R2



call printout(nl//'>>>>>> DUST potential test beginning >>>>>>'//nl)
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

! output:
call prms%CreateStringOption('basename','oputput basename','./')
call prms%CreateStringOption('basename_debug','oputput basename for debug','./')
call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output','0')

call prms%CreateIntOption('n_particles', 'number of wake particles', &
                                                                  '10000')
call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
     &particles bounding box', '(/-10.0, -10.0, -10.0/)')
call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
     &particles bounding box', '(/10.0, 10.0, 10.0/)')

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

!== Octree stuff == 
call prms%CreateRealOption('BoxLength','length of the octree box')
call prms%CreateIntArrayOption('NBox','number of boxes in each direction','(/1,1,1/)')
call prms%CreateRealArrayOption( 'OctreeOrigin', "rigid wake velocity" )
call prms%CreateIntOption('NOctreeLevels','number of octree levels')
call prms%CreateIntOption('MinOctreePart','minimum number of octree particles')
call prms%CreateIntOption('MultipoleDegree','')


! get the parameters and print them out
call prms%read_options(input_file_name, printout_val=.true.)

debug_level = getint(prms, 'debug_level')
n_parts = getint(prms, 'n_particles')
part_box_min = getrealarray(prms, 'particles_box_min',3)
part_box_max = getrealarray(prms, 'particles_box_max',3)

basename = getstr(prms,'basename')
basename_debug = getstr(prms,'basename_debug')

ff_ratio_dou  = getreal(prms, 'FarFieldRatioDoublet')
ff_ratio_sou  = getreal(prms, 'FarFieldRatioSource')
eps_dou   = getreal(prms, 'DoubletThreshold')
r_Rankine = getreal(prms, 'RankineRad')
r_cutoff  = getreal(prms, 'CutoffRad')

BoxLength = getreal(prms, 'BoxLength')
NBox = getintarray(prms, 'NBox',3)
OctreeOrigin = getrealarray(prms, 'OctreeOrigin',3)
NOctreeLevels = getint(prms, 'NOctreeLevels')
MinOctreePart = getint(prms, 'MinOctreePart')
MultipoleDegree = getint(prms, 'MultipoleDegree')


!Octree stuff

!-- Parameters Initializations --
call initialize_doublet(ff_ratio_dou, eps_dou, r_Rankine, r_cutoff);
call initialize_vortline(r_Rankine, r_cutoff);
call initialize_vortpart(r_Rankine, r_cutoff);
call initialize_surfpan(ff_ratio_sou);

nout = 0


!---- Simulation parameters ----
sim_param%debug_level = debug_level
sim_param%basename = basename

call ignoredParameters(prms)
call finalizeParameters(prms)


!===========EXPERIMENTAL PART, OCTREE========
call printout(nl//'====== Initializing Octree ======')
t0 = dust_time()
call initialize_octree(BoxLength, NBox, OctreeOrigin, &
                       NOctreeLevels, MinOctreePart, MultipoleDegree, r_Rankine, octree)
t1 = dust_time()
write(message,'(A,F9.3,A)') 'Initialized octree in: ' , t1 - t0,' s.'
call printout(message)
!============================================


!==========INITIALIZE THE PARTICLES==========
call printout(nl//'====== Initializing Particles ======')
t0 = dust_time()
allocate(parts(n_parts,n_parts,n_parts))
allocate(prt_mag(n_parts,n_parts,n_parts))
allocate(part_p(n_parts*n_parts*n_parts))
allocate(velocity_bf(3,n_parts,n_parts,n_parts))
allocate(velocity_fm(3,n_parts,n_parts,n_parts))
allocate(velocity_fm_long(3,n_parts*n_parts*n_parts))
ip = 0
partmag = 1.0_wp/real(n_parts**3,wp)
do k=1,n_parts; do j=1,n_parts; do i=1,n_parts 
  parts(i,j,k)%mag => prt_mag(i,j,k)
  parts(i,j,k)%mag = partmag
  parts(i,j,k)%dir = (/1.0_wp, 0.0_wp, 0.0_wp/)
  !parts(i,j,k)%cen(1) = -BoxLength/2.0_wp+real(i,wp)*BoxLength/real(n_parts,wp)
  !parts(i,j,k)%cen(2) = -BoxLength/2.0_wp+real(j,wp)*BoxLength/real(n_parts,wp)
  !parts(i,j,k)%cen(3) = -BoxLength/2.0_wp+real(k,wp)*BoxLength/real(n_parts,wp)
  parts(i,j,k)%cen(1) = OctreeOrigin(1)+real(i,wp)*BoxLength/real(n_parts,wp)
  parts(i,j,k)%cen(2) = OctreeOrigin(2)+real(j,wp)*BoxLength/real(n_parts,wp)
  parts(i,j,k)%cen(3) = OctreeOrigin(3)+real(k,wp)*BoxLength/real(n_parts,wp)
  ip = ip + 1
  part_p(ip)%p => parts(i,j,k)
  parts(i,j,k)%vel => velocity_fm(:,i,j,k)

enddo; enddo; enddo;
npart_tot = n_parts**3

t1 = dust_time()
write(message,'(A,F9.3,A)') 'Initialized particles in: ' , t1 - t0,' s.'
call printout(message)
!============================================


!==========BRUTE FORCE INTERACTION===========
call printout(nl//'====== Performing brute force interactions ======')
t0 = dust_time()
velocity_bf = 0.0_wp
do kq=1,n_parts; do jq=1,n_parts; do iq=1,n_parts 
  do k=1,n_parts; do j=1,n_parts; do i=1,n_parts 
    if(.not.all( (/i-iq, j-jq, k-kq/) .eq. 0)) then
    R2 = sum((parts(iq,jq,kq)%cen-parts(i,j,k)%cen)**2)+r_Rankine**2
    velocity_bf(:,iq,jq,kq) =  velocity_bf(:,iq,jq,kq) - &
    parts(i,j,k)%mag/(4.0_wp*pi*sqrt(R2)**3) * &
      cross((-parts(i,j,k)%cen+parts(iq,jq,kq)%cen),parts(i,j,k)%dir )

    endif
  enddo; enddo; enddo;
enddo; enddo; enddo;
t1 = dust_time()
write(message,'(A,F9.3,A)') 'Brute force interactions in: ' , t1 - t0,' s.'
call printout(message)
!============================================


!=========OCTREE INTERACTION=================

call printout(nl//'====== Performing octree interactions ======')
t0 = dust_time()
call sort_particles(part_p, octree)
call perform_multipole(part_p, octree)


t1 = dust_time()
write(message,'(A,F9.3,A)') 'Octree interactions in: ' , t1 - t0,' s.'
call printout(message)
!============================================


!=========CHECKS AND OUTPUT==================
!ip = 0
!do kq=1,n_parts; do jq=1,n_parts; do iq=1,n_parts 
!  ip = ip + 1
!  velocity_fm(:,iq,jq,kq) =  velocity_fm_long(:,ip)
!enddo; enddo; enddo;

err = norm2(velocity_fm-velocity_bf)
rel_err = err/norm2(velocity_bf)
write(*,*) 'error',err
write(*,*) 'rel_err',rel_err

  call write_basic(velocity_fm(2,:,:,2),'fast_multipole.dat')
  call write_basic(velocity_bf(2,:,:,2),'brute_force.dat')

!============================================








t22 = dust_time()
write(message,'(A,F9.3,A)') 'Completed all computations in ',t22-t00,' s'
call printout(message)
call printout(nl//'<<<<<< DUST potential test end <<<<<<'//nl)



!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
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

!------------------------------------------------------------------------------

!> Perform preliminary procedures each timestep, mainly chech if it is time
!! to perform output or not
subroutine init_timestep(t)
 real(wp), intent(in) :: t

  if (real(t-t_last_out) .ge. real(dt_out)) then
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
  if((it .eq. nstep) .and. (dt_out .le. tend)) then
    time_2_out = .true.
    time_2_debug_out = .true.
  endif
 
  !if requested, print the output also at the first step (t0)
  if ((t.eq.tstart) .and. output_start) then
    t_last_out = t
    t_last_debug_out = t
    time_2_out = .true.
    time_2_debug_out = .true.
  endif

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


  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  ! surf_vel and vel_phi
  allocate( surf_vel(3,size(elems)), vel_phi(3,size(elems)) )
  surf_vel = -666.6 ; vel_phi = -666.6 
 
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

 ! surf_vel and vel_phi
 real(wp), allocatable :: surf_vel(:,:), vel_phi(:,:) 
 integer :: i_e


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




