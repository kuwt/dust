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
!!          Alberto Savino
!!=========================================================================


!> Module to treat the whole wake
module mod_wake

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_math, only: &
  cross, infinite_plate_spline

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, t_tedge, calc_geo_data_pan, calc_node_vel



use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_liftlin, only: &
  t_liftlin

use mod_vortline, only: &
  t_vortline

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_actuatordisk, only: &
  t_actdisk

use mod_hdf5_io, only: &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  write_hdf5_attr, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5

use mod_octree, only: &
  t_octree, sort_particles, calculate_multipole, apply_multipole, &
  apply_multipole_panels

!$ use omp_lib, only: &
!$   omp_get_thread_num, omp_get_num_threads

use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: t_wake, initialize_wake, update_wake, &
          prepare_wake, complete_wake, load_wake, destroy_wake

private


!> Type containing wake panels information
type :: t_wake

  !> Number of maximum streamwise panels
  integer :: nmax_pan

  !> Number of maximum "rows" of disks
  integer :: nmax_rin

  !> Number of actual streamwise panels
  integer :: pan_wake_len

  !> Actual number of rings in the ring wake
  integer :: rin_wake_len

  !> Number of wake stripes ("spanwise" panels)
  integer :: n_pan_stripes

  !> Number of wake points in the "spanwise" direction
  integer :: n_pan_points

  !> Number of generating disks
  integer :: ndisks

  !> Number of points for each "row"
  integer :: np_row

  !> Pointer and index of the 2 generating elements of the panel wake
  !! (2 x n_pan_stripes)
  type(t_pot_elem_p), allocatable :: pan_gen_elems(:,:)
  integer, allocatable :: pan_gen_elems_id(:,:)

  !> Generating actuator disk elements
  type(t_pot_elem_p), allocatable :: rin_gen_elems(:)

  !> Index of the 2 generating points of each wake point
  !! (2 x n_pan_points)
  integer, allocatable :: pan_gen_points(:,:)

  !> Direction of the wake at the trailing edge
  !! (3xn_pan_points)
  real(wp), allocatable :: pan_gen_dir(:,:)

  !> Reference frame of the generating points (id)
  !! (n_pan_points)
  integer, allocatable :: pan_gen_ref(:)

  !> Component of the generating points
  integer, allocatable :: pan_gen_icomp(:)

  !> Individual scaling of the firs element of the wake
  real(wp), allocatable :: pan_gen_scaling(:)

  !> Panels neighbours in wake numbering
  integer, allocatable :: pan_neigh(:,:)

  !> Relative orientation of neighbours
  integer, allocatable :: pan_neigh_o(:,:)

  !> Index in the stripes of end elements in panel wake
  integer, allocatable :: pan_i_ends(:)

  !> Index of the joined trailing edges
  integer, allocatable :: joined_tes(:,:,:)

  !> Wake starting points: calculated from the 2 (or 1) starting points of the
  !! geometry possibly moved
  real(wp), allocatable :: w_start_points(:,:)

  !> Index of the 2 wake starting points for each wake stripe
  !! (2 x n_pan_stripes)
  integer, allocatable :: i_start_points(:,:)

  !> Points of the wake, in a structured way
  !! (3 x n_pan_points x npan+1)
  real(wp), allocatable :: pan_w_points(:,:,:)

  !> Velocities of the wake panels
  !! (3 x n_pan_points x npan+1) !!!! now used for rotational effects on pressure !!!!
  real(wp), allocatable :: pan_w_vel(:,:,:)

  !> Relative velocity ( u_inf - ub ) of the nodes at the TE
  !! (3 x n_pan_points)
  !! used to determine the first prescribed panel of the wake
  real(wp), allocatable :: w_vel_te(:,:)

  !> Velocity at the nodes of the wake. For output only
  !! (3 x n_pan_points x npan+1)
  real(wp), allocatable :: w_vel(:,:,:)

  !> elements of the panels
  type(t_vortlatt), allocatable :: wake_panels(:,:)

  !> Ring elements
  type(t_actdisk), allocatable :: wake_rings(:,:)

  !> end vortices
  type(t_vortline), allocatable :: end_vorts(:)

  !> doublets intensities
  real(wp), allocatable :: pan_idou(:,:)

  !> vortex intensities
  real(wp), allocatable :: rin_idou(:,:)

  !> pointer to the wake elements to be passed to the linsys
  !! solver
  type(t_pot_elem_p), allocatable :: pan_p(:)

  !> pointer to the wake elements to be passed to the linsys
  !! solver
  type(t_pot_elem_p), allocatable :: rin_p(:)


  !! Particles data

  !> Maximum number of particles
  integer :: nmax_prt

  !> Actual number of particles
  integer :: n_prt

  !> Wake particles
  type(t_vortpart), allocatable :: wake_parts(:)

  !> Magnitude of particles vorticity
  real(wp), allocatable :: prt_ivort(:)

  !> Wake particles pointer
  type(t_vortpart_p), allocatable :: part_p(:)

  !> Bounding box
  real(wp) :: part_box_min(3), part_box_max(3)

  type(t_vort_elem_p), allocatable :: vort_p(:)

  !> Last vortex intensity from removed panels
  real(wp), allocatable :: last_pan_idou(:)

  !> Last vortex intensity from removed panels
  real(wp), allocatable :: end_pan_idou(:)

  !> Are the panels full? (and so need to produce particles...)
  logical :: full_panels=.false.

  !> Are the rings full? (and so need to produce particles...)
  logical :: full_rings=.false.
  
  !> Employ wake refinement
  logical :: refine_wake

  !> Wake refinement factor
  integer :: k_refine
  
  !> Employ wake refinement
  logical :: interpolate_wake

#if USE_PRECICE  
  !> In order to correctly build the wake we need to esplicitly store
  ! the end points of the last step...
  real(wp), allocatable :: old_second_row(:,:)
  !> ... but not within precice iterations
  logical :: update_old_second_row
#endif
end type

!module variables to share among the different subroutines
  real(wp), allocatable :: points_end(:,:)
  real(wp), allocatable :: points_end_ring(:,:)

!> Class to change methods from different wake implementations
type, abstract :: c_wake_mov
  contains
  procedure(i_get_vel), deferred, pass(this) :: get_vel
end type

abstract interface
  subroutine i_get_vel(this, elems, wake, pos, vel_hcas, vel)
    import                                :: c_wake_mov, wp, t_pot_elem_p, t_wake
    class(c_wake_mov)                     :: this
    type(t_pot_elem_p), intent(in)        :: elems(:)
    type(t_wake), intent(in)              :: wake
    real(wp), intent(in)                  :: pos(3)
    real(wp), intent(in)                  :: vel_hcas(3)
    real(wp), intent(out)                 :: vel(3)
  end subroutine
end interface

type, extends(c_wake_mov) :: t_free_wake
contains
  procedure, pass(this) :: get_vel => get_vel_free
end type

type, extends(c_wake_mov) :: t_rigid_wake
contains
  procedure, pass(this) :: get_vel => get_vel_rigid
end type

class(c_wake_mov), allocatable  :: wake_movement
character(len=max_char_len)     :: msg
real(t_realtime)                :: t1 , t0
character(len=*), parameter     :: this_mod_name='mod_wake'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the panel wake
subroutine initialize_wake(wake, geo, te,  npan, nrings, nparts)
  type(t_wake), intent(out),target     :: wake
  type(t_geo), intent(in), target      :: geo
  type(t_tedge), intent(inout)         :: te
  integer, intent(in)                  :: npan
  integer, intent(in)                  :: nrings
  integer, intent(in)                  :: nparts

  integer                              :: iw, ip, nsides
  integer                              :: ic, nad, ie, npt, id, ir, nend
  integer                              :: p1, p2
  real(wp)                             :: dist(3) , vel_te(3), wind(3)

  if (sim_param%rigid_wake) then
    allocate(t_rigid_wake::wake_movement)
  else
    allocate(t_free_wake::wake_movement)
  endif

  !panel wake: set and allocate all the relevant variables
  wake%nmax_pan = npan
  wake%n_pan_stripes = size(te%e,2)
  wake%n_pan_points  = size(te%i,2)

  allocate(wake%pan_gen_elems(2,wake%n_pan_stripes))
  allocate(wake%pan_gen_elems_id(2,wake%n_pan_stripes))
  allocate(wake%pan_gen_points(2,wake%n_pan_points))
  allocate(wake%pan_gen_dir(3,wake%n_pan_points))
  allocate(wake%pan_gen_ref(wake%n_pan_points))
  allocate(wake%pan_gen_icomp(wake%n_pan_points))
  allocate(wake%pan_gen_scaling(wake%n_pan_points))
  allocate(wake%w_start_points(3,wake%n_pan_points))
  allocate(wake%i_start_points(2,wake%n_pan_stripes))
  allocate(wake%pan_neigh(2,wake%n_pan_stripes))
  allocate(wake%pan_neigh_o(2,wake%n_pan_stripes))
  allocate(wake%pan_w_points(3,wake%n_pan_points,npan+1))
  allocate(wake%pan_w_vel(   3,wake%n_pan_points,npan+1))
  allocate(wake%w_vel(3,wake%n_pan_points,npan+1))
  allocate(wake%wake_panels(wake%n_pan_stripes,npan))
  allocate(wake%end_vorts(wake%n_pan_stripes))
  allocate(wake%pan_idou(wake%n_pan_stripes,npan))

  !ring wake: count the number of actuator disks and points
  nad = 0; npt = 0;
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      nad = nad + geo%components(ic)%nelems
      do ie = 1, size(geo%components(ic)%el)
        npt = npt + geo%components(ic)%el(ie)%n_ver
      enddo
    endif
  enddo

  !ring wake: set and allocate all the relevant variables
  wake%nmax_rin = nrings
  wake%ndisks = nad
  wake%np_row = npt
  allocate(wake%rin_gen_elems(wake%ndisks))
  allocate(wake%wake_rings(wake%ndisks,wake%nmax_rin))
  allocate(wake%rin_idou(wake%ndisks,wake%nmax_rin))

  !panel wake: Associate for each panel the relevant intensity
  ! and allocate all the relevant fields
  nsides = 4
  do ip = 1,npan
    do iw=1,wake%n_pan_stripes
      wake%wake_panels(iw,ip)%mag => wake%pan_idou(iw,ip)
      wake%wake_panels(iw,ip)%n_ver = nsides
      allocate(wake%wake_panels(iw,ip)%ver(3,nsides))
      allocate(wake%wake_panels(iw,ip)%edge_vec(3,nsides))
      allocate(wake%wake_panels(iw,ip)%edge_len(nsides))
      allocate(wake%wake_panels(iw,ip)%edge_uni(3,nsides))
    enddo
  enddo

  !ring wake: associate the generating elements
  id = 1
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      do ie = 1,size(geo%components(ic)%el)
        wake%rin_gen_elems(id)%p => geo%components(ic)%el(ie)
        id = id+1
      enddo
    endif
  enddo

  !ring wake: allocate the relevant fields of each ring
  do id = 1,wake%ndisks
    do ir = 1,wake%nmax_rin
      wake%wake_rings(id,ir)%mag => wake%rin_idou(id,ir)
      nsides = wake%rin_gen_elems(id)%p%n_ver
      wake%wake_rings(id,ir)%n_ver = nsides
      allocate(wake%wake_rings(id,ir)%ver(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_vec(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_len(nsides))
      allocate(wake%wake_rings(id,ir)%edge_uni(3,nsides))
    enddo
    wake%wake_rings(id,:)%moving = wake%rin_gen_elems(id)%p%moving
  enddo
  !Starting length of the ring wake is
  wake%rin_wake_len = 0
  allocate(wake%rin_p(0))

  !panel wake: set the generating elements and points from the trailing edge
  wake%pan_gen_elems  = te%e
  !do iw=1,size(te%e,2)
  !!TODO: consider shifting the whole t.e. to c_pot_elem to avoid this nightmare
  !select type(el => te%e(1,iw)%p); class is(c_pot_elem)
  !  wake%gen_elems(1,iw)%p => el
  !end select
  !select type(el => te%e(2,iw)%p); class is(c_pot_elem)
  !  wake%gen_elems(2,iw)%p => el
  !end select
  !enddo

  ! FIX FOR TE WITH HINGES

  wake%pan_gen_points = te%i
  wake%pan_gen_dir = te%t_hinged
  wake%pan_gen_ref = te%ref
  wake%pan_gen_icomp = te%icomp
  wake%pan_gen_scaling = te%scaling
  wake%i_start_points = te%ii
  wake%pan_neigh = te%neigh
  wake%pan_neigh_o = te%o
  nend = 0

  do iw = 1,wake%n_pan_stripes
    wake%pan_gen_elems_id(1,iw) = wake%pan_gen_elems(1,iw)%p%id
    if(associated(wake%pan_gen_elems(2,iw)%p)) then
      wake%pan_gen_elems_id(2,iw) = wake%pan_gen_elems(2,iw)%p%id
    else
      wake%pan_gen_elems_id(2,iw) = 0
    endif
    if(any(wake%pan_neigh(:,iw).eq.0)) nend = nend+1
  enddo

  if (sim_param%join_te) then
    !Store the indices of the ends of the panel wakes
    allocate(wake%pan_i_ends(nend))
    allocate(wake%joined_tes(2,nend,npan+1)); wake%joined_tes = 0
    ie = 1
    do iw=1,wake%n_pan_stripes
      if(any(wake%pan_neigh(:,iw).eq.0)) then
        wake%pan_i_ends(ie) = iw
        ie = ie+1
      endif
    enddo
  endif

  wake%pan_idou = 0.0_wp
  wake%rin_idou = 0.0_wp

  !panel wake: The whole first line of (implicit) panels is initialized
  !here to be used during linear system initialization

  !the first line of points is calculated from the mesh points
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%pan_gen_points(1,:)) + &
                                  geo%points(:,wake%pan_gen_points(2,:)))
  wake%pan_w_points(:,:,1) = wake%w_start_points
  wake%pan_w_vel(   :,:,:) = 0.0_wp

  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1 , wake%n_pan_points
  !dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))
    call calc_node_vel( wake%w_start_points(:,ip), &
            geo%refs(wake%pan_gen_ref(ip))%G_g, &
            geo%refs(wake%pan_gen_ref(ip))%f_g, &
            vel_te )

    wind = variable_wind(wake%w_start_points(:,ip), sim_param%time)

    if ( norm2(wind-vel_te) .gt. sim_param%min_vel_at_te ) then
      
      dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))

      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*wake%pan_gen_scaling(ip)* &
                  norm2(wind-vel_te)* &
                  sim_param%dt*real(sim_param%ndt_update_wake,wp) / norm2(dist)
  ! normalisation occurs here! --------------------------------------^

    else

      dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))

      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*wake%pan_gen_scaling(ip) * & ! next line may be commented
                  sim_param%min_vel_at_te* &
                  sim_param%dt*real(sim_param%ndt_update_wake,wp)
    end if

  enddo

  !Starting length of the panel wake is 1
  wake%pan_wake_len = 1

  !Associate the first implicit row
  allocate(wake%pan_p(wake%n_pan_stripes))
  do iw = 1,wake%n_pan_stripes
    wake%wake_panels(iw,:)%moving = geo%refs(wake%pan_gen_ref(iw))%moving
    wake%pan_p(iw)%p => wake%wake_panels(iw,1)
  enddo

  !set the end vortexes to null
  do iw = 1,wake%n_pan_stripes
    wake%end_vorts(iw)%mag => null()
  enddo

  ! Calculate geometrical quantities of first row
  do iw = 1,wake%n_pan_stripes
    p1 = wake%i_start_points(1,iw)
    p2 = wake%i_start_points(2,iw)
    call wake%wake_panels(iw,1)%calc_geo_data( &
        reshape((/wake%pan_w_points(:,p1,1),   wake%pan_w_points(:,p2,1), &
                  wake%pan_w_points(:,p2,1+1), wake%pan_w_points(:,p1,1+1)/),&
                                                                  (/3,4/)))
  enddo

  if (sim_param%join_te) call join_first_panels(wake,sim_param%join_te_factor)


  !Particles
  wake%nmax_prt = nparts
  allocate(wake%wake_parts(wake%nmax_prt))
  allocate(wake%prt_ivort(wake%nmax_prt))
  wake%n_prt = 0
  allocate(wake%part_p(0))
  do ip = 1,wake%nmax_prt
    wake%wake_parts(ip)%mag => wake%prt_ivort(ip)
  enddo
  wake%part_box_min = sim_param%particles_box_min
  wake%part_box_max = sim_param%particles_box_max

  allocate(wake%vort_p(0))
  allocate(wake%last_pan_idou(wake%n_pan_stripes))
  allocate(wake%end_pan_idou(wake%n_pan_stripes))
  wake%last_pan_idou = 0.0_wp

  wake%full_panels = .false.
  wake%refine_wake = sim_param%refine_wake
  wake%k_refine = sim_param%k_refine
  wake%interpolate_wake = sim_param%interpolate_wake
#if USE_PRECICE
  allocate(wake%old_second_row(3,wake%n_pan_points)) ! check for memory leak
  wake%update_old_second_row = .true. ! set to false inside precice iters
#endif

end subroutine initialize_wake

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_wake(wake)
  type(t_wake), intent(out) :: wake

  !dummy to avoid compiler warnings
  wake%nmax_pan = -1
  if(allocated(points_end)) deallocate(points_end)
  if(allocated(points_end_ring)) deallocate(points_end_ring)

end subroutine

!----------------------------------------------------------------------

!> Load the wake panels solution from a previous result
subroutine load_wake(filename, wake, elems)
  character(len=*), intent(in)          :: filename
  type(t_wake), intent(inout), target   :: wake
  type(t_pot_elem_p), intent(inout)     :: elems(:)

  integer(h5loc)                        :: floc, gloc
  real(wp), allocatable                 :: wpoints(:,:,:), wvels(:,:,:), wvort(:,:)
  real(wp), allocatable                 :: vppoints(:,:), vpvort(:,:) , vpvels(:,:), vprad(:)
  integer, allocatable                  :: start_points(:,:)
  integer, allocatable                  :: conn_pe(:)
  integer                               :: ipan, iw, p1, p2, ipt
  integer                               :: id, ir, ip, np
  character(len=*), parameter           :: this_sub_name = 'load_wake'

  call open_hdf5_file(filename, floc)

  !=== Panels ===
  !Read the past results
  call open_hdf5_group(floc, 'PanelWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wvels  ,'WakeVels'  ,gloc)
  call read_hdf5_al(start_points,'StartPoints',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  !call close_hdf5_file(floc)

  !Perform a few checks to be sure that the correct size solution was loaded
  if(.not. all(start_points .eq. wake%i_start_points)) call error( &
    this_sub_name, this_mod_name, 'Different wake trailing edge connectivity&
    & between the loded and built geometry')

  !the wake length is the loaded one, or less if imposed so
  wake%pan_wake_len = min(wake%nmax_pan, size(wvort,2))

  !store points position and doublets intensity
  wake%pan_w_points(:,:,1:wake%pan_wake_len+1) = wpoints(:,:,1:wake%pan_wake_len+1)
  wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1) = wvels(  :,:,1:wake%pan_wake_len+1)
  wake%pan_idou(:,1:wake%pan_wake_len) = wvort(:,1:wake%pan_wake_len)

  deallocate(wake%pan_p)
  allocate(wake%pan_p(wake%n_pan_stripes*wake%pan_wake_len))
  ! Update the panels geometrical quantities of all the panels
  ipt = 1
  do ipan = 1,wake%pan_wake_len
    do iw = 1,wake%n_pan_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%wake_panels(iw,ipan)%calc_geo_data( &
            reshape((/wake%pan_w_points(:,p1,ipan),   wake%pan_w_points(:,p2,ipan), &
                    wake%pan_w_points(:,p2,ipan+1), wake%pan_w_points(:,p1,ipan+1)/),&
                                                                      (/3,4/)))
      wake%pan_p(ipt)%p => wake%wake_panels(iw,ipan)
      ipt = ipt + 1
    enddo
  enddo

  deallocate(wpoints, wvels, wvort, start_points)

  !=== Rings ===
  !Read the past results
  !call open_hdf5_file(filename, floc)
  call open_hdf5_group(floc, 'RingWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(conn_pe,'Conn_pe',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  !call close_hdf5_file(floc)

  !check the consistency
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    if (.not. all(conn_pe(ip:ip+np-1) .eq. id)) call error( &
    this_sub_name, this_mod_name, 'Different wake disk connectivity&
    & between the loded and built geometry')
    ip = ip+np
  enddo

  !the wake length is the loaded one, or less if imposed so
  wake%rin_wake_len = min(wake%nmax_rin, size(wvort,2))

  !Load the old results
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    do ir = 1,wake%rin_wake_len
      wake%wake_rings(id,ir)%ver(:,:) = wpoints(:,ip:ip+np-1,ir)
    enddo
    ip = ip+np
  enddo
  wake%rin_idou(:,1:wake%rin_wake_len) = wvort(:,1:wake%rin_wake_len)

  deallocate(wake%rin_p); allocate(wake%rin_p(wake%ndisks*wake%rin_wake_len))

  !Update the geometrical quantities
  ipt = 1
  do ir = 1,wake%rin_wake_len
    do id = 1,wake%ndisks
      call wake%wake_rings(id,ir)%calc_geo_data( &
                    wake%wake_rings(id,ir)%ver)
      wake%rin_p(ipt)%p => wake%wake_rings(id,ir)
      ipt = ipt + 1
    enddo
  enddo

  !=== Particles ===
  call open_hdf5_group(floc, 'ParticleWake', gloc)
  call read_hdf5_al(vppoints,'WakePoints',gloc)
  if ( check_dset_hdf5('WakeVels',gloc) ) then
    call read_hdf5_al(vpvels,'WakeVels',gloc) ! <<<< restart with Bernoulli integral equation
  end if
  call read_hdf5_al(vpvort,'WakeVort',gloc)
  call read_hdf5(wake%last_pan_idou,'LastPanIdou',gloc)
  call read_hdf5_al(vprad,'VortexRad',gloc)
  call close_hdf5_group(gloc)

  if(size(vpvort,2) .gt. wake%nmax_prt) call error(this_sub_name, &
    this_mod_name, 'Loading a number of particles &
    & greater than the maximum allowed for this run: cannot truncate the &
    & particles in a meaningful way. Consider restarting the run with a &
    & higher amount of maximum particles')

  wake%n_prt = size(vpvort,2)

  deallocate(wake%part_p)
  allocate(wake%part_p(wake%n_prt))
  if(wake%n_prt .gt. 0) then
    deallocate(wake%vort_p)
    if(sim_param%use_fmm_pan) then ! particles treated in FMM
      allocate(wake%vort_p(wake%n_pan_stripes))
    else ! particles velocity calculated alongside the end vortices
      allocate(wake%vort_p(wake%n_prt+wake%n_pan_stripes))
    endif

    do iw = 1, wake%n_pan_stripes
      if(sim_param%use_fmm_pan) then
        wake%vort_p(iw)%p => wake%end_vorts(iw)
      else
        wake%vort_p(wake%n_prt+iw)%p => wake%end_vorts(iw)
      endif
    enddo
  endif

  do ip = 1,wake%n_prt
    wake%wake_parts(ip)%cen = vppoints(:,ip)
    wake%wake_parts(ip)%vel = vpvels(:,ip)
    wake%wake_parts(ip)%mag = norm2(vpvort(:,ip))
    wake%wake_parts(ip)%r_Vortex = vprad(ip)
    if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
      wake%wake_parts(ip)%dir = vpvort(:,ip)/wake%wake_parts(ip)%mag
    else
      wake%wake_parts(ip)%dir = vpvort(:,ip)
    endif
    wake%wake_parts(ip)%free = .false.
    wake%part_p(ip)%p => wake%wake_parts(ip)
    if(.not.sim_param%use_fmm_pan) wake%vort_p(ip)%p => wake%wake_parts(ip)
  enddo

  deallocate(vppoints, vpvort, vpvels, vprad)

  call close_hdf5_file(floc)


  if(wake%pan_wake_len .eq. wake%nmax_pan) wake%full_panels=.true.
  !If the wake is full, attach the end vortex
  if (wake%full_panels) then
    do iw = 1,wake%n_pan_stripes
      wake%end_vorts(iw)%mag => wake%last_pan_idou(iw)
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%end_vorts(iw)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,wake%pan_wake_len+1),  &
                    wake%pan_w_points(:,p2,wake%pan_wake_len+1)/), (/3,2/)))
    enddo
  endif

!!!  !If there are particles and the multipole is employed, i
!!!  !need to pre-compute the particles induced velocity now
!!!  if(sim_param%use_fmm .and. wake%n_prt .gt. 0) then
!!!!!$omp parallel do private(ie, ip, vel)
!!!    do ie = 1, size(elems)
!!!      elems(ie)%p%uvort = 0.0
!!!      do ip = 1,wake%n_prt
!!!        call wake%part_p(ip)%p%compute_vel(elems(ie)%p%cen, &
!!!                                               sim_param%u_inf, vel)
!!!         elems(ie)%p%uvort = elems(ie)%p%uvort + vel/(4.0_wp*pi)
!!!      enddo
!!!    enddo
!!!!!$omp end parallel do
!!!  endif


end subroutine load_wake

!----------------------------------------------------------------------

!> Prepare the wake before the timestep
!!
!! Mainly prepare all the structures for the octree
subroutine prepare_wake(wake, elems, octree)
  type(t_wake), intent(inout), target   :: wake
  type(t_pot_elem_p), intent(inout)     :: elems(:)
  type(t_octree), intent(inout)         :: octree

  integer                               :: k, ip, ir, iw, ie, n_end_vort

  !reset all the vorticity induced velocity
  do ie = 1,size(elems)
    elems(ie)%p%uvort = 0.0_wp
  enddo

  if (sim_param%use_fmm) then
    call sort_particles(wake%wake_parts, wake%n_prt, elems, octree)
    call calculate_multipole(wake%part_p, octree)
    call apply_multipole_panels(octree, elems)
  endif

  !==>Recreate sturctures and pointers, if particles are present
  if(wake%full_panels .or. wake%full_rings .or. (wake%n_prt.gt.0) ) then

    !Recreate the pointer vector
    if(allocated(wake%part_p)) then 
      deallocate(wake%part_p)
    endif

    allocate(wake%part_p(wake%n_prt))
    deallocate(wake%vort_p)
    
    ! to add or not line vortices at the (only when ring or panel wakes are full )
    n_end_vort = 0
    if ( wake%full_panels .or. wake%full_rings ) then
      n_end_vort = wake%n_pan_stripes
    end if

    if(sim_param%use_fmm_pan) then ! particles vel treated in FMM
      allocate(wake%vort_p( n_end_vort))
    else !particles velocity treated alongside the end vortices
      allocate(wake%vort_p(wake%n_prt + n_end_vort))
    endif

    !TODO: consider inverting these two cycles
    k = 1
    do ip = 1, wake%n_prt
      do ir=k,wake%nmax_prt
        if(.not. wake%wake_parts(ir)%free) then
          k = ir+1
          wake%part_p(ip)%p => wake%wake_parts(ir)
          if(.not.sim_param%use_fmm_pan) wake%vort_p(ip)%p => wake%wake_parts(ir)
          exit
        endif
      enddo
    enddo

    !Add the end vortex to the votical elements pointer
    do iw = 1, n_end_vort
      if(sim_param%use_fmm_pan) then
        wake%vort_p(iw)%p => wake%end_vorts(iw)
      else
        wake%vort_p(wake%n_prt+iw)%p => wake%end_vorts(iw)
      endif
    enddo
  endif

end subroutine prepare_wake

!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels 
!  Brings them to the next time step
!  Only updates the "old" panels, ie not the first two, and existing particles
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll)
!! elements
subroutine update_wake(wake, elems, octree)
  type(t_wake), intent(inout), target :: wake
  type(t_pot_elem_p), intent(in)      :: elems(:)
  type(t_octree), intent(inout)       :: octree

  integer                             :: iw, ipan, ie, ip, np, iq
  integer                             :: id, ir
  real(wp)                            :: pos_p(3), vel_p(3)
  real(wp)                            :: str(3), stretch(3)
  real(wp)                            :: ru(3), rotu(3)
  real(wp)                            :: df(3), diff(3)
  real(wp)                            :: hcas_vel(3)
  type(t_pot_elem_p), allocatable     :: pan_p_temp(:)
  real(wp), allocatable               :: point_old(:,:,:)
  real(wp), allocatable               :: points(:,:,:)
  logical                             :: increase_wake
  integer                             :: size_old
  character(len=*), parameter         :: this_sub_name='update_wake'

  wake%w_vel = 0.0_wp
  if(sim_param%HCAS) hcas_vel = get_vel_hcas()

  if(wake%pan_wake_len .eq. wake%nmax_pan) wake%full_panels=.true.

  !==> Panels:  Update the first row of vortex intensities:
  !      it was already calculated (implicitly) in the linear system
  do iw = 1,wake%n_pan_stripes
    
    if      ( associated(wake%pan_gen_elems(2,iw)%p) ) then
      wake%wake_panels(iw,1)%mag  = wake%pan_gen_elems(1,iw)%p%mag - &
                                    wake%pan_gen_elems(2,iw)%p%mag
    else if ( .not. associated(wake%pan_gen_elems(2,iw)%p) ) then
      wake%wake_panels(iw,1)%mag  = wake%pan_gen_elems(1,iw)%p%mag
    end if
  enddo

  !==> Panels: Update wake points position ==

  !Save the old positions for the integration
  allocate(point_old(size(wake%pan_w_points,1),size(wake%pan_w_points,2), &
                                              size(wake%pan_w_points,3)))
  point_old = wake%pan_w_points

#if USE_PRECICE
  ! pan_w_points is overwritten at every precice iteration, the actual old
  ! points have been saved to old_second_row (the old first row is useless)
  point_old(:,:,2) = wake%old_second_row
#endif

  !calculate the velocities at the old positions of the points and then
  !update the positions (from the third row of points: the first is the
  !trailing edge, the second is extrapolated from the trailing edge)
  np = wake%pan_wake_len+1
  if(.not.wake%full_panels) then
    np = np + 1
  endif

!$omp parallel do collapse(2) private(pos_p, vel_p, ie, ipan, iw) schedule(dynamic)
  do ipan = 3,np
    do iw = 1,wake%n_pan_points
      pos_p = point_old(:,iw,ipan-1)
      vel_p = 0.0_wp

      ! for OUTPUT only -----
      call wake_movement%get_vel(elems, wake, pos_p, hcas_vel, vel_p)

      !update the position in time
      wake%pan_w_vel(   :,iw,ipan) = vel_p
      if ( ipan .ne. 3 ) then
        wake%pan_w_points(:,iw,ipan) = point_old(:,iw,ipan-1) + &
                                      vel_p*sim_param%dt*real(sim_param%ndt_update_wake,wp)
      else
        wake%pan_w_points(:,iw,ipan) = point_old(:,iw,ipan-1) + &
                                      vel_p*sim_param%dt
      end if
    enddo
  enddo
!$omp end parallel do

  !if the wake is full, calculate another row of points to generate the
  !particles
  if(wake%full_panels) then
    if(.not.allocated(points_end)) allocate(points_end(3,wake%n_pan_points))

  ! create another row of points
!$omp parallel do private(iw, pos_p, vel_p) schedule(dynamic)
    do iw = 1,wake%n_pan_points
      pos_p = point_old(:,iw,wake%pan_wake_len+1)
      vel_p = 0.0_wp

      call wake_movement%get_vel(elems, wake, pos_p, hcas_vel, vel_p)

      !update the position in time
      points_end(:,iw) = pos_p + vel_p*sim_param%dt*real(sim_param%ndt_update_wake,wp)
    enddo
!$omp end parallel do
  endif


  deallocate(point_old)

  !==> Rings: Update wake points position ==

  !Save the old positions for the integration
  increase_wake = .false.
  if(wake%rin_wake_len .lt. wake%nmax_rin) then
    wake%rin_wake_len = wake%rin_wake_len+1
    increase_wake = .true.
  else
    wake%full_rings = .true.
  endif
  allocate(points(3,wake%np_row,wake%rin_wake_len))

  !Store at the beginning the disk points
  ip=1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    points(:,ip:ip+np-1,1) = wake%rin_gen_elems(id)%p%ver
    ip = ip+np
  enddo

  !Then store the old points of the rest of the wake (the shift forward
  ! of the points in the array is happening here)
  do ir = 1,wake%rin_wake_len-1
    ip=1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      points(:,ip:ip+np-1,ir+1) = wake%wake_rings(id,ir)%ver(:,:)
      ip = ip+np
    enddo
  enddo

  !calculate the velocities at the old positions of the points and then
  !update the positions

!$omp parallel do private(pos_p, vel_p, ip, ir)
  do ip = 1,size(points,2)
    do ir = 1,size(points,3)
      pos_p = points(:,ip,ir)
      vel_p = 0.0_wp

      call wake_movement%get_vel(elems, wake, pos_p, hcas_vel, vel_p)

      !update the position in time
      points(:,ip,ir) = points(:,ip,ir) + &
                        vel_p*sim_param%dt*real(sim_param%ndt_update_wake,wp)
    enddo !ir
  enddo !ip
!$omp end parallel do

  !if the wake is full, calculate another row of points to generate the
  !particles
  if(wake%full_rings) then
    if(.not.allocated(points_end_ring)) allocate(points_end_ring(3,wake%np_row))
    !Store the end points
    ip=1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      points_end_ring(:,ip:ip+np-1) = wake%wake_rings(id,ir)%ver(:,:)
      ip = ip+np
    enddo
    ! calculate velocity and evolve them
    do ip = 1,wake%np_row
      pos_p = points_end_ring(:,ip)
      vel_p = 0.0_wp

      call wake_movement%get_vel(elems, wake, pos_p, hcas_vel,  vel_p)

      !update the position in time
      points_end_ring(:,ip) = pos_p + &
                              vel_p*sim_param%dt*real(sim_param%ndt_update_wake,wp)
    enddo
  endif

  !redistribute the points
  do ir = 1,wake%rin_wake_len
    ip = 1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      wake%wake_rings(id,ir)%ver = points(:,ip:ip+np-1,ir)
      ip = ip + np
    enddo
  enddo

  deallocate(points)


  !==>    Particles: evolve the position in time

  !calculate the velocities at the points
!$omp parallel do private(pos_p, vel_p, ip, iq,  stretch, diff, df, str, ru, rotu)
  do ip = 1, wake%n_prt
    wake%part_p(ip)%p%vel_old = wake%part_p(ip)%p%vel
    wake%part_p(ip)%p%stretch_old = wake%part_p(ip)%p%stretch
    wake%part_p(ip)%p%stretch = 0.0_wp
    wake%part_p(ip)%p%rotu = 0.0_wp

    !If not using the fast multipole, update particles position now
    if (.not.sim_param%use_fmm) then
      pos_p = wake%part_p(ip)%p%cen

      call wake_movement%get_vel(elems, wake, pos_p, hcas_vel, vel_p)

      wake%part_p(ip)%p%vel =  vel_p

      !if using vortex stretching, calculate it now
      if(sim_param%use_vs) then
        stretch = 0.0_wp
        rotu = 0.0_wp
        do iq = 1, wake%n_prt
          if (ip.ne.iq) then
            call wake%part_p(iq)%p%compute_stretch(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, str)
            ! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            stretch = stretch + str/(4.0_wp*pi)

            if(sim_param%use_divfilt) then
              call wake%part_p(iq)%p%compute_rotu(wake%part_p(ip)%p%cen, &
                    wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, ru)
              rotu = rotu + ru/(4.0_wp*pi)

            endif
          endif
        enddo
      
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + stretch
        
        if(sim_param%use_divfilt) then 
          wake%part_p(ip)%p%rotu = wake%part_p(ip)%p%rotu + rotu
        endif
      
      endif !use_vs

      !if using the vortex diffusion, calculate it now
      if(sim_param%use_vd) then
        diff = 0.0_wp

        do iq = 1, wake%n_prt

          if (ip.ne.iq) then
            call wake%part_p(iq)%p%compute_diffusion(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, &
                  wake%part_p(ip)%p%r_Vortex, df)
            diff = diff + df*sim_param%nu_inf
          endif

        enddo !iq
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + diff
      endif !use_vd

    end if !use_fmm

  enddo
!$omp end parallel do

  if (sim_param%use_fmm) then
    t0 = dust_time()
    call apply_multipole(wake%part_p, octree, elems, wake%pan_p, wake%rin_p, &
                          wake%end_vorts)
    if (sim_param%HCAS) then
      do ip = 1, wake%n_prt
        wake%part_p(ip)%p%vel = wake%part_p(ip)%p%vel + hcas_vel
      enddo
    endif
    t1 = dust_time()
    write(msg,'(A,F9.3,A)') 'Multipoles calculation: ' , t1 - t0,' s.'
    if(sim_param%debug_level.ge.3) call printout(msg)
    write(msg,'(A,I0)') 'Number of particles: ' , wake%n_prt
    if(sim_param%debug_level.ge.5) call printout(msg)
  endif

  !==> Panels:  Increase the length of the wake, if it is necessary
  if (.not.wake%full_panels) then
      wake%pan_wake_len = wake%pan_wake_len + 1
      allocate(pan_p_temp(wake%n_pan_stripes*wake%pan_wake_len))

      do ip = 1,size(wake%pan_p)
        pan_p_temp(ip) = wake%pan_p(ip)
      enddo
      do iw = 1,wake%n_pan_stripes
        pan_p_temp(wake%n_pan_stripes*(wake%pan_wake_len-1)+iw)%p &
                                        => wake%wake_panels(iw,wake%pan_wake_len)
      enddo

      if(allocated(wake%pan_p)) deallocate(wake%pan_p)
      allocate(wake%pan_p(size(pan_p_temp)))

      do ip = 1,size(wake%pan_p)
        wake%pan_p(ip) = pan_p_temp(ip)
      enddo

      deallocate(pan_p_temp)

  endif


  !==> Rings:  Increase the length of the wake, if it is necessary
  if (increase_wake) then
    allocate(pan_p_temp(wake%ndisks*wake%rin_wake_len))
    size_old = 0; if(allocated(wake%rin_p)) size_old = size(wake%rin_p)

    do ip = 1,size_old
      pan_p_temp(ip) = wake%rin_p(ip)
    enddo
    do id = 1,wake%ndisks
      pan_p_temp(size_old+id)%p => wake%wake_rings(id,wake%rin_wake_len)
    enddo

    !manual move alloc
    if(allocated(wake%rin_p)) deallocate(wake%rin_p)
    allocate(wake%rin_p(size(pan_p_temp)))
    do ip = 1,size(wake%rin_p)
      wake%rin_p(ip) = pan_p_temp(ip)
    enddo
    deallocate(pan_p_temp)
  endif

  !==> Panels:  Update the intensities of the panels
  !       From the back, all the vortex intensities come from
  !       the previous panel

  !QUICK AND DIRTY
  do iw = 1,wake%n_pan_stripes
    wake%end_pan_idou(iw) = wake%wake_panels(iw,wake%pan_wake_len)%mag
  enddo
  do ipan = wake%pan_wake_len,2,-1
    do iw = 1,wake%n_pan_stripes
      wake%wake_panels(iw,ipan)%mag = wake%wake_panels(iw,ipan-1)%mag
    enddo
  enddo

  !==> Rings: Update the intensities of the rings
  !       From the back, all the vortex intensities come from the
  !       previous panel
  do ir = wake%rin_wake_len,2,-1
    do id = 1,wake%ndisks
      wake%wake_rings(id,ir)%mag = wake%wake_rings(id,ir-1)%mag
    enddo
  enddo
  do id = 1,wake%ndisks
    wake%wake_rings(id,1)%mag  = wake%rin_gen_elems(id)%p%mag
  enddo

  !==> Rings: Update rings geometrical quantities
  do ir = 1,wake%rin_wake_len
    do id = 1,wake%ndisks
      call wake%wake_rings(id,ir)%calc_geo_data( &
                    wake%wake_rings(id,ir)%ver)
    enddo
  enddo

  ! The geometrical quantities of the panels will be all updated in prepare
  ! wake after the geometry have been updated

end subroutine update_wake

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!! Completes the updating to the next time step begun in update_wake
!! first and second row are updated to the next step and new particles are
!! created if necessary; they will appear at the save_date in the next time step
subroutine complete_wake(wake, geo, elems, te)
  type(t_wake), target, intent(inout)   :: wake
  type(t_geo), intent(in)               :: geo
  type(t_pot_elem_p), intent(in)        :: elems(:)
  type(t_tedge), intent(inout)          :: te

  integer                               :: p1, p2
  integer                               :: ip, iw, ipan, id, is, nprev
  real(wp)                              :: dist(3) , vel_te(3), pos_p(3)
  real(wp)                              :: dir(3), partvec(3), ave, alpha_p(3), alpha_p_n
  integer                               :: k, n_part
  real(wp)                              :: vel_in(3), vel_out(3), wind(3)
  real(wp)                              :: area
  real(wp)                              :: min_side, max_side, chord_side_len, span_side_len, step_span, step_chord
  integer                               :: n_span, n_chord, ic
  
  integer                               :: iwc, isp, n_sbprt, n_max_pan_comp, n_max_sbprt_comp, pan_count, sbprt_count
  real(wp), allocatable                 :: cen_parent(:,:), dir_parent(:,:), mag_parent(:)
  real(wp), allocatable                 :: cen_sbprt(:,:), area_sbprt(:), mag_sbprt(:), dir_sbprt(:,:)
  real(wp), allocatable                 :: W(:,:), w_i(:)
  ! flow separation variables
  integer                                :: i_comp , i_elem , n_elem

  character(len=max_char_len)            :: msg
  character(len=*), parameter            :: this_sub_name='prepare_wake'

#if USE_PRECICE
  ! first and second row of the wake were already taken care of by update_near_field_wake
#else
  !==> Panels: update the first rows of panels

  !first row  of new points comes from geometry
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%pan_gen_points(1,:)) + &
                                  geo%points(:,wake%pan_gen_points(2,:)))
  wake%pan_w_points(:,:,1) = wake%w_start_points
  wake%pan_gen_dir = te%t_hinged

  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1,wake%n_pan_points
    dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))
    call calc_node_vel( wake%w_start_points(:,ip), &
            geo%refs(wake%pan_gen_ref(ip))%G_g, &
            geo%refs(wake%pan_gen_ref(ip))%f_g, &
            vel_te )
    wind = variable_wind(wake%w_start_points(:,ip), sim_param%time)

    if ( norm2(wind-vel_te) .gt. sim_param%min_vel_at_te ) then
      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) + &
                          dist*wake%pan_gen_scaling(ip)* &
                          norm2(wind-vel_te)* &
                  sim_param%dt*real(sim_param%ndt_update_wake,wp) / norm2(dist)
  ! normalisation occurs here! -------------------------------------------^

    else
      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*wake%pan_gen_scaling(ip)* & ! next line may be commented
                  sim_param%min_vel_at_te* &
                  sim_param%dt*real(sim_param%ndt_update_wake,wp) / norm2(dist)
    end if

  enddo

  !==> Check if the panels need to be joined
  if(sim_param%join_te) then
    wake%joined_tes(:,:,2:wake%nmax_pan+1) = &
          wake%joined_tes(:,:,1:wake%nmax_pan)
    wake%joined_tes(:,:,1) = 0
    call join_first_panels(wake,sim_param%join_te_factor)
  endif

  !==> Panels:  update the panels geometrical quantities of all the panels,
  !      the first two row of points have just been updated, the other
  !      rows of points were updated at the end of the last iteration
  do ipan = 1,wake%pan_wake_len
    do iw = 1,wake%n_pan_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%wake_panels(iw,ipan)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,ipan),   wake%pan_w_points(:,p2,ipan), &
                    wake%pan_w_points(:,p2,ipan+1), wake%pan_w_points(:,p1,ipan+1)/),&
                                                                      (/3,4/)))
    enddo
  enddo
#endif

  !==> Particles: update the position and intensity in time, avoid penetration
  !               and chech if remain into the boundaries
  n_part = wake%n_prt
!$omp parallel do schedule(dynamic,4) private(ip,pos_p,alpha_p,alpha_p_n,vel_in,vel_out)
  do ip = 1, n_part

    if(sim_param%use_pa) then
      vel_in = wake%part_p(ip)%p%vel
      call avoid_collision(elems, wake, &
                        wake%part_p(ip)%p, vel_in, vel_out)
      wake%part_p(ip)%p%vel = vel_out
    endif

    if(.not. wake%part_p(ip)%p%free) then
      pos_p = wake%part_p(ip)%p%cen + wake%part_p(ip)%p%vel* &
              sim_param%dt*real(sim_param%ndt_update_wake,wp)

      if(all(pos_p .ge. wake%part_box_min) .and. &
          all(pos_p .le. wake%part_box_max)) then
        wake%part_p(ip)%p%cen = pos_p

        if(sim_param%use_vs .or. sim_param%use_vd) then

          !add filtering
          if(sim_param%use_divfilt) then
            wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch - &
              sim_param%filt_eta/real(sim_param%ndt_update_wake,wp)*( wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag - &
              wake%part_p(ip)%p%rotu*wake%part_p(ip)%p%mag/norm2(wake%part_p(ip)%p%rotu))
          endif

          !Explicit Euler
          alpha_p = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag + &
                          wake%part_p(ip)%p%stretch* &
                          sim_param%dt*real(sim_param%ndt_update_wake,wp)
          alpha_p_n = norm2(alpha_p)

! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
          if(alpha_p_n .ne. 0.0_wp) &
              wake%part_p(ip)%p%dir = alpha_p/alpha_p_n
        endif
      else
        wake%part_p(ip)%p%free = .true.
!$omp atomic update
        wake%n_prt = wake%n_prt -1
!$omp end atomic
      endif
    endif
  enddo
!$omp end parallel do

  !==> Particles: if the panel wake is at the end, create a particle
  if(wake%full_panels) then
    if(wake%interpolate_wake) then ! TODO crosscheck interp_parts .and. refine_wake
    ! WAKE INTERPOLATION
    ! general idea:
    ! - for each wake panel collect data on cen and mag
    ! - for each wake panel, determine division in subparts
    ! - RBF interpolation
    ! - particle creation and insertion
    !
    ! The previous wake row is also considered
    ! Ghost panels are also considered, ie points outside the wake where there is no vorticity (approx) TODO check validity
    !
    ! The process is done on a component based loop, with the index iwc resetting after each component
    
    !TODO consider moving to a separate subroutine and/or make dedicated sub for left+right+end side computations
      k = 1
      n_sbprt = 0 ! number of subparticles to insert
      n_max_pan_comp = 0 ! number of wake panels in the biggest component, for allocation
      n_max_sbprt_comp = 0 ! number of subparts from the biggest component, for allocation
      pan_count = 0 ! aux counter
      sbprt_count = 0 ! aux counter
      ! cycle all panels to find n_sbprt
      do iw = 1,wake%n_pan_stripes
        p1 = wake%i_start_points(1,iw)
        p2 = wake%i_start_points(2,iw)
        !Left side
        dir = wake%pan_w_points(:,p1,wake%nmax_pan+1)-points_end(:,p1)
        chord_side_len = norm2(dir)
        !End side
        dir = points_end(:,p1) - points_end(:,p2)
        span_side_len = norm2(dir)
        
        max_side = max(chord_side_len, span_side_len)
        min_side = min(chord_side_len, span_side_len)
        n_chord = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
        n_span = wake%k_refine 
        n_sbprt = n_sbprt + n_chord*n_span
        
        if (wake%pan_neigh(1,iw) .gt. 0) then ! increase counters
          pan_count = pan_count+1
          sbprt_count = sbprt_count + n_chord*n_span
        else ! increase, check max and reset
          pan_count = pan_count+1
          sbprt_count = sbprt_count + n_chord*n_span
          n_max_pan_comp = max(n_max_pan_comp, pan_count)
          n_max_sbprt_comp = max(n_max_sbprt_comp, sbprt_count)
          pan_count = 0
          sbprt_count = 0
        endif
      enddo
      
      ! the following quantities are accumulated only within one component, then particles are generated
      ! and the cycle repeats for the next component. Since we don't want to keep track of components, we
      ! allocate for the max size, ie for the biggest component
      allocate(cen_sbprt(3,n_max_sbprt_comp)) ! TODO take max
      allocate(area_sbprt(n_max_sbprt_comp))
      allocate(mag_sbprt(n_max_sbprt_comp))
      allocate(dir_sbprt(3,n_max_sbprt_comp))
      ! +2 accounts for ghost panels
      ! *2 because we also take the previous row of panels
      allocate(cen_parent(3,(n_max_pan_comp+2)*2)) 
      allocate(mag_parent((n_max_pan_comp+2)*2))
      allocate(dir_parent(3,(n_max_pan_comp+2)*2))
      
      ! cycle again all panels to compute ave and dir, cen and cen_ref
      iw = 1 ! global parent panel index
      
      do while (iw .lt. wake%n_pan_stripes)
        iwc = 1 ! local component parent panel index
        isp = 0 ! subparts index
        
        ! ghost and first panel
        p1 = wake%i_start_points(1,iw)
        p2 = wake%i_start_points(2,iw)     
        partvec = 0.0_wp
        !Left side
        dir = wake%pan_w_points(:,p1,wake%nmax_pan+1)-points_end(:,p1)
        chord_side_len = norm2(dir)
        if (wake%pan_neigh(1,iw) .gt. 0) then
          ave = wake%end_pan_idou(iw) - &
                real(wake%pan_neigh_o(1,iw),wp)* &
                wake%end_pan_idou(wake%pan_neigh(1,iw))
          ave = ave/2.0_wp
        else !has no fixed neighbour
          if(sim_param%join_te) then
            ave = get_joined_intensity(wake, iw, 1)
          else
            ave = wake%end_pan_idou(iw)
          endif
        endif
        partvec = partvec + dir*ave
  
        !Right side
        dir = -wake%pan_w_points(:,p2,wake%nmax_pan+1)+points_end(:,p2)
        if (wake%pan_neigh(2,iw) .gt. 0) then
          ave = wake%end_pan_idou(iw) - &
                real(wake%pan_neigh_o(2,iw),wp)* &
                wake%end_pan_idou(wake%pan_neigh(2,iw))
          ave = ave/2.0_wp
        else
          if(sim_param%join_te) then
            ave = get_joined_intensity(wake, iw, 2)
          else
            ave = wake%end_pan_idou(iw)
          endif
        endif
        partvec = partvec + dir*ave
  
        !End side
        dir = points_end(:,p1) - points_end(:,p2)
        span_side_len = norm2(dir)
        ave = wake%end_pan_idou(iw)-wake%last_pan_idou(iw)
        wake%last_pan_idou(iw) = wake%end_pan_idou(iw)
        partvec = partvec + dir*ave

        ! centre of parent panel
        pos_p = (points_end(:,p1)+points_end(:,p2)+ &
                wake%pan_w_points(:,p1,wake%nmax_pan+1) + &
                wake%pan_w_points(:,p2,wake%nmax_pan+1) )/4.0_wp
                
        ! ghost panel
        cen_parent(:,iwc) = pos_p-(wake%pan_w_points(:,p1,wake%nmax_pan+1)-&
                      wake%pan_w_points(:,p2,wake%nmax_pan+1))
        mag_parent(iwc) = 0.0_wp
        dir_parent(:,iwc) = 0.0_wp
        iwc = iwc + 1 ! added ghost panel
     
        ! ghost panel from previous row
        cen_parent(:,iwc) = wake%wake_panels(iw,wake%nmax_pan)%cen-&
               (wake%pan_w_points(:,p1,wake%nmax_pan)-wake%pan_w_points(:,p2,wake%nmax_pan)) 
        mag_parent(iwc) = 0.0_wp
        dir_parent(:,iwc) = 0.0_wp
        iwc = iwc + 1 ! added ghost panel 
        
        ! parent panel
        cen_parent(:,iwc) = pos_p
        mag_parent(iwc) = norm2(partvec)
        if(mag_parent(iwc) .gt. 1.0e-13_wp) then
          dir_parent(:,iwc) = partvec/mag_parent(iwc)
        else
          dir_parent(:,iwc) = partvec
        endif
        iwc = iwc +1 ! added parent panel
       
        ! parent panel from previous row
        !Left side
        dir = wake%pan_w_points(:,p1,wake%nmax_pan)-wake%pan_w_points(:,p1,wake%nmax_pan+1)
        if (wake%pan_neigh(1,iw) .gt. 0) then
          ave = wake%wake_panels(iw,wake%nmax_pan)%mag- &
                real(wake%pan_neigh_o(1,iw),wp)* &
                wake%wake_panels(wake%pan_neigh(1,iw),wake%nmax_pan)%mag
          ave = ave/2.0_wp
        else !has no fixed neighbour
            ! TODO join_te?
            ave = wake%wake_panels(iw,wake%nmax_pan)%mag
        endif
        partvec = partvec + dir*ave
  
        !Right side
        dir = -wake%pan_w_points(:,p2,wake%nmax_pan)+wake%pan_w_points(:,p2,wake%nmax_pan+1)
        if (wake%pan_neigh(2,iw) .gt. 0) then
          ave = wake%wake_panels(iw,wake%nmax_pan)%mag- &
                real(wake%pan_neigh_o(2,iw),wp)* &
                wake%wake_panels(wake%pan_neigh(2,iw),wake%nmax_pan)%mag
          ave = ave/2.0_wp
        else
            ave = wake%wake_panels(iw,wake%nmax_pan)%mag
        endif
        partvec = partvec + dir*ave
  
        !End side
        dir = wake%pan_w_points(:,p2,wake%nmax_pan+1)-wake%pan_w_points(:,p2,wake%nmax_pan+1)
        ave = wake%wake_panels(iw,wake%nmax_pan)%mag-wake%end_pan_idou(iw) ! TODO check if it's not identically 0 
        partvec = partvec + dir*ave

        cen_parent(:,iwc) = wake%wake_panels(iw,wake%nmax_pan)%cen
        mag_parent(iwc) = norm2(partvec)
        if(mag_parent(iwc) .gt. 1.0e-13_wp) then
          dir_parent(:,iwc) = partvec/mag_parent(iwc)
        else
          dir_parent(:,iwc) = partvec
        endif
        
        ! values for subparts
        if ( chord_side_len .ge. span_side_len) then
            ! chord side is longer
            max_side = chord_side_len
            min_side = span_side_len
            n_chord = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
            n_span = wake%k_refine
            step_chord = max_side/real(n_chord,wp)
            step_span = min_side/real(n_span,wp)
        else
            ! span side is longer
            max_side = span_side_len
            min_side = chord_side_len
            n_chord = wake%k_refine
            n_span = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
            step_chord = min_side/real(n_chord,wp)
            step_span = max_side/real(n_span,wp)
        end if
        
        do ic = 1, n_chord
          do is = 1, n_span   
            cen_sbprt(:,isp+(ic-1)*n_span+is) = wake%pan_w_points(:,p1,wake%nmax_pan+1)+& ! top left corner
                    (real(ic,wp)-0.5_wp)/real(n_chord,wp)*& ! move chordwise
                    (points_end(:,p1)-wake%pan_w_points(:,p1,wake%nmax_pan+1))+& 
                    (real(is,wp)-0.5_wp)/real(n_span,wp)*& ! move spanwise
                    (wake%pan_w_points(:,p2,wake%nmax_pan+1)-wake%pan_w_points(:,p1,wake%nmax_pan+1)) 
          enddo
        enddo
        area_sbprt(1:n_chord*n_span) = step_chord*step_span ! TODO refine it
        !mag_sbprt(1:n_chord*n_span) = mag_parent(iwc)
        do ic = 1,3
          dir_sbprt(ic,1:n_chord*n_span) = dir_parent(ic,iwc)
        enddo
        isp = isp + n_chord*n_span ! added subparticles

        iwc = iwc+1 ! added parent panel
               
        ! loop over all the other panels until component ends
        do while (wake%pan_neigh(1,iw) .gt. 0)
          iw = iw + 1 ! move to next panel
          
          p1 = wake%i_start_points(1,iw)
          p2 = wake%i_start_points(2,iw)
          partvec = 0.0_wp
          !Left side
          dir = wake%pan_w_points(:,p1,wake%nmax_pan+1)-points_end(:,p1)
          chord_side_len = norm2(dir)
          if (wake%pan_neigh(1,iw) .gt. 0) then
            ave = wake%end_pan_idou(iw) - &
                  real(wake%pan_neigh_o(1,iw),wp)* &
                  wake%end_pan_idou(wake%pan_neigh(1,iw))
            ave = ave/2.0_wp
          else !has no fixed neighbour
            if(sim_param%join_te) then
              ave = get_joined_intensity(wake, iw, 1)
            else
              ave = wake%end_pan_idou(iw)
            endif
          endif
          partvec = partvec + dir*ave
    
          !Right side
          dir = -wake%pan_w_points(:,p2,wake%nmax_pan+1)+points_end(:,p2)
          if (wake%pan_neigh(2,iw) .gt. 0) then
            ave = wake%end_pan_idou(iw) - &
                  real(wake%pan_neigh_o(2,iw),wp)* &
                  wake%end_pan_idou(wake%pan_neigh(2,iw))
            ave = ave/2.0_wp
          else
            if(sim_param%join_te) then
              ave = get_joined_intensity(wake, iw, 2)
            else
              ave = wake%end_pan_idou(iw)
            endif
          endif
          partvec = partvec + dir*ave
    
          !End side
          dir = points_end(:,p1) - points_end(:,p2)
          span_side_len = norm2(dir)
          ave = wake%end_pan_idou(iw)-wake%last_pan_idou(iw)
          wake%last_pan_idou(iw) = wake%end_pan_idou(iw)
          partvec = partvec + dir*ave
  
          ! centre of parent panel
          pos_p = (points_end(:,p1)+points_end(:,p2)+ &
                  wake%pan_w_points(:,p1,wake%nmax_pan+1) + &
                  wake%pan_w_points(:,p2,wake%nmax_pan+1) )/4.0_wp
                           
          ! parent panel
          cen_parent(:,iwc) = pos_p
          mag_parent(iwc) = norm2(partvec)
          if(mag_parent(iwc) .gt. 1.0e-13_wp) then
            dir_parent(:,iwc) = partvec/mag_parent(iwc)
          else
            dir_parent(:,iwc) = partvec
          endif
                    
          ! values for subparts
          if ( chord_side_len .ge. span_side_len) then
              ! chord side is longer
              max_side = chord_side_len
              min_side = span_side_len
              n_chord = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
              n_span = wake%k_refine
              step_chord = max_side/real(n_chord,wp)
              step_span = min_side/real(n_span,wp)
          else
              ! span side is longer
              max_side = span_side_len
              min_side = chord_side_len
              n_chord = wake%k_refine
              n_span = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
              step_chord = min_side/real(n_chord,wp)
              step_span = max_side/real(n_span,wp)
          end if
       
          do ic = 1, n_chord
            do is = 1, n_span   
              cen_sbprt(:,isp+(ic-1)*n_span+is) = wake%pan_w_points(:,p1,wake%nmax_pan+1)+& ! top left corner
                      (real(ic,wp)-0.5_wp)/real(n_chord,wp)*& ! move chordwise
                      (points_end(:,p1)-wake%pan_w_points(:,p1,wake%nmax_pan+1))+& 
                      (real(is,wp)-0.5_wp)/real(n_span,wp)*& ! move spanwise
                      (wake%pan_w_points(:,p2,wake%nmax_pan+1)-wake%pan_w_points(:,p1,wake%nmax_pan+1)) 
            enddo
          enddo
          ! in the following, note that +1 is because isp started at 0 (cfr isp+ic above)
          area_sbprt(isp+1:isp+n_chord*n_span) = step_chord*step_span ! TODO refine it
          !mag_sbprt(isp+1:isp+n_chord*n_span) = mag_parent(iwc)
          do ic = 1,3
            dir_sbprt(ic,isp+1:isp+n_chord*n_span) = dir_parent(ic,iwc)
          enddo
          isp = isp + n_chord*n_span ! added subparticles
          
          iwc = iwc+1 ! added parent panel
           
          ! parent panel from previous row
          !Left side
          dir = wake%pan_w_points(:,p1,wake%nmax_pan)-wake%pan_w_points(:,p1,wake%nmax_pan+1)
          if (wake%pan_neigh(1,iw) .gt. 0) then
            ave = wake%wake_panels(iw,wake%nmax_pan)%mag- &
                  real(wake%pan_neigh_o(1,iw),wp)* &
                  wake%wake_panels(wake%pan_neigh(1,iw),wake%nmax_pan)%mag
            ave = ave/2.0_wp
          else !has no fixed neighbour
              ! TODO join_te?
              ave = wake%wake_panels(iw,wake%nmax_pan)%mag
          endif
          partvec = partvec + dir*ave
    
          !Right side
          dir = -wake%pan_w_points(:,p2,wake%nmax_pan)+wake%pan_w_points(:,p2,wake%nmax_pan+1)
          if (wake%pan_neigh(2,iw) .gt. 0) then
            ave = wake%wake_panels(iw,wake%nmax_pan)%mag- &
                  real(wake%pan_neigh_o(2,iw),wp)* &
                  wake%wake_panels(wake%pan_neigh(2,iw),wake%nmax_pan)%mag
            ave = ave/2.0_wp
          else
              ave = wake%wake_panels(iw,wake%nmax_pan)%mag
          endif
          partvec = partvec + dir*ave
    
          !End side
          dir = wake%pan_w_points(:,p2,wake%nmax_pan+1)-wake%pan_w_points(:,p2,wake%nmax_pan+1)
          ave = wake%wake_panels(iw,wake%nmax_pan)%mag-wake%end_pan_idou(iw) ! TODO check if it's not identically 0 
          partvec = partvec + dir*ave
  
          cen_parent(:,iwc) = wake%wake_panels(iw,wake%nmax_pan)%cen
          mag_parent(iwc) = norm2(partvec)
          if(mag_parent(iwc) .gt. 1.0e-13_wp) then
            dir_parent(:,iwc) = partvec/mag_parent(iwc)
          else
            dir_parent(:,iwc) = partvec
          endif
          iwc = iwc +1 ! added parent panel
           
        enddo ! loop over panels for this component
        
        ! add final ghost panel
        p1 = wake%i_start_points(1,iw-1)
        p2 = wake%i_start_points(2,iw-1)

        cen_parent(:,iwc) = pos_p+(wake%pan_w_points(:,p1,wake%nmax_pan+1)-&
                      wake%pan_w_points(:,p2,wake%nmax_pan+1))
        mag_parent(iwc) = 0.0_wp
        dir_parent(:,iwc) = 0.0_wp      
        iwc = iwc + 1 ! added ghost panel 
    
        ! ghost panel from previous row
        cen_parent(:,iwc) = wake%wake_panels(iw,wake%nmax_pan)%cen+&
               (wake%pan_w_points(:,p1,wake%nmax_pan)-wake%pan_w_points(:,p2,wake%nmax_pan)) 
        mag_parent(iwc) = 0.0_wp
        dir_parent(:,iwc) = 0.0_wp

        ! perform interpolation, only passing actual number of parent panels and subparts
        call infinite_plate_spline(cen_sbprt(:,1:isp), cen_parent(:,1:iwc), W)
        
        ! array of weights
        allocate(w_i(isp))

        ! mag
        w_i = matmul(W,mag_parent)/maxval(mag_parent)
        mag_sbprt = w_i*sum(mag_parent)/sum(w_i)
        
        ! dir ! TODO check normalization
        ! will be used below in the insertion loop
        !w_i = matmul(W,norm2(dir_sbprt,1))/max(norm2(dir_sbprt,1))
        
        deallocate(W, w_i)        
        
        ! actually insert particles
        do ic = 1,isp
        
          pos_p = cen_sbprt(:,ic)
          
          if(all(pos_p .ge. wake%part_box_min) .and. &
              all(pos_p .le. wake%part_box_max)) then
    
            do ip = k, size(wake%wake_parts)
              if (wake%wake_parts(ip)%free) then
                wake%wake_parts(ip)%free = .false.
                k = ip+1
                wake%n_prt = wake%n_prt+1
                wake%wake_parts(ip)%mag = mag_sbprt(ic)
                wake%wake_parts(ip)%dir = dir_sbprt(:,ic)  
                wake%wake_parts(ip)%cen = pos_p
              if (sim_param%KVortexRad .ge. 1e-10_wp) then ! Variable vortex rad
                wake%wake_parts(ip)%r_Vortex = sim_param%KVortexRad*&
                              sqrt(2.0_wp*area_sbprt(ic)) ! k*radius of the circumscribed circle
                else ! revert to sim_param vortex rad
                  wake%wake_parts(ip)%r_Vortex = sim_param%VortexRad
                end if
                wake%wake_parts(ip)%r_cutoff  = sim_param%CutoffRad
                wake%wake_parts(ip)%vel = 0.5_wp * &
                                  ( wake%pan_w_vel(:,p1,wake%nmax_pan+1) + &
                                    wake%pan_w_vel(:,p2,wake%nmax_pan+1) )
                exit
              endif
            enddo
      
            if (ip .gt. wake%nmax_prt) then
              write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
                wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
                &restarting with a higher number of maximum wake particles'
      
              call error(this_sub_name, this_mod_name, trim(msg))
            endif !max number of particles
          endif !inside the box
        enddo ! insert particles
      
      iw = iw + 1 ! move to next panel (first of next component)
      
      enddo ! iw, finished all wake
      
      deallocate(cen_sbprt, area_sbprt, mag_sbprt, dir_sbprt)
      deallocate(cen_parent, mag_parent, dir_parent)
    
    else ! old behaviour, possibly with refined wake
      k = 1
      do iw = 1,wake%n_pan_stripes
        p1 = wake%i_start_points(1,iw)
        p2 = wake%i_start_points(2,iw)
        partvec = 0.0_wp
        !Left side
        dir = wake%pan_w_points(:,p1,wake%nmax_pan+1)-points_end(:,p1)
        chord_side_len = norm2(dir)
        if (wake%pan_neigh(1,iw) .gt. 0) then
          ave = wake%end_pan_idou(iw) - &
                real(wake%pan_neigh_o(1,iw),wp)* &
                wake%end_pan_idou(wake%pan_neigh(1,iw))
          ave = ave/2.0_wp
        else !has no fixed neighbour
          if(sim_param%join_te) then
            ave = get_joined_intensity(wake, iw, 1)
          else
            ave = wake%end_pan_idou(iw)
          endif
        endif
        partvec = partvec + dir*ave
  
        !Right side
        dir = -wake%pan_w_points(:,p2,wake%nmax_pan+1)+points_end(:,p2)
        if (wake%pan_neigh(2,iw) .gt. 0) then
          ave = wake%end_pan_idou(iw) - &
                real(wake%pan_neigh_o(2,iw),wp)* &
                wake%end_pan_idou(wake%pan_neigh(2,iw))
          ave = ave/2.0_wp
        else
          if(sim_param%join_te) then
            ave = get_joined_intensity(wake, iw, 2)
          else
            ave = wake%end_pan_idou(iw)
          endif
        endif
        partvec = partvec + dir*ave
  
        !End side
        dir = points_end(:,p1) - points_end(:,p2)
        span_side_len = norm2(dir)
        ave = wake%end_pan_idou(iw)-wake%last_pan_idou(iw)
        wake%last_pan_idou(iw) = wake%end_pan_idou(iw)
        partvec = partvec + dir*ave
  
      if (wake%refine_wake) then
        
        if ( chord_side_len .ge. span_side_len) then
            ! chord side is longer
            max_side = chord_side_len
            min_side = span_side_len
            n_chord = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
            n_span = wake%k_refine
            step_chord = max_side/real(n_chord,wp)
            step_span = min_side/real(n_span,wp)
  
        else
            ! span side is longer
            max_side = span_side_len
            min_side = chord_side_len
            n_chord = wake%k_refine
            n_span = ceiling(max_side/(min_side/real(wake%k_refine,wp)))
            step_chord = min_side/real(n_chord,wp)
            step_span = max_side/real(n_span,wp)
            
        end if
  
        do ic = 1, n_chord
          do is = 1, n_span
            ! find centre of the subpanel by starting from top left of the panel and
            ! move half steps chord- and span-wise
            pos_p = wake%pan_w_points(:,p1,wake%nmax_pan+1)+& ! top left corner
                    (real(ic,wp)-0.5_wp)/real(n_chord,wp)*& ! move chordwise
                    (points_end(:,p1)-wake%pan_w_points(:,p1,wake%nmax_pan+1))+& 
                    (real(is,wp)-0.5_wp)/real(n_span,wp)*& ! move spanwise
                    (wake%pan_w_points(:,p2,wake%nmax_pan+1)-wake%pan_w_points(:,p1,wake%nmax_pan+1)) 
  
            area = step_chord*step_span ! TODO refine it
  
            !Add the particle (if it is in the box)
            if(all(pos_p .ge. wake%part_box_min) .and. &
                all(pos_p .le. wake%part_box_max)) then
      
              do ip = k, size(wake%wake_parts)
                if (wake%wake_parts(ip)%free) then
                  wake%wake_parts(ip)%free = .false.
                  k = ip+1
                  wake%n_prt = wake%n_prt+1
                  ! mag of the particle is 1/n_subpart * mag_panel
                  wake%wake_parts(ip)%mag = norm2(partvec)/real(n_span*n_chord,wp) ! TODO generalize n_subpart
      
                  if(norm2(partvec) .gt. 1.0e-13_wp) then
                    wake%wake_parts(ip)%dir = partvec/norm2(partvec)
                  else
                    wake%wake_parts(ip)%dir = partvec
                  endif
      
                  wake%wake_parts(ip)%cen = pos_p
                if (sim_param%KVortexRad .ge. 1e-10_wp) then ! Variable vortex rad
                  wake%wake_parts(ip)%r_Vortex = sim_param%KVortexRad*sqrt(2.0_wp*area) ! k*radius of the circumscribed circle
                  else ! revert to sim_param vortex rad
                    wake%wake_parts(ip)%r_Vortex = sim_param%VortexRad
                  end if
                  wake%wake_parts(ip)%r_cutoff  = sim_param%CutoffRad
                  wake%wake_parts(ip)%vel = 0.5_wp * &
                                    ( wake%pan_w_vel(:,p1,wake%nmax_pan+1) + &
                                      wake%pan_w_vel(:,p2,wake%nmax_pan+1) )
                  exit
                endif
              enddo
        
              if (ip .gt. wake%nmax_prt) then
                write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
                  wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
                  &restarting with a higher number of maximum wake particles'
        
                call error(this_sub_name, this_mod_name, trim(msg))
              endif !max number of particles
            endif !inside the box
          enddo !is
        enddo !ic
        
      else !classical behaviour
        !Calculate the center
        pos_p = (points_end(:,p1)+points_end(:,p2)+ &
                wake%pan_w_points(:,p1,wake%nmax_pan+1) + &
                wake%pan_w_points(:,p2,wake%nmax_pan+1) )/4.0_wp
        ! A = 1/2 [()        
        area = norm2(cross(points_end(:,p1)- wake%pan_w_points(:,p2,wake%nmax_pan+1),&
                     points_end(:,p2)-wake%pan_w_points(:,p2,wake%nmax_pan+1)))
        !Add the particle (if it is in the box)
        if(all(pos_p .ge. wake%part_box_min) .and. &
            all(pos_p .le. wake%part_box_max)) then
  
          do ip = k, size(wake%wake_parts)
            if (wake%wake_parts(ip)%free) then
              wake%wake_parts(ip)%free = .false.
              k = ip+1
              wake%n_prt = wake%n_prt+1
              wake%wake_parts(ip)%mag = norm2(partvec)
  
              if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
                wake%wake_parts(ip)%dir = partvec/wake%wake_parts(ip)%mag
              else
                wake%wake_parts(ip)%dir = partvec
              endif
  
              wake%wake_parts(ip)%cen = pos_p
            if (sim_param%KVortexRad .ge. 1e-10_wp) then ! Variable vortex rad
              wake%wake_parts(ip)%r_Vortex = sim_param%KVortexRad*sqrt(2.0_wp*area) ! k*radius of the circumscribed circle
              else ! revert to sim_param vortex rad
                wake%wake_parts(ip)%r_Vortex = sim_param%VortexRad
              end if
              wake%wake_parts(ip)%r_cutoff  = sim_param%CutoffRad
              wake%wake_parts(ip)%vel = 0.5_wp * &
                                ( wake%pan_w_vel(:,p1,wake%nmax_pan+1) + &
                                  wake%pan_w_vel(:,p2,wake%nmax_pan+1) )
              exit
            endif
          enddo
    
          if (ip .gt. wake%nmax_prt) then
            write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
              wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
              &restarting with a higher number of maximum wake particles'
    
            call error(this_sub_name, this_mod_name, trim(msg))
          endif !max number of particles
        endif !inside the box
        
      endif !refine_wake
      enddo !iw
    endif !interp_parts
  endif !full wake

  if(wake%full_rings) then
    k = 1
    nprev = 0
    do id = 1,wake%ndisks
      do is = 1,wake%wake_rings(id,wake%rin_wake_len)%n_ver
        p1 = is + nprev
        p2 = nprev + 1 + mod(is,wake%wake_rings(id,wake%rin_wake_len)%n_ver)
        partvec = 0.0_wp

        dir = points_end_ring(:,p2)-points_end_ring(:,p1)
        ave = wake%wake_rings(id,wake%rin_wake_len)%mag
        partvec =dir*ave

        !Calculate the center
        pos_p = (points_end_ring(:,p1)+points_end_ring(:,p2))/2.0_wp

        !Add the particle
        if(all(pos_p .ge. wake%part_box_min) .and. &
            all(pos_p .le. wake%part_box_max)) then

          do ip = k, size(wake%wake_parts)

            if (wake%wake_parts(ip)%free) then

              wake%wake_parts(ip)%free = .false.
              k = ip+1
              wake%n_prt = wake%n_prt+1
              wake%wake_parts(ip)%mag = norm2(partvec)

              if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
                wake%wake_parts(ip)%dir = partvec/wake%wake_parts(ip)%mag
              else
                wake%wake_parts(ip)%dir = partvec
              endif

              wake%wake_parts(ip)%cen = pos_p

              exit
            endif
          enddo  !ip

          if (ip .gt. wake%nmax_prt) then
            write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
              wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
              &restarting with a higher number of maximum wake particles'
          call error(this_sub_name, this_mod_name, trim(msg))
          endif
        endif
      enddo !is
      nprev = nprev + wake%wake_rings(id,wake%rin_wake_len)%n_ver
    enddo !id
  endif


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Treat flow separations
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(sim_param%use_ve) then
    do i_comp = 1 , size( geo%components)! ***** loop #1 over components *****

      if ( k .lt. 1 ) k = 1

      ! flow separation allowed only for surfpan elements -----
        if ( trim( geo%components(i_comp)%comp_el_type ) .eq. 'p' ) then

          n_elem = size( geo%components(i_comp)%el )

          do i_elem = 1 , n_elem     ! ***** loop #2 over elements   *****

            select type( el => geo%components(i_comp)%el(i_elem) )
            type is (t_surfpan)

            ! flow separation
            if ( el % al_free .gt. 0.0_wp ) then
              pos_p = el%cen + el%nor * el%h_bl + &
                       el % surf_vel * sim_param%dt*real(sim_param%ndt_update_wake,wp)

              if(all(pos_p .ge. wake%part_box_min) .and. &
                all(pos_p .le. wake%part_box_max)) then
                !Add the particle
                do ip = k, size(wake%wake_parts)

                  if (wake%wake_parts(ip)%free) then

                    wake%wake_parts(ip)%free = .false.
                    k = ip+1
                    wake%n_prt = wake%n_prt+1
                    wake%wake_parts(ip)%mag = norm2(el%free_vort)
                    if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
                      wake%wake_parts(ip)%dir = el%free_vort/wake%wake_parts(ip)%mag
                    else
                      wake%wake_parts(ip)%dir = el%free_vort
                    endif
                    !wake%wake_parts(ip)%cen = el%cen + el%nor * el%h_bl + &
                    !       el % surf_vel * sim_param%dt
                    wake%wake_parts(ip)%cen = pos_p
                    wake%wake_parts(ip)%vel = el % surf_vel
                    exit
                  endif
                enddo
                if (ip .gt. wake%nmax_prt) then
                  write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
                    wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
                    &restarting with a higher number of maximum wake particles'
                  call error(this_sub_name, this_mod_name, trim(msg))
                endif

              endif !particle in box
            end if !generate the particle
          end select!select panels
        end do ! ***** loop #2 over elements   *****
      end if !if is a 3dp
    end do ! ***** loop #1 over components *****
  endif !calculate the viscous particles detachment


  !If the wake is full, attach the end vortex
  if (wake%full_panels) then
    do iw = 1,wake%n_pan_stripes
      !wake%end_vorts(iw)%mag => wake%wake_panels(iw,wake%pan_wake_len)%mag
      wake%end_vorts(iw)%mag => wake%last_pan_idou(iw)
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%end_vorts(iw)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,wake%pan_wake_len+1),  &
                    wake%pan_w_points(:,p2,wake%pan_wake_len+1)/), (/3,2/)))
      wake%end_vorts(iw)%ver_vel(:,1) = wake%pan_w_vel(:,p1,wake%pan_wake_len+1)
      wake%end_vorts(iw)%ver_vel(:,2) = wake%pan_w_vel(:,p2,wake%pan_wake_len+1)
    enddo
  endif

end subroutine complete_wake

!----------------------------------------------------------------------

subroutine compute_vel_from_all(elems, wake, pos, vel)
  type(t_pot_elem_p), intent(in)  :: elems(:)
  type(t_wake), intent(in)        :: wake
  real(wp), intent(in)            :: pos(3)
  real(wp), intent(out)           :: vel(3)

  integer                         :: ie
  real(wp)                        :: v(3)

  vel = 0.0_wp

  !calculate the influence of the solid bodies
  do ie=1,size(elems)
    call elems(ie)%p%compute_vel(pos, v)
    vel = vel + v/(4*pi)
  enddo

  ! calculate the influence of the wake panels
  do ie=1,size(wake%pan_p)
    call wake%pan_p(ie)%p%compute_vel(pos, v)
    vel = vel + v/(4*pi)
  enddo

  ! calculate the influence of the wake rings
  do ie=1,size(wake%rin_p)
    call wake%rin_p(ie)%p%compute_vel(pos, v)
    vel = vel + v/(4*pi)
  enddo

  !calculate the influence of the end vortex
  do ie=1,size(wake%end_vorts)
    call wake%end_vorts(ie)%compute_vel(pos, v)
    vel = vel + v/(4*pi)
  enddo

  !calculate the influence of particles
  do ie=1,size(wake%part_p)
    call wake%part_p(ie)%p%compute_vel(pos, v)
    vel = vel+ v/(4*pi)
  enddo

end subroutine compute_vel_from_all

!----------------------------------------------------------------------

subroutine get_vel_free(this, elems, wake, pos, vel_hcas, vel)
  class(t_free_wake)                    :: this
  type(t_pot_elem_p), intent(in)        :: elems(:)
  type(t_wake), intent(in)              :: wake
  real(wp), intent(in)                  :: pos(3)
  real(wp), intent(in)                  :: vel_hcas(3)
  real(wp), intent(out)                 :: vel(3)

  call compute_vel_from_all(elems, wake, pos, vel)

  !vel = vel + sim_param%u_inf
  vel = vel + variable_wind(pos, sim_param%time)

  if (sim_param%HCAS) vel = vel + vel_hcas

end subroutine get_vel_free

!----------------------------------------------------------------------

subroutine get_vel_rigid(this, elems, wake, pos, vel_hcas, vel)
  class(t_rigid_wake)                      :: this
  type(t_pot_elem_p), intent(in)           :: elems(:)
  type(t_wake), intent(in)                 :: wake
  real(wp), intent(in)                     :: pos(3)
  real(wp), intent(in)                     :: vel_hcas(3)
  real(wp), intent(out)                    :: vel(3)

  vel = sim_param%rigid_wake_vel

end subroutine get_vel_rigid

!----------------------------------------------------------------------

function get_vel_hcas() result(vel_hcas)
  real(wp)                  :: vel_hcas(3)

  real(wp)                  :: hcas_reltime

  hcas_reltime = (sim_param%time-sim_param%t0)/sim_param%hcas_time
  vel_hcas = sim_param%hcas_vel*max(0.0_wp, (1.0_wp-hcas_reltime)) !linear reduction

end function get_vel_hcas

!----------------------------------------------------------------------

subroutine join_first_panels(wake, te_fact)
  type(t_wake), intent(inout) :: wake
  real(wp)    , intent(in)    :: te_fact

  real(wp)                    :: lmin
  integer                     :: i, j, iw, ir , i_case
  real(wp)                    :: sp1(3), sp2(3), sp1_1(3), sp2_1(3), l1, l2, tol, pos_p(3)

!integer :: n_joined_te

  do i = 1,size(wake%pan_i_ends)
    iw = wake%pan_i_ends(i)
    sp1 = wake%w_start_points(:,wake%i_start_points(1,iw))
    sp2 = wake%w_start_points(:,wake%i_start_points(2,iw))
    l1 = norm2(sp1-sp2)
    do j = i+1,size(wake%pan_i_ends)
      ir = wake%pan_i_ends(j)
      sp1_1 = wake%w_start_points(:,wake%i_start_points(1,ir))
      sp2_1 = wake%w_start_points(:,wake%i_start_points(2,ir))
      l2 = norm2(sp1_1-sp2_1)
      tol = te_fact*min(l1,l2)

      !> check which is the pair of nodes to be joined, very rude approach
      i_case = 1
      lmin = norm2(sp1_1-sp1)
      if ( norm2(sp2_1-sp1) .lt. lmin ) then
        lmin = norm2(sp2_1-sp1) ; i_case = 2
      end if
      if ( norm2(sp1_1-sp2) .lt. lmin ) then
        lmin = norm2(sp1_1-sp2) ; i_case = 3
      end if
      if ( norm2(sp2_1-sp2) .lt. lmin ) then
        lmin = norm2(sp2_1-sp2) ; i_case = 4
      end if

      !checking only one side is enough, it is the opposite side for the other
      !end
      if ( i_case .eq. 1 ) then
        if( norm2(sp1-sp1_1).lt.tol) then

          pos_p = (sp1+sp1_1)/2.0_wp
          
          wake%pan_w_points(:,wake%i_start_points(1,iw),1) = pos_p
          wake%pan_w_points(:,wake%i_start_points(1,ir),1) = pos_p

          pos_p = (wake%pan_w_points(:,wake%i_start_points(1,iw),2) + &
                  wake%pan_w_points(:,wake%i_start_points(1,ir),2))/2.0_wp
          
          wake%pan_w_points(:,wake%i_start_points(1,iw),2) = pos_p
          wake%pan_w_points(:,wake%i_start_points(1,ir),2) = pos_p

          wake%joined_tes(1,i,1) = ir
          wake%joined_tes(1,j,1) = iw

        end if
      elseif ( i_case .eq. 2 ) then
        if( norm2(sp1-sp2_1).lt.tol) then

          pos_p = (sp1+sp2_1)/2.0_wp

          wake%pan_w_points(:,wake%i_start_points(1,iw),1) = pos_p
          wake%pan_w_points(:,wake%i_start_points(2,ir),1) = pos_p

          pos_p = (wake%pan_w_points(:,wake%i_start_points(1,iw),2) + &
                  wake%pan_w_points(:,wake%i_start_points(2,ir),2))/2.0_wp

          wake%pan_w_points(:,wake%i_start_points(1,iw),2) = pos_p
          wake%pan_w_points(:,wake%i_start_points(2,ir),2) = pos_p

          wake%joined_tes(1,i,1) = ir
          wake%joined_tes(2,j,1) = iw

        end if

      elseif ( i_case .eq. 3 ) then
        if( norm2(sp2-sp1_1).lt.tol) then

          pos_p = (sp2+sp1_1)/2.0_wp
          wake%pan_w_points(:,wake%i_start_points(2,iw),1) = pos_p
          wake%pan_w_points(:,wake%i_start_points(1,ir),1) = pos_p

          pos_p = (wake%pan_w_points(:,wake%i_start_points(2,iw),2) + &
                  wake%pan_w_points(:,wake%i_start_points(1,ir),2))/2.0_wp
          
          wake%pan_w_points(:,wake%i_start_points(2,iw),2) = pos_p
          wake%pan_w_points(:,wake%i_start_points(1,ir),2) = pos_p

          wake%joined_tes(2,i,1) = ir
          wake%joined_tes(1,j,1) = iw

        end if
      elseif ( i_case .eq. 4 ) then
        if( norm2(sp2-sp2_1).lt.tol) then

          pos_p = (sp2+sp2_1)/2.0_wp
          wake%pan_w_points(:,wake%i_start_points(2,iw),1) = pos_p
          wake%pan_w_points(:,wake%i_start_points(2,ir),1) = pos_p

          pos_p = (wake%pan_w_points(:,wake%i_start_points(2,iw),2) + &
                  wake%pan_w_points(:,wake%i_start_points(2,ir),2))/2.0_wp

          wake%pan_w_points(:,wake%i_start_points(2,iw),2) = pos_p
          wake%pan_w_points(:,wake%i_start_points(2,ir),2) = pos_p

          wake%joined_tes(2,i,1) = ir
          wake%joined_tes(2,j,1) = iw

        end if
      endif
    enddo
  enddo

end subroutine join_first_panels

!----------------------------------------------------------------------

function get_joined_intensity(wake, iw, side) result(ave)
  type(t_wake), intent(in)            :: wake
  integer, intent(in)                 :: iw
  integer, intent(in)                 :: side
  real(wp)                            :: ave

  integer                              :: iend, ineigh
  real(wp)                             :: orient

  
  do iend = 1,size(wake%pan_i_ends)
    if(wake%pan_i_ends(iend) .eq. iw) exit
  enddo


  ineigh = wake%joined_tes(side,iend,wake%nmax_pan+1)
  if (ineigh .ne. 0) then
    !check orientation
    if(sum(wake%wake_panels(iw,1)%nor * wake%wake_panels(ineigh,1)%nor) &
        .ge. 0.0_wp) then
      orient = 1.0_wp
    else
      orient = -1.0_wp
    endif

    ave = wake%wake_panels(iw,wake%pan_wake_len)%mag - &
          orient * wake%wake_panels(ineigh,wake%pan_wake_len)%mag
    ave = ave/2.0_wp
  else
    ave = wake%wake_panels(iw,wake%pan_wake_len)%mag
  endif

end function get_joined_intensity


subroutine avoid_collision(elems, wake, part, vel_in, vel_out)
  type(t_pot_elem_p), intent(in)      :: elems(:)
  type(t_wake), intent(inout)         :: wake
  type(t_vortpart), intent(inout)     :: part
  real(wp), intent(in)                :: vel_in(3)
  real(wp), intent(out)               :: vel_out(3)

  real(wp)                            :: vel(3)

  integer    :: ie
  real(wp)   :: dist1(3), dist2(3), n(3)
  real(wp)   :: pos1(3), pos2(3), relvel(3), newvel(3), nveldiff
  real(wp)   :: distn, dist1_nor, dist1_tan, dist2_nor, dist2_tan
  real(wp)   :: normvel, normvel_corr, tanvel(3), dt_part, tanvel2
  real(wp)   :: check_radius, rad_mult, elrad_mult, r2d2, elrad, blthick
  real(wp)   :: dampf

  !Multiplication factor for the check radius
  rad_mult = sim_param%pa_rad_mult

  !Multiplication factor for the element radius
  elrad_mult = sim_param%pa_elrad_mult

  ! Thickness of the "surface blocks"
  ! ( TODO: evaluate the 'phisical' boundary layer thickness )
  blthick = part%r_Vortex

  r2d2 = sqrt(2.0_wp)/2.0_wp

  !starting position
  pos1 = part%cen
  vel = vel_in

  !Cycle on all the elements
  do ie=1,size(elems)

    select type( elem => elems(ie)%p )

      class is (t_surfpan)
        !Get the position of the particle with respect to the element
        dist1 = pos1-(elem%cen)
        elrad = maxval(elem%edge_len)*r2d2
        check_radius = sim_param%dt*sim_param%u_ref*rad_mult + elrad
        distn = norm2(dist1)

        !if it is in the check radius perform calculations
        if ((distn .lt. check_radius) ) then

          !use the relative velocity to take into account also the element
          !movement
          relvel = vel - elem%ub
          n = elem%nor
          dist1_nor = sum(dist1 * n)
          dist1_tan = norm2(dist1-(n*dist1_nor))
          normvel = sum(relvel*n)

          !If it is on the opposite side of the element, or it is moving away,
          !do not perform any modification
          if (dist1_nor .lt. 0.0_wp .or. normvel.ge.0.0_wp) cycle

            !get the time at which the particle hits the surface
            !make sure not to go back in time
            dt_part = min(max((dist1_nor)/(-normvel),0.0_wp), sim_param%dt)
            !if it is lower than the timestep (i.e. the particles hits the surface
            !within next step) start the correction
            !Get the position at the time in which the particle hits the
            !surface plane
            pos2  = part%cen+relvel*dt_part
            dist2 = pos2 - ( elem%cen )
            dist2_nor = sum(dist2 * n)
            
          if(dist2_nor .le. blthick) then
            dist2_tan = norm2(dist2-(n*dist2_nor))

            if (dist2_tan .lt. elrad_mult*elrad) then

              !correct the normal velocity to avoid penetration
              if (dt_part .lt. sim_param%dt) then
                normvel_corr = -dist1_nor/ &
                                (sim_param%dt*real(sim_param%ndt_update_wake,wp))
              else
                dampf = (dist2_nor/blthick)
                normvel_corr = normvel + (blthick-dist2_nor)/ &
                              (sim_param%dt*real(sim_param%ndt_update_wake,wp)) * dampf
              endif

              ! should be
              ! vel = relvel + (normvel_corr-normvel)*n + elem%ub
              ! but simplifying
              newvel = vel + (normvel_corr-normvel)*n

              !fix the tangential velocity to keep the magnitude constant
              normvel_corr = sum(newvel*n)
              normvel  = sum(vel*n)
              tanvel   = vel - normvel*n
              tanvel2  = sum(tanvel**2)
              nveldiff = normvel**2-normvel_corr**2+tanvel2
              if(tanvel2.gt.0.0_wp .and. normvel_corr*normvel.ge.0.0_wp .and. &
                                                      nveldiff.ge.0.0_wp) then
                tanvel = tanvel * (sqrt(nveldiff)/sqrt(tanvel2) )
              endif
              !reconstruct the velocity
              vel = tanvel + normvel_corr*n
            endif
          endif

        endif

      class is ( t_vortlatt )
        !Get the position of the particle with respect to the element
        dist1 = pos1-(elem%cen)
        elrad = maxval(elem%edge_len)*r2d2
        check_radius = sim_param%dt*sim_param%u_ref*rad_mult + elrad
        distn = norm2(dist1)

        !if it is in the check radius perform calculations
        if ((distn .lt. check_radius) ) then

          !use the relative velocity to take into account also the element
          !movement
          relvel = vel - elem%ub
          n = elem%nor
          dist1_nor = norm2(dist1 * n)
          dist1_tan = norm2(dist1-(n*dist1_nor))
          normvel = sum(relvel*n)

          !If it is moving away, do not perform any modification
          if (sum(dist1 * n)*normvel .ge. 0.0_wp) cycle

            !get the time at which the particle hits the surface
            !make sure not to go back in time
            dt_part = min(max(abs(dist1_nor)/(abs(normvel)),0.0_wp), sim_param%dt)
            !if it is lower than the timestep (i.e. the particles hits the surface
            !within next step) start the correction
            !Get the position at the time in which the particle hits the
            !surface plane
            pos2  = part%cen+relvel*dt_part
            dist2 = pos2 - ( elem%cen )
            dist2_nor = norm2(dist2 * n)
            blthick = 2.0_wp*blthick
          if(dist2_nor .le. blthick) then
            dist2_tan = norm2(dist2-(n*dist2_nor))

            if (dist2_tan .lt. elrad_mult*elrad) then

              !correct the normal velocity to avoid penetration
              if (dt_part .lt. sim_param%dt) then
                normvel_corr = -dist1_nor/ &
                                (sim_param%dt*real(sim_param%ndt_update_wake,wp))
              else
                dampf = (dist2_nor/blthick)
                normvel_corr = normvel + (blthick-dist2_nor)/ &
                              (sim_param%dt*real(sim_param%ndt_update_wake,wp)) * dampf
              endif

              newvel = vel + (normvel_corr-normvel)*n

              !fix the tangential velocity to keep the magnitude constant
              normvel_corr = sum(newvel*n)
              normvel  = sum(vel*n)
              tanvel   = vel - normvel*n
              tanvel2  = sum(tanvel**2)
              nveldiff = normvel**2-normvel_corr**2+tanvel2
              if(tanvel2.gt.0.0_wp .and. normvel_corr*normvel.ge.0.0_wp .and. &
                                                      nveldiff.ge.0.0_wp) then
                tanvel = tanvel * (sqrt(nveldiff)/sqrt(tanvel2) )
              endif
              !reconstruct the velocity
              vel = tanvel + normvel_corr*n
            endif
          endif

        endif

      class is ( t_liftlin )

      class default

    end select

  enddo

  vel_out = vel

end subroutine avoid_collision

end module mod_wake
